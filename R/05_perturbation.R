# ============================================================
# 05_perturbation.R
# Test wrażliwości przez perturbację przestrzenno-temporalną
#
# Wymaga: powiaty_spdf w pamięci (opcjonalnie — dla lambda_s)
# Wczytuje: output/{PERIOD}/stvc_gam/gam_results.rds
# Wymaga: 00_config.R
# ============================================================

if (!exists("LANCET") || !exists("YEAR_MAX"))
  source(file.path(getwd(), "R", "00_config.R"))

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr)
  library(stringr); library(tibble); library(ggplot2)
  library(patchwork); library(forcats)
  library(mgcv); library(scales)
})

set.seed(20240525)

GAM_RDS    <- file.path(OUT_BASE, "stvc_gam",
                         "gam_results.rds")
OUT_DIR    <- file.path(OUT_BASE, "perturbation")
TAB_DIR    <- file.path(OUT_DIR, "tables")
FIG_DIR    <- file.path(OUT_DIR, "figures")

B_PERTURB      <- 50
NOISE_SCALE    <- 0.50
N_GRID_PERT    <- 50
RUN_MULTISCALE <- FALSE
MULTISCALE_B   <- 30
NOISE_SCALES   <- c(0.25, 0.50, 1.00)

dir.create(TAB_DIR, recursive = TRUE,
           showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE,
           showWarnings = FALSE)

cat("=== Perturbation ===\n")
cat("Period:", PERIOD, "|", YEAR_LABEL, "\n")
cat("B =", B_PERTURB,
    "| noise_scale =", NOISE_SCALE, "\n")

if (!file.exists(GAM_RDS))
  stop("Brak: ", GAM_RDS)

# ============================================================
# 1. WCZYTAJ MODEL
# ============================================================

res          <- readRDS(GAM_RDS)
m            <- res$m
dat          <- res$dat_all
between_terms <- res$between_terms
within_terms  <- res$within_terms
all_terms     <- res$all_terms
VAR_LABELS_P  <- res$VAR_LABELS
OUTCOME_P     <- res$settings$OUTCOME
LAG_N         <- res$settings$WITHIN_LAG_N %||% 1L
GAMMA_P       <- res$settings$GAMMA %||% GAMMA

if (!"code_pow" %in% names(dat))
  dat <- dat %>%
    mutate(code_pow = make_code_pow(code))

clean_var_p <- function(t) {
  t %>% str_remove("_between_z$") %>%
    str_remove(paste0("_within_lag", LAG_N, "_z$"))
}

comp_from_term_p <- function(t) {
  case_when(
    str_detect(t, "_between_z$") ~ "between",
    str_detect(t, paste0("_within_lag",
                          LAG_N, "_z$")) ~ "within",
    TRUE ~ "other"
  )
}

label_var_p <- function(x) {
  out <- VAR_LABELS_P[x]
  out[is.na(out)] <- x[is.na(out)]
  as.character(out)
}

# ============================================================
# 2. PARAMETRY SZUMU Z RESIDUÓW
# ============================================================

cat("Szacowanie parametrów szumu...\n")

dat_r <- dat %>%
  mutate(
    fitted = as.numeric(predict(m, newdata=dat,
                                 type="response")),
    resid  = .data[[OUTCOME_P]] - fitted
  ) %>%
  arrange(code, year)

sigma_resid <- sd(dat_r$resid, na.rm = TRUE)

rho_t <- dat_r %>%
  group_by(code) %>%
  summarise(
    rho = if (n() > 2)
      suppressWarnings(
        cor(resid[-n()], resid[-1],
            use="complete.obs"))
    else NA_real_,
    .groups = "drop"
  ) %>%
  summarise(rho = median(rho, na.rm=TRUE)) %>%
  pull(rho)

if (!is.finite(rho_t)) rho_t <- 0

lambda_s <- 0

if (exists("powiaty_spdf") &&
    requireNamespace("spdep", quietly=TRUE)) {

  pow_ord <- powiaty_spdf %>%
    mutate(
      JPT_KOD_JE = as.character(JPT_KOD_JE),
      JPT_KOD_JE = ifelse(nchar(JPT_KOD_JE)==3,
                          paste0("0",JPT_KOD_JE),
                          JPT_KOD_JE)
    ) %>%
    arrange(JPT_KOD_JE)

  nb_p <- tryCatch(
    spdep::poly2nb(pow_ord, queen=TRUE),
    error=function(e) NULL
  )

  if (!is.null(nb_p)) {
    lw_p <- spdep::nb2listw(nb_p, style="W",
                              zero.policy=TRUE)
    lambda_by_yr <- map_dbl(
      sort(unique(dat_r$year)), function(yr) {
        d_yr <- dat_r %>% filter(year==yr) %>%
          arrange(code_pow)
        idx  <- match(pow_ord$JPT_KOD_JE,
                       d_yr$code_pow)
        r    <- d_yr$resid[idx]
        w_r  <- tryCatch(
          spdep::lag.listw(lw_p, r,
                            zero.policy=TRUE,
                            na.action=na.pass),
          error=function(e) rep(NA_real_,length(r))
        )
        valid <- !is.na(r) & !is.na(w_r)
        if (sum(valid) < 5) return(NA_real_)
        tryCatch(
          coef(lm(r[valid] ~ w_r[valid]))[2],
          error=function(e) NA_real_
        )
      }
    )
    lambda_s <- median(lambda_by_yr, na.rm=TRUE)
    if (!is.finite(lambda_s)) lambda_s <- 0
  }
}

cat(sprintf(
  "  sigma=%.3f | rho_t=%.3f | lambda_s=%.3f\n",
  sigma_resid, rho_t, lambda_s
))

# ============================================================
# 3. GENERATOR SZUMU
# ============================================================

gen_noise <- function(dat_in, rho_t, lambda_s,
                       sigma, scale = NOISE_SCALE,
                       W_mat = NULL,
                       pow_codes = NULL) {

  powiaty_v <- sort(unique(dat_in$code_pow))
  years_v   <- sort(unique(dat_in$year))
  n_p <- length(powiaty_v)
  n_y <- length(years_v)
  sig <- sigma * scale

  nu <- matrix(rnorm(n_p * n_y, 0, sig),
               nrow=n_p, ncol=n_y,
               dimnames=list(powiaty_v, years_v))

  eps <- matrix(0, n_p, n_y,
                dimnames=list(powiaty_v, years_v))
  eps[, 1] <- nu[, 1]
  for (t in seq(2, n_y))
    eps[, t] <- rho_t * eps[, t-1] + nu[, t]

  if (!is.null(W_mat) && abs(lambda_s) > 0.01) {
    w_codes <- rownames(W_mat)
    for (t in seq_len(n_y)) {
      e_t   <- eps[, t]
      idx_w <- match(powiaty_v, w_codes)
      valid <- !is.na(idx_w)
      if (sum(valid) > 0) {
        e_w <- rep(0, nrow(W_mat))
        e_w[idx_w[valid]] <- e_t[valid]
        w_e <- as.numeric(W_mat %*% e_w)
        e_t[valid] <- e_t[valid] +
          lambda_s * w_e[idx_w[valid]]
        eps[, t] <- e_t
      }
    }
  }

  as.data.frame(eps) %>%
    tibble::rownames_to_column("code_pow") %>%
    pivot_longer(-code_pow,
                 names_to  = "year",
                 values_to = "noise") %>%
    mutate(year = as.integer(year))
}

# Macierz W
W_mat <- NULL
if (exists("pow_ord") && exists("lw_p")) {
  W_mat <- tryCatch(
    spdep::listw2mat(lw_p),
    error = function(e) NULL
  )
}

# ============================================================
# 4. KRZYWA ŁĄCZNA Z MODELU
# ============================================================

compute_curve <- function(mod, data_in, var,
                           n_grid = N_GRID_PERT,
                           v_ref  = 0) {

  model_vars   <- names(mod$var.summary)
  numeric_cols <- names(data_in)[
    sapply(data_in, is.numeric)
  ]
  vars_to_agg <- intersect(model_vars, numeric_cols)

  pow_base <- data_in %>%
    group_by(code, name) %>%
    summarise(
      x    = first(x), y = first(y),
      year_c = median(year_c, na.rm=TRUE),
      covid_period = {
        lv <- levels(data_in$covid_period)
        if (is.null(lv))
          lv <- c("pre_covid","covid_2020_2021",
                  "post_covid")
        factor(lv[1], levels=lv)
      },
      across(all_of(vars_to_agg),
             ~ mean(.x, na.rm=TRUE)),
      .groups = "drop"
    ) %>%
    mutate(code_pow = make_code_pow(code))

  vseq <- seq(
    quantile(data_in[[var]], 0.05, na.rm=TRUE),
    quantile(data_in[[var]], 0.95, na.rm=TRUE),
    length.out = n_grid
  )

  nd_ref <- pow_base
  for (v in all_terms)
    if (v %in% names(nd_ref)) nd_ref[[v]] <- 0
  nd_ref[[var]] <- v_ref

  p_ref <- tryCatch(
    as.numeric(predict(mod, newdata=nd_ref,
                       type="response")),
    error=function(e) rep(NA_real_, nrow(nd_ref))
  )

  map_dbl(vseq, function(vv) {
    nd_v        <- nd_ref
    nd_v[[var]] <- vv
    p_v <- tryCatch(
      as.numeric(predict(mod, newdata=nd_v,
                         type="response")),
      error=function(e) rep(NA_real_, nrow(nd_v))
    )
    median(p_v - p_ref, na.rm=TRUE)
  }) -> effects

  tibble(var_value = vseq, effect = effects,
         var       = var,
         variable  = clean_var_p(var),
         var_label = label_var_p(clean_var_p(var)),
         component = comp_from_term_p(var))
}

# ============================================================
# 5. KRZYWE BAZOWE
# ============================================================

cat("Krzywe bazowe...\n")
curves_base <- map_dfr(all_terms, function(var) {
  tryCatch(
    compute_curve(m, dat, var) %>%
      mutate(run=0L, type="baseline"),
    error=function(e){warning(e$message);NULL}
  )
})
cat("  wierszy:", nrow(curves_base), "\n")

# ============================================================
# 6. PERTURBACJE
# ============================================================

run_single <- function(b, scale_val) {
  noise_df <- gen_noise(dat, rho_t, lambda_s,
                         sigma_resid, scale_val,
                         W_mat)
  dat_p <- dat %>%
    left_join(noise_df,
              by=c("code_pow","year")) %>%
    mutate(
      noise = replace_na(noise, 0),
      !!sym(OUTCOME_P) :=
        .data[[OUTCOME_P]] + noise
    )

  m_p <- tryCatch(
    mgcv::bam(
      formula(m),
      data     = dat_p,
      family   = m$family,
      method   = BAM_METHOD,
      gamma    = GAMMA_P,
      discrete = BAM_DISCRETE,
      nthreads = BAM_NTHREADS,
      control  = gam.control(trace=FALSE)
    ),
    error=function(e){warning(e$message);NULL}
  )
  if (is.null(m_p)) return(NULL)

  map_dfr(all_terms, function(var) {
    tryCatch(
      compute_curve(m_p, dat_p, var) %>%
        mutate(run=b, type="perturbed",
               noise_scale=scale_val),
      error=function(e) NULL
    )
  })
}

cat("Perturbacje B =", B_PERTURB, "...\n")
curves_pert <- map_dfr(seq_len(B_PERTURB),
                        function(b) {
  if (b %% 10 == 0)
    cat("  b =", b, "/", B_PERTURB, "\n")
  run_single(b, NOISE_SCALE)
}) %>%
  mutate(noise_scale = NOISE_SCALE)

cat("Udane perturbacje:",
    n_distinct(curves_pert$run), "\n")

write.csv(curves_pert,
          file.path(TAB_DIR,
                    "curves_perturbed.csv"),
          row.names=FALSE)

# ============================================================
# 7. SNR TABLE
# ============================================================

pert_summary <- curves_pert %>%
  group_by(var, variable, var_label,
           component, var_value) %>%
  summarise(
    effect_median = median(effect, na.rm=TRUE),
    effect_q25    = quantile(effect,.25,na.rm=TRUE),
    effect_q75    = quantile(effect,.75,na.rm=TRUE),
    effect_p10    = quantile(effect,.10,na.rm=TRUE),
    effect_p90    = quantile(effect,.90,na.rm=TRUE),
    effect_sd     = sd(effect, na.rm=TRUE),
    n_runs        = n_distinct(run),
    .groups="drop"
  )

snr_table <- curves_base %>%
  group_by(var, variable, var_label, component) %>%
  summarise(
    signal_range = diff(range(effect, na.rm=TRUE)),
    aps_proxy = {
      eff <- effect[order(var_value)]
      vv  <- var_value[order(var_value)]
      vr  <- diff(range(vv, na.rm=TRUE))
      if (is.finite(vr) && vr > 0)
        mean(diff(eff)/diff(vv), na.rm=TRUE)
      else NA_real_
    },
    .groups="drop"
  ) %>%
  left_join(
    pert_summary %>%
      group_by(var, variable, var_label,
               component) %>%
      summarise(
        noise_sd = mean(effect_sd, na.rm=TRUE),
        .groups="drop"
      ),
    by=c("var","variable","var_label","component")
  ) %>%
  left_join(
    # r_shape na zagregowanej krzywej (mediana perturbacji)
    curves_base %>%
      select(var, variable, component,
             var_value, effect_base=effect) %>%
      inner_join(
        pert_summary %>%
          select(var, variable, component,
                 var_value, effect_median),
        by=c("var","variable","component",
             "var_value")
      ) %>%
      group_by(var, variable, component) %>%
      summarise(
        r_shape = if (n()>2)
          suppressWarnings(
            cor(effect_base, effect_median,
                use="complete.obs"))
        else NA_real_,
        .groups="drop"
      ),
    by=c("var","variable","component")
  ) %>%
  left_join(
    # r_shape per pojedyncza perturbacja — mean i SD po replikacjach
    curves_base %>%
      select(var, variable, component,
             var_value, effect_base=effect) %>%
      inner_join(
        curves_pert %>%
          select(var, variable, component,
                 var_value, run, effect_pert=effect),
        by=c("var","variable","component","var_value")
      ) %>%
      group_by(var, variable, component, run) %>%
      summarise(
        r_shape_run = if (n()>2)
          suppressWarnings(
            cor(effect_base, effect_pert,
                use="complete.obs"))
        else NA_real_,
        .groups="drop"
      ) %>%
      group_by(var, variable, component) %>%
      summarise(
        mean_r_shape = mean(r_shape_run, na.rm=TRUE),
        sd_r_shape   = sd(r_shape_run,   na.rm=TRUE),
        pct_positive_r = mean(r_shape_run > 0, na.rm=TRUE),
        .groups="drop"
      ),
    by=c("var","variable","component")
  ) %>%
  mutate(
    snr       = signal_range / noise_sd,
    snr_vs_aps = abs(aps_proxy) / noise_sd,
    snr_class = case_when(
      snr > 4 ~ "strong",
      snr > 2 ~ "moderate",
      TRUE    ~ "weak"
    ),
    meaningful   = abs(aps_proxy) >= 1.0,
    final_class  = case_when(
      !meaningful ~ "negligible_effect",
      snr > 2     ~ "reliable",
      TRUE        ~ "unreliable"
    ),
    # robust: SNR > 2 AND kształt zachowany zarówno na medianie
    # jak i średnio per-replikację
    robust       = snr > 2 & r_shape > 0.8 & mean_r_shape > 0.8,
    # dodatkowa flaga: czy kształt jest spójny między replikacjami
    shape_consistent = sd_r_shape < 0.2 & pct_positive_r > 0.9,
    period       = PERIOD
  ) %>%
  arrange(component, desc(snr))

write.csv(snr_table,
          file.path(TAB_DIR, "snr_table.csv"),
          row.names=FALSE)

cat("\n=== SNR table ===\n")
print(snr_table %>%
        select(var_label, component,
               snr, aps_proxy,
               r_shape, mean_r_shape, sd_r_shape,
               pct_positive_r, final_class),
      n=60)


# ============================================================
# 7b. TABELE I BOXPLOTY PERTURBACJI
# ============================================================

perturbation_run_summary <- curves_pert %>%
  group_by(run, variable, var_label, component) %>%
  summarise(
    effect_range = diff(range(effect, na.rm=TRUE)),
    effect_mean  = mean(effect, na.rm=TRUE),
    effect_sd    = sd(effect, na.rm=TRUE),
    .groups = "drop"
  )

perturbation_box_stats <- perturbation_run_summary %>%
  group_by(variable, var_label, component) %>%
  summarise(
    n_runs = n_distinct(run),
    effect_range_q25 = quantile(effect_range, .25, na.rm=TRUE),
    effect_range_median = median(effect_range, na.rm=TRUE),
    effect_range_q75 = quantile(effect_range, .75, na.rm=TRUE),
    effect_sd_q25 = quantile(effect_sd, .25, na.rm=TRUE),
    effect_sd_median = median(effect_sd, na.rm=TRUE),
    effect_sd_q75 = quantile(effect_sd, .75, na.rm=TRUE),
    .groups = "drop"
  ) %>%
  arrange(component, var_label)

perturbation_snr_box_stats <- snr_table %>%
  group_by(component, final_class) %>%
  summarise(
    n_terms = n(),
    snr_q25 = quantile(snr, .25, na.rm=TRUE),
    snr_median = median(snr, na.rm=TRUE),
    snr_q75 = quantile(snr, .75, na.rm=TRUE),
    aps_abs_median = median(abs(aps_proxy), na.rm=TRUE),
    .groups = "drop"
  )

write.csv(perturbation_run_summary,
          file.path(TAB_DIR, "perturbation_run_summary.csv"),
          row.names=FALSE)
write.csv(perturbation_box_stats,
          file.path(TAB_DIR, "perturbation_box_stats.csv"),
          row.names=FALSE)
write.csv(perturbation_snr_box_stats,
          file.path(TAB_DIR, "perturbation_snr_box_stats.csv"),
          row.names=FALSE)

p_box_pert_range <- perturbation_run_summary %>%
  mutate(label = paste0(var_label, " (", component, ")"),
         label = fct_reorder(label, effect_range, .fun=median, na.rm=TRUE)) %>%
  ggplot(aes(x=effect_range, y=label)) +
  geom_boxplot(outlier.alpha=.25, outlier.size=.5, fill=LANCET$grey_light,
               colour=LANCET$grey_dark) +
  labs(title="Perturbation effect range by variable",
       subtitle=paste0(YEAR_LABEL, " | B=", B_PERTURB,
                       " | noise_scale=", NOISE_SCALE),
       x="Range of perturbation curve [ASMR]", y=NULL) +
  theme_lancet(base_size=9)

p_box_pert_sd <- perturbation_run_summary %>%
  mutate(label = paste0(var_label, " (", component, ")"),
         label = fct_reorder(label, effect_sd, .fun=median, na.rm=TRUE)) %>%
  ggplot(aes(x=effect_sd, y=label)) +
  geom_boxplot(outlier.alpha=.25, outlier.size=.5, fill=LANCET$grey_light,
               colour=LANCET$grey_dark) +
  labs(title="Perturbation within-curve SD by variable",
       subtitle=paste0(YEAR_LABEL, " | B=", B_PERTURB,
                       " | noise_scale=", NOISE_SCALE),
       x="SD of perturbation curve [ASMR]", y=NULL) +
  theme_lancet(base_size=9)

final_colours_box <- c(
  reliable          = LANCET$navy,
  unreliable        = LANCET$red,
  negligible_effect = LANCET$grey_mid
)

p_box_snr_class <- snr_table %>%
  mutate(final_class = factor(final_class, levels=names(final_colours_box))) %>%
  ggplot(aes(x=final_class, y=snr, fill=final_class)) +
  geom_hline(yintercept=2, linetype="dashed", colour=LANCET$grey_mid) +
  geom_boxplot(outlier.alpha=.35, outlier.size=.7, colour=LANCET$grey_dark) +
  scale_fill_manual(values=final_colours_box, name=NULL) +
  labs(title="SNR distribution by perturbation robustness class",
       subtitle=YEAR_LABEL, x=NULL, y="SNR") +
  theme_lancet(base_size=10) +
  theme(legend.position="none")

save_hp_plot(file.path(FIG_DIR, "boxplot_perturbation_effect_range.png"),
             p_box_pert_range, width=10, height=7)
save_hp_plot(file.path(FIG_DIR, "boxplot_perturbation_effect_sd.png"),
             p_box_pert_sd, width=10, height=7)
save_hp_plot(file.path(FIG_DIR, "boxplot_snr_by_class.png"),
             p_box_snr_class, width=8, height=5)

# ============================================================
# 8. ANALIZA WIELOSKALOWA (opcjonalna)
# ============================================================

snr_multiscale <- tibble()

if (RUN_MULTISCALE) {
  cat("\nAnaliza wieloskalowa",
      length(NOISE_SCALES), "skal...\n")

  # signal_range z krzywych bazowych (staly)
  signal_base <- curves_base %>%
    group_by(var, variable, var_label, component) %>%
    summarise(
      signal_range = diff(range(effect, na.rm=TRUE)),
      aps_proxy    = {
        eff <- effect[order(var_value)]
        vv  <- var_value[order(var_value)]
        vr  <- diff(range(vv, na.rm=TRUE))
        if (is.finite(vr) && vr > 0)
          mean(diff(eff)/diff(vv), na.rm=TRUE)
        else NA_real_
      },
      .groups="drop"
    )

  snr_multiscale <- map_dfr(
    NOISE_SCALES, function(ns) {
      cat("  scale =", ns, "\n")
      curves_ns <- map_dfr(
        seq_len(MULTISCALE_B), function(b)
          run_single(b, ns)
      ) %>%
        mutate(noise_scale = ns)

      if (nrow(curves_ns) == 0) return(NULL)

      curves_ns %>%
        group_by(var, variable, var_label,
                 component, var_value) %>%
        summarise(sd_pert=sd(effect,na.rm=TRUE),
                  .groups="drop") %>%
        group_by(var, variable, var_label,
                 component) %>%
        summarise(noise_sd=mean(sd_pert,na.rm=TRUE),
                  .groups="drop") %>%
        left_join(signal_base,
                  by=c("var","variable",
                       "var_label","component")) %>%
        mutate(
          noise_scale  = ns,
          snr          = signal_range / noise_sd,
          snr_vs_aps   = abs(aps_proxy) / noise_sd,
          meaningful   = abs(aps_proxy) >= 1.0,
          final_class  = case_when(
            !meaningful ~ "negligible_effect",
            snr > 2     ~ "reliable",
            TRUE        ~ "unreliable"
          )
        )
    }
  )

  write.csv(snr_multiscale,
            file.path(TAB_DIR,
                      "snr_multiscale.csv"),
            row.names=FALSE)
}

# ============================================================
# 9. WYKRESY
# ============================================================

final_colours <- c(
  reliable          = LANCET$navy,
  unreliable        = LANCET$red,
  negligible_effect = LANCET$grey_mid
)

# Panel A: SNR
p_snr_a <- snr_table %>%
  mutate(
    var_label   = fct_reorder(var_label, snr),
    final_class = factor(final_class,
                          levels=names(final_colours))
  ) %>%
  ggplot(aes(x=snr, y=var_label,
             colour=final_class,
             shape=component)) +
  geom_vline(xintercept=2, linetype="dashed",
             colour=LANCET$grey_mid,
             linewidth=0.4) +
  geom_point(size=3,
             position=position_dodge(0.5)) +
  scale_colour_manual(values=final_colours,
                      name=NULL) +
  scale_shape_manual(
    values=c(between=16, within=17),
    name=NULL) +
  labs(title="SNR",
       subtitle="signal_range(base) / mean_SD(perturbations)",
       x="SNR", y=NULL) +
  theme_lancet(base_size=10) +
  theme(legend.position="none")

# Panel B: |APS proxy|
p_snr_b <- snr_table %>%
  mutate(
    var_label   = fct_reorder(var_label, snr),
    final_class = factor(final_class,
                          levels=names(final_colours))
  ) %>%
  ggplot(aes(x=abs(aps_proxy), y=var_label,
             colour=final_class,
             shape=component)) +
  geom_vline(xintercept=1, linetype="dashed",
             colour=LANCET$grey_mid,
             linewidth=0.4) +
  geom_point(size=3,
             position=position_dodge(0.5)) +
  scale_colour_manual(values=final_colours,
                      name=NULL) +
  scale_shape_manual(
    values=c(between=16, within=17),
    name=NULL) +
  labs(title="|APS proxy|",
       subtitle="average partial slope of baseline curve",
       x="|APS proxy| [ASMR/SD]", y=NULL) +
  theme_lancet(base_size=10) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank(),
        legend.position="none")

# Panel C: scatter SNR vs APS
p_snr_c <- snr_table %>%
  mutate(
    final_class = factor(final_class,
                          levels=names(final_colours)),
    label = paste0(var_label, " (", component, ")")
  ) %>%
  ggplot(aes(x=abs(aps_proxy), y=snr,
             colour=final_class,
             shape=component,
             label=label)) +
  geom_hline(yintercept=2, linetype="dashed",
             colour=LANCET$grey_mid, linewidth=0.35) +
  geom_vline(xintercept=1, linetype="dashed",
             colour=LANCET$grey_mid, linewidth=0.35) +
  annotate("rect", xmin=1, xmax=Inf,
           ymin=2, ymax=Inf,
           fill=LANCET$navy, alpha=0.04) +
  geom_point(size=3) +
  {
    if (requireNamespace("ggrepel", quietly=TRUE))
      ggrepel::geom_text_repel(
        size=2.5, max.overlaps=20,
        show.legend=FALSE) else
      ggplot2::geom_text(
        hjust=-0.1, size=2.5,
        show.legend=FALSE)
  } +
  scale_colour_manual(values=final_colours,
                      name=NULL) +
  scale_shape_manual(
    values=c(between=16, within=17),
    name=NULL) +
  labs(title="SNR vs |APS proxy|",
       subtitle="Top-right = reliable & meaningful",
       x="|APS proxy|", y="SNR") +
  theme_lancet(base_size=10)

p_snr_combined <- (p_snr_a | p_snr_b | p_snr_c) +
  plot_layout(widths=c(1,0.8,1.2),
              guides="collect") +
  plot_annotation(
    title=paste0(
      "Signal-to-Noise Ratio | B=", B_PERTURB,
      " | noise_scale=", NOISE_SCALE,
      " | ", YEAR_LABEL
    ),
    subtitle=paste0(
      "SNR = signal_range(baseline) / mean_SD(perturbations).  ",
      "|APS proxy| = average partial slope of baseline curve.\n",
      "Reliable = SNR>2 AND |APS|\u22651."
    ),
    caption=paste0(
      "Noise: rho_t=", round(rho_t,3),
      "  lambda_s=", round(lambda_s,3),
      "  sigma_scale=", NOISE_SCALE
    ),
    theme=theme(
      plot.title=element_text(
        size=11,face="bold",hjust=0,
        colour=LANCET$black),
      plot.subtitle=element_text(
        size=8,colour=LANCET$grey_mid),
      plot.caption=element_text(
        size=7,colour=LANCET$grey_mid),
      legend.position="bottom"
    )
  ) & theme(legend.position="bottom")

save_hp_plot(file.path(FIG_DIR, "snr_panel.png"),
             p_snr_combined,
             width=18, height=9)
cat("Zapisano: snr_panel.png\n")

# Krzywe perturbacji per zmienna
all_vars_p <- snr_table %>%
  distinct(variable, var_label, component) %>%
  arrange(component, var_label)

for (comp_n in c("between","within")) {
  vars_c <- all_vars_p %>%
    filter(component == comp_n)
  if (nrow(vars_c) == 0) next

  plots_p <- map(seq_len(nrow(vars_c)),
                  function(i) {
    vclean <- vars_c$variable[i]
    vlabel <- vars_c$var_label[i]
    col    <- if (comp_n=="between")
      LANCET$navy else LANCET$red

    db <- curves_base %>%
      filter(variable==vclean,
             component==comp_n)
    dp <- pert_summary %>%
      filter(variable==vclean,
             component==comp_n)

    snr_v <- snr_table %>%
      filter(variable==vclean,
             component==comp_n) %>%
      pull(snr) %>% first()
    fc_v  <- snr_table %>%
      filter(variable==vclean,
             component==comp_n) %>%
      pull(final_class) %>% first()

    snr_col <- if (fc_v=="reliable")
      LANCET$navy
    else if (fc_v=="negligible_effect")
      LANCET$grey_mid
    else LANCET$red

    if (nrow(db)==0 || nrow(dp)==0)
      return(ggplot()+theme_void())

    ggplot() +
      geom_hline(yintercept=0,
                 linetype="dashed",
                 colour=LANCET$grey_mid,
                 linewidth=0.3) +
      geom_ribbon(data=dp,
                  aes(x=var_value,
                      ymin=effect_p10,
                      ymax=effect_p90),
                  fill=col,alpha=0.10,
                  colour=NA) +
      geom_ribbon(data=dp,
                  aes(x=var_value,
                      ymin=effect_q25,
                      ymax=effect_q75),
                  fill=col,alpha=0.20,
                  colour=NA) +
      geom_line(data=dp,
                aes(x=var_value,
                    y=effect_median),
                colour=col,
                linewidth=0.5,
                linetype="dashed",
                alpha=0.7) +
      geom_line(data=db,
                aes(x=var_value,
                    y=effect),
                colour=col,
                linewidth=1.0) +
      annotate("text", x=Inf, y=Inf,
               label=sprintf("SNR=%.1f", snr_v),
               hjust=1.1, vjust=1.5,
               size=3, colour=snr_col,
               fontface="bold") +
      labs(title=paste0(vlabel),
           x="[SD]", y="\u0394 ASMR") +
      theme_lancet(base_size=8) +
      theme(plot.title=element_text(size=8))
  })

  n_cols <- min(4, length(plots_p))
  n_rows <- ceiling(length(plots_p) / n_cols)

  p_curves <- wrap_plots(plots_p, ncol=n_cols) +
    plot_annotation(
      title=paste0(
        "Perturbation curves | ", comp_n,
        " | ", YEAR_LABEL
      ),
      subtitle=paste0(
        "Solid = baseline | Dashed+ribbon = ",
        B_PERTURB, " perturbations (IQR+P10/P90) | ",
        "noise_scale=", NOISE_SCALE
      ),
      theme=theme(
        plot.title=element_text(
          face="bold",size=11),
        plot.subtitle=element_text(
          size=8,colour=LANCET$grey_mid)
      )
    )

  save_hp_plot(
    file.path(FIG_DIR,
              paste0("perturbation_curves_",
                     comp_n, ".png")),
    p_curves,
    width=4*n_cols+1,
    height=3.5*n_rows+1.5
  )
  cat("Zapisano: perturbation_curves_",
      comp_n, ".png\n", sep="")
}

# Wieloskalowy SNR
if (nrow(snr_multiscale) > 0) {
  p_multi <- snr_multiscale %>%
    mutate(
      var_label   = fct_reorder(var_label, snr),
      final_class = factor(final_class,
                            levels=names(final_colours))
    ) %>%
    ggplot(aes(x=snr, y=var_label,
               colour=final_class,
               shape=component)) +
    geom_vline(xintercept=2, linetype="dashed",
               colour=LANCET$grey_mid,
               linewidth=0.35) +
    geom_point(size=2.4,
               position=position_dodge(0.5)) +
    scale_colour_manual(values=final_colours,
                        name=NULL) +
    scale_shape_manual(
      values=c(between=16,within=17),
      name=NULL) +
    facet_wrap(~paste0("scale=",noise_scale),
               nrow=1) +
    labs(
      title="SNR at different noise scales",
      subtitle=paste0(
        "signal_range from baseline model (constant).  ",
        "SNR should decrease with noise scale."
      ),
      x="SNR", y=NULL
    ) +
    theme_lancet(base_size=10) +
    theme(strip.background=element_blank(),
          strip.text=element_text(face="bold"),
          legend.position="bottom")

  save_hp_plot(file.path(FIG_DIR, "snr_multiscale.png"),
               p_multi, width=16, height=8)
  cat("Zapisano: snr_multiscale.png\n")
}

# ============================================================
# 10. ZAPIS
# ============================================================

saveRDS(
  list(
    curves_base    = curves_base,
    curves_pert    = curves_pert,
    pert_summary   = pert_summary,
    snr_table      = snr_table,
    snr_multiscale = snr_multiscale,
    noise_params   = list(
      rho_t       = rho_t,
      lambda_s    = lambda_s,
      sigma_resid = sigma_resid,
      noise_scale = NOISE_SCALE,
      B_perturb   = B_PERTURB
    ),
    perturbation_run_summary   = perturbation_run_summary,
    perturbation_box_stats     = perturbation_box_stats,
    perturbation_snr_box_stats = perturbation_snr_box_stats,
    # Dokumentacja metryk r_shape:
    # r_shape       = cor(baseline_curve, median_perturbed_curve) po gridzie
    # mean_r_shape  = mean po replikacjach cor(baseline, single_pert)
    # sd_r_shape    = SD po replikacjach — miara stabilnosci samego ksztaltu
    # pct_positive_r = odsetek replikacji z r > 0 (kierunek)
    period         = PERIOD,
    year_label     = YEAR_LABEL
  ),
  file.path(OUT_DIR, "perturbation_results.rds")
)

cat("\n=== Perturbation DONE ===\n")
cat("Period:", PERIOD, "|", YEAR_LABEL, "\n")
cat("Output:", OUT_DIR, "\n")
