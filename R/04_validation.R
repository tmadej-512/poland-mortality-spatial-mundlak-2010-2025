# ============================================================
# 04_validation.R
# Walidacja przestrzennego Mundlak GAM
#
# Wymaga: powiaty_spdf w pamięci
# Wczytuje: output/{PERIOD}/stvc_gam/gam_results.rds
# Wymaga: 00_config.R (lub PERIOD ustawionego w env)
# ============================================================

if (!exists("LANCET") || !exists("YEAR_MAX"))
  source(file.path(getwd(), "R", "00_config.R"))

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr)
  library(stringr); library(tibble); library(ggplot2)
  library(sf); library(mgcv); library(scales)
  library(patchwork)
})

GAM_RDS <- file.path(OUT_BASE, "stvc_gam",
                      "gam_results.rds")
OUT_DIR <- file.path(OUT_BASE, "validation")
TAB_DIR <- file.path(OUT_DIR, "tables")
FIG_DIR <- file.path(OUT_DIR, "figures")

dir.create(TAB_DIR, recursive = TRUE,
           showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE,
           showWarnings = FALSE)

cat("=== Validation ===\n")
cat("Period:", PERIOD, "|", YEAR_LABEL, "\n")

if (!file.exists(GAM_RDS))
  stop("Brak: ", GAM_RDS,
       "\nUruchom najpierw 02_stvc_gam.R")

# ============================================================
# 1. WCZYTAJ MODEL
# ============================================================

res      <- readRDS(GAM_RDS)
m        <- res$m
dat      <- res$dat_all
OUTCOME_V <- res$settings$OUTCOME
LAG_N    <- res$settings$WITHIN_LAG_N %||% 1L
GAMMA_V  <- res$settings$GAMMA %||% GAMMA

cat("n =", nrow(dat),
    "| R-sq =", round(summary(m)$r.sq, 3), "\n")

# ============================================================
# 2. PREDYKCJE I RESIDUA
# ============================================================

cat("Obliczanie residuów...\n")

dat_val <- dat %>%
  mutate(
    fitted     = as.numeric(predict(m, newdata = dat,
                                     type = "response")),
    observed   = .data[[OUTCOME_V]],
    resid      = observed - fitted,
    resid_std  = as.numeric(scale(resid)),
    code_pow   = make_code_pow(code)
  )

# ---- Metryki globalne ----
rmse_fn <- function(x) sqrt(mean(x^2, na.rm = TRUE))
mae_fn  <- function(x) mean(abs(x),   na.rm = TRUE)

val_global <- dat_val %>%
  summarise(
    n           = n(),
    n_powiat    = n_distinct(code),
    year_min    = min(year, na.rm = TRUE),
    year_max    = max(year, na.rm = TRUE),
    bias        = mean(resid,   na.rm = TRUE),
    rmse        = rmse_fn(resid),
    mae         = mae_fn(resid),
    resid_sd    = sd(resid,     na.rm = TRUE),
    r2_insample = 1 - sum(resid^2, na.rm = TRUE) /
      sum((observed - mean(observed, na.rm=TRUE))^2,
           na.rm = TRUE),
    period      = PERIOD
  )

val_by_year <- dat_val %>%
  group_by(year) %>%
  summarise(
    n        = n(),
    bias     = mean(resid,   na.rm = TRUE),
    rmse     = rmse_fn(resid),
    mae      = mae_fn(resid),
    r2       = suppressWarnings(
      cor(observed, fitted, use="complete.obs")^2),
    .groups = "drop"
  )

val_by_powiat <- dat_val %>%
  group_by(code, code_pow, name) %>%
  summarise(
    n_years    = n(),
    obs_mean   = mean(observed, na.rm = TRUE),
    fit_mean   = mean(fitted,   na.rm = TRUE),
    bias       = mean(resid,    na.rm = TRUE),
    rmse       = rmse_fn(resid),
    resid_sd   = sd(resid,      na.rm = TRUE),
    .groups    = "drop"
  ) %>%
  arrange(desc(abs(bias)))

write.csv(val_global,
          file.path(TAB_DIR, "val_global.csv"),
          row.names = FALSE)
write.csv(val_by_year,
          file.path(TAB_DIR, "val_by_year.csv"),
          row.names = FALSE)
write.csv(val_by_powiat,
          file.path(TAB_DIR, "val_by_powiat.csv"),
          row.names = FALSE)

cat("\n=== Metryki globalne ===\n")
print(val_global)

# ============================================================
# 3. GAM.CHECK — k-indeks
# ============================================================

cat("\ngam.check...\n")

k_check <- tryCatch(
  as_tibble(mgcv::k.check(m), rownames = "term"),
  error = function(e) {
    warning("k.check: ", e$message); tibble()
  }
)

concurv <- tryCatch({
  cc <- mgcv::concurvity(m, full = TRUE)
  as.data.frame(t(cc)) %>%
    rownames_to_column("term") %>%
    as_tibble()
}, error = function(e) {
  warning("concurvity: ", e$message); tibble()
})

if (nrow(k_check) > 0) {
  write.csv(k_check,
            file.path(TAB_DIR, "k_check.csv"),
            row.names = FALSE)
  cat("\nk-check (uwaga: k-index < 1 przy p<0.05):\n")
  print(k_check)
}

if (nrow(concurv) > 0)
  write.csv(concurv,
            file.path(TAB_DIR, "concurvity.csv"),
            row.names = FALSE)

# ============================================================
# 4. AUTOKORELACJA TEMPORALNA RESIDUÓW (ACF1)
# ============================================================

temporal_acf <- dat_val %>%
  arrange(code, year) %>%
  group_by(code) %>%
  summarise(
    n    = n(),
    acf1 = if (n() > 2)
      suppressWarnings(
        cor(resid[-n()], resid[-1],
            use = "complete.obs"))
    else NA_real_,
    .groups = "drop"
  ) %>%
  filter(is.finite(acf1))

acf_summary <- temporal_acf %>%
  summarise(
    mean_acf1    = mean(acf1,   na.rm = TRUE),
    median_acf1  = median(acf1, na.rm = TRUE),
    sd_acf1      = sd(acf1,     na.rm = TRUE),
    pct_gt_02    = mean(abs(acf1) > 0.2,
                         na.rm = TRUE) * 100
  )

write.csv(temporal_acf,
          file.path(TAB_DIR, "acf_by_powiat.csv"),
          row.names = FALSE)
write.csv(acf_summary,
          file.path(TAB_DIR, "acf_summary.csv"),
          row.names = FALSE)

cat("\nACF(1) residuów:\n")
print(acf_summary)

# ============================================================
# 5. MORAN I RESIDUÓW PER ROK
# ============================================================

moran_by_year    <- tibble()
moran_map_metrics <- tibble()

if (exists("powiaty_spdf")) {

  if (!requireNamespace("spdep", quietly = TRUE)) {
    warning("Brak pakietu spdep — pomijam Moran I.")
  } else {

    cat("Moran I...\n")

    pow_ord <- powiaty_spdf %>%
      mutate(
        JPT_KOD_JE = as.character(JPT_KOD_JE),
        JPT_KOD_JE = ifelse(nchar(JPT_KOD_JE) == 3,
                            paste0("0", JPT_KOD_JE),
                            JPT_KOD_JE)
      ) %>%
      arrange(JPT_KOD_JE)

    nb <- spdep::poly2nb(pow_ord, queen = TRUE)
    lw <- spdep::nb2listw(nb, style = "W",
                           zero.policy = TRUE)

    safe_moran <- function(x, lw) {
      tryCatch({
        mt <- spdep::moran.test(
          x, lw, zero.policy = TRUE,
          na.action = na.exclude
        )
        tibble(
          moran_i  = unname(mt$estimate[1]),
          expected = unname(mt$estimate[2]),
          p_value  = mt$p.value,
          sig      = mt$p.value < 0.05
        )
      }, error = function(e)
        tibble(moran_i=NA_real_, expected=NA_real_,
               p_value=NA_real_, sig=NA)
      )
    }

    moran_by_year <- map_dfr(
      sort(unique(dat_val$year)), function(yr) {
        d_yr <- dat_val %>% filter(year == yr)
        idx  <- match(pow_ord$JPT_KOD_JE,
                       d_yr$code_pow)
        r    <- d_yr$resid[idx]
        safe_moran(r, lw) %>%
          mutate(year = yr, .before = 1)
      }
    )

    # Moran dla metryk powiatowych
    bias_ord <- val_by_powiat$bias[
      match(pow_ord$JPT_KOD_JE,
             val_by_powiat$code_pow)
    ]
    rmse_ord <- val_by_powiat$rmse[
      match(pow_ord$JPT_KOD_JE,
             val_by_powiat$code_pow)
    ]

    moran_map_metrics <- bind_rows(
      safe_moran(bias_ord, lw) %>%
        mutate(metric = "mean_bias", .before=1),
      safe_moran(rmse_ord, lw) %>%
        mutate(metric = "rmse", .before=1)
    )

    write.csv(moran_by_year,
              file.path(TAB_DIR,
                        "moran_by_year.csv"),
              row.names = FALSE)
    write.csv(moran_map_metrics,
              file.path(TAB_DIR,
                        "moran_map_metrics.csv"),
              row.names = FALSE)

    cat("\nMoran I per rok:\n")
    print(moran_by_year %>%
            select(year, moran_i, p_value, sig))
  }
}

# ============================================================
# 6. LOO-YEAR CV
# ============================================================

cat("\nLOO-year CV...\n")

fml_main <- formula(m)
years_cv <- sort(unique(dat$year))

loo_results <- map_dfr(years_cv, function(yr) {
  cat("  LOO year:", yr, "\n")

  train <- dat %>% filter(year != yr)
  test  <- dat %>% filter(year == yr)

  m_loo <- tryCatch(
    mgcv::bam(
      fml_main,
      data     = train,
      family   = m$family,
      method   = BAM_METHOD,
      gamma    = GAMMA_V,
      discrete = BAM_DISCRETE,
      nthreads = BAM_NTHREADS,
      control  = gam.control(trace = FALSE)
    ),
    error = function(e) {
      warning("LOO ", yr, ": ", e$message); NULL
    }
  )

  if (is.null(m_loo))
    return(tibble(year=yr, rmse=NA_real_,
                  mae=NA_real_, r2=NA_real_,
                  bias=NA_real_))

  pred_loo <- tryCatch(
    as.numeric(predict(m_loo, newdata=test,
                       type="response")),
    error=function(e) rep(NA_real_, nrow(test))
  )
  obs <- test[[OUTCOME_V]]

  tibble(
    year = yr,
    rmse = rmse_fn(obs - pred_loo),
    mae  = mae_fn(obs  - pred_loo),
    r2   = 1 - sum((obs-pred_loo)^2, na.rm=TRUE) /
      sum((obs - mean(obs, na.rm=TRUE))^2,
           na.rm=TRUE),
    bias = mean(pred_loo - obs, na.rm=TRUE),
    period = PERIOD
  )
})

loo_summary <- loo_results %>%
  summarise(
    rmse_mean = mean(rmse, na.rm=TRUE),
    rmse_sd   = sd(rmse,   na.rm=TRUE),
    r2_mean   = mean(r2,   na.rm=TRUE),
    bias_mean = mean(bias, na.rm=TRUE)
  )

write.csv(loo_results,
          file.path(TAB_DIR, "loo_year_cv.csv"),
          row.names = FALSE)
write.csv(loo_summary,
          file.path(TAB_DIR, "loo_summary.csv"),
          row.names = FALSE)

cat("\nLOO-year CV:\n")
print(loo_results %>% select(year, rmse, r2, bias))
cat("\nLOO summary:\n")
print(loo_summary)

# ============================================================
# 7. WYKRESY
# ============================================================

# Observed vs fitted
p1 <- ggplot(dat_val,
             aes(x=observed, y=fitted)) +
  geom_point(alpha=0.3, size=0.7,
             colour=LANCET$navy) +
  geom_abline(slope=1, intercept=0,
              linetype="dashed",
              colour=LANCET$red) +
  coord_equal() +
  labs(title="Observed vs fitted",
       subtitle=paste0(PERIOD, " | R\u00b2 = ",
                       round(val_global$r2_insample, 3)),
       x="Observed ASMR",
       y="Fitted ASMR") +
  theme_lancet()

# Residuals vs fitted
p2 <- ggplot(dat_val,
             aes(x=fitted, y=resid_std)) +
  geom_hline(yintercept=0, linetype="dashed",
             colour=LANCET$grey_mid) +
  geom_hline(yintercept=c(-2,2),
             linetype="dotted",
             colour=LANCET$grey_mid) +
  geom_point(alpha=0.3, size=0.7,
             colour=LANCET$navy) +
  geom_smooth(method="loess", formula=y~x,
              se=TRUE, colour=LANCET$red,
              linewidth=0.7) +
  labs(title="Residuals vs fitted",
       x="Fitted ASMR",
       y="Standardised residual") +
  theme_lancet()

# QQ plot
p3 <- ggplot(dat_val,
             aes(sample=resid_std)) +
  stat_qq(alpha=0.4, size=0.7,
          colour=LANCET$navy) +
  stat_qq_line(colour=LANCET$red,
               linewidth=0.7) +
  labs(title="Q-Q plot",
       x="Theoretical quantiles",
       y="Sample quantiles") +
  theme_lancet()

# Bias by year
p4 <- val_by_year %>%
  ggplot(aes(x=year, y=bias)) +
  geom_hline(yintercept=0, linetype="dashed",
             colour=LANCET$grey_mid) +
  geom_line(colour=LANCET$navy, linewidth=0.8) +
  geom_point(colour=LANCET$navy, size=2) +
  labs(title="Mean residual bias by year",
       x="Year", y="Mean residual") +
  theme_lancet()

# ACF histogram
p5 <- ggplot(temporal_acf, aes(x=acf1)) +
  geom_histogram(bins=35, fill=LANCET$navy,
                 alpha=0.7, colour="white") +
  geom_vline(xintercept=c(-0.2, 0, 0.2),
             linetype=c("dotted","dashed","dotted"),
             colour=LANCET$red) +
  labs(title="AR(1) residual autocorrelation",
       subtitle=paste0(
         "Mean = ", round(acf_summary$mean_acf1, 3),
         "  |  % |ACF1|>0.2 = ",
         round(acf_summary$pct_gt_02, 1)
       ),
       x="ACF(1)", y="Count") +
  theme_lancet()

# LOO-year CV
p6 <- loo_results %>%
  pivot_longer(c(rmse, r2),
               names_to="metric",
               values_to="value") %>%
  mutate(metric=recode(metric, rmse="RMSE",
                        r2="R\u00b2 OOS")) %>%
  ggplot(aes(x=year, y=value,
             colour=metric)) +
  geom_line(linewidth=0.8) +
  geom_point(size=2) +
  scale_colour_manual(
    values=c("RMSE"=LANCET$navy,
             "R\u00b2 OOS"=LANCET$red),
    name=NULL
  ) +
  facet_wrap(~metric, scales="free_y") +
  labs(title="LOO-year cross-validation",
       x="Held-out year", y=NULL) +
  theme_lancet() +
  theme(strip.background=element_blank(),
        strip.text=element_text(face="bold"))

# Moran I per rok
if (nrow(moran_by_year) > 0) {
  p7 <- moran_by_year %>%
    ggplot(aes(x=year, y=moran_i,
               colour=sig)) +
    geom_hline(yintercept=0, linetype="dashed",
               colour=LANCET$grey_mid) +
    geom_point(size=3) +
    geom_line(colour=LANCET$grey_mid,
              linewidth=0.4) +
    scale_colour_manual(
      values=c("TRUE"=LANCET$red,
               "FALSE"=LANCET$navy),
      labels=c("TRUE"="p<0.05",
               "FALSE"="p\u22650.05"),
      name=NULL
    ) +
    labs(title="Moran's I of residuals by year",
         subtitle=paste0(
           "Red = significant spatial autocorrelation"
         ),
         x="Year", y="Moran's I") +
    theme_lancet()

  save_hp_plot(file.path(FIG_DIR, "moran_by_year.png"),
               p7, width=8, height=5)
}

# Panel diagnostyczny główny
p_diag <- (p1 | p2) / (p3 | p4) / (p5 | p6) +
  plot_annotation(
    title=paste0(
      "Model diagnostics | ",
      "Spatial Mundlak GAM | ", YEAR_LABEL
    ),
    theme=theme(
      plot.title=element_text(
        size=12, face="bold", hjust=0.5,
        colour=LANCET$black)
    )
  )

save_hp_plot(file.path(FIG_DIR, "diagnostics_panel.png"),
             p_diag, width=16, height=18)

# Mapa residuów przestrzennych
if (exists("powiaty_spdf") && nrow(val_by_powiat) > 0) {

  powiaty_geom <- powiaty_spdf %>%
    mutate(
      JPT_KOD_JE = as.character(JPT_KOD_JE),
      JPT_KOD_JE = ifelse(nchar(JPT_KOD_JE)==3,
                          paste0("0",JPT_KOD_JE),
                          JPT_KOD_JE)
    ) %>%
    dplyr::rename(code_pow = JPT_KOD_JE) %>%
    dplyr::select(code_pow, geometry)

  map_df <- val_by_powiat %>%
    left_join(powiaty_geom, by="code_pow") %>%
    st_as_sf()

  lim_b <- max(abs(map_df$bias), na.rm=TRUE)

  p_map_bias <- ggplot(map_df) +
    geom_sf(aes(fill=bias), color=NA) +
    scale_fill_gradient2(
      low = HAP_COLS$blue, mid = HAP_COLS$white, high = HAP_COLS$vermillion,
      midpoint = 0, limits=c(-lim_b, lim_b),
      name="Mean residual", oob=squish
    ) +
    labs(title="Spatial pattern of mean residuals",
         subtitle=paste0(
           "Moran's I = ",
           if (nrow(moran_map_metrics)>0)
             round(moran_map_metrics$moran_i[
               moran_map_metrics$metric=="mean_bias"],3)
           else "n/a"
         )) +
    theme_lancet_map()

  p_map_rmse <- ggplot(map_df) +
    geom_sf(aes(fill=rmse), color=NA) +
    scale_fill_gradient(
      low = HAP_COLS$grey_panel, high = HAP_COLS$blue,
      name="RMSE", oob=squish
    ) +
    labs(title="Spatial pattern of RMSE") +
    theme_lancet_map()

  p_maps <- (p_map_bias | p_map_rmse) +
    plot_annotation(
      title="Spatial residual diagnostics",
      theme=theme(plot.title=element_text(
        size=11, face="bold", hjust=0.5))
    )

  save_hp_plot(file.path(FIG_DIR, "maps_residuals.png"),
               p_maps, width=14, height=7)
}


# ============================================================
# 7b. TABELE I BOXPLOTY WALIDACYJNE
# ============================================================

validation_box_stats_year <- dat_val %>%
  group_by(year) %>%
  summarise(
    n = n(),
    resid_min = min(resid, na.rm=TRUE),
    resid_q25 = quantile(resid, .25, na.rm=TRUE),
    resid_median = median(resid, na.rm=TRUE),
    resid_q75 = quantile(resid, .75, na.rm=TRUE),
    resid_max = max(resid, na.rm=TRUE),
    abs_resid_median = median(abs(resid), na.rm=TRUE),
    rmse = rmse_fn(resid),
    .groups = "drop"
  )

validation_box_stats_powiat <- val_by_powiat %>%
  summarise(
    n_powiat = n(),
    bias_q25 = quantile(bias, .25, na.rm=TRUE),
    bias_median = median(bias, na.rm=TRUE),
    bias_q75 = quantile(bias, .75, na.rm=TRUE),
    rmse_q25 = quantile(rmse, .25, na.rm=TRUE),
    rmse_median = median(rmse, na.rm=TRUE),
    rmse_q75 = quantile(rmse, .75, na.rm=TRUE),
    resid_sd_median = median(resid_sd, na.rm=TRUE)
  )

write.csv(validation_box_stats_year,
          file.path(TAB_DIR, "validation_box_stats_by_year.csv"),
          row.names = FALSE)
write.csv(validation_box_stats_powiat,
          file.path(TAB_DIR, "validation_box_stats_by_powiat.csv"),
          row.names = FALSE)

p_box_resid_year <- dat_val %>%
  mutate(year = factor(year)) %>%
  ggplot(aes(x=year, y=resid)) +
  geom_hline(yintercept=0, linetype="dashed", colour=LANCET$grey_mid) +
  geom_boxplot(outlier.alpha=.25, outlier.size=.5, fill=LANCET$grey_light,
               colour=LANCET$grey_dark) +
  labs(title="Validation residuals by year",
       subtitle=YEAR_LABEL, x="Year", y="Observed - fitted ASMR") +
  theme_lancet()

p_box_abs_year <- dat_val %>%
  mutate(year = factor(year)) %>%
  ggplot(aes(x=year, y=abs(resid))) +
  geom_boxplot(outlier.alpha=.25, outlier.size=.5, fill=LANCET$grey_light,
               colour=LANCET$grey_dark) +
  labs(title="Absolute validation error by year",
       subtitle=YEAR_LABEL, x="Year", y="|Observed - fitted| ASMR") +
  theme_lancet()

p_box_powiat <- val_by_powiat %>%
  pivot_longer(c(bias, rmse), names_to="metric", values_to="value") %>%
  mutate(metric = recode(metric, bias="Mean residual bias", rmse="RMSE")) %>%
  ggplot(aes(x=metric, y=value)) +
  geom_boxplot(outlier.alpha=.25, outlier.size=.5, fill=LANCET$grey_light,
               colour=LANCET$grey_dark) +
  facet_wrap(~metric, scales="free_y") +
  labs(title="Validation metrics across powiaty",
       subtitle=YEAR_LABEL, x=NULL, y="ASMR") +
  theme_lancet() +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        strip.background=element_blank(), strip.text=element_text(face="bold"))

save_hp_plot(file.path(FIG_DIR, "boxplot_residuals_by_year.png"),
             p_box_resid_year, width=10, height=5)
save_hp_plot(file.path(FIG_DIR, "boxplot_abs_error_by_year.png"),
             p_box_abs_year, width=10, height=5)
save_hp_plot(file.path(FIG_DIR, "boxplot_powiat_validation_metrics.png"),
             p_box_powiat, width=8, height=5)

# ============================================================
# 8. ZAPIS
# ============================================================

saveRDS(
  list(
    val_global        = val_global,
    val_by_year       = val_by_year,
    val_by_powiat     = val_by_powiat,
    temporal_acf      = temporal_acf,
    acf_summary       = acf_summary,
    moran_by_year     = moran_by_year,
    moran_map_metrics = moran_map_metrics,
    loo_results       = loo_results,
    loo_summary       = loo_summary,
    validation_box_stats_year = validation_box_stats_year,
    validation_box_stats_powiat = validation_box_stats_powiat,
    k_check           = k_check,
    concurvity        = concurv,
    period            = PERIOD,
    year_label        = YEAR_LABEL
  ),
  file.path(OUT_DIR, "validation_results.rds")
)

cat("\n=== Validation DONE ===\n")
cat("Period:", PERIOD, "|", YEAR_LABEL, "\n")
cat("Output:", OUT_DIR, "\n")
