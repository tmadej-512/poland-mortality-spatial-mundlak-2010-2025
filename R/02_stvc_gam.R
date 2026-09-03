# ============================================================
# 02_stvc_gam.R
# Przestrzenny Mundlak GAM — estymacja + APS + krzywe + mapy
#
# Wymaga: dat, powiaty_spdf w pamięci
# Wymaga: 00_config.R
# ============================================================

if (!exists("LANCET")) source(file.path(getwd(), "R", "00_config.R"))

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr)
  library(stringr); library(tibble); library(mgcv)
  library(sf); library(scales)
})

OUT_DIR <- file.path(OUT_BASE, "stvc_gam")
TAB_DIR <- file.path(OUT_DIR, "tables")
MOD_DIR <- file.path(OUT_DIR, "models")
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(MOD_DIR, recursive = TRUE, showWarnings = FALSE)

cat("=== Spatial Mundlak GAM ===\n")
cat("Period:", PERIOD, "|", YEAR_LABEL, "\n")
validate_period_data(dat)
cat("K_GLOBAL=", K_GLOBAL, "K_SVC=", K_SVC,
    "K_BASELINE=", K_BASELINE, "gamma=", GAMMA, "\n")

# ============================================================
# 1. CENTROIDY
# ============================================================

if (!all(c("x", "y") %in% names(dat))) {
  cat("Dodaję centroidy...\n")
  centroids_xy <- powiaty_spdf %>%
    st_transform(2180) %>%
    mutate(
      JPT_KOD_JE = as.character(JPT_KOD_JE),
      JPT_KOD_JE = ifelse(nchar(JPT_KOD_JE) == 3,
                          paste0("0", JPT_KOD_JE),
                          JPT_KOD_JE)
    ) %>%
    mutate(cent = st_centroid(geometry),
           x    = st_coordinates(cent)[, 1],
           y    = st_coordinates(cent)[, 2]) %>%
    st_drop_geometry() %>%
    select(code_pow = JPT_KOD_JE, x, y)
  dat <- dat %>%
    mutate(
      code     = as.character(code),
      code_pow = make_code_pow(code),
      code_pow = ifelse(nchar(code_pow) == 3,
                        paste0("0", code_pow), code_pow)
    ) %>%
    select(-any_of(c("x", "y"))) %>%
    left_join(centroids_xy, by = "code_pow")
  cat("  x/y dodane.\n")
}

# ============================================================
# 2. PANEL MUNDLAK
# ============================================================

all_vars    <- unique(c(MODEL_VARS, ONLY_BETWEEN_VARS))
all_vars    <- intersect(all_vars, names(dat))
within_vars <- setdiff(all_vars,
                        intersect(ONLY_BETWEEN_VARS,
                                  all_vars))

dat_all <- dat %>%
  filter(year >= YEAR_MIN, year <= YEAR_MAX,
         !is.na(.data[[OUTCOME]]),
         !is.na(x), !is.na(y)) %>%
  mutate(
    code     = as.character(code),
    year     = as.integer(year),
    year_c   = year - mean(year, na.rm = TRUE),
    code_pow = make_code_pow(code),
    covid_period = droplevels(covid_period_factor(year))
  ) %>%
  select(code, code_pow, name, year, year_c,
         covid_period, x, y,
         all_of(OUTCOME), all_of(all_vars)) %>%
  mutate(across(all_of(all_vars), as.numeric)) %>%
  arrange(code, year) %>%
  group_by(code) %>%
  mutate(
    across(all_of(all_vars),
           ~ mean(.x, na.rm = TRUE),
           .names = "{.col}_between"),
    across(all_of(within_vars),
           ~ .x - mean(.x, na.rm = TRUE),
           .names = "{.col}_within")
  ) %>%
  mutate(
    across(
      paste0(within_vars, "_within"),
      ~ dplyr::lag(.x, n = WITHIN_LAG_N),
      .names = paste0("{.col}_lag", WITHIN_LAG_N)
    )
  ) %>%
  ungroup()

# Standaryzuj
terms_to_scale <- c(
  paste0(all_vars, "_between"),
  paste0(within_vars, "_within_lag", WITHIN_LAG_N)
)
dat_all <- dat_all %>%
  mutate(across(
    all_of(intersect(terms_to_scale, names(.))),
    scale_safe,
    .names = "{.col}_z"
  ))

between_terms <- intersect(
  c(paste0(MODEL_VARS, "_between_z"),
    paste0(ONLY_BETWEEN_VARS, "_between_z")),
  names(dat_all)
)
within_terms <- intersect(
  paste0(MODEL_VARS, "_within_lag", WITHIN_LAG_N, "_z"),
  names(dat_all)
)
all_terms <- c(between_terms, within_terms)

cat("Between:", length(between_terms),
    "| Within:", length(within_terms),
    "| n =", nrow(dat_all), "\n")

# ============================================================
# 3. FORMUŁA
# ============================================================

build_formula <- function(bt, wt) {
  ctrl <- sprintf("s(year_c,bs='%s',m=%d,k=%d)",
                  BS_YEAR, M_YEAR, K_YEAR)

  # Dla pełnego okresu dodaj covid_period tylko gdy ma >=2 obserwowane poziomy.
  covid_ctrl <- if (use_covid_control(dat_all))
    "covid_period" else character(0)

  baseline <- sprintf("s(x,y,bs='%s',k=%d)",
                      BS_BASELINE, K_BASELINE)

  smooths <- function(terms) unlist(map(terms, ~ c(
    sprintf("s(%s,bs='%s',m=%d,k=%d)",
            .x, BS_GLOBAL, M_GLOBAL, K_GLOBAL),
    sprintf("s(x,y,by=%s,bs='%s',m=%d,k=%d)",
            .x, BS_SVC, M_SVC, K_SVC)
  )))

  rhs <- paste(
    c(ctrl, covid_ctrl, baseline,
      smooths(bt), smooths(wt)),
    collapse = " +\n  "
  )
  as.formula(paste(OUTCOME, "~\n ", rhs))
}

cat("covid_period w danych GAM:\n")
print(table(dat_all$covid_period, useNA = "ifany"))
if (PERIOD == "full" && !use_covid_control(dat_all)) {
  warning("PERIOD='full', ale covid_period ma mniej niż 2 obserwowane poziomy; kontrola covid_period zostanie pominięta.")
}

fml <- build_formula(between_terms, within_terms)
cat("Formuła zbudowana.\n")

# ============================================================
# 4. ESTYMACJA
# ============================================================

cat("Start bam()...\n")
m <- mgcv::bam(
  fml,
  data     = dat_all,
  family   = FAMILY,
  method   = BAM_METHOD,
  discrete = BAM_DISCRETE,
  gamma    = GAMMA,
  nthreads = BAM_NTHREADS,
  control  = gam.control(trace = FALSE)
)

cat("R-sq =", round(summary(m)$r.sq, 3),
    "| dev.expl =",
    round(summary(m)$dev.expl, 3), "\n")

saveRDS(m, file.path(MOD_DIR, "gam_model.rds"))

# EDF tabela
edf_tbl <- as_tibble(summary(m)$s.table,
                      rownames = "term") %>%
  mutate(
    term_type = case_when(
      str_detect(term, "s\\(x,y\\)") &
        str_detect(term, ":")     ~ "SVC",
      str_detect(term, "s\\(x,y\\)") ~ "baseline",
      str_detect(term, "year_c")      ~ "year",
      TRUE                            ~ "global"
    ),
    raw_var   = str_extract(term, "(?<=s\\()[^,)]+"),
    component = comp_from_term(raw_var),
    variable  = clean_var(raw_var),
    var_label = label_var(variable),
    edf       = round(edf, 2),
    `p-value` = round(`p-value`, 4),
    sig = case_when(
      `p-value` < 0.001 ~ "***",
      `p-value` < 0.01  ~ "**",
      `p-value` < 0.05  ~ "*",
      `p-value` < 0.10  ~ ".",
      TRUE              ~ ""
    )
  )

write.csv(edf_tbl,
          file.path(TAB_DIR, "edf_table.csv"),
          row.names = FALSE)

cat("\nEDF (global + SVC):\n")
print(edf_tbl %>%
        filter(term_type %in% c("global", "SVC")) %>%
        select(var_label, component, term_type,
               edf, `p-value`, sig),
      n = 80)

# ============================================================
# 5. HELPERY PREDYKCJI
# ============================================================

build_pow_base <- function(var) {
  model_vars   <- names(m$var.summary)
  numeric_cols <- names(dat_all)[sapply(dat_all, is.numeric)]
  vars_to_agg  <- intersect(model_vars, numeric_cols)

  dat_all %>%
    group_by(code, name) %>%
    summarise(
      x = first(x),
      y = first(y),
      year_c = median(year_c, na.rm = TRUE),
      covid_period = covid_factor_ref(),
      v_actual = mean(.data[[var]], na.rm = TRUE),
      across(all_of(vars_to_agg), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    mutate(code_pow = make_code_pow(code))
}

delta_pred <- function(nd0, nd1) {
  p0 <- tryCatch(
    as.numeric(predict(m, newdata = nd0,
                       type = "response")),
    error = function(e) rep(NA_real_, nrow(nd0))
  )
  p1 <- tryCatch(
    as.numeric(predict(m, newdata = nd1,
                       type = "response")),
    error = function(e) rep(NA_real_, nrow(nd1))
  )
  p1 - p0
}

# ============================================================
# 6. AVERAGE PARTIAL SLOPES (APS)
# ============================================================
# APS zastępuje wcześniejsze AME. Liczymy ją w tym samym
# reżimie referencyjnym co krzywe i mapy:
#   - zmieniamy tylko focal determinant,
#   - wszystkie pozostałe determinant terms ustawiamy na 0,
#   - zachowujemy lokalizacje powiatów,
#   - predykcja obejmuje global smooth + retained SVC.
# Wynik to średnia po gridzie P2–P98 county-specific finite-
# difference slope; summary w tabeli to mediana i IQR po powiatach.

zero_det_terms <- function(nd) {
  for (v in all_terms) {
    if (v %in% names(nd)) nd[[v]] <- 0
  }
  nd
}

# Mediany standaryzowanych terminów — punkt referencyjny
# dla krzywych i map (zamiast zera = sredniej)
det_medians <- tryCatch({
  dat_all %>%
    dplyr::summarise(dplyr::across(
      dplyr::all_of(all_terms),
      ~ median(.x, na.rm = TRUE)
    )) %>% as.list()
}, error = function(e) {
  warning("det_medians: ", e$message)
  setNames(as.list(rep(0, length(all_terms))), all_terms)
})

median_det_terms <- function(nd) {
  for (v in all_terms) {
    if (v %in% names(nd) && v %in% names(det_medians))
      nd[[v]] <- det_medians[[v]]
  }
  nd
}

compute_aps <- function(var, eps = AME_EPS,
                        n_grid = N_GRID) {
  pow_base <- build_pow_base(var) %>%
    filter(!is.na(.data[[var]]), !is.na(x), !is.na(y))

  if (nrow(pow_base) == 0) return(NULL)

  q_lo <- quantile(dat_all[[var]], 0.02, na.rm = TRUE)
  q_hi <- quantile(dat_all[[var]], 0.98, na.rm = TRUE)
  if (!is.finite(q_lo) || !is.finite(q_hi) || q_lo == q_hi)
    return(NULL)

  vseq <- seq(q_lo, q_hi, length.out = n_grid)

  slope_grid_county <- map_dfr(vseq, function(vg) {
    nd0 <- zero_det_terms(pow_base)
    nd1 <- zero_det_terms(pow_base)

    nd0[[var]] <- vg
    nd1[[var]] <- vg + eps

    slope <- delta_pred(nd0, nd1) / eps

    tibble(
      code          = pow_base$code,
      code_pow      = pow_base$code_pow,
      name          = pow_base$name,
      var           = var,
      variable      = clean_var(var),
      var_label     = label_var(clean_var(var)),
      component     = comp_from_term(var),
      grid_value    = vg,
      partial_slope = slope
    )
  })

  county_slopes <- slope_grid_county %>%
    group_by(code, code_pow, name, var, variable,
             var_label, component) %>%
    summarise(
      aps_county = mean(partial_slope, na.rm = TRUE),
      aps_county_median = median(partial_slope, na.rm = TRUE),
      aps_county_q25 = quantile(partial_slope, 0.25, na.rm = TRUE),
      aps_county_q75 = quantile(partial_slope, 0.75, na.rm = TRUE),
      n_grid = sum(!is.na(partial_slope)),
      .groups = "drop"
    )

  aps_summary <- county_slopes %>%
    group_by(var, variable, var_label, component) %>%
    summarise(
      aps        = median(aps_county, na.rm = TRUE),
      aps_mean   = mean(aps_county, na.rm = TRUE),
      aps_q25    = quantile(aps_county, 0.25, na.rm = TRUE),
      aps_q75    = quantile(aps_county, 0.75, na.rm = TRUE),
      aps_p10    = quantile(aps_county, 0.10, na.rm = TRUE),
      aps_p90    = quantile(aps_county, 0.90, na.rm = TRUE),
      n_counties = sum(!is.na(aps_county)),
      .groups = "drop"
    )

  list(
    aps_summary = aps_summary,
    county_slopes = county_slopes,
    slope_grid_county = slope_grid_county
  )
}

cat("Average partial slopes (APS)...\n")
aps_results <- map(all_terms, function(var) {
  cat("  APS:", clean_var(var), comp_from_term(var), "\n")
  tryCatch(
    compute_aps(var),
    error = function(e) {
      warning(e$message); NULL
    }
  )
})

aps_all <- map_dfr(aps_results, "aps_summary") %>%
  mutate(
    direction = case_when(
      aps_q75 < 0 ~ "lower ASMR",
      aps_q25 > 0 ~ "higher ASMR",
      TRUE        ~ "uncertain"
    ),
    period = PERIOD
  )

aps_county_all <- map_dfr(aps_results, "county_slopes") %>%
  mutate(period = PERIOD)

aps_grid_county_all <- map_dfr(aps_results, "slope_grid_county") %>%
  mutate(period = PERIOD)

write.csv(aps_all,
          file.path(TAB_DIR, "aps_all.csv"),
          row.names = FALSE)
write.csv(aps_county_all,
          file.path(TAB_DIR, "aps_county_all.csv"),
          row.names = FALSE)
write.csv(aps_grid_county_all,
          file.path(TAB_DIR, "aps_grid_county_all.csv"),
          row.names = FALSE)

# ============================================================
# 7. KRZYWE IQR
# ============================================================

compute_curve_iqr <- function(var, v_ref = 0) {
  pow_base <- build_pow_base(var)
  vseq <- seq(
    quantile(dat_all[[var]], 0.02, na.rm = TRUE),
    quantile(dat_all[[var]], 0.98, na.rm = TRUE),
    length.out = N_GRID
  )
  nd_ref <- median_det_terms(pow_base)
  nd_ref[[var]] <- v_ref

  p_ref <- tryCatch(
    as.numeric(predict(m, newdata = nd_ref,
                       type = "response")),
    error = function(e) rep(NA_real_, nrow(nd_ref))
  )

  curve_raw <- map_dfr(seq_along(vseq), function(j) {
    nd_v        <- nd_ref
    nd_v[[var]] <- vseq[j]
    p_v <- tryCatch(
      as.numeric(predict(m, newdata = nd_v,
                         type = "response")),
      error = function(e) rep(NA_real_, nrow(nd_v))
    )
    tibble(var_value = vseq[j],
           code      = pow_base$code,
           effect    = p_v - p_ref)
  })

  curve_raw %>%
    group_by(var_value) %>%
    summarise(
      effect_median = median(effect, na.rm = TRUE),
      effect_q25    = quantile(effect, 0.25, na.rm = TRUE),
      effect_q75    = quantile(effect, 0.75, na.rm = TRUE),
      effect_p10    = quantile(effect, 0.10, na.rm = TRUE),
      effect_p90    = quantile(effect, 0.90, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      var       = var,
      variable  = clean_var(var),
      var_label = label_var(clean_var(var)),
      component = comp_from_term(var),
      period    = PERIOD
    )
}

cat("Krzywe...\n")
curves_all <- map_dfr(all_terms, function(var) {
  cat("  curve:", clean_var(var), "\n")
  tryCatch(
    compute_curve_iqr(var),
    error = function(e) {
      warning(e$message); NULL
    }
  )
})

write.csv(curves_all,
          file.path(TAB_DIR, "curves_iqr.csv"),
          row.names = FALSE)

# ============================================================
# 8. MAPY EFEKTÓW: POCHODNA + DELTA
# ============================================================

compute_effect_maps <- function(var, eps = AME_EPS) {
  component_use <- comp_from_term(var)
  is_within     <- component_use == "within"

  cat("  maps:", clean_var(var), component_use, "
")

  # ==========================================================
  # BETWEEN:
  #   marginal map:
  #     d(ASMR)/d(var) przy typowym poziomie strukturalnym powiatu
  #
  #   contribution map:
  #     predict(v_county) - predict(P2)
  # ==========================================================

  if (!is_within) {
    pow_base <- build_pow_base(var)

    meta <- pow_base %>%
      select(code, code_pow, name, x, y, v_actual)

    # Mapa 1: pochodna przy v_actual
    nd0 <- zero_det_terms(pow_base)
    nd0[[var]] <- pow_base$v_actual

    nd1 <- nd0
    nd1[[var]] <- pow_base$v_actual + eps

    deriv <- delta_pred(nd0, nd1) / eps

    # Mapa 2: contribution od P2 do v_actual
    v_ref <- quantile(dat_all[[var]], 0.02, na.rm = TRUE)

    nd_ref <- median_det_terms(pow_base)
    nd_ref[[var]] <- v_ref

    nd_act <- median_det_terms(pow_base)
    nd_act[[var]] <- pow_base$v_actual

    delta <- delta_pred(nd_ref, nd_act)

    cat(sprintf(
      "    between ref=P2 %.2f | deriv [%.2f, %.2f]  contribution [%.2f, %.2f]
",
      v_ref,
      min(deriv, na.rm = TRUE), max(deriv, na.rm = TRUE),
      min(delta, na.rm = TRUE), max(delta, na.rm = TRUE)
    ))

    return(
      bind_rows(
        meta %>% mutate(
          value     = deriv,
          map_type  = "deriv",
          map_label = "Marginal effect",
          subtitle  = "d(ASMR)/d(var) at county structural level | rest=0"
        ),
        meta %>% mutate(
          value     = delta,
          map_type  = "delta",
          map_label = "Contribution",
          subtitle  = "predict(v_county) - predict(P2) | rest=national median"
        )
      ) %>%
        mutate(
          var       = var,
          variable  = clean_var(var),
          var_label = label_var(clean_var(var)),
          component = component_use,
          v_ref     = v_ref,
          period    = PERIOD
        )
    )
  }

  # ==========================================================
  # WITHIN:
  #   marginal map:
  #     mean_t [ d(ASMR)/d(var_it) ]
  #
  #   contribution map:
  #     mean_t [ predict(v_it) - predict(0) ]
  #
  # gdzie v_it = obserwowane odchylenie powiatu od własnej średniej
  # w roku t, z lagiem użytym w modelu.
  # ==========================================================

  model_vars   <- names(m$var.summary)
  numeric_cols <- names(dat_all)[sapply(dat_all, is.numeric)]
  vars_keep    <- intersect(model_vars, numeric_cols)

  nd_year <- dat_all %>%
    select(any_of(c(
      "code", "code_pow", "name", "year",
      "x", "y", "year_c", "covid_period",
      vars_keep
    ))) %>%
    filter(
      !is.na(.data[[var]]),
      !is.na(x),
      !is.na(y)
    ) %>%
    mutate(
      code_pow = make_code_pow(code),
      covid_period = covid_factor_ref()
    )

  if (nrow(nd_year) == 0) {
    warning("Brak danych do map within dla: ", var)
    return(NULL)
  }

  # Mapa 1: roczna pochodna przy obserwowanym v_it
  nd0 <- median_det_terms(nd_year)
  nd0[[var]] <- nd_year[[var]]

  nd1 <- nd0
  nd1[[var]] <- nd_year[[var]] + eps

  p0_deriv <- tryCatch(
    as.numeric(predict(m, newdata = nd0, type = "response")),
    error = function(e) rep(NA_real_, nrow(nd0))
  )

  p1_deriv <- tryCatch(
    as.numeric(predict(m, newdata = nd1, type = "response")),
    error = function(e) rep(NA_real_, nrow(nd1))
  )

  deriv_it <- (p1_deriv - p0_deriv) / eps

  # Mapa 2: roczna contribution względem 0
  nd_ref <- median_det_terms(nd_year)
  nd_ref[[var]] <- 0

  nd_act <- median_det_terms(nd_year)
  nd_act[[var]] <- nd_year[[var]]

  p_ref <- tryCatch(
    as.numeric(predict(m, newdata = nd_ref, type = "response")),
    error = function(e) rep(NA_real_, nrow(nd_ref))
  )

  p_act <- tryCatch(
    as.numeric(predict(m, newdata = nd_act, type = "response")),
    error = function(e) rep(NA_real_, nrow(nd_act))
  )

  delta_it <- p_act - p_ref

  yearly_effects <- nd_year %>%
    transmute(
      code,
      code_pow,
      name,
      year,
      x,
      y,
      v_it     = .data[[var]],
      deriv_it = deriv_it,
      delta_it = delta_it
    )

  pow_effects <- yearly_effects %>%
    group_by(code, code_pow, name) %>%
    summarise(
      x = first(x),
      y = first(y),
      v_actual = mean(v_it, na.rm = TRUE),
      n_years = sum(!is.na(v_it)),

      # Finalna mapa marginalna within:
      # średnia po latach z lokalnych pochodnych.
      deriv = mean(deriv_it, na.rm = TRUE),

      # Finalna mapa contribution within:
      # średnia po latach z predict(v_it) - predict(0).
      delta = mean(delta_it, na.rm = TRUE),

      mean_v_it = mean(v_it, na.rm = TRUE),
      median_v_it = median(v_it, na.rm = TRUE),
      share_v_above_0 = mean(v_it > 0, na.rm = TRUE),
      share_v_below_0 = mean(v_it < 0, na.rm = TRUE),

      .groups = "drop"
    )

  deriv <- pow_effects$deriv
  delta <- pow_effects$delta

  cat(sprintf(
    "    within ref=0 averaged over years | deriv [%.2f, %.2f]  contribution [%.2f, %.2f]
",
    min(deriv, na.rm = TRUE), max(deriv, na.rm = TRUE),
    min(delta, na.rm = TRUE), max(delta, na.rm = TRUE)
  ))

  bind_rows(
    pow_effects %>%
      transmute(
        code, code_pow, name, x, y,
        v_actual, n_years,
        mean_v_it, median_v_it,
        share_v_above_0, share_v_below_0,
        value     = deriv,
        map_type  = "deriv",
        map_label = "Average partial slope",
        subtitle  = "mean over years of local partial slope | rest=national median"
      ),
    pow_effects %>%
      transmute(
        code, code_pow, name, x, y,
        v_actual, n_years,
        mean_v_it, median_v_it,
        share_v_above_0, share_v_below_0,
        value     = delta,
        map_type  = "delta",
        map_label = "Average within contribution",
        subtitle  = "mean over years of predict(v_it) - predict(0) | rest=national median"
      )
  ) %>%
    mutate(
      var       = var,
      variable  = clean_var(var),
      var_label = label_var(clean_var(var)),
      component = component_use,
      v_ref     = 0,
      period    = PERIOD
    )
}

cat("Mapy efektów...
")
effect_maps <- map_dfr(all_terms, function(var) {
  tryCatch(
    compute_effect_maps(var),
    error = function(e) {
      warning(e$message); NULL
    }
  )
})

write.csv(effect_maps,
          file.path(TAB_DIR, "effect_maps.csv"),
          row.names = FALSE)

# ============================================================
# 9. G_I
# ============================================================

tn_all <- tryCatch(
  colnames(predict(m, newdata = dat_all[1, ],
                   type = "terms")),
  error = function(e) character(0)
)

compute_gi <- function(var) {
  st <- tn_all[
    str_detect(tn_all, "s\\(x,\\s*y\\)") &
      str_detect(tn_all, fixed(var))
  ]
  if (length(st) == 0) return(NULL)
  pow_base      <- build_pow_base(var)
  nd_svc        <- pow_base
  nd_svc[[var]] <- 1
  pr <- tryCatch(
    predict(m, newdata = nd_svc,
            type = "terms", terms = st[1]),
    error = function(e) NULL
  )
  if (is.null(pr)) return(NULL)
  pow_base %>%
    select(code, code_pow, name, x, y) %>%
    mutate(
      g_i       = as.numeric(as.matrix(pr)[, 1]),
      var       = var,
      variable  = clean_var(var),
      var_label = label_var(clean_var(var)),
      component = comp_from_term(var),
      period    = PERIOD
    )
}

cat("g_i...\n")
gi_all <- map_dfr(all_terms, function(var) {
  tryCatch(compute_gi(var),
           error = function(e) {
             warning(e$message); NULL
           })
})

write.csv(gi_all,
          file.path(TAB_DIR, "gi_all.csv"),
          row.names = FALSE)


# ============================================================
# 10. PREDYKCJE ŁĄCZNE DLA POWIATÓW
# ============================================================

cat("Predykcje powiatowe...\n")

pred_response <- tryCatch(
  as.numeric(predict(m, newdata = dat_all, type = "response")),
  error = function(e) {
    warning("Nie udało się policzyć predykcji powiatowych: ", e$message)
    rep(NA_real_, nrow(dat_all))
  }
)

pred_yearly <- dat_all %>%
  transmute(
    period       = PERIOD,
    code         = code,
    code_pow     = code_pow,
    name         = name,
    year         = year,
    observed     = .data[[OUTCOME]],
    predicted    = pred_response,
    residual     = observed - predicted,
    abs_error    = abs(residual),
    covid_period = covid_period
  )

pred_powiat <- pred_yearly %>%
  group_by(period, code, code_pow, name) %>%
  summarise(
    year_min       = min(year, na.rm = TRUE),
    year_max       = max(year, na.rm = TRUE),
    n_years        = n_distinct(year),
    observed_mean  = mean(observed, na.rm = TRUE),
    predicted_mean = mean(predicted, na.rm = TRUE),
    residual_mean  = mean(residual, na.rm = TRUE),
    abs_error_mean = mean(abs_error, na.rm = TRUE),
    observed_total  = sum(observed, na.rm = TRUE),
    predicted_total = sum(predicted, na.rm = TRUE),
    residual_total  = sum(residual, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    year_label = YEAR_LABEL,
    outcome    = OUTCOME
  )

write.csv(pred_yearly,
          file.path(TAB_DIR, "powiat_predictions_yearly.csv"),
          row.names = FALSE)
write.csv(pred_powiat,
          file.path(TAB_DIR, "powiat_predictions_summary.csv"),
          row.names = FALSE)

cat("  zapisano: powiat_predictions_yearly.csv, powiat_predictions_summary.csv\n")

# ============================================================
# 10. ZAPIS
# ============================================================

saveRDS(
  list(
    m             = m,
    fml           = fml,
    dat_all       = dat_all,
    between_terms = between_terms,
    within_terms  = within_terms,
    all_terms     = all_terms,
    edf_tbl       = edf_tbl,
    aps_all       = aps_all,
    aps_county_all = aps_county_all,
    aps_grid_county_all = aps_grid_county_all,
    curves_all    = curves_all,
    effect_maps   = effect_maps,
    gi_all        = gi_all,
    pred_yearly   = pred_yearly,
    pred_powiat   = pred_powiat,
    VAR_LABELS    = VAR_LABELS,
    settings      = list(
      OUTCOME           = OUTCOME,
      YEAR_MIN          = YEAR_MIN,
      YEAR_MAX          = YEAR_MAX,
      PERIOD            = PERIOD,
      YEAR_LABEL        = YEAR_LABEL,
      WITHIN_LAG_N      = WITHIN_LAG_N,
      MODEL_VARS        = MODEL_VARS,
      ONLY_BETWEEN_VARS = ONLY_BETWEEN_VARS,
      GAMMA             = GAMMA,
      FAMILY            = deparse(FAMILY),
      K_GLOBAL          = K_GLOBAL,
      K_SVC             = K_SVC,
      K_BASELINE        = K_BASELINE,
      K_YEAR            = K_YEAR,
      N_GRID            = N_GRID,
      AME_EPS           = AME_EPS,
      APS_DEFINITION    = "Average finite-difference slope of the reference partial fitted surface over P2-P98; non-focal determinant terms set to zero, county locations retained, global smooth plus retained SVC included"
    )
  ),
  file.path(OUT_DIR, "gam_results.rds")
)

cat("\n=== GAM DONE ===\n")
cat("Period:", PERIOD, "|", YEAR_LABEL, "\n")
cat("Output:", OUT_DIR, "\n")
