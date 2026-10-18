# ============================================================
# 01_mundlak_ridge.R — linear Mundlak ridge benchmark
# v40.5: M1-reference scaling + recentered dependent multiplier CI
# ============================================================

if (!exists("SCRIPTS_DIR")) SCRIPTS_DIR <- file.path(getwd(), "R")
if (!exists("PERIOD_SPECS")) source(file.path(SCRIPTS_DIR, "00_config.R"), local = FALSE)
if (!exists("build_mundlak_panel")) source(file.path(SCRIPTS_DIR, "00_helpers.R"), local = FALSE)

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(tibble)
  library(stringr); library(glmnet); library(readr)
})

set.seed(20260908)

OUT_DIR <- file.path(OUT_BASE, "mundlak_ridge")
TAB_DIR <- file.path(OUT_DIR, "tables")
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)

cat("=== MUNDLAK RIDGE v40.5 ===\n")
cat(MODEL_LABEL, "|", YEAR_LABEL, "| controls:",
    ifelse(length(CONTROL_VARS), paste(CONTROL_VARS, collapse = ", "), "none"), "\n")
cat("Final scaling:", SCALING_MODE, "| reference:", SCALING_REFERENCE_PERIOD, "\n")

prep <- build_mundlak_panel(dat, period_id = PERIOD, scaling_mode = SCALING_MODE)
dat_m <- prep$data
between_terms <- prep$between_terms
within_terms <- prep$within_terms
x_vars <- prep$all_terms

readr::write_csv(prep$scaler_table, file.path(TAB_DIR, "scaling_parameters.csv"))

d_fit <- dat_m %>%
  model_complete_cases(OUTCOME, x_vars, CONTROL_VARS, require_xy = FALSE) %>%
  arrange(code, year)

if (nrow(d_fit) == 0L) stop("Brak kompletnych obserwacji do Mundlak ridge.")
cat("n model =", nrow(d_fit), "| counties =", n_distinct(d_fit$code),
    "| years =", min(d_fit$year), "–", max(d_fit$year), "\n")

control_formula <- function() {
  rhs <- c("year_c", CONTROL_VARS)
  as.formula(paste("~", paste(rhs, collapse = " + ")))
}

residualize <- function(d, y_var, x_vars) {
  Z <- model.matrix(control_formula(), data = d)
  y <- as.numeric(d[[y_var]])
  X <- as.matrix(d[, x_vars, drop = FALSE])
  qrZ <- qr(Z)
  list(
    y_res = as.numeric(y - qr.fitted(qrZ, y)),
    X_res = X - qr.fitted(qrZ, X),
    y = y, X = X, Z = Z
  )
}

choose_alpha <- function(X, y, lambda, grid = seq(0.1, 1.0, 0.1)) {
  n <- nrow(X); p <- ncol(X)
  S <- crossprod(X) / n
  M <- solve(S + lambda * diag(p))
  beta_r <- as.numeric(M %*% (crossprod(X, y) / n))
  loo_mse <- vapply(grid, function(a) {
    beta_d <- beta_r + a * lambda * as.numeric(M %*% beta_r)
    resid <- as.numeric(y - X %*% beta_d)
    H <- rowSums((X %*% M) * X) / n
    mean((resid / (1 - pmin(H, 0.99)))^2, na.rm = TRUE)
  }, numeric(1))
  grid[which.min(loo_mse)]
}

ridge_solve <- function(X, y, lambda, alpha = NULL) {
  n <- nrow(X); p <- ncol(X)
  S <- crossprod(X) / n
  M <- solve(S + lambda * diag(p))
  beta_r <- as.numeric(M %*% (crossprod(X, y) / n))
  names(beta_r) <- colnames(X)
  if (is.null(alpha)) alpha <- choose_alpha(X, y, lambda)
  Cmat <- diag(p) + alpha * lambda * M
  beta_d <- as.numeric(Cmat %*% beta_r)
  names(beta_d) <- colnames(X)
  list(beta_ridge = beta_r, beta_debiased = beta_d, alpha = alpha, M = M, Cmat = Cmat)
}

rz <- residualize(d_fit, OUTCOME, x_vars)

# Foldy po powiatach: wszystkie lata tego samego powiatu są w tym samym foldzie.
codes <- sort(unique(d_fit$code))
fold_map <- tibble(code = codes,
                   fold = sample(rep(seq_len(CV_FOLDS), length.out = length(codes))))
foldid <- fold_map$fold[match(d_fit$code, fold_map$code)]

cvfit <- glmnet::cv.glmnet(
  rz$X_res, rz$y_res,
  alpha = 0,
  family = "gaussian",
  standardize = FALSE,
  intercept = FALSE,
  foldid = foldid
)
lambda <- cvfit[[LAMBDA_RULE]]
fit <- ridge_solve(rz$X_res, rz$y_res, lambda)
cat(sprintf("lambda = %.6f | correction alpha = %.2f\n", lambda, fit$alpha))

beta_hat <- fit$beta_debiased[x_vars]

# Kontrolna część pełnego fitted — tylko do diagnostyki modelu.
gamma_hat <- qr.coef(qr(rz$Z), rz$y - as.numeric(rz$X %*% beta_hat))
gamma_hat[!is.finite(gamma_hat)] <- 0
fitted_vals <- as.numeric(rz$Z %*% gamma_hat + rz$X %*% beta_hat)
resids <- rz$y - fitted_vals

# ------------------------------------------------------------
# Dependent multiplier bootstrap dla liniowego operatora ridge
# ------------------------------------------------------------
# Poprzednia wersja generowała y* z beta_debiased i ponownie nakładała
# shrinkage, przez co rozkład bootstrap był przesunięty względem estymaty.
# Tutaj perturbujemy score w przestrzeni zresidualizowanej, przy stałych
# lambda i alpha. To daje bootstrap różnic beta* - beta_hat, który następnie
# jawnie recentrujemy przed wyznaczeniem przedziałów.
gen_dw <- function(n, rho) {
  if (n <= 0L) return(numeric(0))
  w <- numeric(n)
  w[1] <- stats::rnorm(1)
  if (n >= 2L) {
    for (tt in 2:n) {
      w[tt] <- rho * w[tt - 1] + sqrt(max(0, 1 - rho^2)) * stats::rnorm(1)
    }
  }
  w
}

n_obs <- nrow(rz$X_res)
resid_res <- as.numeric(rz$y_res - rz$X_res %*% beta_hat)

cat("Bootstrap B =", BOOT_B, "| type = recentered dependent multiplier ...\n")
boot_delta <- purrr::map_dfr(seq_len(BOOT_B), function(b) {
  if (b %% 100L == 0L) cat("  bootstrap", b, "/", BOOT_B, "\n")
  wt <- d_fit %>%
    transmute(.row_id = dplyr::row_number(), code, year) %>%
    group_by(code) %>%
    mutate(w = gen_dw(n(), BOOT_RHO)) %>%
    ungroup() %>%
    arrange(.row_id)
  score_b <- as.numeric(crossprod(rz$X_res, resid_res * wt$w) / n_obs)
  delta_r <- as.numeric(fit$M %*% score_b)
  delta_d <- as.numeric(fit$Cmat %*% delta_r)
  tibble(term = x_vars, delta = delta_d, boot_id = b)
})

if (n_distinct(boot_delta$boot_id) < max(20L, floor(0.9 * BOOT_B))) {
  stop("Za mało poprawnych replikacji bootstrap.")
}

boot_ci <- boot_delta %>%
  group_by(term) %>%
  mutate(delta_centered = delta - mean(delta, na.rm = TRUE)) %>%
  summarise(
    boot_sd = sd(delta_centered, na.rm = TRUE),
    q_low = quantile(delta_centered, 0.025, na.rm = TRUE),
    q_high = quantile(delta_centered, 0.975, na.rm = TRUE),
    bootstrap_center_shift = mean(delta, na.rm = TRUE),
    n_boot = n_distinct(boot_id),
    .groups = "drop"
  )

coef_final <- tibble(
  term = x_vars,
  beta_debiased = as.numeric(beta_hat),
  lambda = lambda,
  alpha_used = fit$alpha
) %>%
  left_join(boot_ci, by = "term") %>%
  mutate(
    conf.low = beta_debiased + q_low,
    conf.high = beta_debiased + q_high,
    component = comp_from_term(term),
    variable = clean_var(term),
    var_label = label_var(variable),
    direction = case_when(
      conf.high < 0 ~ "lower outcome",
      conf.low > 0 ~ "higher outcome",
      TRUE ~ "uncertain"
    ),
    period = PERIOD,
    year_label = YEAR_LABEL,
    scaling_mode = prep$scaling_mode,
    scaling_source = prep$scaler_source,
    bootstrap_type = BOOT_TYPE
  ) %>%
  select(-q_low, -q_high)

# Sanity check: po recentrowaniu estimate powinno leżeć w centralnym 95% CI
# poza ekstremalnie asymetrycznymi realizacjami. Nie przerywamy, ale zapisujemy flagę.
coef_final <- coef_final %>%
  mutate(estimate_inside_ci = beta_debiased >= conf.low & beta_debiased <= conf.high)
if (any(!coef_final$estimate_inside_ci, na.rm = TRUE)) {
  warning("Co najmniej jedna estymata leży poza multiplier CI; sprawdź coef_final.csv.")
}

readr::write_csv(coef_final, file.path(TAB_DIR, "coef_final.csv"))
readr::write_csv(fold_map, file.path(TAB_DIR, "ridge_cv_county_folds.csv"))
readr::write_csv(boot_delta, file.path(TAB_DIR, "bootstrap_delta_draws.csv"))

fit_summary <- tibble(
  period = PERIOD,
  year_label = YEAR_LABEL,
  n_model = nrow(d_fit),
  n_counties = n_distinct(d_fit$code),
  year_min_model = min(d_fit$year),
  year_max_model = max(d_fit$year),
  lambda = lambda,
  alpha_used = fit$alpha,
  bootstrap_B = BOOT_B,
  bootstrap_type = BOOT_TYPE,
  bootstrap_rho = BOOT_RHO,
  scaling_mode = prep$scaling_mode,
  scaling_source = prep$scaler_source,
  controls = paste(CONTROL_VARS, collapse = ";")
)
readr::write_csv(fit_summary, file.path(TAB_DIR, "ridge_fit_summary.csv"))

saveRDS(list(
  coef_final = coef_final,
  boot_delta = boot_delta,
  lambda = lambda,
  alpha = fit$alpha,
  x_vars = x_vars,
  dat_m = dat_m,
  dat_fit = d_fit,
  prep = prep,
  period = PERIOD,
  year_label = YEAR_LABEL,
  settings = list(
    OUTCOME = OUTCOME,
    CONTROL_VARS = CONTROL_VARS,
    WITHIN_LAG_N = WITHIN_LAG_N,
    BOOT_B = BOOT_B,
    BOOT_TYPE = BOOT_TYPE,
    SCALING_MODE = SCALING_MODE,
    SCALING_REFERENCE_PERIOD = SCALING_REFERENCE_PERIOD
  )
), file.path(OUT_DIR, "mundlak_results.rds"))

cat("\nMundlak ridge DONE:", OUT_DIR, "\n")
print(coef_final %>% select(var_label, component, beta_debiased, conf.low, conf.high,
                            estimate_inside_ci), n = 50)
