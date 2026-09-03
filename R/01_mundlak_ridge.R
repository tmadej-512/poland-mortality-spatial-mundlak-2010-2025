# ============================================================
# 01_mundlak_ridge.R
# Mundlak ridge + debiased + dependent wild bootstrap
#
# Wymaga: dat, powiaty_spdf w pamięci
# Wymaga: 00_config.R w tym samym katalogu
#
# Uruchom dwa razy:
#   PERIOD <- "pre_covid"  # w 00_config.R
#   PERIOD <- "full"
# ============================================================

# config jest sourcowany przez RUN_ALL przed kazdym krokiem
# przy recznym uruchomieniu: source("R/00_config.R")
if (!exists("LANCET") || !exists("YEAR_MAX")) {
  source(file.path(getwd(), "R", "00_config.R"))
}

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr)
  library(stringr); library(ggplot2); library(tibble)
  library(glmnet); library(forcats); library(scales)
})

set.seed(20240525)

OUT_DIR <- file.path(OUT_BASE, "mundlak_ridge")
FIG_DIR <- file.path(OUT_DIR, "figures")
TAB_DIR <- file.path(OUT_DIR, "tables")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)

cat("=== Mundlak ridge ===\n")
cat("Period:", PERIOD, "|", YEAR_LABEL, "\n")
validate_period_data(dat)

# ============================================================
# 1. DANE — panel Mundlak
# ============================================================

all_vars    <- unique(c(MODEL_VARS, ONLY_BETWEEN_VARS))
all_vars    <- intersect(all_vars, names(dat))
within_vars <- setdiff(all_vars,
                        intersect(ONLY_BETWEEN_VARS,
                                  all_vars))

dat_m <- dat %>%
  filter(year >= YEAR_MIN, year <= YEAR_MAX,
         !is.na(.data[[OUTCOME]])) %>%
  mutate(code   = as.character(code),
         year   = as.integer(year),
         year_c = year - mean(year, na.rm = TRUE)) %>%
  select(code, name, year, year_c,
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
    across(paste0(within_vars, "_within"),
           ~ dplyr::lag(.x, n = WITHIN_LAG_N),
           .names = paste0("{.col}_lag", WITHIN_LAG_N))
  ) %>%
  ungroup()

between_raw <- paste0(all_vars, "_between")
within_raw  <- paste0(within_vars,
                       "_within_lag", WITHIN_LAG_N)
terms_scale <- c(between_raw, within_raw)

dat_m <- dat_m %>%
  mutate(across(
    all_of(intersect(terms_scale, names(.))),
    scale_safe,
    .names = "{.col}_z"
  )) %>%
  filter(!is.na(.data[[OUTCOME]]))

between_terms_m <- paste0(all_vars, "_between_z")
within_terms_m  <- paste0(within_vars,
                            "_within_lag", WITHIN_LAG_N, "_z")
x_vars <- intersect(c(between_terms_m, within_terms_m),
                     names(dat_m))

cat("n =", nrow(dat_m),
    "| powiaty =", n_distinct(dat_m$code),
    "| x_vars =", length(x_vars), "\n")

# ============================================================
# 2. RIDGE + DEBIASING
# ============================================================

choose_alpha <- function(X, y, lambda,
                          grid = seq(0.1, 1.0, 0.1)) {
  n <- nrow(X); p <- ncol(X)
  S <- crossprod(X) / n
  M <- solve(S + lambda * diag(p))
  beta_r <- as.numeric(M %*% (crossprod(X, y) / n))
  loo_mse <- map_dbl(grid, function(a) {
    beta_d <- beta_r + a * lambda *
      as.numeric(M %*% beta_r)
    resid <- y - X %*% beta_d
    H <- rowSums((X %*% M) * X) / n
    mean((resid / (1 - pmin(H, 0.99)))^2,
         na.rm = TRUE)
  })
  grid[which.min(loo_mse)]
}

ridge_solve <- function(X, y, lambda, alpha = NULL) {
  n <- nrow(X); p <- ncol(X)
  S <- crossprod(X) / n
  M <- solve(S + lambda * diag(p))
  beta_r <- as.numeric(M %*% (crossprod(X, y) / n))
  names(beta_r) <- colnames(X)
  if (is.null(alpha))
    alpha <- choose_alpha(X, y, lambda)
  beta_d <- beta_r + alpha * lambda *
    as.numeric(M %*% beta_r)
  names(beta_d) <- colnames(X)
  list(beta_ridge    = beta_r,
       beta_debiased = beta_d,
       alpha         = alpha,
       S             = S)
}

residualize <- function(d, y_var, x_vars) {
  # Kontrola roku + covid_period tylko wtedy, gdy faktor ma >=2 obserwowane poziomy.
  # To zapobiega błędowi: "kontrasty można zastosować jedynie do czynników...".
  ctrl_vars <- if (use_covid_control(d)) {
    "~ year_c + covid_period"
  } else {
    "~ year_c"
  }
  fml <- as.formula(ctrl_vars)
  Z   <- model.matrix(fml, data = d)
  y   <- d[[y_var]]
  X   <- as.matrix(d %>% select(all_of(x_vars)))
  qrZ <- qr(Z)
  list(y_res = as.numeric(y - qr.fitted(qrZ, y)),
       X_res = X - qr.fitted(qrZ, X),
       y = y, X = X)
}

# ---- Fit ----
d_fit <- dat_m %>%
  mutate(
    covid_period = covid_period_factor(year)
  ) %>%
  mutate(covid_period = droplevels(covid_period)) %>%
  select(code, name, year, year_c, covid_period,
         all_of(OUTCOME), all_of(x_vars)) %>%
  mutate(across(all_of(x_vars), impute_median)) %>%
  drop_na(all_of(c(OUTCOME, "year_c")))

cat("covid_period w danych modelowych:
")
print(table(d_fit$covid_period, useNA = "ifany"))
if (PERIOD == "full" && !use_covid_control(d_fit)) {
  warning("PERIOD='full', ale covid_period ma mniej niż 2 obserwowane poziomy; kontrola covid_period zostanie pominięta.")
}

rz <- residualize(d_fit, OUTCOME, x_vars)

cvfit <- cv.glmnet(
  rz$X_res, rz$y_res,
  alpha      = 0,
  family     = "gaussian",
  standardize = FALSE,
  intercept  = FALSE,
  nfolds     = CV_FOLDS
)

lambda <- cvfit[[LAMBDA_RULE]]
cat(sprintf("lambda = %.4f\n", lambda))

fit <- ridge_solve(rz$X_res, rz$y_res, lambda)
cat(sprintf("alpha = %.2f\n", fit$alpha))

# ============================================================
# 3. BOOTSTRAP
# ============================================================

fitted_vals <- as.numeric(
  lm(as.formula(paste0(
    OUTCOME, " ~ year_c",
    if (use_covid_control(d_fit)) " + covid_period" else ""
  )), data = d_fit)$fitted.values
) + as.numeric(rz$X %*% fit$beta_debiased)

resids  <- d_fit[[OUTCOME]] - fitted_vals
codes_u <- unique(d_fit$code)

gen_dw <- function(n, rho) {
  w <- numeric(n)
  w[1] <- sample(c(-1, 1), 1)
  for (t in seq(2, n))
    w[t] <- rho * w[t-1] +
      sqrt(1 - rho^2) * sample(c(-1, 1), 1)
  w
}

cat("Bootstrap B =", BOOT_B, "...\n")

boot_draws <- map_dfr(seq_len(BOOT_B), function(b) {
  if (b %% 100 == 0)
    cat("  bootstrap", b, "/", BOOT_B, "\n")
  wt <- map_dfr(codes_u, function(cc) {
    rows <- d_fit[d_fit$code == cc, ]
    tibble(code = rows$code, year = rows$year,
           w    = gen_dw(nrow(rows), BOOT_RHO))
  })
  w_vec <- d_fit %>%
    left_join(wt, by = c("code", "year")) %>%
    pull(w)
  y_star <- fitted_vals + resids * w_vec
  rz2 <- residualize(
    d_fit %>% mutate(!!sym(OUTCOME) := y_star),
    OUTCOME, x_vars
  )
  fit2 <- tryCatch(
    ridge_solve(rz2$X_res, rz2$y_res,
                lambda, fit$alpha),
    error = function(e) NULL
  )
  if (is.null(fit2)) return(NULL)
  tibble(term    = x_vars,
         beta    = fit2$beta_debiased[x_vars],
         boot_id = b)
})

boot_ci <- boot_draws %>%
  left_join(
    tibble(term     = x_vars,
           beta_hat = fit$beta_debiased[x_vars]),
    by = "term"
  ) %>%
  group_by(term) %>%
  summarise(
    beta_hat  = first(beta_hat),
    boot_sd   = sd(beta, na.rm = TRUE),
    conf.low  = beta_hat - 1.96 * boot_sd,
    conf.high = beta_hat + 1.96 * boot_sd,
    .groups   = "drop"
  )

coef_final <- tibble(
  term          = x_vars,
  beta_debiased = fit$beta_debiased[x_vars],
  lambda        = lambda,
  alpha_used    = fit$alpha
) %>%
  left_join(boot_ci, by = "term") %>%
  mutate(
    component = case_when(
      str_detect(term, "_between_z$") ~ "between",
      TRUE                            ~ "within"
    ),
    variable  = clean_var(term),
    var_label = label_var(variable),
    direction = case_when(
      conf.high < 0 ~ "lower ASMR",
      conf.low  > 0 ~ "higher ASMR",
      TRUE         ~ "uncertain"
    ),
    period    = PERIOD
  )

write.csv(coef_final,
          file.path(TAB_DIR, "coef_final.csv"),
          row.names = FALSE)

cat("\nWyniki Mundlak ridge:\n")
print(coef_final %>%
        select(var_label, component,
               beta_debiased, conf.low,
               conf.high, direction),
      n = 60)

# ============================================================
# 4. FOREST PLOT
# ============================================================

var_order <- coef_final %>%
  filter(component == "between") %>%
  arrange(beta_debiased) %>%
  pull(var_label) %>%
  unique()

all_labels <- unique(coef_final$var_label)
var_order  <- c(var_order,
                setdiff(all_labels, var_order))

p_forest <- coef_final %>%
  mutate(
    var_label = factor(var_label, levels = var_order),
    pt_col = case_when(
      direction == "lower ASMR"  ~ LANCET$navy,
      direction == "higher ASMR" ~ LANCET$red,
      TRUE                        ~ LANCET$grey_mid
    ),
    pt_alpha = if_else(direction == "uncertain",
                        0.4, 1.0)
  ) %>%
  ggplot(aes(x = beta_debiased, y = var_label)) +
  geom_vline(xintercept = 0, linetype = "dashed",
             colour = LANCET$grey_mid, linewidth = 0.35) +
  geom_linerange(
    aes(xmin = conf.low, xmax = conf.high,
        alpha = pt_alpha),
    colour = LANCET$navy, linewidth = 0.85
  ) +
  geom_point(
    aes(colour = pt_col, alpha = pt_alpha),
    shape = 16, size = 2.2
  ) +
  scale_colour_identity() +
  scale_alpha_identity() +
  scale_y_discrete(limits = var_order, drop = FALSE) +
  facet_wrap(~ component, scales = "free_x") +
  labs(
    title    = paste0("Linear Mundlak β | ", YEAR_LABEL),
    subtitle = paste0("Debiased β; 95% CI from the dependent wild bootstrap (B=", BOOT_B, ")"),
    x        = "\u0394 ASMR per 100k at +1 SD",
    y        = NULL
  ) +
  theme_lancet(base_size = 10) +
  theme(
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold")
  )

save_hp_plot(file.path(FIG_DIR, "forest_mundlak.png"),
             p_forest, width = 12, height = 8)

# ============================================================
# 5. ZAPIS
# ============================================================

saveRDS(
  list(coef_final  = coef_final,
       boot_draws  = boot_draws,
       lambda      = lambda,
       x_vars      = x_vars,
       dat_m       = dat_m,
       period      = PERIOD,
       year_label  = YEAR_LABEL),
  file.path(OUT_DIR, "mundlak_results.rds")
)

cat("\n=== Mundlak ridge DONE ===\n")
cat("Period:", PERIOD, "|", YEAR_LABEL, "\n")
cat("Output:", OUT_DIR, "\n")
