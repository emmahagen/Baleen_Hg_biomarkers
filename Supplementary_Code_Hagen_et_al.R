# Supplementary Code (Hagen et al.)

# ———————————————————————————————————————————————————————————————
#  Load libraries and data file
# ———————————————————————————————————————————————————————————————

library(nlme)
library(ggplot2)
library(dplyr)
library(broom)
library(mgcv)
library(itsadug)  
library(gratia)
library(performance)


df <- read.csv("data/Supplementary_Data_Hagen_et_al.csv", fileEncoding = "UTF-8-BOM")


# ———————————————————————————————————————————————————————————————
#  Make sure data is in proper format
# ———————————————————————————————————————————————————————————————

df$ID <- as.factor(df$ID)
df$cm <- as.numeric(df$cm)
df$negs_to_zero_THg <- as.numeric(df$THg_negs_to_zero_extract_eff_corrected)


# Organize the data frame to group by whale, ordering their plates by gumline distance (cm)

df <- df %>% arrange(ID, cm)


df <- df %>%
  group_by(ID) %>%
  mutate(
    # Subtracting from max(time) sets each individual's final timepoint to 0
    cm_countdown = max(cm) - cm
  ) %>%
  ungroup()

# Methods 2.5.1--------------------------------------------------
# ———————————————————————————————————————————————————————————————
#  Longitudinal THg patterns
# ———————————————————————————————————————————————————————————————

# raw negative-corrected THg data

df$negs_zero_log_THg <- log(df$negs_to_zero_THg)

negs_zero_model <- lme(
  negs_to_zero_THg ~ cm,
  random = ~1 | ID,
  data = df)


# log1p-transformed negative-corrected THg data

df$log1p_negs_to_zero_THg <- log1p(df$negs_to_zero_THg)

ntz_log1p_model <- lme(
  log1p_negs_to_zero_THg ~ cm,
  random = ~1 | ID,
  data = df)


# ———————————————————————————————————————————————————————————————
# Normality tests (See Methods 2.4)
# ———————————————————————————————————————————————————————————————

# Q-Q plots

qqnorm(negs_zero_model)
qqnorm(ntz_log1p_model)


# Shapiro-Wilk tests

shapiro.test(residuals(negs_zero_model))
shapiro.test(residuals(ntz_log1p_model))


# Results: residuals normal across both models, but especially log1p:
  # negs_zero_model (p = 0.155); ntz_log1p_model (p = 0.9905)
  # Thus proceeded with log1p transformed data
# ———————————————————————————————————————————————————————————————————————
# Check models for autocorrelation
# ———————————————————————————————————————————————————————————————————————

negs_zero_model_res <- resid(negs_zero_model, type = "normalized")
acf(negs_zero_model_res, plot = FALSE) ## check the acf values
acf(negs_zero_model_res, main="log(THg) ~ cm + (cm|ID)") ## check the acf plot

ntz_log1p_model_res <- resid(ntz_log1p_model, type = "normalized")
acf(ntz_log1p_model_res, plot = FALSE) ## check the acf values
acf(ntz_log1p_model_res, main="log1p(THg) ~ cm + (cm|ID)") ## check the acf plot

# ———————————————————————————————————————————————————————————————————————
# AR1 model for autocorrelation
# ———————————————————————————————————————————————————————————————————————


# log1p transformed THg + AR1 (see Methods 2.4.1)

AR1_ntz_log1p_model <- lme(
  log1p_negs_to_zero_THg ~ cm_countdown,
  random = ~1 | ID,
  correlation = corAR1(form = ~ cm | ID),
  data = df)

summary(AR1_ntz_log1p_model)
anova(AR1_ntz_log1p_model)
print(AR1_ntz_log1p_model)
AR1_ntz_log1p_model_res <- resid(AR1_ntz_log1p_model, type = "normalized")
acf(AR1_ntz_log1p_model_res, plot = FALSE) ## check the acf values
acf(AR1_ntz_log1p_model_res, main="log1p(THg) ~ cm + (cm|ID) + AR(1)") ## check the acf plot
r2(AR1_ntz_log1p_model)

# Results: Phi = 0.694 & THg & cm association still p <.0001 (See Results 3.1)


qqnorm(resid(AR1_ntz_log1p_model)) 
shapiro.test(resid(AR1_ntz_log1p_model))

# Results: passed normality test (p = 0.556)   

# ———————————————————————————————————————————————————————————————
#  Population- & individual-level THg trends
# ———————————————————————————————————————————————————————————————

# Population trend

summary(ntz_log1p_model)

# Results: p < 0.00001 (See Results 3.1, Table 1)


# Per-idividual trends

log1p_ntz_whale_results <- df %>%
  group_by(ID) %>%
  do({
    fit <- lm(log1p_negs_to_zero_THg ~ cm_countdown, data = .)
    slope <- tidy(fit) %>% filter(term=="cm_countdown")
    data.frame(
      intercept = coef(fit)[1],
      std.Err = summary(fit)$coefficients[1, "Std. Error"],
      beta = slope$estimate,
      std.Err = summary(fit)$coefficients[2, "Std. Error"],
      p_value = slope$p.value,
      r2 = summary(fit)$r.squared)})


print(log1p_ntz_whale_results)

# Results: all whales p < 0.001, r2 > 0.5 (See Results 3.1, Table 1)


# Methods 2.5.2----------------------------------------------------------
# ———————————————————————————————————————————————————————————————————————
# THg associations with stable isotope and hormone biomarkers
# ———————————————————————————————————————————————————————————————————————

## Test AR(1)
# we first test the full model against a model that accounts for autocorrelation

full_model <- lme(
  log1p_negs_to_zero_THg ~ cm + dN15 + dC13 + log_T3 +
    log_Cortisol,
  random = ~1 | ID,
  method = "ML",
  data = df
)

summary(full_model)
anova(full_model)
shapiro.test(residuals(full_model))
res <- resid(full_model, type = "normalized")
acf(res, plot = FALSE) ## check the acf values
acf(res, main="Full model without AR1") ## check the acf plot

# Add AR1

full_model_AR1 <- lme(
  log1p_negs_to_zero_THg ~ cm + dN15 + dC13 + log_T3 +
    log_Cortisol,
  random = ~1 | ID,
  correlation = corAR1(form = ~ cm | ID, value = 0.231),
  method = "ML",
  data = df
)

summary(full_model_AR1)
anova(full_model_AR1)
shapiro.test(residuals(full_model_AR1))
res_ar1 <- resid(full_model_AR1, type = "normalized")
acf(res_ar1, main="Full model with AR1") ## check that including the acf has improved the fit of the model

## compare the AIC of the two models as another way to check that it has improved
## model fit. Model with the AR(1) has a lower AIC and reduced autocorrelation
## USE AR(1) AS THE BASE MODEL

## Full Model 

full_model <- lme(
  log1p_negs_to_zero_THg ~ cm + dN15 + dC13 + log_T3 +
    log_Cortisol,
  random = ~1 | ID,
  method = "ML",
  correlation = corAR1(form = ~ cm | ID, value = 0.231),
  data = df)


# Candidate predictors
predictors <- c(
  "cm",
  "dN15",
  "dC13",
  "log_T3",
  "log_Cortisol")


# To be safe:
# Keep only observations with complete data for all candidate variables
aic_df <- df[complete.cases(df[, c("log1p_negs_to_zero_THg", predictors, "ID")]), ]

# Check sample size
nrow(aic_df)
length(unique(aic_df$ID))


# Fit all possible models


# Create all possible combinations of predictors
model_list <- list()

for (k in 0:length(predictors)) {
  combinations <- combn(predictors, k, simplify = FALSE)
  for (vars in combinations) {
    if (length(vars) == 0) {
      formula <- log1p_negs_to_zero_THg ~ 1
    } else {
      formula <- as.formula(
        paste(
          "log1p_negs_to_zero_THg ~",
          paste(vars, collapse = " + ")))}
    
    model_name <- if (length(vars) == 0) {
      "null"
    } else {
      paste(vars, collapse = " + ")}
    
    model_list[[model_name]] <- lme(
      formula,
      random = ~1 | ID,
      correlation = corAR1(form = ~ cm | ID, value = 0.231),
      method = "ML",
      data = aic_df)}}



# Compare AIC

AIC_comparison <- data.frame(
  model = names(model_list),
  df     = sapply(model_list, function(m) attr(logLik(m), "df")),
  logLik = sapply(model_list, function(m) as.numeric(logLik(m))),
  AIC = sapply(model_list, AIC)
)


AIC_comparison <- AIC_comparison[order(AIC_comparison$AIC), ]

# Calculate delta AIC
AIC_comparison$Delta_AIC <- AIC_comparison$AIC -
  min(AIC_comparison$AIC)

# Add AIC weights
AIC_comparison$AIC_weight <- exp(-0.5 * AIC_comparison$Delta_AIC) /
  sum(exp(-0.5 * AIC_comparison$Delta_AIC))

# View results
print(AIC_comparison, row.names = FALSE)

## the model with the best fit was:

final_model <- lme(
  log1p_negs_to_zero_THg ~ cm_countdown + dN15,
  random = ~1 | ID,
  correlation = corAR1(form = ~ cm | ID, value = 0.231),
  method = "REML",
  data = df
)

summary(final_model)
anova(final_model)
shapiro.test(residuals(final_model))
final_model_res <- resid(final_model, type = "normalized")
acf(final_model_res, plot = FALSE) ## check the acf values
acf(final_model_res, main="log1p(THg) ~ cm + δ15N + (cm|ID) + AR(1)") ## check the acf plot
r2(final_model)


# Methods 2.5.3----------------------------------------------------------
# ———————————————————————————————————————————————————————————————————————
# THg and cortisol inflection point analysis
# ———————————————————————————————————————————————————————————————————————
# ———————————————————————————————————————————————————————————————————————
# 1. Data assumptions
# ———————————————————————————————————————————————————————————————————————
# df has columns:
#   THg           : numeric response (variable 1), MODELED ON THE LOG1P SCALE
#   cortisol      : numeric response (variable 2), MODELED ON THE LOG1P SCALE
#   cm            : numeric predictor (e.g. body length / position)
#   ID            : factor, individual identifier (random effect grouping)
#
df$ID <- factor(df$ID)

# ———————————————————————————————————————————————————————————————————————
# 1a. cm's chronological direction: 0 = present, larger = further in the
# PAST.
# ———————————————————————————————————————————————————————————————————————
# Things in the script that handle the fact that cm runs
# BACKWARD through time:
#   1. AR1 chronological order: since forward-in-time means DESCENDING
#      cm (past/large -> present/0), run_change_analysis() below sorts
#      each individual by `desc(cm)` (not ascending cm) before building
#      AR.start, so bam()'s rho correction is applied in true
#      chronological order.
#   2. Display: every plot with cm on the x-axis uses scale_x_reverse(),
#      so larger cm (past) appears on the LEFT and cm = 0 (present)
#      appears on the RIGHT
#   3. Derivative sign: rate-of-change is computed with respect to cm,
#      but moving forward through TIME is the opposite
#      direction (decreasing cm). So immediately after each derivative is
#      computed, its sign is flipped so that a
#      positive value means "increasing as time moves forward" — the
#      ordinary reading — while the plotted cm axis itself keeps its
#      original values and is just displayed in reverse.

# ———————————————————————————————————————————————————————————————————————
# 2. run_change_analysis(): fit model + diagnostics + derivative plots
#    for a single response variable
# ———————————————————————————————————————————————————————————————————————
# `transform` tells the function what was done to the response BEFORE
# modeling, so it knows how to invert it for plotting:
#   "log1p" -> response column holds log1p(original); back-transform via expm1
#   "none"  -> response column already on the original scale; no change
#
# This is a function that wraps the analysis, model checking and plotting in a 
# single function so that all we need to do is change the variable and we can 
# run the analysis all together
run_change_analysis <- function(response, data, n_grid_pop = 200, n_grid_id = 100,
                                eps = 1e-7, transform = c("none", "log1p", "log1p")) {
  
  transform <- match.arg(transform)
  inv_link <- switch(transform,
                     log1p = expm1,
                     none  = identity)
  deriv_scale_fun <- switch(transform,
                            log1p = exp,   
                            none  = function(f) 1)
  
  ## 2a. Drop rows with NA/NaN/Inf in the columns 
  needed_cols <- c(response, "cm", "ID")
  bad <- !stats::complete.cases(data[, needed_cols]) |
    !is.finite(data[[response]]) | !is.finite(data$cm)
  n_dropped <- sum(bad)
  if (n_dropped > 0) {
    warning(sprintf(
      "run_change_analysis('%s'): dropping %d of %d rows with NA/NaN/Inf in %s or cm",
      response, n_dropped, nrow(data), response
    ))
  }
  data <- data[!bad, , drop = FALSE]
  data <- droplevels(data)
  
  ## Sort in TRUE chronological (forward-time) order
  data <- data %>% arrange(ID, desc(cm))
  data$start_event <- with(data, ID != dplyr::lag(ID, default = dplyr::first(ID)))
  data$start_event[1] <- TRUE
  
  ## original-scale response, computed once here so it's available for
  ## standardizing fitted-curve comparisons later (and for the raw-data
  ## overlay on individual trajectory plots)
  data$.response_orig <- inv_link(data[[response]])
  mean_response_orig <- mean(data$.response_orig, na.rm = TRUE)
  sd_response_orig   <- sd(data$.response_orig, na.rm = TRUE)
  
  form_uncorrected <- as.formula(
    paste0(response, " ~ s(cm, k = 20, bs = 'tp') + s(cm, ID, bs = 'fs', k = 5, m = 2)")
  )
  
  # ———————————————————————————————————————————————————————————————————————
  # 2b. Estimate rho, then fit with AR1 correction
  # ———————————————————————————————————————————————————————————————————————
  
  m0 <- bam(form_uncorrected, data = data, method = "fREML")
  
  rho_est <- start_value_rho(m0, plot = FALSE)   # itsadug helper: AR.start-aware
  
  m_alt <- bam(form_uncorrected,
               data     = data,
               method   = "fREML",
               rho      = rho_est,
               AR.start = data$start_event)
  
  cat("\n============================================================\n")
  cat("Response variable:", response, " (transform:", transform, ")\n")
  cat("Estimated rho:", round(rho_est, 3), "\n")
  cat("============================================================\n")
  print(summary(m_alt))
  gam.check(m_alt)   # basis dimension check (k.check) for both smooth terms
  
  # ———————————————————————————————————————————————————————————————————————
  # 2c. Diagnostics: autocorrelation + normality + variance
  # ———————————————————————————————————————————————————————————————————————
  acf_resid(m_alt, split_pred = list(data$ID),
            main = paste("ACF of residuals (AR-corrected) —", response))
  acf_resid(m0, split_pred = list(data$ID),
            main = paste("ACF of residuals (uncorrected) —", response))
  
  data$.res <- resid(m_alt)
  id_acf1 <- data %>%
    group_by(ID) %>%
    group_map(~ acf(.x$.res, plot = FALSE)$acf[2]) %>%
    unlist()
  cat("Per-individual lag-1 autocorrelation (should hover near 0):\n")
  print(id_acf1)
  
  res <- resid(m_alt)
  qqnorm(res, main = paste("QQ plot of residuals —", response))
  qqline(res, col = "red")
  hist(res, breaks = 30, main = paste("Histogram of residuals —", response),
       xlab = "Residual")
  print(shapiro.test(res))   # formal test; treat QQ plot as primary evidence,
  ## this as a secondary/confirmatory check 
  
  plot(fitted(m_alt), res,
       xlab = "Fitted values", ylab = "Residuals",
       main = paste("Residuals vs fitted —", response))
  abline(h = 0, lty = 2, col = "red")
  
  print(gratia::appraise(m_alt))
  
  # ———————————————————————————————————————————————————————————————————————
  # 3. Population-level derivative
  # ———————————————————————————————————————————————————————————————————————
  sd_response <- sd(data[[response]], na.rm = TRUE)
  
  d <- derivatives(m_alt, select = "s(cm)", n = n_grid_pop, eps = eps,
                   level = 0.95, interval = "confidence") %>%
    mutate(
      .derivative = -.derivative,
      .lower_ci   = .derivative - 1.96 * .se,
      .upper_ci   = .derivative + 1.96 * .se,
      .derivative_std = .derivative / sd_response,
      .lower_ci_std   = .lower_ci   / sd_response,
      .upper_ci_std   = .upper_ci   / sd_response
    )
  
  deriv_grid_pred <- data.frame(cm = d$cm, ID = levels(data$ID)[1])
  fit_log_at_deriv <- predict(m_alt, newdata = deriv_grid_pred,
                              type = "iterms", terms = "s(cm)")[, "s(cm)"]
  
  scale_factor <- deriv_scale_fun(fit_log_at_deriv)   # = 1 if transform == "none"
  
  d <- d %>%
    mutate(
      .derivative_orig = .derivative * scale_factor,
      .se_orig          = .se * scale_factor,
      .lower_ci_orig    = .derivative_orig - 1.96 * .se_orig,
      .upper_ci_orig    = .derivative_orig + 1.96 * .se_orig,
      signif    = !(.lower_ci <= 0 & .upper_ci >= 0),
      direction = case_when(
        signif & .derivative > 0 ~ "increasing",
        signif & .derivative < 0 ~ "decreasing",
        TRUE                     ~ "no change"
      )
    )
  
  # ———————————————————————————————————————————————————————————————————————
  # 4. Population-level plots, ORIGINAL SCALE
  # ———————————————————————————————————————————————————————————————————————
  p_deriv_pop <- ggplot(d, aes(cm, .derivative_orig)) +
    geom_ribbon(aes(ymin = .lower_ci_orig, ymax = .upper_ci_orig), alpha = 0.2) +
    geom_line() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_point(data = filter(d, direction == "increasing"), colour = "blue", size = 1) +
    geom_point(data = filter(d, direction == "decreasing"), colour = "red", size = 1) +
    scale_x_reverse() +
    labs(x = "Distance along growth axis (cm)") +
    theme_bw()
  print(p_deriv_pop)
  
  
  # ———————————————————————————————————————————————————————————————————————
  # 5. Per-individual derivative + fitted-trajectory plots, ORIGINAL
  #    SCALE (population + individual deviation combined
  # ———————————————————————————————————————————————————————————————————————
  
  id_grids <- data %>%
    group_by(ID) %>%
    summarise(cm_min = min(cm), cm_max = max(cm), .groups = "drop") %>%
    group_by(ID) %>%
    reframe(cm = seq(cm_min, cm_max, length.out = n_grid_id))
  
  X0 <- predict(m_alt, newdata = id_grids, type = "lpmatrix")
  grid_eps <- id_grids
  grid_eps$cm <- grid_eps$cm + eps
  X1 <- predict(m_alt, newdata = grid_eps, type = "lpmatrix")
  Xp <- (X1 - X0) / eps
  
  term_cols <- grep("^s\\(cm\\)|^s\\(cm,ID\\)", colnames(X0))
  Xi <- matrix(0, nrow(Xp), ncol(Xp), dimnames = list(NULL, colnames(Xp)))
  Xi[, term_cols] <- Xp[, term_cols]
  
  beta <- coef(m_alt)
  Vp   <- vcov(m_alt)   
  
  id_grids <- id_grids %>%
    mutate(
      .derivative = -as.numeric(Xi %*% beta),
      .se         = sqrt(rowSums((Xi %*% Vp) * Xi)),
      .lower_ci   = .derivative - 1.96 * .se,
      .upper_ci   = .derivative + 1.96 * .se,
      signif      = !(.lower_ci <= 0 & .upper_ci >= 0),
      direction   = case_when(
        signif & .derivative > 0 ~ "increasing",
        signif & .derivative < 0 ~ "decreasing",
        TRUE                     ~ "no change"
      ),
      .derivative_std = .derivative / sd_response,
      .lower_ci_std   = .lower_ci   / sd_response,
      .upper_ci_std   = .upper_ci   / sd_response
    )
  
  fv_id <- predict(m_alt, newdata = id_grids, type = "response", se.fit = TRUE)
  id_grids <- id_grids %>%
    mutate(fit_log = fv_id$fit, se_log = fv_id$se.fit,
           upper_log = fit_log + 1.96 * se_log, lower_log = fit_log - 1.96 * se_log,
           fit   = inv_link(fit_log),
           upper = inv_link(upper_log),
           lower = inv_link(lower_log),
           fit_std   = (fit   - mean_response_orig) / sd_response_orig,
           upper_std = (upper - mean_response_orig) / sd_response_orig,
           lower_std = (lower - mean_response_orig) / sd_response_orig)
  
  id_scale_factor <- deriv_scale_fun(id_grids$fit_log)
  id_grids <- id_grids %>%
    mutate(
      .derivative_orig = .derivative * id_scale_factor,
      .se_orig          = .se * id_scale_factor,
      .lower_ci_orig    = .derivative_orig - 1.96 * .se_orig,
      .upper_ci_orig    = .derivative_orig + 1.96 * .se_orig
    )
  
  p_deriv_id <- ggplot(id_grids, aes(cm, .derivative_orig)) +
    geom_ribbon(aes(ymin = .lower_ci_orig, ymax = .upper_ci_orig), alpha = 0.2) +
    geom_line() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_point(data = filter(id_grids, direction == "increasing"), colour = "blue", size = 0.8) +
    geom_point(data = filter(id_grids, direction == "decreasing"), colour = "red", size = 0.8) +
    facet_wrap(~ ID, scales = "free_y", nrow = 1) +
    scale_x_reverse() +
    labs(x = "Distance along growth axis (cm)") +
    theme_bw()
  print(p_deriv_id)
  
  
  ## return everything for later reuse (e.g. cross-hormone comparison)
  list(
    m0 = m0, m_alt = m_alt, rho = rho_est, sd_response = sd_response,
    mean_response_orig = mean_response_orig, sd_response_orig = sd_response_orig,
    transform = transform,
    pop_derivatives = d,
    id_derivatives  = id_grids,
    plots = list(deriv_pop = p_deriv_pop,
                 deriv_id = p_deriv_id)
  )
}

# ———————————————————————————————————————————————————————————————————————
# SECTION A — THg ~ cm  (modeled as log1p(THg))
# ———————————————————————————————————————————————————————————————————————
THg_results <- run_change_analysis("log1p_negs_to_zero_THg", df, transform = "log1p")

# relabels the plots
THg_results$plots$deriv_pop +
  labs(y = "Rate of change of THg (µg/g) per cm") 
THg_results$plots$deriv_id +
  labs(y = "Rate of change of \nTHg (µg/g) per cm")

# ———————————————————————————————————————————————————————————————————————
# SECTION B — cortisol ~ cm  (modeled as log1p(cortisol))
# ———————————————————————————————————————————————————————————————————————
cortisol_results <- run_change_analysis("log_Cortisol", df, transform = "log1p")

# relabels the plots
cortisol_results$plots$deriv_pop +
  labs(y = "Rate of change of Cortisol (ng/g) per cm")
cortisol_results$plots$deriv_id +
  labs(y = "Rate of change of \nCortisol (ng/g) per cm")

# ———————————————————————————————————————————————————————————————————————
# Remove Er5
# ———————————————————————————————————————————————————————————————————————
## to test if the model changes without the inclusion of Er5

clean_data <- subset(df, ID != "Er_5")

# ———————————————————————————————————————————————————————————————————————
# SECTION C — THg ~ cm  (modeled as log1p(THg))
# ———————————————————————————————————————————————————————————————————————
clean_THg_results <- run_change_analysis("log1p_negs_to_zero_THg", clean_data, transform = "log1p")

# relabels the plots
clean_THg_results$plots$deriv_pop +
  labs(y = "Rate of change of THg (µg/g) per cm")
clean_THg_results$plots$deriv_id +
  labs(y = "Rate of change of \nTHg (µg/g) per cm")

# ———————————————————————————————————————————————————————————————————————
# SECTION D — cortisol ~ cm  (modeled as log1p(cortisol))
# ———————————————————————————————————————————————————————————————————————
clean_cortisol_results <- run_change_analysis("log_Cortisol", clean_data, transform = "log1p")

# relabels the plots
clean_cortisol_results$plots$deriv_pop +
  labs(y = "Rate of change of Cortisol (ng/g) per cm")
clean_cortisol_results$plots$deriv_id +
  labs(y = "Rate of change of \nCortisol (ng/g) per cm")


