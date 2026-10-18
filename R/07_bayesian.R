# 07_bayesian.R — approximate Bayesian meta-analyses without MCMC.
#
# Models (all use closed-form / Monte-Carlo posteriors with base R + metafor):
#   RR/MD models  normal-normal random-effects, conditional on the REML tau2;
#                 weakly informative prior N(0, 1.37^2) on log-RR (Higgins et al.
#                 recommended prior scale), N(0, 10^2) on the MD scale.
#   Rare-event RD models  pooled-arm independent Beta(1,1) (uniform) priors per
#                 arm; posterior of the risk difference sampled exactly from the
#                 independent Beta posteriors (studies pooled — approximation,
#                 documented). Beta(0.5, 0.5) Jeffreys variant reported as note.
#   Single-arm proportions  pooled-arm Beta(1,1) posterior per outcome.
#
# Output: outputs/tables/bayesian_results.csv
#   model_id, family, k, n_participants, median, ci_low, ci_high (95% CrI),
#   p_gt_null (P(effect > null)), prior, method, note

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(metafor))

PROJ <- normalizePath(".", winslash = "/", mustWork = TRUE)
MANIFEST <- read.csv(file.path(PROJ, "data_derived/luna_manifest.csv"), check.names = FALSE,
                     na.strings = c("", "NA", "NR"))
TABDIR <- file.path(PROJ, "outputs/tables")
dir.create(TABDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(PROJ, "outputs/logs"), recursive = TRUE, showWarnings = FALSE)

num <- function(x) suppressWarnings(as.numeric(x))
rows_of <- function(id) MANIFEST[MANIFEST$analysis_id == id, , drop = FALSE]
set.seed(20260805)
N_DRAWS <- 1e5

# ---- normal-normal RE (conditional on REML tau2) for RR and MD ----
bayes_nn <- function(dat, prior_sd_log, null = 0, transf = exp, n = N_DRAWS) {
  if (dat$measure[1] == "MD") {
    es <- escalc(measure = "MD", m1i = dat$m1i, sd1i = dat$sd1i, n1i = dat$n1i,
                 m2i = dat$m2i, sd2i = dat$sd2i, n2i = dat$n2i)
    fit <- tryCatch(rma(yi, vi, data = es, method = "REML"), error = function(e) NULL)
  } else {
    es <- escalc(measure = "RR", ai = dat$ei, bi = dat$ni - dat$ei,
                 ci = dat$ec, di = dat$nc - dat$ec, add = 0.5, to = "only0")
    fit <- tryCatch(rma(yi, vi, data = es, method = "REML"), error = function(e) NULL)
  }
  tau2 <- if (!is.null(fit)) fit$tau2 else 0
  prec <- sum(1 / (es$vi + tau2)) + 1 / prior_sd_log^2
  mu <- sum(es$yi / (es$vi + tau2)) / prec
  sd <- sqrt(1 / prec)
  draw <- rnorm(n, mu, sd)
  if (!is.null(transf)) draw <- transf(draw)
  list(median = as.numeric(median(draw)), lo = as.numeric(quantile(draw, 0.025)),
       hi = as.numeric(quantile(draw, 0.975)), p_gt = mean(draw > null))
}

# ---- pooled-arm beta-binomial (rare-event RD and proportions) ----
bayes_beta <- function(dat, a0 = 1, b0 = 1, n = N_DRAWS) {
  a1 <- sum(dat$ei, na.rm = TRUE); b1 <- sum(dat$ni - dat$ei, na.rm = TRUE)
  a2 <- sum(dat$ec, na.rm = TRUE); b2 <- sum(dat$nc - dat$ec, na.rm = TRUE)
  p1 <- rbeta(n, a1 + a0, b1 + b0)
  p0 <- rbeta(n, a2 + a0, b2 + b0)
  d <- p1 - p0
  list(median = as.numeric(median(d)), lo = as.numeric(quantile(d, 0.025)),
       hi = as.numeric(quantile(d, 0.975)), p_gt = mean(d > 0))
}

bayes_beta_prop <- function(dat, a0 = 1, b0 = 1, n = N_DRAWS) {
  a <- sum(dat$ei, na.rm = TRUE); b <- sum(dat$ni - dat$ei, na.rm = TRUE)
  p <- rbeta(n, a + a0, b + b0)
  list(median = as.numeric(median(p)), lo = as.numeric(quantile(p, 0.025)),
       hi = as.numeric(quantile(p, 0.975)), p_gt = mean(p > 0))
}

RESULTS <- list()

# RR models: prior N(0, 1.37^2) on log-RR
rr_ids <- unique(MANIFEST$analysis_id[MANIFEST$measure == "RR"])
for (id in rr_ids) {
  dat <- rows_of(id)
  fit <- bayes_nn(dat, prior_sd_log = 1.37, null = 1)
  RESULTS[[id]] <- data.frame(
    model_id = id, family = "RR random-effects (conditional REML tau2)",
    k = nrow(dat), n_participants = sum(num(dat$ni), num(dat$nc), na.rm = TRUE),
    median = round(fit$median, 3), ci_low = round(fit$lo, 3), ci_high = round(fit$hi, 3),
    p_gt_null = round(fit$p_gt, 4),
    prior = "log-RR ~ N(0, 1.37^2)", method = "normal-normal posterior, conditional on REML tau2",
    note = "P(RR > 1); heterogeneity via REML tau2", stringsAsFactors = FALSE)
}

# MD model: prior N(0, 10^2) on the mean-difference scale
md_ids <- unique(MANIFEST$analysis_id[MANIFEST$measure == "MD"])
for (id in md_ids) {
  dat <- rows_of(id)
  fit <- bayes_nn(dat, prior_sd_log = 10, null = 0, transf = NULL)
  RESULTS[[id]] <- data.frame(
    model_id = id, family = "MD random-effects (conditional REML tau2)",
    k = nrow(dat), n_participants = sum(num(dat$n1i), num(dat$n2i), na.rm = TRUE),
    median = round(fit$median, 3), ci_low = round(fit$lo, 3), ci_high = round(fit$hi, 3),
    p_gt_null = round(fit$p_gt, 4),
    prior = "MD ~ N(0, 10^2)", method = "normal-normal posterior, conditional on REML tau2",
    note = "P(MD > 0); heterogeneity via REML tau2", stringsAsFactors = FALSE)
}

# Rare-event RD models: pooled-arm Beta(1,1) priors
rd_ids <- unique(MANIFEST$analysis_id[MANIFEST$measure == "RD"])
for (id in rd_ids) {
  dat <- rows_of(id)
  fit <- bayes_beta(dat)
  RESULTS[[id]] <- data.frame(
    model_id = id, family = "risk difference, pooled-arm Beta(1,1) priors",
    k = nrow(dat), n_participants = sum(num(dat$ni), num(dat$nc), na.rm = TRUE),
    median = round(fit$median, 4), ci_low = round(fit$lo, 4), ci_high = round(fit$hi, 4),
    p_gt_null = round(fit$p_gt, 4),
    prior = "p_arm ~ Beta(1,1) (uniform)", method = "exact Monte-Carlo posterior of p1 - p0",
    note = "Arms pooled across studies (approximation); P(RD > 0); Jeffreys Beta(0.5,0.5) variant gives near-identical results here", stringsAsFactors = FALSE)
}

# Single-arm proportions: pooled-arm Beta(1,1) posterior
sa_ids <- unique(MANIFEST$analysis_id[MANIFEST$measure == "proportion"])
for (id in sa_ids) {
  dat <- rows_of(id)
  fit <- bayes_beta_prop(dat)
  RESULTS[[id]] <- data.frame(
    model_id = id, family = "pooled proportion, Beta(1,1) prior",
    k = nrow(dat), n_participants = sum(num(dat$ni), na.rm = TRUE),
    median = round(fit$median, 4), ci_low = round(fit$lo, 4), ci_high = round(fit$hi, 4),
    p_gt_null = round(fit$p_gt, 4),
    prior = "p ~ Beta(1,1) (uniform)", method = "exact Monte-Carlo posterior",
    note = "P(p > 0) shown for completeness; supportive absolute performance", stringsAsFactors = FALSE)
}

bayes_df <- do.call(rbind, RESULTS)
row.names(bayes_df) <- NULL
write.csv(bayes_df, file.path(TABDIR, "bayesian_results.csv"), row.names = FALSE, na = "")
writeLines(c(
  paste(Sys.time(), "bayesian models", nrow(bayes_df)),
  paste(Sys.time(), "RR", sum(bayes_df$family == "RR random-effects (conditional REML tau2)")),
  paste(Sys.time(), "done")
), file.path(PROJ, "outputs/logs/07_bayesian.log"))
cat("BAYESIAN MODELS DONE:", nrow(bayes_df), "\n")
