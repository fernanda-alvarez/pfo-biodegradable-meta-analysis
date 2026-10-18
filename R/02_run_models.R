# Run the corrected Luna model registry from the explicit manifest.

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(metafor))

PROJ <- normalizePath(".", winslash = "/", mustWork = TRUE)
MANIFEST <- read.csv(file.path(PROJ, "data_derived/luna_manifest.csv"), check.names = FALSE,
                     na.strings = c("", "NA", "NR"))
dir.create(file.path(PROJ, "outputs/tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(PROJ, "outputs/logs"), recursive = TRUE, showWarnings = FALSE)

num <- function(x) suppressWarnings(as.numeric(x))
logfile <- file.path(PROJ, "outputs/logs/02_run_models.log")
writeLines(paste(Sys.time(), "start"), logfile)

empty_result <- function(model_id, spec, note, method) {
  data.frame(model_id = model_id, outcome_family = spec$outcome_family,
             window = spec$window, population_stratum = spec$population_stratum,
             role = spec$role, measure = spec$measure, k = 0, n_participants = NA,
             estimate = NA, ci_low = NA, ci_high = NA, pi_low = NA, pi_high = NA,
             tau2 = NA, i2 = NA, q = NA, q_p = NA, included = "", method = method,
             note = note, stringsAsFactors = FALSE)
}

fit_comparative <- function(dat, measure) {
  if (measure == "MD") {
    dat <- dat[!is.na(dat$m1i) & !is.na(dat$sd1i) & !is.na(dat$n1i) &
                 !is.na(dat$m2i) & !is.na(dat$sd2i) & !is.na(dat$n2i), , drop = FALSE]
    if (!nrow(dat)) return(list(fit = NULL, estimate = NA, ci_low = NA, ci_high = NA, method = "MD REML+HK"))
    es <- escalc(measure = "MD", m1i = dat$m1i, sd1i = dat$sd1i, n1i = dat$n1i,
                 m2i = dat$m2i, sd2i = dat$sd2i, n2i = dat$n2i)
    fit <- rma(yi, vi, data = es, method = "REML", test = "knha")
    pr <- predict(fit)
    return(list(fit = fit, estimate = as.numeric(pr$pred), ci_low = as.numeric(pr$ci.lb),
                ci_high = as.numeric(pr$ci.ub), pi_low = as.numeric(pr$pi.lb),
                pi_high = as.numeric(pr$pi.ub), method = "MD REML+HK"))
  }

  dat <- dat[!is.na(dat$ei) & !is.na(dat$ni) & !is.na(dat$ec) & !is.na(dat$nc), , drop = FALSE]
  if (!nrow(dat)) return(list(fit = NULL, estimate = NA, ci_low = NA, ci_high = NA, method = measure))
  ai <- dat$ei; bi <- dat$ni - dat$ei; ci <- dat$ec; di <- dat$nc - dat$ec

  if (measure == "RD") {
    # No universal continuity correction. MH RD retains zero-event rows as far as
    # the estimator permits; this is preferable to the v4 targeted 0.5 correction.
    fit <- suppressWarnings(rma.mh(ai = ai, bi = bi, ci = ci, di = di, measure = "RD",
                                   add = 0, to = "none", correct = FALSE))
    pr <- predict(fit)
    return(list(fit = fit, estimate = as.numeric(pr$pred), ci_low = as.numeric(pr$ci.lb),
                ci_high = as.numeric(pr$ci.ub), pi_low = NA_real_, pi_high = NA_real_,
                method = "MH RD; no continuity correction"))
  }

  es <- escalc(measure = measure, ai = ai, bi = bi, ci = ci, di = di,
               add = 0.5, to = "only0")
  fit <- rma(yi, vi, data = es, method = "REML", test = "knha")
  pr <- predict(fit, transf = exp)
  list(fit = fit, estimate = as.numeric(pr$pred), ci_low = as.numeric(pr$ci.lb),
       ci_high = as.numeric(pr$ci.ub), pi_low = as.numeric(pr$pi.lb),
       pi_high = as.numeric(pr$pi.ub), method = paste(measure, "REML+HK"))
}

fit_proportion <- function(dat) {
  dat <- dat[!is.na(dat$ei) & !is.na(dat$ni), , drop = FALSE]
  if (!nrow(dat)) return(list(fit = NULL, estimate = NA, ci_low = NA, ci_high = NA, method = "proportion"))
  # Use the logit-normal inverse-variance implementation rather than
  # metafor::rma.glmm(), whose Windows R 4.6.1 cleanup path can terminate the
  # process with native heap corruption after returning a valid estimate.
  es <- escalc(measure = "PLO", xi = dat$ei, ni = dat$ni,
               add = 0.5, to = "only0")
  fit <- rma(yi, vi, data = es, method = "REML", test = "knha")
  pr <- predict(fit, transf = plogis)
  list(fit = fit, estimate = as.numeric(pr$pred), ci_low = as.numeric(pr$ci.lb),
       ci_high = as.numeric(pr$ci.ub), pi_low = as.numeric(pr$pi.lb),
       pi_high = as.numeric(pr$pi.ub), method = "logit-normal REML+HK")
}

run_one <- function(model_id, dat) {
  spec <- list(outcome_family = dat$outcome_family[1], window = dat$window[1],
               population_stratum = dat$population_stratum[1], role = dat$role[1],
               measure = dat$measure[1])
  if (spec$measure == "proportion") fit <- fit_proportion(dat) else fit <- fit_comparative(dat, spec$measure)
  n_total <- if (spec$measure == "proportion") {
    sum(num(dat$ni), na.rm = TRUE)
  } else if (spec$measure == "MD") {
    sum(num(dat$n1i), num(dat$n2i), na.rm = TRUE)
  } else {
    sum(num(dat$ni), num(dat$nc), na.rm = TRUE)
  }
  f <- fit$fit
  data.frame(
    model_id = model_id, outcome_family = spec$outcome_family, window = spec$window,
    population_stratum = spec$population_stratum, role = spec$role, measure = spec$measure,
    k = nrow(dat), n_participants = n_total,
    estimate = fit$estimate, ci_low = fit$ci_low, ci_high = fit$ci_high,
    pi_low = if (is.null(fit$pi_low)) NA_real_ else fit$pi_low,
    pi_high = if (is.null(fit$pi_high)) NA_real_ else fit$pi_high,
    tau2 = if (!is.null(f) && !is.null(f$tau2)) as.numeric(f$tau2) else NA_real_,
    i2 = if (!is.null(f) && !is.null(f$I2)) as.numeric(f$I2) else NA_real_,
    q = if (!is.null(f) && !is.null(f$QE)) as.numeric(f$QE) else NA_real_,
    q_p = if (!is.null(f) && !is.null(f$QEp)) as.numeric(f$QEp) else NA_real_,
    included = paste(dat$study_id, collapse = "; "), method = fit$method,
    note = paste(unique(dat$note), collapse = " "), stringsAsFactors = FALSE
  )
}

model_ids <- unique(MANIFEST$analysis_id)
results <- do.call(rbind, lapply(model_ids, function(id) run_one(id, MANIFEST[MANIFEST$analysis_id == id, , drop = FALSE])))
registry <- results

write.csv(results[results$measure != "proportion", ], file.path(PROJ, "outputs/tables/comparative_results.csv"), row.names = FALSE, na = "")
write.csv(results[results$measure == "proportion", ], file.path(PROJ, "outputs/tables/single_arm_results.csv"), row.names = FALSE, na = "")
write.csv(registry, file.path(PROJ, "outputs/tables/model_registry.csv"), row.names = FALSE, na = "")
write.csv(MANIFEST, file.path(PROJ, "outputs/tables/model_input_manifest.csv"), row.names = FALSE, na = "")

writeLines(c(
  paste(Sys.time(), "models", nrow(results)),
  paste(Sys.time(), "primary", results$estimate[results$model_id == "ecl6_core"]),
  paste(Sys.time(), "rare methods", paste(unique(results$method[results$model_id %in% c("drt_any", "embolization", "pericardial")]), collapse = "; ")),
  paste(Sys.time(), "done")
), logfile)
cat("MODELS WRITTEN", nrow(results), "\n")
