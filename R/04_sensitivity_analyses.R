# 04_sensitivity_analyses.R — prespecified sensitivity and subgroup analyses
# (protocol section 7; v4 plan section 7) reimplemented on the Luna manifest.
#
# Families implemented here (refit-based):
#   S1  Mantel-Haenszel fixed-effect RR for every RR model
#   S2  leave-one-study-out influence (k>=3 bases)
#   S3  exclude critical/high risk-of-bias studies (rob_joined.csv; Critical only,
#       plus a stricter Critical+Serious variant for the primary model)
#   S4  device-generation subgroups: fully biodegradable only / BioSTAR only
#   S5  exclude derived/reconstructed/flagged rows
#   S6  alternative continuity correction (targeted add=0.5) for rare-event RD
#       models + exact per-study Clopper-Pearson CIs for single-arm DRT
#   S8  exploratory design subgroups (randomized vs observational), labeled
#       exploratory, no meta-regression
# Families S7 (alternative follow-up selection), S9 (strict vs effective),
# S10 (adjusted vs crude) and S11 (overlap exclusion) are documented in the
# outcome-tables document; S7 models live in the model registry
# (ecl6_expanded_early, ecl12_gen_with_luani) and S9-S11 are restatements.
#
# Output: outputs/tables/sensitivity_results.csv

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(metafor))

PROJ <- normalizePath(".", winslash = "/", mustWork = TRUE)
dir.create(file.path(PROJ, "outputs/tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(PROJ, "outputs/logs"), recursive = TRUE, showWarnings = FALSE)

MANIFEST <- read.csv(file.path(PROJ, "data_derived/luna_manifest.csv"), check.names = FALSE,
                     na.strings = c("", "NA", "NR"))
rob <- read.csv(file.path(PROJ, "data_derived/rob_joined.csv"), check.names = FALSE)

num <- function(x) suppressWarnings(as.numeric(x))
rows_of <- function(id) MANIFEST[MANIFEST$analysis_id == id, , drop = FALSE]
is_derived <- function(dat) grepl("^DERIVED:", dat$row_id) | dat$status == "derived"
crit_studies <- rob$study_id[rob$overall_judgment %in% c("Critical", "High")]
serious_studies <- rob$study_id[rob$overall_judgment %in% c("Serious", "Critical", "High")]

SENS <- list()
record <- function(model_id, base, family, type, dat, fit, measure, method, note = "") {
  SENS[[model_id]] <<- data.frame(
    model_id = model_id, base_analysis = base, sensitivity_family = family,
    sensitivity_type = type,
    k = nrow(dat),
    k_used = if (!is.null(fit$k_used)) fit$k_used else if (!is.null(fit$fit)) fit$fit$k else NA_integer_,
    n_participants = sum(num(dat$ni), num(dat$nc), na.rm = TRUE),
    measure = measure, method = method,
    estimate = fit$estimate, ci_low = fit$ci_low, ci_high = fit$ci_high,
    tau2 = if (!is.null(fit$tau2)) fit$tau2 else NA_real_,
    i2 = if (!is.null(fit$i2)) fit$i2 else NA_real_,
    q = if (!is.null(fit$q)) fit$q else NA_real_,
    q_p = if (!is.null(fit$q_p)) fit$q_p else NA_real_,
    included = paste(dat$study_id, collapse = "; "),
    note = paste(unique(c(note, dat$note)), collapse = " "),
    stringsAsFactors = FALSE)
}

# ---------------------------------------------------------------- fit helpers
fit_rr_re <- function(dat) {
  ai <- dat$ei; bi <- dat$ni - dat$ei; ci <- dat$ec; di <- dat$nc - dat$ec
  es <- escalc(measure = "RR", ai = ai, bi = bi, ci = ci, di = di, add = 0.5, to = "only0")
  fit <- rma(yi, vi, data = es, method = "REML", test = "knha")
  pr <- predict(fit, transf = exp)
  list(fit = fit, estimate = as.numeric(pr$pred), ci_low = as.numeric(pr$ci.lb),
       ci_high = as.numeric(pr$ci.ub), tau2 = as.numeric(fit$tau2), i2 = as.numeric(fit$I2),
       q = as.numeric(fit$QE), q_p = as.numeric(fit$QEp), method = "RR REML+HK")
}

fit_rr_fixed <- function(dat) {
  ai <- dat$ei; bi <- dat$ni - dat$ei; ci <- dat$ec; di <- dat$nc - dat$ec
  fit <- tryCatch(suppressWarnings(rma.mh(ai = ai, bi = bi, ci = ci, di = di, measure = "RR")),
                  error = function(e) e)
  if (inherits(fit, "error")) {
    return(list(estimate = NA_real_, ci_low = NA_real_, ci_high = NA_real_, k_used = 0L,
                tau2 = NA_real_, i2 = NA_real_, q = NA_real_, q_p = NA_real_,
                method = "MH fixed (not estimable)", note = conditionMessage(fit)))
  }
  pr <- predict(fit, transf = exp)
  list(estimate = as.numeric(pr$pred), ci_low = as.numeric(pr$ci.lb), ci_high = as.numeric(pr$ci.ub),
       k_used = fit$k, tau2 = NA_real_, i2 = NA_real_, q = as.numeric(fit$QE), q_p = as.numeric(fit$QEp),
       method = "Mantel-Haenszel fixed effect (RR)", note = "")
}

fit_rd_mh <- function(dat, add) {
  ai <- dat$ei; bi <- dat$ni - dat$ei; ci <- dat$ec; di <- dat$nc - dat$ec
  fit <- suppressWarnings(rma.mh(ai = ai, bi = bi, ci = ci, di = di, measure = "RD",
                                 add = add, to = "only0", correct = FALSE))
  pr <- predict(fit)
  list(estimate = as.numeric(pr$pred), ci_low = as.numeric(pr$ci.lb), ci_high = as.numeric(pr$ci.ub),
       k_used = fit$k, tau2 = NA_real_, i2 = NA_real_, q = as.numeric(fit$QE), q_p = as.numeric(fit$QEp),
       method = paste0("MH RD; targeted continuity correction (add=", add, ")"))
}

desc_fit <- function(dat, method) {
  list(estimate = NA_real_, ci_low = NA_real_, ci_high = NA_real_, k_used = nrow(dat),
       tau2 = NA_real_, i2 = NA_real_, q = NA_real_, q_p = NA_real_, method = method)
}

# ------------------------------------------------- S1. MH fixed-effect (all RR) ----
rr_ids <- unique(MANIFEST$analysis_id[MANIFEST$measure == "RR"])
for (id in rr_ids) {
  dat <- rows_of(id)
  fit <- fit_rr_fixed(dat)
  record(paste0(id, "__mh_fixed"), id, "S1 MH fixed effect",
         "Mantel-Haenszel fixed effect (RR)", dat, fit, "RR", fit$method,
         if (fit$k_used < nrow(dat)) "MH excludes double-zero rows; k_used < k." else "")
}

# ------------------------------------------------- S2. leave-one-study-out ----
loo_bases <- c("ecl6_core", "strict_early", "cve_late_core", "cve_late_with_luani",
               "af_late_full", "af_late_full_with_luani", "af_any_late",
               "af_any_late_with_darsaklis", "ecl6_expanded_early",
               "drt_any", "embolization", "pericardial")
for (base in loo_bases) {
  all <- rows_of(base)
  if (nrow(all) < 3) next
  meas <- all$measure[1]
  for (sid in unique(all$study_id)) {
    dat <- all[all$study_id != sid, , drop = FALSE]
    fit <- if (nrow(dat) < 2) {
      desc_fit(dat, paste0(meas, " (k=1 after exclusion)"))
    } else if (meas == "RR") {
      fit_rr_re(dat)
    } else {
      fit_rd_mh(dat, add = 0)
    }
    record(paste0(base, "__loo_", gsub("[^A-Za-z0-9]", "_", sid)), base, "S2 leave-one-out",
           paste("leave-one-study-out (excluded:", sid, ")"), dat, fit, meas, fit$method,
           if (nrow(dat) < 2) "single study remains; descriptive only, not pooled." else "")
  }
}

# ------------------------------------------------ S3. risk-of-bias exclusion ----
rob_bases <- c("ecl6_core", "ecl24", "cve_late_core", "cve_late_with_luani",
               "af_late_full", "af_late_full_with_luani", "af_any_late",
               "af_any_late_with_darsaklis", "ecl12_gen_with_luani", "af_late_with_luani")
# S3a: exclude critical/high risk-of-bias studies (permissive ROBINS-I tier 2
# consensus: no comparative study is critical, so these are documented no-ops)
for (base in rob_bases) {
  all <- rows_of(base)
  kept <- all[!(all$study_id %in% crit_studies), , drop = FALSE]
  if (nrow(kept) == nrow(all)) {
    record(paste0(base, "__no_critical_rob"), base, "S3 risk-of-bias exclusion",
           "exclude critical/high risk-of-bias studies", kept, desc_fit(kept, "no-op"),
           kept$measure[1], "no-op",
           "no critical/high-RoB studies under the permissive ROBINS-I consensus (tier 2); result identical to base model.")
    next
  }
  fit <- if (nrow(kept) < 2) {
    desc_fit(kept, paste0(kept$measure[1], " (k=1 after exclusion)"))
  } else if (kept$measure[1] == "RR") {
    fit_rr_re(kept)
  } else {
    fit_rd_mh(kept, add = 0)
  }
  record(paste0(base, "__no_critical_rob"), base, "S3 risk-of-bias exclusion",
         paste("exclude critical/high risk-of-bias studies (excluded:",
               paste(setdiff(all$study_id, kept$study_id), collapse = ", "), ")"),
         kept, fit, kept$measure[1], fit$method,
         if (nrow(kept) < 2) "single study remains; descriptive only, not pooled." else "")
}
# S3b: exclude critical+serious risk-of-bias studies (informative under the
# permissive consensus: Luani 2015 and Cernic Suligoj 2020 are serious)
for (base in rob_bases) {
  all <- rows_of(base)
  kept <- all[!(all$study_id %in% serious_studies), , drop = FALSE]
  if (nrow(kept) == nrow(all)) {
    record(paste0(base, "__no_serious_rob"), base, "S3 risk-of-bias exclusion",
           "exclude critical+serious risk-of-bias studies", kept, desc_fit(kept, "no-op"),
           kept$measure[1], "no-op",
           "no serious/critical-RoB studies in this base; result identical to base model.")
    next
  }
  fit <- if (nrow(kept) < 2) {
    desc_fit(kept, paste0(kept$measure[1], " (k=1 after exclusion)"))
  } else if (kept$measure[1] == "RR") {
    fit_rr_re(kept)
  } else {
    fit_rd_mh(kept, add = 0)
  }
  record(paste0(base, "__no_serious_rob"), base, "S3 risk-of-bias exclusion",
         paste("exclude critical+serious risk-of-bias studies (excluded:",
               paste(setdiff(all$study_id, kept$study_id), collapse = ", "), ")"),
         kept, fit, kept$measure[1], fit$method,
         if (nrow(kept) < 2) "single study remains; descriptive only, not pooled." else "")
}

# ------------------------------------------ S4. device-generation subgroups ----
base <- "ecl6_core"
all <- rows_of(base)
fb <- all[all$study_id != "Cernic Suligoj 2020", , drop = FALSE]
record(paste0(base, "__fully_biodegradable_only"), base, "S4 device generation",
       "device generation: fully biodegradable only (exclude partially bioabsorbable BioSTAR)",
       fb, fit_rr_re(fb), "RR", "RR REML+HK")
bs <- all[all$study_id == "Cernic Suligoj 2020", , drop = FALSE]
record(paste0(base, "__biostar_only"), base, "S4 device generation",
       "device generation: BioSTAR (partially bioabsorbable) only, k=1 descriptive",
       bs, desc_fit(bs, "RR (k=1)"), "RR", "RR (k=1 descriptive)",
       "single study remains; descriptive only, not pooled.")
base <- "ecl24"
all <- rows_of(base)
zh <- all[all$study_id == "Zhang 2026", , drop = FALSE]
record(paste0(base, "__zhang_only_fullybio"), base, "S4 device generation",
       "device generation: long-term durability, fully biodegradable only (Zhang 2026, k=1 descriptive)",
       zh, desc_fit(zh, "RR (k=1)"), "RR", "RR (k=1 descriptive)",
       "single study remains; descriptive only, not pooled.")

# ------------------------------------ S5. exclude derived/reconstructed rows ----
s5_specs <- list(
  list(base = "ecl12_mig", keep = function(dat) dat[dat$study_id != "Li 2025", , drop = FALSE],
       type = "exclude inferred Li 2025 rows (counts derived from printed percentages)"),
  list(base = "close12_mig_strict", keep = function(dat) dat[dat$study_id != "Li 2025", , drop = FALSE],
       type = "exclude inferred Li 2025 rows (counts derived from printed percentages)"),
  list(base = "shunt_g6_general", keep = function(dat) dat[!is_derived(dat), , drop = FALSE],
       type = "exclude derived rows (Li 2026 grade 2-3 combined)"),
  list(base = "af_any_late", keep = function(dat) dat[!is_derived(dat), , drop = FALSE],
       type = "exclude derived rows (Li 2026 AF+SVT summed row)")
)
for (sp in s5_specs) {
  all <- rows_of(sp$base)
  kept <- sp$keep(all)
  fit <- if (nrow(kept) < 2) {
    desc_fit(kept, paste0(kept$measure[1], " (k=1 after exclusion)"))
  } else {
    fit_rr_re(kept)
  }
  record(paste0(sp$base, "__direct_only"), sp$base, "S5 exclude derived/conditional rows",
         sp$type, kept, fit, kept$measure[1], fit$method,
         if (nrow(kept) < 2) "single study remains; descriptive only, not pooled." else "")
}
# restatements: bases that are already all-direct
for (base in c("ecl6_core", "ecl12_gen_core")) {
  all <- rows_of(base)
  record(paste0(base, "__direct_only"), base, "S5 exclude derived/conditional rows",
         "all rows direct; base model already excludes derived/reconstructed rows (restated)",
         all, desc_fit(all, "no-op"), all$measure[1], "no-op",
         "no derived/conditional rows in base; result identical to base model.")
}

# ----------------------- S6. alternative continuity correction / exact CIs ----
for (base in c("drt_any", "embolization", "pericardial")) {
  dat <- rows_of(base)
  fit <- fit_rd_mh(dat, add = 0.5)
  record(paste0(base, "__cc05"), base, "S6 alternative continuity correction",
         "MH RD with targeted 0.5 continuity correction (v4 approach); primary uses no universal correction",
         dat, fit, "RD", fit$method)
}
record("sa_drt__glmm_note", "sa_drt", "S6 alternative continuity correction",
       "binomial GLMM (rma.glmm) unavailable: documented native heap-corruption crash on Windows R 4.6.1 (see tests/diagnose_fits.R); exact per-study Clopper-Pearson CIs reported instead",
       rows_of("sa_drt"), desc_fit(rows_of("sa_drt"), "exact (per study)"), "proportion", "exact (per study)")
for (sid in c("He 2025", "Du 2022")) {
  dat <- rows_of("sa_drt")
  dat <- dat[dat$study_id == sid, , drop = FALSE]
  bt <- binom.test(dat$ei, dat$ni)
  record(paste0("sa_drt__exact_", gsub("[^A-Za-z0-9]", "_", sid)), "sa_drt",
         "S6 alternative continuity correction",
         paste("per-study exact Clopper-Pearson CI:", sid), dat,
         list(estimate = as.numeric(bt$estimate), ci_low = as.numeric(bt$conf.int[1]),
              ci_high = as.numeric(bt$conf.int[2]), k_used = 1L,
              tau2 = NA_real_, i2 = NA_real_, q = NA_real_, q_p = NA_real_,
              method = "exact Clopper-Pearson"),
         "proportion", "exact Clopper-Pearson")
}

# -------------------------- S8. exploratory design subgroups (labeled) ----
base <- "ecl6_core"
all <- rows_of(base)
rct <- all[all$study_id == "Zhang 2026", , drop = FALSE]
record(paste0(base, "__subgrp_rct"), base, "S8 exploratory subgroups",
       "exploratory design subgroup: randomized (Zhang 2026), k=1 descriptive; labeled exploratory, no meta-regression",
       rct, desc_fit(rct, "RR (k=1)"), "RR", "RR (k=1 descriptive)",
       "single study remains; descriptive only, not pooled.")
obs <- all[all$study_id != "Zhang 2026", , drop = FALSE]
record(paste0(base, "__subgrp_obs"), base, "S8 exploratory subgroups",
       "exploratory design subgroup: observational (Li 2026, Cernic Suligoj 2020); labeled exploratory, no meta-regression",
       obs, fit_rr_re(obs), "RR", "RR REML+HK")

# ---------------------------------------------------------------- output ----
sens_df <- do.call(rbind, SENS)
row.names(sens_df) <- NULL
write.csv(sens_df, file.path(PROJ, "outputs/tables/sensitivity_results.csv"), row.names = FALSE, na = "")
writeLines(c(
  paste(Sys.time(), "sensitivity rows", nrow(sens_df)),
  paste(Sys.time(), "S1 MH fixed", sum(grepl("S1", sens_df$sensitivity_family))),
  paste(Sys.time(), "S2 leave-one-out", sum(grepl("S2", sens_df$sensitivity_family))),
  paste(Sys.time(), "done")
), file.path(PROJ, "outputs/logs/04_sensitivity_analyses.log"))
cat("SENSITIVITY MODELS DONE:", nrow(sens_df), "\n")
