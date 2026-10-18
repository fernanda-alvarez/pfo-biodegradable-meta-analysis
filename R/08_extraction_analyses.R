# 08_extraction_analyses.R — analyses enabled by targeted re-extraction from
# source PDFs and previously non-pooled extracted values.
#
# Part 1 — Migraine relief (outcome family 9).
#   Sun 2025 contributes arm-level 12-month relief counts (39/56 vs 54/77).
#   Li 2025 reports Kaplan-Meier curves only (figure raster, no axis values in
#   the text layer) with log-rank P = 0.53; reliable digitization requires OCR
#   and imaging libraries that are not available in this environment, so the
#   plan-specified fallback applies: quantified narrative row, no pooled model.
#   The KM figure is exported to outputs/figures/li2025_km_figure.png.
#
# Part 2 — Person-time incidence-rate ratios for recurrent CVE.
#   Assumption-based person-time (documented per study):
#     Karagianni 2011  mean FU midpoint of reported 18-36 months = 27 months
#     Darsaklis 2014   overall mean follow-up 20.8 months (per paper)
#     Cernic Suligoj 2020  long-term TOE "beyond 2 years" -> 2.5 y midpoint
#                         (sensitivity at 2.0 and 3.0 years)
#     Luani 2015       per-group median follow-up (BioSTAR 3.3 y); total 1,815
#                      patient-years scaled to the n x median allocation
#   IRR models mirror the RR sets (core and with-Luani).
#
# Part 3 — MIDAS change-from-baseline mean difference.
#   Baseline and 12-month means/SDs already extracted (v4 tidy C8:46/49/72/73,
#   verified against PDFs). Within-arm change SD derived assuming correlation
#   r = 0.5 (sensitivity 0.25 / 0.75); pooled MD of arm-level changes (k=2).
#
# Outputs: outputs/tables/migraine_relief_results.csv, irr_results.csv,
#          midas_change_results.csv

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(metafor))

PROJ <- normalizePath(".", winslash = "/", mustWork = TRUE)
TABDIR <- file.path(PROJ, "outputs/tables")
dir.create(TABDIR, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------- Part 1: migraine relief
migraine_relief <- data.frame(
  study_id = c("Sun 2025", "Li 2025"),
  definition = c(">=50% reduction in monthly migraine attacks or complete elimination (12 months)",
                 ">=50% reduction in MIDAS score (12 months); Kaplan-Meier curves only"),
  events_bio = c(39, NA), n_bio = c(56, NA),
  events_comp = c(54, NA), n_comp = c(77, NA),
  pct_bio = c("69.6%", "not reported"),
  pct_comp = c("70.1%", "not reported"),
  quantitative = c("Arm-level counts; single-study descriptive (RR 0.99, 95% CI 0.80-1.23)",
                   "No arm-level counts; log-rank P = 0.53 (no significant between-group difference)"),
  stringsAsFactors = FALSE)
write.csv(migraine_relief, file.path(TABDIR, "migraine_relief_results.csv"), row.names = FALSE, na = "")

# ---------------------------------------------------- Part 2: person-time IRRs
# per-arm person-time (person-years) and events; mirrors the CVE RR row sets
pt_data <- data.frame(
  study_id = c("Karagianni 2011", "Darsaklis 2014", "Cernic Suligoj 2020", "Luani 2015"),
  fu_bio_y = c(2.25, 20.8 / 12, 2.5, 3.3),
  fu_comp_y = c(2.25, 20.8 / 12, 2.5, NA),
  n_bio = c(30, 32, 44, 68),
  n_comp = c(29, 63, 40, 267),
  ev_bio = c(4, 1, 3, 1),
  ev_comp = c(2, 1, 1, 14),
  fu_note = c("mean FU midpoint of reported 18-36 months (assumption)",
              "overall mean follow-up 20.8 months per paper, applied to both arms (assumption; group FU unbalanced)",
              "long-term TOE 'beyond 2 years'; midpoint 2.5 y assumed (sensitivity 2.0/3.0)",
              "per-group median follow-up (BioSTAR 3.3 y); comparator person-time = 1,815 total patient-years minus BioSTAR allocation (n x median scaled to the stated total)"),
  stringsAsFactors = FALSE)
# Luani comparator person-time: allocate total 1815 PY proportionally to n x median
luani_alloc <- 68 * 3.3
tot_alloc <- luani_alloc + 109 * 5 + 104 * 9 + 54 * 4.6
scale <- 1815 / tot_alloc
pt_data$pt_bio <- pt_data$n_bio * pt_data$fu_bio_y
pt_data$pt_comp <- ifelse(is.na(pt_data$fu_comp_y),
                          (1815 - pt_data$n_bio * pt_data$fu_bio_y * scale),
                          pt_data$n_comp * pt_data$fu_comp_y)
pt_data$pt_bio[4] <- pt_data$n_bio[4] * pt_data$fu_bio_y[4] * scale

fit_irr <- function(dat) {
  es <- escalc(measure = "IRR", x1i = dat$ev_bio, t1i = dat$pt_bio,
               x2i = dat$ev_comp, t2i = dat$pt_comp)
  fit <- rma(yi, vi, data = es, method = "REML", test = "knha")
  pr <- predict(fit, transf = exp)
  list(fit = fit, est = as.numeric(pr$pred), lo = as.numeric(pr$ci.lb), hi = as.numeric(pr$ci.ub),
       tau2 = as.numeric(fit$tau2), i2 = as.numeric(fit$I2), q = as.numeric(fit$QE), q_p = as.numeric(fit$QEp))
}
per_study_irr <- function(dat) {
  sapply(seq_len(nrow(dat)), function(i) {
    irr <- (dat$ev_bio[i] / dat$pt_bio[i]) / (dat$ev_comp[i] / dat$pt_comp[i])
    se <- sqrt(1 / dat$ev_bio[i] + 1 / dat$ev_comp[i])
    c(irr = irr, lo = exp(log(irr) - 1.96 * se), hi = exp(log(irr) + 1.96 * se))
  })
}

irr_sets <- list(
  core = pt_data[1:3, , drop = FALSE],
  with_luani = pt_data
)
IRR_OUT <- list()
for (nm in names(irr_sets)) {
  dat <- irr_sets[[nm]]
  psi <- per_study_irr(dat)
  pool <- fit_irr(dat)
  IRR_OUT[[nm]] <- data.frame(
    model_id = paste0("cve_irr_", nm), k = nrow(dat),
    irr = round(pool$est, 3), ci_low = round(pool$lo, 3), ci_high = round(pool$hi, 3),
    tau2 = round(pool$tau2, 5), i2 = round(pool$i2, 1), q = round(pool$q, 3), q_p = round(pool$q_p, 4),
    per_study = paste(sprintf("%s IRR %.2f (%.2f-%.2f)", dat$study_id, psi[1, ], psi[2, ], psi[3, ]), collapse = "; "),
    fu_assumptions = paste(dat$fu_note, collapse = " | "),
    stringsAsFactors = FALSE)
}
# sensitivity to the Cernic FU midpoint (2.0 and 3.0 years)
for (fu in c(2.0, 3.0)) {
  dat <- pt_data[1:3, , drop = FALSE]
  dat$fu_bio_y[3] <- fu; dat$fu_comp_y[3] <- fu
  dat$pt_bio <- dat$n_bio * dat$fu_bio_y
  dat$pt_comp <- dat$n_comp * dat$fu_comp_y
  pool <- fit_irr(dat)
  IRR_OUT[[paste0("cve_irr_core_cernic_", fu, "y")]] <- data.frame(
    model_id = paste0("cve_irr_core_cernic_", fu, "y"), k = 3,
    irr = round(pool$est, 3), ci_low = round(pool$lo, 3), ci_high = round(pool$hi, 3),
    tau2 = round(pool$tau2, 5), i2 = round(pool$i2, 1), q = round(pool$q, 3), q_p = round(pool$q_p, 4),
    per_study = paste(sprintf("%s IRR %.2f", dat$study_id, per_study_irr(dat)[1, ]), collapse = "; "),
    fu_assumptions = paste0("Cernic follow-up set to ", fu, " years (sensitivity)"),
    stringsAsFactors = FALSE)
}
irr_df <- do.call(rbind, IRR_OUT)
write.csv(irr_df, file.path(TABDIR, "irr_results.csv"), row.names = FALSE, na = "")

# ------------------------------------------- Part 3: MIDAS change-from-baseline
midas_dat <- data.frame(
  study_id = c("Sun 2025", "Li 2025"),
  m0_bio = c(37.80, 41.54), sd0_bio = c(7.30, 17.18), n0_bio = c(56, 81),
  m12_bio = c(10.45, 10.07), sd12_bio = c(9.19, 8.81), n12_bio = c(56, 71),
  m0_comp = c(40.05, 38.09), sd0_comp = c(7.75, 20.26), n0_comp = c(77, 77),
  m12_comp = c(11.32, 10.74), sd12_comp = c(9.62, 8.81), n12_comp = c(77, 68),
  stringsAsFactors = FALSE)

change_sd <- function(sd0, sd12, r) sqrt(sd0^2 + sd12^2 - 2 * r * sd0 * sd12)

fit_midas_change <- function(r) {
  change <- data.frame(
    study_id = midas_dat$study_id,
    change_bio = midas_dat$m12_bio - midas_dat$m0_bio,
    sd_change_bio = change_sd(midas_dat$sd0_bio, midas_dat$sd12_bio, r),
    n_bio = midas_dat$n12_bio,
    change_comp = midas_dat$m12_comp - midas_dat$m0_comp,
    sd_change_comp = change_sd(midas_dat$sd0_comp, midas_dat$sd12_comp, r),
    n_comp = midas_dat$n12_comp,
    stringsAsFactors = FALSE)
  es <- escalc(measure = "MD", m1i = change$change_bio, sd1i = change$sd_change_bio, n1i = change$n_bio,
               m2i = change$change_comp, sd2i = change$sd_change_comp, n2i = change$n_comp)
  fit <- tryCatch(rma(yi, vi, data = es, method = "REML", test = "knha"), error = function(e) NULL)
  pr <- if (!is.null(fit)) predict(fit) else NULL
  data.frame(
    model_id = paste0("midas_change_r", r), correlation = r, k = nrow(change),
    md = if (!is.null(pr)) round(pr$pred, 3) else NA,
    ci_low = if (!is.null(pr)) round(pr$ci.lb, 3) else NA,
    ci_high = if (!is.null(pr)) round(pr$ci.ub, 3) else NA,
    per_study = paste(sprintf("%s MD %.2f", change$study_id, es$yi), collapse = "; "),
    note = "Within-arm change = 12-month minus baseline MIDAS; change SD derived with assumed correlation r (sensitivity 0.25/0.75); pooled MD of changes, k=2",
    stringsAsFactors = FALSE)
}
midas_change_df <- do.call(rbind, lapply(c(0.25, 0.5, 0.75), fit_midas_change))
write.csv(midas_change_df, file.path(TABDIR, "midas_change_results.csv"), row.names = FALSE, na = "")

writeLines(c(
  paste(Sys.time(), "migraine relief rows", nrow(migraine_relief)),
  paste(Sys.time(), "IRR models", nrow(irr_df)),
  paste(Sys.time(), "MIDAS change models", nrow(midas_change_df)),
  paste(Sys.time(), "done")
), file.path(PROJ, "outputs/logs/08_extraction_analyses.log"))
cat("EXTRACTION ANALYSES DONE:", nrow(migraine_relief), "relief;", nrow(irr_df),
    "IRR;", nrow(midas_change_df), "MIDAS-change\n")
