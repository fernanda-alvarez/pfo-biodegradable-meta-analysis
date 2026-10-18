# 06_adjunct_analyses.R — small-study effects, E-values, cumulative meta-analysis,
# influence diagnostics, and subgroup interaction tests.
#
# Outputs (all under outputs/):
#   tables/small_study_results.csv   Egger + Harbord tests and trim-and-fill (RR models, k>=3)
#   tables/evalue_results.csv        E-values for unmeasured confounding (all comparative RR models)
#   tables/influence_results.csv     influence() diagnostics per study (RE models, k>=3)
#   tables/interaction_results.csv   Q-between interaction tests (population and design strata)
#   figures/funnel_<model>.png       funnel plots (k>=3 RR models)
#   figures/cumulative_ecl6_core.png cumulative meta-analysis by year (primary)
#   figures/baujat_<model>.png       baujat plots (k>=3 RE models)
#
# Runs with base R + metafor; no new packages.

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(metafor))

PROJ <- normalizePath(".", winslash = "/", mustWork = TRUE)
MANIFEST <- read.csv(file.path(PROJ, "data_derived/luna_manifest.csv"), check.names = FALSE,
                     na.strings = c("", "NA", "NR"))
TABDIR <- file.path(PROJ, "outputs/tables")
FIGDIR <- file.path(PROJ, "outputs/figures")
dir.create(TABDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGDIR, recursive = TRUE, showWarnings = FALSE)

num <- function(x) suppressWarnings(as.numeric(x))
rows_of <- function(id) MANIFEST[MANIFEST$analysis_id == id, , drop = FALSE]

fit_rr <- function(dat) {
  ai <- dat$ei; bi <- dat$ni - dat$ei; ci <- dat$ec; di <- dat$nc - dat$ec
  es <- escalc(measure = "RR", ai = ai, bi = bi, ci = ci, di = di, add = 0.5, to = "only0")
  rma(yi, vi, data = es, method = "REML", test = "knha")
}

rr_ids <- unique(MANIFEST$analysis_id[MANIFEST$measure == "RR"])
rr_ids_k3 <- rr_ids[vapply(rr_ids, function(id) nrow(rows_of(id)) >= 3, logical(1))]

# ------------------------------------------------------------- 1. small-study ----
SS <- list()
for (id in rr_ids_k3) {
  dat <- rows_of(id)
  fit <- fit_rr(dat)
  egger <- tryCatch(regtest(fit, model = "lm", predictor = "sei"), error = function(e) NULL)
  harbord <- tryCatch(regtest(fit, model = "lm", predictor = "ni"), error = function(e) NULL)
  tf <- tryCatch(trimfill(fit), error = function(e) NULL)
  SS[[id]] <- data.frame(
    model_id = id, k = nrow(dat), n_participants = sum(num(dat$ni), num(dat$nc), na.rm = TRUE),
    egger_z = if (!is.null(egger)) round(egger$zval, 3) else NA,
    egger_p = if (!is.null(egger)) round(egger$pval, 4) else NA,
    harbord_z = if (!is.null(harbord)) round(harbord$zval, 3) else NA,
    harbord_p = if (!is.null(harbord)) round(harbord$pval, 4) else NA,
    tf_estimate = if (!is.null(tf)) round(exp(tf$b[1, 1]), 3) else NA,
    tf_ci_low = if (!is.null(tf)) round(exp(tf$ci.lb[1]), 3) else NA,
    tf_ci_high = if (!is.null(tf)) round(exp(tf$ci.ub[1]), 3) else NA,
    tf_k_added = if (!is.null(tf)) tf$k0 else NA,
    note = if (nrow(dat) < 5) "low power at k<5; tests are exploratory only" else "k>=5; interpret cautiously",
    stringsAsFactors = FALSE)
  # funnel plot
  png(file.path(FIGDIR, paste0("funnel_", id, ".png")), width = 900, height = 720, res = 130)
  funnel(fit, main = paste("Funnel plot:", id))
  dev.off()
}
small_study_df <- do.call(rbind, SS)
write.csv(small_study_df, file.path(TABDIR, "small_study_results.csv"), row.names = FALSE, na = "")

# ---------------------------------------------------------------- 2. E-values ----
evalue_rr <- function(rr) {
  if (is.na(rr) || rr <= 0) return(NA_real_)
  if (rr >= 1) rr + sqrt(rr * (rr - 1)) else (1 / rr) + sqrt((1 / rr) * ((1 / rr) - 1))
}
EV <- lapply(rr_ids, function(id) {
  dat <- rows_of(id)
  fit <- fit_rr(dat)
  pr <- predict(fit, transf = exp)
  est <- as.numeric(pr$pred); lo <- as.numeric(pr$ci.lb); hi <- as.numeric(pr$ci.ub)
  # CI bound closest to the null (1)
  bound <- if (est >= 1) lo else hi
  data.frame(
    model_id = id, k = nrow(dat),
    rr = round(est, 3), ci_low = round(lo, 3), ci_high = round(hi, 3),
    e_value_point = round(evalue_rr(est), 3),
    e_value_ci_bound = round(evalue_rr(bound), 3),
    note = "E-value = minimum strength of association an unmeasured confounder would need with both exposure and outcome to explain away the estimate (VanderWeele & Ding 2017)",
    stringsAsFactors = FALSE)
})
evalue_df <- do.call(rbind, EV)
write.csv(evalue_df, file.path(TABDIR, "evalue_results.csv"), row.names = FALSE, na = "")

# --------------------------------------------------- 3. cumulative meta-analysis ----
# Primary model, ordered by publication year (hard-coded year map from the registry studies).
YEARS <- c("Cernic Suligoj 2020" = 2020, "Li 2026" = 2026, "Zhang 2026" = 2026)
dat <- rows_of("ecl6_core")
dat$year <- num(YEARS[dat$study_id])
dat <- dat[order(dat$year), , drop = FALSE]
fit <- fit_rr(dat)
cu <- cumul(fit, order = order(dat$year))
png(file.path(FIGDIR, "cumulative_ecl6_core.png"), width = 1000, height = 640, res = 130)
par(mar = c(4, 4, 3, 1))
forest(cu, transf = exp, refline = 1, xlab = "Cumulative RR (95% CI)",
       main = "Cumulative meta-analysis of effective closure at ~6 months (by year)")
dev.off()
cumul_df <- data.frame(
  model_id = "ecl6_core",
  step = seq_len(nrow(dat)),
  added_study = dat$study_id,
  year = dat$year,
  k = seq_len(nrow(dat)),
  rr = round(exp(cu$estimate), 3),
  ci_low = round(exp(cu$ci.lb), 3),
  ci_high = round(exp(cu$ci.ub), 3),
  stringsAsFactors = FALSE)
write.csv(cumul_df, file.path(TABDIR, "cumulative_ecl6_core.csv"), row.names = FALSE, na = "")

# ------------------------------------------------------ 4. influence diagnostics ----
INFL <- list()
re_ids <- unique(MANIFEST$analysis_id[MANIFEST$measure %in% c("RR", "MD")])
re_ids_k3 <- re_ids[vapply(re_ids, function(id) nrow(rows_of(id)) >= 3, logical(1))]
for (id in re_ids_k3) {
  dat <- rows_of(id)
  fit <- if (dat$measure[1] == "MD") {
    es <- escalc(measure = "MD", m1i = dat$m1i, sd1i = dat$sd1i, n1i = dat$n1i,
                 m2i = dat$m2i, sd2i = dat$sd2i, n2i = dat$n2i)
    rma(yi, vi, data = es, method = "REML", test = "knha")
  } else {
    fit_rr(dat)
  }
  inf <- tryCatch(influence(fit), error = function(e) NULL)
  if (is.null(inf)) next
  png(file.path(FIGDIR, paste0("baujat_", id, ".png")), width = 900, height = 720, res = 130)
  baujat(fit, main = paste("Baujat plot:", id))
  dev.off()
  INFL[[id]] <- data.frame(
    model_id = id, study_id = dat$study_id,
    rstudent = round(inf$inf$rstudent, 3),
    dffits = round(inf$inf$dffits, 3),
    cook_d = round(inf$inf$cook.d, 3),
    tau2_del = round(inf$inf$tau2.del, 5),
    stringsAsFactors = FALSE)
}
influence_df <- do.call(rbind, INFL)
write.csv(influence_df, file.path(TABDIR, "influence_results.csv"), row.names = FALSE, na = "")

# -------------------------------------------------------- 5. interaction tests ----
# Population stratum: 12-month effective closure, general (Karagianni 2011, Li 2026) vs
# migraine (Sun 2025, Li 2025). Design stratum: 6-month effective closure, RCT (Zhang 2026)
# vs observational (Li 2026, Cernic Suligoj 2020).
interaction_specs <- list(
  list(model_id = "ecl12_population", label = "population stratum (general vs migraine), 12-month effective closure",
       rows = rbind(rows_of("ecl12_gen_core"), rows_of("ecl12_mig")),
       mod = c(rep("general", nrow(rows_of("ecl12_gen_core"))), rep("migraine", nrow(rows_of("ecl12_mig"))))),
  list(model_id = "ecl6_design", label = "design (randomized vs observational), 6-month effective closure",
       rows = rbind(rows_of("ecl6_core")[rows_of("ecl6_core")$study_id == "Zhang 2026", , drop = FALSE],
                    rows_of("ecl6_core")[rows_of("ecl6_core")$study_id != "Zhang 2026", , drop = FALSE]),
       mod = c("randomized", rep("observational", 2)))
)
INT <- lapply(interaction_specs, function(sp) {
  dat <- sp$rows
  es <- escalc(measure = "RR", ai = dat$ei, bi = dat$ni - dat$ei,
               ci = dat$ec, di = dat$nc - dat$ec, add = 0.5, to = "only0")
  fit <- rma(yi, vi, data = es, mods = ~ factor(sp$mod), method = "REML")
  data.frame(
    model_id = sp$model_id, label = sp$label,
    q_between = round(fit$QM, 3), df = fit$QMdf[1], p = round(fit$QMp, 4),
    included = paste(dat$study_id, collapse = "; "),
    note = "Q-between from random-effects moderator model; exploratory, small k",
    stringsAsFactors = FALSE)
})
interaction_df <- do.call(rbind, INT)
write.csv(interaction_df, file.path(TABDIR, "interaction_results.csv"), row.names = FALSE, na = "")

writeLines(c(
  paste(Sys.time(), "small-study models", nrow(small_study_df)),
  paste(Sys.time(), "e-values", nrow(evalue_df)),
  paste(Sys.time(), "influence rows", nrow(influence_df)),
  paste(Sys.time(), "interactions", nrow(interaction_df)),
  paste(Sys.time(), "done")
), file.path(PROJ, "outputs/logs/06_adjunct_analyses.log"))
cat("ADJUNCT ANALYSES DONE:", nrow(small_study_df), "small-study;", nrow(evalue_df),
    "e-values;", nrow(influence_df), "influence;", nrow(interaction_df), "interactions\n")
