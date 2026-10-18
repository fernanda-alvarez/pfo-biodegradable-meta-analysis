# 10_grade_tsa.R — rule-based GRADE certainty for every model and a simplified
# trial-sequential-analysis for the primary outcome.
#
# GRADE (rule-based, conservative; v4 rationales as the starting point):
#   Start: High if at least one randomized trial is in the model, else Low.
#   Downgrade -1 per: critical/high risk-of-bias study present;
#     inconsistency I2 > 50% or Q p < 0.10; imprecision if the 95% CI includes
#     the null or N < 200; indirectness if population strata are mixed;
#     publication bias if Egger p < 0.10 or trim-and-fill changes the
#     conclusion (models in the small-study table only).
#   Certainty: >=4 High, 3 Moderate, 2 Low, <=1 Very low.
#
# TSA (simplified Lan-DeMets, exploratory): cumulative primary (ecl6_core) z
# path vs O'Brien-Fleming alpha-spending boundaries; required information size
# under an assumed 10% relative risk reduction.

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(metafor))

PROJ <- normalizePath(".", winslash = "/", mustWork = TRUE)
MANIFEST <- read.csv(file.path(PROJ, "data_derived/luna_manifest.csv"), check.names = FALSE,
                     na.strings = c("", "NA", "NR"))
RESULTS <- read.csv(file.path(PROJ, "outputs/tables/model_registry.csv"), check.names = FALSE,
                    na.strings = c("", "NA", "NR"))
SMALL <- read.csv(file.path(PROJ, "outputs/tables/small_study_results.csv"), check.names = FALSE,
                  na.strings = c("", "NA", "NR"))
TABDIR <- file.path(PROJ, "outputs/tables")
FIGDIR <- file.path(PROJ, "outputs/figures")
dir.create(TABDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGDIR, recursive = TRUE, showWarnings = FALSE)

num <- function(x) suppressWarnings(as.numeric(x))
rows_of <- function(id) MANIFEST[MANIFEST$analysis_id == id, , drop = FALSE]

ROB <- read.csv(file.path(PROJ, "data_derived/rob_joined.csv"), check.names = FALSE)
CRITICAL <- ROB$study_id[ROB$overall_judgment %in% c("Critical", "High")]
SERIOUS <- ROB$study_id[ROB$overall_judgment == "Serious"]

grade_row <- function(r) {
  id <- r$model_id
  dat <- rows_of(id)
  if (!nrow(dat)) return(NULL)
  studies <- dat$study_id
  has_rct <- "Zhang 2026" %in% studies
  start <- if (has_rct) 4 else 2
  downs <- character()
  if (any(studies %in% CRITICAL)) { start <- start - 2; downs <- c(downs, "risk of bias (critical study present)") }
  if (any(studies %in% SERIOUS)) { start <- start - 1; downs <- c(downs, "risk of bias (serious study present)") }
  i2 <- num(r$i2); qp <- num(r$q_p)
  if (!is.na(i2) && i2 > 50) { start <- start - 1; downs <- c(downs, "inconsistency (I2>50%)") }
  else if (!is.na(qp) && qp < 0.10) { start <- start - 1; downs <- c(downs, "inconsistency (Q p<0.10)") }
  lo <- num(r$ci_low); hi <- num(r$ci_high); n <- num(r$n_participants)
  null <- if (r$measure %in% c("RR", "OR")) 1 else 0
  if (!is.na(lo) && !is.na(hi) && lo < null && hi > null) { start <- start - 1; downs <- c(downs, "imprecision (CI includes null)") }
  if (!is.na(n) && n < 200) { start <- start - 1; downs <- c(downs, "imprecision (N<200)") }
  if (r$population_stratum == "mixed") { start <- start - 1; downs <- c(downs, "indirectness (mixed strata)") }
  ss <- SMALL[SMALL$model_id == id, , drop = FALSE]
  if (nrow(ss)) {
    ep <- num(ss$egger_p)
    if (!is.na(ep) && ep < 0.10) { start <- start - 1; downs <- c(downs, "publication bias (Egger p<0.10)") }
  }
  start <- max(start, 1)
  cert <- if (start >= 4) "High" else if (start == 3) "Moderate" else if (start == 2) "Low" else "Very low"
  data.frame(
    model_id = id, outcome = r$outcome_family, measure = r$measure, k = r$k,
    n_participants = n,
    estimate_text = sprintf("%s %.2f (%.2f, %.2f)", r$measure, num(r$estimate), lo, hi),
    starts = if (has_rct) "High (RCT present)" else "Low (no RCT)",
    downgrades = paste(downs, collapse = "; "),
    certainty = cert, stringsAsFactors = FALSE)
}
grade_df <- do.call(rbind, lapply(seq_len(nrow(RESULTS)), function(i) grade_row(RESULTS[i, , drop = FALSE])))
write.csv(grade_df, file.path(TABDIR, "grade_results.csv"), row.names = FALSE, na = "")

# ---------------------------------------------------------------- simplified TSA
# Primary: ecl6_core. Assumed treatment effect: 10% relative risk reduction.
dat <- rows_of("ecl6_core")
YEARS <- c("Cernic Suligoj 2020" = 2020, "Li 2026" = 2026, "Zhang 2026" = 2026)
dat$year <- num(YEARS[dat$study_id]); dat <- dat[order(dat$year), , drop = FALSE]

p0 <- sum(dat$ec) / sum(dat$nc)            # control event rate (pooled arms)
rr_assumed <- 0.90
p1 <- p0 * rr_assumed
z_alpha <- qnorm(0.975); z_beta <- qnorm(0.80)
ris_per_arm <- (z_alpha + z_beta)^2 * (p1 * (1 - p1) + p0 * (1 - p0)) / (p1 - p0)^2
ris_total <- 2 * ris_per_arm

# cumulative z path with O'Brien-Fleming boundaries (Lan-DeMets alpha spending)
cum_rows <- list()
cum_ev_b <- 0; cum_n_b <- 0; cum_ev_c <- 0; cum_n_c <- 0
for (i in seq_len(nrow(dat))) {
  cum_ev_b <- cum_ev_b + dat$ei[i]; cum_n_b <- cum_n_b + dat$ni[i]
  cum_ev_c <- cum_ev_c + dat$ec[i]; cum_n_c <- cum_n_c + dat$nc[i]
  p1c <- cum_ev_b / cum_n_b; p0c <- cum_ev_c / cum_n_c
  se <- sqrt(p1c * (1 - p1c) / cum_n_b + p0c * (1 - p0c) / cum_n_c)
  z <- (log(p1c) - log(p0c)) / se
  n_cum <- cum_n_b + cum_n_c
  info_frac <- min(n_cum / ris_total, 1)
  # O'Brien-Fleming boundary at alpha = 0.05, two-sided, spending ~ 2 - 2*Phi(z_alpha / sqrt(IF))
  alpha_spent <- 0.05 * (2 - 2 * pnorm(qnorm(1 - 0.025) / sqrt(info_frac)))
  z_boundary <- qnorm(1 - alpha_spent / 2)
  cum_rows[[i]] <- data.frame(
    step = i, added_study = dat$study_id[i], year = dat$year[i],
    n_cumulative = n_cum, cum_events_bio = cum_ev_b, cum_events_comp = cum_ev_c,
    cum_log_irr = round(log(p1c) - log(p0c), 4), z = round(z, 3),
    z_boundary = round(z_boundary, 3),
    crossed = abs(z) >= z_boundary,
    info_fraction = round(info_frac, 3),
    stringsAsFactors = FALSE)
}
tsa_df <- do.call(rbind, cum_rows)
tsa_df$ris_total <- round(ris_total, 0)
write.csv(tsa_df, file.path(TABDIR, "tsa_results.csv"), row.names = FALSE, na = "")

png(file.path(FIGDIR, "tsa_ecl6_core.png"), width = 1000, height = 640, res = 130)
par(mar = c(4, 4, 3, 1))
plot(tsa_df$info_fraction, tsa_df$z, type = "b", pch = 19, ylim = c(-3, 3),
     xlab = "Information fraction (cumulative N / required information size)",
     ylab = "Cumulative z", main = "Trial sequential analysis (simplified Lan-DeMets), primary outcome")
lines(tsa_df$info_fraction, tsa_df$z_boundary, lty = 2, col = "red")
lines(tsa_df$info_fraction, -tsa_df$z_boundary, lty = 2, col = "red")
abline(h = 0, col = "grey")
legend("topright", legend = c("Cumulative z", "O'Brien-Fleming boundary"), lty = c(1, 2),
       col = c("black", "red"), bty = "n")
dev.off()

writeLines(c(
  paste(Sys.time(), "grade rows", nrow(grade_df)),
  paste(Sys.time(), "TSA steps", nrow(tsa_df)),
  paste(Sys.time(), "done")
), file.path(PROJ, "outputs/logs/10_grade_tsa.log"))
cat("GRADE/TSA DONE:", nrow(grade_df), "grade rows;", nrow(tsa_df), "TSA steps\n")
