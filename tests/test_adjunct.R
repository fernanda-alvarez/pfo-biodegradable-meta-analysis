# Assertions for the adjunct analyses (small-study effects, E-values, influence,
# cumulative meta-analysis, interaction tests).
adjunct_files <- c(
  "small_study_results.csv", "evalue_results.csv", "influence_results.csv",
  "interaction_results.csv", "cumulative_ecl6_core.csv"
)
for (f in adjunct_files) {
  if (!file.exists(file.path("outputs/tables", f))) stop("Missing adjunct table: ", f)
}
ss <- read.csv("outputs/tables/small_study_results.csv", check.names = FALSE)
ev <- read.csv("outputs/tables/evalue_results.csv", check.names = FALSE)
inf <- read.csv("outputs/tables/influence_results.csv", check.names = FALSE)
it <- read.csv("outputs/tables/interaction_results.csv", check.names = FALSE)
cu <- read.csv("outputs/tables/cumulative_ecl6_core.csv", check.names = FALSE)

if (nrow(ss) < 10L) stop("Small-study table must cover all k>=3 RR models")
if (any(!is.na(ss$egger_p) & ss$egger_p < 0) || any(!is.na(ss$egger_p) & ss$egger_p > 1)) stop("Invalid Egger p-values")
if (!"ecl6_core" %in% ss$model_id) stop("Primary model missing from small-study table")

if (nrow(ev) < 17L) stop("E-value table must cover all comparative RR models")
if (any(!is.na(ev$e_value_point) & ev$e_value_point < 1)) stop("E-values below 1 are invalid")
cve <- ev[ev$model_id == "cve_late_core", , drop = FALSE]
if (nrow(cve) != 1L || abs(cve$e_value_point - 3.697) > 0.01) stop("cve_late_core E-value changed unexpectedly")

if (nrow(inf) < 30L) stop("Influence table too small")
if (any(!is.finite(inf$rstudent)) || any(!is.finite(inf$dffits))) stop("Influence values must be finite")
# Rare-event models with continuity correction at k=3 legitimately produce large
# studentized residuals (e.g. cve_late_core / Cernic Suligoj 2020); pin the known value.
cernic <- inf[inf$model_id == "cve_late_core" & inf$study_id == "Cernic Suligoj 2020", , drop = FALSE]
if (nrow(cernic) != 1L || abs(cernic$rstudent - 22.674) > 0.01) stop("cve_late_core influence changed unexpectedly")

if (nrow(it) != 2L) stop("Interaction tests must have exactly 2 rows")
if (any(it$p < 0 | it$p > 1)) stop("Invalid interaction p-values")

if (nrow(cu) != 3L) stop("Cumulative meta-analysis must have 3 steps")
if (cu$k[3] != 3L || abs(cu$rr[3] - 1.01) > 0.01) stop("Cumulative final step must equal the primary RR")

for (f in c("funnel_ecl6_core.png", "cumulative_ecl6_core.png", "baujat_ecl6_core.png")) {
  if (!file.exists(file.path("outputs/figures", f))) stop("Missing adjunct figure: ", f)
}
cat("ADJUNCT TESTS PASS\n")
