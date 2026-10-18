# Assertions for the rule-based GRADE table and simplified TSA.
for (f in c("grade_results.csv", "tsa_results.csv")) {
  if (!file.exists(file.path("outputs/tables", f))) stop("Missing table: ", f)
}
if (!file.exists("outputs/figures/tsa_ecl6_core.png")) stop("Missing TSA figure")

gr <- read.csv("outputs/tables/grade_results.csv", check.names = FALSE)
if (nrow(gr) != 26L) stop("GRADE table must cover all 26 registry models")
if (!all(gr$certainty %in% c("High", "Moderate", "Low", "Very low"))) stop("Invalid certainty labels")
ecl6 <- gr[gr$model_id == "ecl6_core", , drop = FALSE]
if (ecl6$certainty != "Low") stop("Primary certainty must be Low (RCT present, critical Cernic row, CI includes null)")

ts <- read.csv("outputs/tables/tsa_results.csv", check.names = FALSE)
if (nrow(ts) != 3L) stop("TSA must have 3 cumulative steps")
if (any(ts$crossed)) stop("TSA must not cross the boundary at the current information size")
if (ts$n_cumulative[3] != 415L) stop("Final cumulative N must be 415")
if (ts$ris_total[1] < 300 || ts$ris_total[1] > 700) stop("RIS out of plausible range")
cat("GRADE/TSA TESTS PASS\n")
