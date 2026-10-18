# Assertions for the approximate Bayesian analyses.
path <- "outputs/tables/bayesian_results.csv"
if (!file.exists(path)) stop("Missing bayesian_results.csv")
bay <- read.csv(path, check.names = FALSE)

if (nrow(bay) < 25L) stop("Bayesian table must cover all 26 models")
if (any(bay$p_gt_null < 0 | bay$p_gt_null > 1)) stop("Invalid posterior probabilities")
if (any(bay$ci_low > bay$ci_high, na.rm = TRUE)) stop("Inverted credible intervals")

ecl6 <- bay[bay$model_id == "ecl6_core", , drop = FALSE]
if (nrow(ecl6) != 1L) stop("Primary model missing")
if (abs(ecl6$median - 1.01) > 0.01) stop("Primary Bayesian median changed unexpectedly")
if (ecl6$p_gt_null < 0.5 || ecl6$p_gt_null > 0.75) stop("Primary P(RR>1) out of expected range (0.5-0.75)")
if (ecl6$ci_low > 1 || ecl6$ci_high < 1) stop("Primary credible interval must include the null")

cve <- bay[bay$model_id == "cve_late_core", , drop = FALSE]
if (nrow(cve) != 1L || cve$p_gt_null < 0.75) stop("cve_late_core P(RR>1) should exceed 0.75")

drt <- bay[bay$model_id == "drt_any", , drop = FALSE]
if (nrow(drt) != 1L || drt$median <= 0) stop("drt_any posterior RD median must be positive")

for (id in c("RR random-effects (conditional REML tau2)", "risk difference, pooled-arm Beta(1,1) priors",
             "pooled proportion, Beta(1,1) prior")) {
  if (!id %in% bay$family) stop("Missing Bayesian family: ", id)
}
cat("BAYESIAN TESTS PASS\n")
