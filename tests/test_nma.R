# Assertions for the device-family network meta-analysis.
if (!file.exists("outputs/tables/nma_results.csv")) stop("Missing nma_results.csv")
nma <- read.csv("outputs/tables/nma_results.csv", check.names = FALSE)

if (nrow(nma) != 12L) stop("NMA table must have 12 contrasts (4 networks x 3)")
nets <- unique(nma$model_id)
if (!setequal(nets, c("nma_6m_core", "nma_12m_core", "nma_12m_with_luani", "nma_6m_expanded"))) {
  stop("Unexpected NMA network ids")
}
for (nm in nets) {
  sub <- nma[nma$model_id == nm, , drop = FALSE]
  if (nrow(sub) != 3L) stop("Each network must have 3 contrasts")
  if (any(sub$ci_low > sub$ci_high)) stop("Inverted NMA CIs")
}
core6 <- nma[nma$model_id == "nma_6m_core", , drop = FALSE]
be <- core6[core6$contrast == "biodegradable vs metallic", , drop = FALSE]
if (nrow(be) != 1L || abs(be$or - 1.065) > 0.01) stop("6m biodegradable OR changed unexpectedly")
ba <- core6[core6$contrast == "bioabsorbable vs metallic", , drop = FALSE]
if (nrow(ba) != 1L || ba$or < 1) stop("6m bioabsorbable OR should exceed 1 (BioSTAR numerically better closure)")

exp6 <- nma[nma$model_id == "nma_6m_expanded", , drop = FALSE]
bax <- exp6[exp6$contrast == "bioabsorbable vs metallic", , drop = FALSE]
if (!(bax$ci_low > 1)) stop("Expanded 6m bioabsorbable OR CI should exclude 1 (labeled sensitivity)")

for (f in c("nma_forest_nma_6m_core.png", "nma_forest_nma_12m_core.png")) {
  if (!file.exists(file.path("outputs/figures", f))) stop("Missing NMA forest: ", f)
}
cat("NMA TESTS PASS\n")
