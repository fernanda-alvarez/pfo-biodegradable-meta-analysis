# Assertions for the PDF re-extraction analyses (migraine relief, person-time IRR,
# MIDAS change-from-baseline).
for (f in c("migraine_relief_results.csv", "irr_results.csv", "midas_change_results.csv")) {
  if (!file.exists(file.path("outputs/tables", f))) stop("Missing extraction table: ", f)
}
if (!file.exists("outputs/figures/li2025_km_figure.png")) stop("Missing Li 2025 KM figure export")

rel <- read.csv("outputs/tables/migraine_relief_results.csv", check.names = FALSE)
if (nrow(rel) != 2L) stop("Migraine relief must have 2 rows (Sun 2025 quantified, Li 2025 narrative)")
sun <- rel[rel$study_id == "Sun 2025", , drop = FALSE]
if (sun$events_bio != 39 || sun$n_bio != 56 || sun$events_comp != 54 || sun$n_comp != 77) {
  stop("Sun 2025 migraine counts incorrect")
}
li <- rel[rel$study_id == "Li 2025", , drop = FALSE]
if (!is.na(li$events_bio)) stop("Li 2025 must be narrative-only (no imputed counts)")

irr <- read.csv("outputs/tables/irr_results.csv", check.names = FALSE)
if (nrow(irr) != 4L) stop("IRR table must have 4 rows (core, with-Luani, 2 FU sensitivities)")
core <- irr[irr$model_id == "cve_irr_core", , drop = FALSE]
if (nrow(core) != 1L || abs(core$irr - 2.144) > 0.01) stop("Core IRR changed unexpectedly")
withl <- irr[irr$model_id == "cve_irr_with_luani", , drop = FALSE]
if (nrow(withl) != 1L || abs(withl$irr - 1.495) > 0.02) stop("With-Luani IRR changed unexpectedly")
sens <- irr[grepl("cernic_", irr$model_id), , drop = FALSE]
if (nrow(sens) != 2L) stop("Cernic FU sensitivity rows missing")
if (max(abs(sens$irr - core$irr)) > 0.01) stop("Cernic FU sensitivity should be insensitive to the FU assumption")

mid <- read.csv("outputs/tables/midas_change_results.csv", check.names = FALSE)
if (nrow(mid) != 3L) stop("MIDAS change must have 3 rows (r = 0.25, 0.5, 0.75)")
if (any(mid$md > 0 | mid$md < -5)) stop("MIDAS change MD out of plausible range")
if (any(mid$ci_low > 0 | mid$ci_high < 0)) stop("All MIDAS-change CIs must include zero (no between-group difference)")
cat("EXTRACTION TESTS PASS\n")
