# Model-output assertions for the Luna redo.
res <- read.csv("outputs/tables/model_registry.csv", check.names = FALSE)
manifest <- read.csv("data_derived/luna_manifest.csv", check.names = FALSE)

get_one <- function(id) res[res$model_id == id, , drop = FALSE]
primary <- get_one("ecl6_core")
if (nrow(primary) != 1L) stop("Missing primary result")
if (abs(primary$estimate - 1.009903945) > 1e-6) stop("Primary RR changed unexpectedly")
if (primary$k != 3L || primary$n_participants != 415L) stop("Primary k/N incorrect")

midas <- get_one("midas_12m")
if (midas$n_participants != 272L) stop("MIDAS participant total must be 272")

rare <- res[res$model_id %in% c("drt_any", "embolization", "pericardial"), , drop = FALSE]
if (any(!grepl("no continuity correction", rare$method, ignore.case = TRUE))) stop("Rare-event primary uses continuity correction")
if (any(grepl("Darsaklis 2014", rare$included, fixed = TRUE))) stop("Darsaklis entered rare-event results")

sa <- get_one("sa_proc_success")
if (sa$n_participants != 221L) stop("Single-arm procedural-success N must be 221")

# --- v4 coverage-gap models (approved audit additions) ---
new_ids <- c("close12_mig_strict", "shunt_g6_general", "af_late_full", "af_late_full_with_luani",
             "af_any_6m", "af_any_late", "af_any_late_with_darsaklis",
             "ecl6_expanded_early", "ecl12_gen_with_luani")
for (id in new_ids) {
  if (nrow(get_one(id)) != 1L) stop("Missing coverage model: ", id)
}

strict12 <- get_one("close12_mig_strict")
if (strict12$k != 2L || strict12$n_participants != 228L) stop("close12_mig_strict k/N incorrect")
if (!setequal(strsplit(strict12$included, "; ", fixed = TRUE)[[1]], c("Sun 2025", "Li 2025"))) {
  stop("close12_mig_strict study set incorrect")
}

shunt <- get_one("shunt_g6_general")
if (shunt$k != 2L || shunt$n_participants != 191L) stop("shunt_g6_general k/N incorrect")
if (!grepl("Karagianni 2011", shunt$included, fixed = TRUE) || !grepl("Li 2026", shunt$included, fixed = TRUE)) {
  stop("shunt_g6_general study set incorrect")
}

af_full <- get_one("af_late_full")
if (af_full$k != 3L || !grepl("Zhang 2026", af_full$included, fixed = TRUE)) stop("af_late_full must include the Zhang 2026 double-zero row")

af_any <- get_one("af_any_late")
if (af_any$k != 4L || !grepl("Li 2026", af_any$included, fixed = TRUE)) stop("af_any_late k/included incorrect")

expanded <- get_one("ecl6_expanded_early")
if (expanded$k != 3L || !grepl("Darsaklis 2014", expanded$included, fixed = TRUE)) stop("ecl6_expanded_early incorrect")

with_luani12 <- get_one("ecl12_gen_with_luani")
if (with_luani12$k != 3L || !grepl("Luani 2015", with_luani12$included, fixed = TRUE)) stop("ecl12_gen_with_luani incorrect")

# Every model in the registry must trace to manifest rows with usable counts
if (nrow(res) != length(unique(manifest$analysis_id))) stop("Registry/analysis_id count mismatch")

# --- bleeding/vascular complications (new safety outcome) ---
bleed <- get_one("bleeding")
if (nrow(bleed) != 1L) stop("Missing bleeding model")
if (bleed$k != 4L || bleed$n_participants != 430L) stop("bleeding k/N incorrect")
if (!grepl("Van den Branden 2010a", bleed$included, fixed = TRUE) || !grepl("Karagianni 2011", bleed$included, fixed = TRUE)) {
  stop("bleeding study set incorrect")
}
if (!grepl("no continuity correction", bleed$method, ignore.case = TRUE)) stop("bleeding must use MH RD without universal continuity correction")

# --- prediction intervals (k>=3 RE models must carry them; RD models must not) ---
pi_models <- res[res$measure %in% c("RR", "MD", "proportion") & res$k >= 3, , drop = FALSE]
if (any(is.na(pi_models$pi_low)) || any(is.na(pi_models$pi_high))) stop("k>=3 RE models must report prediction intervals")
rd_models <- res[res$measure == "RD", , drop = FALSE]
if (any(!is.na(rd_models$pi_low))) stop("RD (fixed MH) models must not carry prediction intervals")

cat("MODEL TESTS PASS\n")
