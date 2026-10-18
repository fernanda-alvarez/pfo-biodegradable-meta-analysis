# Assertions for the prespecified sensitivity/subgroup ledger.
res <- read.csv("outputs/tables/sensitivity_results.csv", check.names = FALSE)

if (nrow(res) < 80L) stop("Sensitivity ledger too small: ", nrow(res))
if (any(!nzchar(res$sensitivity_family))) stop("Missing sensitivity family tags")

# S1: every RR model must have an MH fixed-effect row
manifest <- read.csv("data_derived/luna_manifest.csv", check.names = FALSE)
rr_ids <- unique(manifest$analysis_id[manifest$measure == "RR"])
s1 <- res[res$sensitivity_family == "S1 MH fixed effect", , drop = FALSE]
missing_s1 <- setdiff(paste0(rr_ids, "__mh_fixed"), s1$model_id)
if (length(missing_s1)) stop("Missing MH fixed-effect rows: ", paste(missing_s1, collapse = ", "))

# S1 primary must be close to the REML primary and use all 3 studies
primary_mh <- s1$estimate[s1$model_id == "ecl6_core__mh_fixed"]
if (length(primary_mh) != 1L || is.na(primary_mh)) stop("Missing primary MH estimate")
if (abs(primary_mh - 1.009787096) > 1e-6) stop("Primary MH fixed estimate changed unexpectedly")

# S2: leave-one-out rows must exist for the primary model (all three studies)
loo_primary <- res[grepl("^ecl6_core__loo_", res$model_id), , drop = FALSE]
if (nrow(loo_primary) != 3L) stop("Primary leave-one-out must have exactly 3 rows")

# S3: under the permissive ROBINS-I consensus no study is critical, so the
# critical-exclusion row is a documented no-op; the serious-exclusion variant
# for the primary must leave Li 2026 + Zhang 2026 (Cernic Suligoj 2020 excluded)
rob3 <- res[res$model_id == "ecl6_core__no_critical_rob", , drop = FALSE]
if (nrow(rob3) != 1L) stop("Missing primary critical-RoB exclusion row")
if (rob3$k != 3L || !grepl("Cernic", rob3$included, fixed = TRUE)) stop("Critical-RoB exclusion must be a no-op (k=3) under the permissive consensus")
if (!grepl("identical to base", rob3$note, fixed = TRUE)) stop("Critical-RoB exclusion note must say 'identical to base'")
ser3 <- res[res$model_id == "ecl6_core__no_serious_rob", , drop = FALSE]
if (nrow(ser3) != 1L) stop("Missing primary serious-RoB exclusion row")
if (ser3$k != 2L || grepl("Cernic", ser3$included, fixed = TRUE)) stop("Serious-RoB exclusion set must be Li 2026 + Zhang 2026")
if (!grepl("Li 2026", ser3$included, fixed = TRUE) || !grepl("Zhang 2026", ser3$included, fixed = TRUE)) stop("Serious-RoB exclusion set incorrect")

# S4: fully-biodegradable-only primary must be Li 2026 + Zhang 2026
fb <- res[res$model_id == "ecl6_core__fully_biodegradable_only", , drop = FALSE]
if (nrow(fb) != 1L || fb$k != 2L) stop("Fully-biodegradable-only primary row incorrect")

# S5: excluding Li 2025 from the 12-month migraine models leaves Sun 2025 (k=1 descriptive)
for (id in c("ecl12_mig__direct_only", "close12_mig_strict__direct_only")) {
  r <- res[res$model_id == id, , drop = FALSE]
  if (nrow(r) != 1L || r$k != 1L || !grepl("Sun 2025", r$included, fixed = TRUE)) {
    stop(id, " must leave exactly Sun 2025")
  }
}

# S6: rare-event bases have an add=0.5 continuity-correction row
for (id in c("drt_any__cc05", "embolization__cc05", "pericardial__cc05")) {
  if (nrow(res[res$model_id == id, , drop = FALSE]) != 1L) stop("Missing continuity-correction row: ", id)
}

# Consistency: pooled rows (k >= 2, non-descriptive) must carry estimates and CIs
pooled <- res[!is.na(res$estimate) & res$k >= 2, , drop = FALSE]
if (any(is.na(pooled$ci_low)) || any(is.na(pooled$ci_high))) stop("Pooled sensitivity rows must have CIs")
if (any(pooled$ci_low > pooled$ci_high, na.rm = TRUE)) stop("Inverted CIs in sensitivity ledger")

cat("SENSITIVITY TESTS PASS\n")
