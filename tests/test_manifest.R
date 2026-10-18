# Test-first validation for the Luna analysis manifest.
# This file is intentionally runnable with base R only.

manifest_path <- file.path("data_derived", "luna_manifest.csv")
if (!file.exists(manifest_path)) {
  stop("RED: manifest is missing: ", manifest_path)
}

manifest <- read.csv(manifest_path, stringsAsFactors = FALSE, check.names = FALSE)
required <- c("analysis_id", "row_id", "study_id", "outcome_family", "window",
              "population_stratum", "status", "ei", "ni", "ec", "nc")
missing_cols <- setdiff(required, names(manifest))
if (length(missing_cols)) stop("Missing manifest columns: ", paste(missing_cols, collapse = ", "))

key <- paste(manifest$analysis_id, manifest$study_id, sep = "::")
if (anyDuplicated(key)) stop("Duplicate study rows within an analysis model")

primary <- subset(manifest, analysis_id == "ecl6_core")
if (nrow(primary) != 3L) stop("Primary ecl6_core must contain exactly 3 studies")
if (!setequal(primary$study_id, c("Li 2026", "Zhang 2026", "Cernic Suligoj 2020"))) {
  stop("Primary ecl6_core study set is incorrect")
}
if (any(primary$row_id == "C8:3")) stop("Countless Zhang duplicate C8:3 entered the primary model")

rare <- subset(manifest, analysis_id %in% c("drt_any", "embolization", "pericardial"))
if (any(rare$study_id == "Darsaklis 2014")) stop("Narrative-only Darsaklis rare-event rows were included")
if (any(rare$status %in% c("narrative", "NR"))) stop("Narrative/NR rare-event rows were included")

numeric_cols <- c("ei", "ni", "ec", "nc")
for (nm in numeric_cols) {
  vals <- manifest[[nm]]
  if (any(!is.na(vals) & vals < 0)) stop("Negative value in ", nm)
}
if (any(!is.na(manifest$ei) & !is.na(manifest$ni) & manifest$ei > manifest$ni)) stop("Intervention events exceed denominator")
if (any(!is.na(manifest$ec) & !is.na(manifest$nc) & manifest$ec > manifest$nc)) stop("Comparator events exceed denominator")

cat("MANIFEST TESTS PASS\n")
