# 05b_rob_join.R — rebuild the risk-of-bias ledger from the selected sources:
#
#   ROBINS-I (non-randomized comparative): 06_Analysis Bias/robvis_ROBINS-I_v2_permissive_tier2.csv
#   JBI Case Series (single-arm):          06_Analysis Bias/JBI Case Series/JBI_Case_Series_7_final_single_arm.csv
#   RoB 2 (randomized):                    06_Analysis Bias/robvis_RoB2_Zhang2026.csv
#
# (Per user instruction 2026-08-05: the PERMISSIVE ROBINS-I consensus is the
# analysis default, replacing the strict consensus; JBI appraisals come from
# the JBI Case Series folder's finalized 7-study file.)
#
# Output: data_derived/rob_joined.csv (v5) — 16 rows, registry study IDs.

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))

PROJ <- normalizePath(".", winslash = "/", mustWork = TRUE)
data_dir <- if (dir.exists(file.path(PROJ, "data"))) file.path(PROJ, "data") else file.path(dirname(PROJ), "06_Analysis Bias")

robins_file <- if (file.exists(file.path(data_dir, "robvis_ROBINS-I_v2_permissive_tier2.csv"))) {
  file.path(data_dir, "robvis_ROBINS-I_v2_permissive_tier2.csv")
} else {
  file.path(dirname(PROJ), "06_Analysis Bias", "robvis_ROBINS-I_v2_permissive_tier2.csv")
}

jbi_file <- if (file.exists(file.path(data_dir, "JBI_Case_Series_7_final_single_arm.csv"))) {
  file.path(data_dir, "JBI_Case_Series_7_final_single_arm.csv")
} else {
  file.path(dirname(PROJ), "06_Analysis Bias", "JBI Case Series", "JBI_Case_Series_7_final_single_arm.csv")
}

rob2_file <- if (file.exists(file.path(data_dir, "robvis_RoB2_Zhang2026.csv"))) {
  file.path(data_dir, "robvis_RoB2_Zhang2026.csv")
} else {
  file.path(dirname(PROJ), "06_Analysis Bias", "robvis_RoB2_Zhang2026.csv")
}

robins <- read_csv(robins_file, show_col_types = FALSE)
jbi <- read_csv(jbi_file, show_col_types = FALSE)
rob2 <- read_csv(rob2_file, show_col_types = FALSE)

ROBINS_MAP <- c(
  "Li (2026)" = "Li 2026", "Li (2025)" = "Li 2025", "Sun (2025)" = "Sun 2025",
  "Van den Branden (2010)" = "Van den Branden 2010a", "Karagianni (2011)" = "Karagianni 2011",
  "Blerim Luani (2015)" = "Luani 2015", "Konstadina Darsaklis (2014)" = "Darsaklis 2014",
  "Natasa Cernic (2020)" = "Cernic Suligoj 2020")
JBI_MAP <- c(
  "Du 2022" = "Du 2022", "He 2025" = "He 2025", "Karagianni 2020" = "Karagianni 2020",
  "Lin 2022" = "Lin 2022", "Ussia 2009" = "Ussia 2009", "Van den Branden 2009" = "Van den Branden 2009",
  "Van den Branden 2010" = "Van den Branden 2010b")
ROB2_MAP <- c("Zhang (2026)" = "Zhang 2026")

# ---- ROBINS-I (permissive tier 2) ----
rb <- robins
rb$study_id <- unname(ROBINS_MAP[rb$Study])
rb$tool <- "ROBINS-I"
rb$overall_judgment <- sub(" risk of bias$", "", rb$Overall)
rb$high_risk_flag <- TRUE
rb$concern_label <- paste0("ROBINS-I: ", rb$Overall)
rb$overall_rationale <- vapply(seq_len(nrow(rb)), function(i) {
  dom <- c("D1" = rb$D1[i], "D2" = rb$D2[i], "D3" = rb$D3[i], "D4" = rb$D4[i],
           "D5" = rb$D5[i], "D6" = rb$D6[i])
  flagged <- names(dom)[dom != "Low risk of bias"]
  if (!length(flagged)) {
    "All six domains at low risk of bias (permissive consensus, tier 2)."
  } else {
    paste0("Permissive ROBINS-I consensus (tier 2): ", rb$Overall[i],
           " overall; flagged domains: ",
           paste0(flagged, " (", unname(dom[flagged]), ")",
                  collapse = "; "), ".")
  }
}, character(1))
rb <- rb[, c("study_id", "tool", "overall_judgment", "high_risk_flag", "concern_label", "overall_rationale")]

# ---- JBI Case Series (finalized 7-study file) ----
jb <- jbi
jb$study_id <- unname(JBI_MAP[jb$Study])
jb$tool <- "JBI"
jb$overall_judgment <- ifelse(jb$Overall == "Low concern", "Low concern",
                              ifelse(jb$Overall == "Unclear", "Unclear", jb$Overall))
jb$high_risk_flag <- TRUE
jb$concern_label <- paste0("JBI: ", jb$Overall)
jb$overall_rationale <- vapply(seq_len(nrow(jb)), function(i) {
  items <- c("Q1", "Q2", "Q3", "Q4", "Q5", "Q6", "Q7", "Q8", "Q9", "Q10")
  unclear <- items[jb[i, items] == "Unclear"]
  if (jb$Overall[i] == "Low concern") {
    "JBI Case Series 10-item appraisal: all items Yes (Low concern)."
  } else {
    paste0("JBI Case Series 10-item appraisal: overall Unclear because ",
           length(unclear), " item(s) Unclear (", paste(unclear, collapse = ", "),
           "); remaining items Yes.")
  }
}, character(1))
jb <- jb[, c("study_id", "tool", "overall_judgment", "high_risk_flag", "concern_label", "overall_rationale")]

# ---- RoB 2 (Zhang 2026) ----
z <- rob2
z$study_id <- unname(ROB2_MAP[z$Study])
z$tool <- "RoB 2"
z$overall_judgment <- sub(" risk of bias$", "", z$Overall)
z$high_risk_flag <- TRUE
z$concern_label <- paste0("RoB 2: ", z$Overall)
z$overall_rationale <- vapply(seq_len(nrow(z)), function(i) {
  dom <- c("D1" = z$D1[i], "D2" = z$D2[i], "D3" = z$D3[i], "D4" = z$D4[i], "D5" = z$D5[i])
  flagged <- names(dom)[dom != "Low risk of bias"]
  paste0("RoB 2 (Zhang 2026): overall ", z$Overall[i],
         "; flagged domains: ", paste0(flagged, " (", unname(dom[flagged]), ")",
                                       collapse = "; "), ".")
}, character(1))
z <- z[, c("study_id", "tool", "overall_judgment", "high_risk_flag", "concern_label", "overall_rationale")]

join <- rbind(rb, jb, z)
if (anyDuplicated(join$study_id)) stop("Duplicate study IDs in the joined RoB ledger")
if (nrow(join) != 16L) stop("Expected 16 studies, got ", nrow(join))

write.csv(join, file.path(PROJ, "data_derived/rob_joined.csv"), row.names = FALSE, na = "")
if (dir.exists(file.path(PROJ, "data"))) {
  write.csv(join, file.path(PROJ, "data/rob_joined.csv"), row.names = FALSE, na = "")
}
dir.create(file.path(PROJ, "outputs/logs"), recursive = TRUE, showWarnings = FALSE)
writeLines(c(paste(Sys.time(), "rob_joined rows", nrow(join)),
             paste(Sys.time(), "critical", sum(join$overall_judgment %in% c("Critical", "High"))),
             paste(Sys.time(), "serious", sum(join$overall_judgment == "Serious")),
             paste(Sys.time(), "done")),
           file.path(PROJ, "outputs/logs/05b_rob_join.log"))
cat("ROB JOIN DONE:", nrow(join), "rows;",
    sum(join$overall_judgment %in% c("Critical", "High")), "critical;",
    sum(join$overall_judgment == "Serious"), "serious\n")
