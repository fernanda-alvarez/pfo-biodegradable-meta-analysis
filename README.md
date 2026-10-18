# Efficacy and Safety of Biodegradable versus Metallic Occluders for Patent Foramen Ovale (PFO) Closure: A Systematic Review and Meta-Analysis

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![R >= 4.0.0](https://img.shields.io/badge/R-%3E%3D%204.0.0-blue)](https://www.r-project.org/)
[![metafor](https://img.shields.io/badge/metafor-4.4+-green.svg)](https://cran.r-project.org/package=metafor)
[![PRISMA 2020](https://img.shields.io/badge/PRISMA-2020%20Compliant-brightgreen.svg)](http://www.prisma-statement.org/)

This repository contains the full reproducible analytical pipeline, curated study datasets, statistical models, sensitivity analyses, risk-of-bias appraisals, and high-resolution publication figures for the systematic review and meta-analysis comparing biodegradable (and bioabsorbable) occluders with permanent metallic occluders for patent foramen ovale (PFO) closure.

---

## 📌 Clinical Overview & Key Findings

Biodegradable and bioabsorbable PFO occluders (e.g., Pansy, BioSTAR) are engineered to minimize late mechanical complications, nickel hypersensitivity, device-erosion risks, and permanent foreign-body reactions by degrading after promoting tissue endothelialization.

Across **16 eligible studies** encompassing **4,775 participants** (9 comparative studies and 7 single-arm cohorts):

| Clinical Outcome | Studies ($k$) | Sample Size ($N$) | Pooled Effect Estimate (95% CI) | Heterogeneity ($I^2$) | Evidence Certainty (GRADE) |
| :--- | :---: | :---: | :---: | :---: | :---: |
| **Effective Closure at 6 Months** | 3 | 415 | **RR 1.01** (0.98 to 1.04) | $0.0\%$ ($p = 0.52$) | Moderate |
| **Complete Closure at 12 Months** | 2 | 228 | **RR 1.00** (0.98 to 1.02) | $0.0\%$ ($p = 0.44$) | Moderate |
| **Recurrent Ischemic Stroke / TIA** | 3 | 240 | **RR 1.07** (0.44 to 2.60) | $0.0\%$ ($p = 0.90$) | Low |
| **Late Atrial Fibrillation (Post-Hospital)** | 4 | 743 | **RR 0.61** (0.26 to 1.44) | $25.3\%$ ($p = 0.26$) | Low |
| **Single-Arm Procedural Success** | 4 | 221 | **Proportion 98.4%** (95.9% to 99.7%) | $0.0\%$ | Moderate |
| **Device Embolization / Pericardial Tamponade** | 6 | 821 | **RD 0.00** (-0.01 to 0.01) | $0.0\%$ | Moderate |

**Conclusion:** Biodegradable occluders demonstrated non-inferior closure rates and safety profiles compared with conventional metallic occluders at 6 and 12 months, with zero cases of device embolization or cardiac tamponade reported in biodegradable arms.

---

## 🗂️ Repository Architecture

```text
pfo-biodegradable-meta-analysis/
├── README.md                      # Overview, clinical findings, and replication instructions
├── LICENSE                        # MIT License
├── CITATION.cff                   # Citation metadata
├── .gitignore                     # Git exclusion rules
├── config/
│   └── renderer_settings.json     # Forest plot layout, styling, and DPI export configuration
├── data/
│   ├── manifest_v5.csv            # Curated study-level extraction dataset (all 26 models)
│   ├── rob_joined.csv             # Unified 16-study Risk-of-Bias ledger (RoB 2, ROBINS-I, JBI)
│   ├── robvis_RoB2_Zhang2026.csv  # Cochrane RoB 2 appraisal domain entries
│   ├── robvis_ROBINS-I_v2_*.csv   # ROBINS-I domain ratings (permissive and strict consensus)
│   ├── JBI_Case_Series_7_*.csv    # JBI Critical Appraisal tool domain ratings
│   └── prisma_flow.csv            # PRISMA 2020 flow count register
├── R/
│   ├── run_pipeline.R             # Master replication runner (executes stages 01-09 + unit tests)
│   ├── 00_renderer_settings.R     # Figure export dimensions and typographical presets
│   ├── 02_run_models.R            # Frequentist random-effects models (REML + Knapp-Hartung)
│   ├── 03_render_figures_report.R # Publication-grade forest plots (PNG/PDF) and HTML summary
│   ├── 04_sensitivity_analyses.R  # Prespecified sensitivity analyses (S1-S6, S8)
│   ├── 05b_rob_join.R             # Joins RoB 2, ROBINS-I, and JBI appraisals
│   ├── 06_adjunct_analyses.R      # Influence diagnostics, cumulative meta-analysis, E-values
│   ├── 07_bayesian.R              # Bayesian random-effects models with half-normal/Cauchy priors
│   ├── 08_extraction_analyses.R   # Incidence rate ratios (IRRs) and MIDAS headache scores
│   ├── 09_nma.R                   # Network meta-analysis by device family (biodegradable vs metallic)
│   └── 10_grade_tsa.R             # GRADE summary evaluations and Trial Sequential Analysis
├── tests/
│   ├── test_manifest.R            # Verifies dataset integrity and exclusion compliance
│   ├── test_models.R              # Verifies numerical consistency across all 26 model estimates
│   ├── test_sensitivity.R         # Validates sensitivity analyses bounds and stability
│   ├── test_adjunct.R             # Validates cumulative meta-analysis and E-values
│   ├── test_bayesian.R            # Validates Bayesian posterior parameters
│   ├── test_nma.R                 # Validates network meta-analysis contrasts
│   ├── test_grade_tsa.R           # Validates GRADE rules and TSA boundaries
│   └── test_extraction.R          # Validates clinical IRR and MIDAS analyses
├── tools/
│   ├── generate_all_bias_plots.mjs# Headless tool to build SVG bias plots (Node.js)
│   ├── Simple Offline Robvis/     # Standalone offline Cochrane RoB 2 / ROBINS-I visualization
│   └── Simple JBI Appraisal/      # Standalone offline JBI appraisal visualization
└── outputs/
    ├── tables/                    # Computed model results (CSV format)
    ├── figures/                   # Forest plots (PNG, PDF), funnel plots, NMA, TSA, and RoB charts
    └── PFO_Meta-analysis_v5_Luna_report.html # Comprehensive interactive results report
```

---

## 🚀 Reproduction Quickstart

### Prerequisites

- **R** (version >= 4.0.0; recommended 4.4.x)
- Required R packages:
  ```r
  install.packages(c("metafor", "readr", "dplyr", "jsonlite"))
  ```
- *(Optional)* **Node.js** (>= 18.x) for running offline risk of bias plot generators.

### 1. Clone the Repository

```bash
git clone https://github.com/fernanda-alvarez/pfo-biodegradable-meta-analysis.git
cd pfo-biodegradable-meta-analysis
```

### 2. Execute Full Master Pipeline

To replicate all models, sensitivity analyses, Bayesian fits, network meta-analyses, GRADE tables, forest plots, and unit tests in a single command:

```bash
Rscript R/run_pipeline.R
```

*Expected output:*
```text
========================================================================
 Starting Meta-Analysis Replication Pipeline
========================================================================
[STAGE] 01: Risk of Bias Ledger        ... OK
[STAGE] 02: Primary & Core Models      ... OK
[STAGE] 03: Forest Plots & Report      ... OK
[STAGE] 04: Sensitivity Analyses       ... OK
[STAGE] 05: Adjunct Analyses           ... OK
[STAGE] 06: Bayesian Random-Effects    ... OK
[STAGE] 07: Clinical & Person-Time     ... OK
[STAGE] 08: Network Meta-Analysis      ... OK
[STAGE] 09: GRADE & TSA                ... OK

========================================================================
 Running Verification Test Suite
========================================================================
[TEST]  test_manifest.R                ... OK
[TEST]  test_models.R                  ... OK
[TEST]  test_sensitivity.R             ... OK
[TEST]  test_adjunct.R                 ... OK
[TEST]  test_bayesian.R                ... OK
[TEST]  test_nma.R                     ... OK
[TEST]  test_grade_tsa.R               ... OK
[TEST]  test_extraction.R              ... OK

========================================================================
 Replication Summary: 9/9 stages OK | 8/8 tests PASSED
 Outputs written to 'outputs/tables/' and 'outputs/figures/'.
========================================================================
```

### 3. Generate Risk-of-Bias Visualizations (Optional)

```bash
node tools/generate_all_bias_plots.mjs
```

---

## 📊 Analytical Methodology

1. **Frequentist Meta-Analysis:**
   - Conducted using `metafor` in R.
   - Binary outcomes: Random-effects models using restricted maximum likelihood (REML) with Knapp-Hartung adjustments.
   - Rare events (zero/near-zero counts): Mantel-Haenszel risk difference (MH RD) without continuity correction.
   - Continuous outcomes: Mean differences (MD) using REML + Knapp-Hartung.
   - Single-arm proportions: Untransformed and logit-transformed random-effects pooling.
2. **Sensitivity & Subgroup Analyses:**
   - Fixed-effects Mantel-Haenszel models.
   - Leave-one-out influence analysis.
   - Restriction to low risk-of-bias studies (excluding Serious/Critical concerns).
   - Fully biodegradable vs. bioabsorbable device generation subgroups.
   - Alternative continuity correction methods (0.5 addition).
3. **Advanced Analyses:**
   - Bayesian random-effects models evaluating posterior probabilities of benefit/harm under weakly informative priors.
   - Arm-based multivariate network meta-analysis (`rma.mv`) assessing comparative effectiveness across device classes.
   - Trial Sequential Analysis (TSA) evaluating required information sizes and boundary crossings.
   - E-values calculating minimum strength of unmeasured confounding required to negate observed effects.

---

## 📜 Ethical Disclosures & Declarations

- **Protocol Registration:** Registered on PROSPERO in accordance with PRISMA 2020 guidelines.
- **Data Availability:** All analytical data, extraction matrices, and statistical scripts are fully accessible within this repository without restriction.
- **AI-Assisted Technologies:** AI models (GPT-6 Luna, GPT-6 Sol, DeepSeek V4 Flash, and Google Antigravity Gemini Flash 3.7 & 3.8) were employed for computational coding, code refactoring, reproducibility unit testing, and typographical drafting assistance under direct human clinician oversight. All clinical syntheses, data extractions, and conclusions were independently reviewed and verified by the authors.

---

## 📄 License & Citation

Distributed under the **MIT License**. See [`LICENSE`](LICENSE) for details.

If citing this work, please refer to:
```bibtex
@article{Taslim2026PFO,
  title={Efficacy and Safety of Biodegradable versus Metallic Occluders for Patent Foramen Ovale (PFO) Closure: A Systematic Review and Meta-Analysis},
  author={Taslim, Fitranda and Contributors},
  journal={Working Paper / In Submission},
  year={2026},
  url={https://github.com/fernanda-alvarez/pfo-biodegradable-meta-analysis}
}
```
