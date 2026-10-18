# 09_nma.R — arm-based network meta-analysis by device family.
#
# Networks (general-PFO effective closure):
#   NMA-6m core   Li 2026 (Pansy vs Amplatzer), Zhang 2026 (biodegradable vs
#                 metallic), Cernic Suligoj 2020 (BioSTAR vs Amplatzer)
#   NMA-12m core  Karagianni 2011 (BioSTAR vs Amplatzer), Li 2026 (Pansy vs
#                 Amplatzer); sensitivity with Luani 2015 (BioSTAR vs metallic
#                 combined)
#   NMA-6m expanded (sensitivity, mixed ~4-6 months): + Darsaklis 2014 (BioSTAR,
#                 GSO, HELEX, APO — 4 arms), Van den Branden 2010 (BioSTAR vs
#                 CardioSEAL), Karagianni 2011 6-month discrepancy row
#
# Device families: biodegradable (Pansy, novel), bioabsorbable (BioSTAR),
# metallic (Amplatzer, CardioSEAL, GSO, HELEX, APO).
#
# Model: arm-based multivariate NMA, logit scale, rma.mv with
# random = ~ arm | study, struct = "CS", REML; fixed family effects.
# Contrasts reported as OR (95% CI) vs metallic reference and between
# biodegradable and bioabsorbable. No closed loops exist (star network), so
# consistency is not testable and the network is described accordingly.
#
# Outputs: outputs/tables/nma_results.csv, outputs/figures/nma_forest_*.png

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(metafor))

PROJ <- normalizePath(".", winslash = "/", mustWork = TRUE)
TABDIR <- file.path(PROJ, "outputs/tables")
FIGDIR <- file.path(PROJ, "outputs/figures")
dir.create(TABDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGDIR, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------ network assembly
arm <- function(study, fam, e, n) data.frame(study_id = study, family = fam, events = e, total = n,
                                             stringsAsFactors = FALSE)

net_6m_core <- rbind(
  arm("Li 2026", "biodegradable", 39, 44), arm("Li 2026", "metallic", 77, 88),
  arm("Zhang 2026", "biodegradable", 87, 96), arm("Zhang 2026", "metallic", 86, 94),
  arm("Cernic Suligoj 2020", "bioabsorbable", 43, 46), arm("Cernic Suligoj 2020", "metallic", 42, 47))

net_12m_core <- rbind(
  arm("Karagianni 2011", "bioabsorbable", 29, 30), arm("Karagianni 2011", "metallic", 23, 29),
  arm("Li 2026", "biodegradable", 42, 44), arm("Li 2026", "metallic", 72, 78))

net_12m_luani <- rbind(net_12m_core,
  arm("Luani 2015", "bioabsorbable", 55, 56), arm("Luani 2015", "metallic", 205, 215))

net_6m_expanded <- rbind(
  net_6m_core,
  arm("Darsaklis 2014", "bioabsorbable", 30, 32), arm("Darsaklis 2014", "metallic", 38, 41),
  arm("Darsaklis 2014", "metallic", 23, 31), arm("Darsaklis 2014", "metallic", 39, 51),
  arm("Van den Branden 2010a", "bioabsorbable", 36, 37), arm("Van den Branden 2010a", "metallic", 43, 44),
  arm("Karagianni 2011", "bioabsorbable", 26, 30), arm("Karagianni 2011", "metallic", 21, 29))

networks <- list(
  nma_6m_core = list(data = net_6m_core, label = "Effective closure ~6 months, general PFO (core)"),
  nma_12m_core = list(data = net_12m_core, label = "Effective closure 12 months, general PFO (core)"),
  nma_12m_with_luani = list(data = net_12m_luani, label = "Effective closure 12 months with conditional Luani (sensitivity)"),
  nma_6m_expanded = list(data = net_6m_expanded, label = "Effective closure ~4-6 months expanded (sensitivity, mixed windows)")
)

fit_network <- function(dat) {
  dat$arm_id <- ave(seq_len(nrow(dat)), dat$study_id, FUN = seq_along)
  dat$family <- relevel(factor(dat$family), ref = "metallic")
  es <- escalc(measure = "PLO", xi = dat$events, ni = dat$total, data = dat)
  fit <- tryCatch(
    rma.mv(yi, vi, data = es, mods = ~ family, random = ~ arm_id | study_id,
           struct = "CS", method = "REML"),
    error = function(e) NULL)
  if (is.null(fit)) {
    # fallback: marginal model ignoring within-study arm correlation
    fit <- rma(yi, vi, data = es, mods = ~ family, method = "REML")
    note <- "rma.mv (CS) did not converge; marginal rma used (arm correlation ignored)"
  } else {
    note <- "arm-based multivariate NMA, rma.mv with CS arm correlation"
  }
  b <- coef(fit); V <- vcov(fit)
  fams <- levels(factor(dat$family)); fams <- fams[fams != "metallic"]
  mk <- function(nm) {
    if (!nm %in% names(b)) return(c(est = NA, lo = NA, hi = NA))
    est <- b[[nm]]; se <- sqrt(V[nm, nm])
    c(est = est, lo = est - 1.96 * se, hi = est + 1.96 * se)
  }
  rows <- lapply(fams, function(f) {
    v <- mk(paste0("family", f))
    data.frame(family = f,
               or_vs_metallic = round(exp(v["est"]), 3),
               ci_low = round(exp(v["lo"]), 3),
               ci_high = round(exp(v["hi"]), 3),
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  # biodegradable vs bioabsorbable contrast (only when both are in the network)
  if (all(c("biodegradable", "bioabsorbable") %in% fams)) {
    nb <- paste0("familybiodegradable"); na <- paste0("familybioabsorbable")
    d <- b[[nb]] - b[[na]]
    sd <- sqrt(V[nb, nb] + V[na, na] - 2 * V[nb, na])
    out <- rbind(out, data.frame(family = "biodegradable_vs_bioabsorbable",
                                 or_vs_metallic = round(exp(d), 3),
                                 ci_low = round(exp(d - 1.96 * sd), 3),
                                 ci_high = round(exp(d + 1.96 * sd), 3),
                                 stringsAsFactors = FALSE))
  }
  list(fit = fit, out = out, note = note, tau2 = fit$tau2)
}

# ------------------------------------------------------- fit all networks ----
NMA <- list()
for (nm in names(networks)) {
  res <- fit_network(networks[[nm]]$data)
  NMA[[nm]] <- res
  # forest vs metallic reference
  png(file.path(FIGDIR, paste0("nma_forest_", nm, ".png")), width = 1000, height = 560, res = 130)
  par(mar = c(4, 6, 3, 1))
  forest(res$out$or_vs_metallic, ci.lb = res$out$ci_low, ci.ub = res$out$ci_high,
         slab = res$out$family, refline = 1, xlab = "OR vs metallic reference (95% CI)",
         main = networks[[nm]]$label, xlim = c(-1, 4))
  dev.off()
}

# ----------------------------------------------------------- contrast assembly
contrast_rows <- list()
for (nm in names(networks)) {
  res <- NMA[[nm]]
  fams <- res$out$family
  for (i in seq_along(fams)) {
    contrast_rows[[paste0(nm, "_", fams[i])]] <- data.frame(
      model_id = nm, network = networks[[nm]]$label,
      contrast = ifelse(fams[i] == "biodegradable_vs_bioabsorbable",
                        "biodegradable vs bioabsorbable",
                        paste(fams[i], "vs metallic")),
      or = res$out$or_vs_metallic[i], ci_low = res$out$ci_low[i], ci_high = res$out$ci_high[i],
      tau2 = round(res$tau2, 5),
      k = length(unique(networks[[nm]]$data$study_id)),
      note = res$note, stringsAsFactors = FALSE)
  }
}
nma_df <- do.call(rbind, contrast_rows)
row.names(nma_df) <- NULL
write.csv(nma_df, file.path(TABDIR, "nma_results.csv"), row.names = FALSE, na = "")
writeLines(c(
  paste(Sys.time(), "NMA networks", length(NMA)),
  paste(Sys.time(), "contrasts", nrow(nma_df)),
  paste(Sys.time(), "done")
), file.path(PROJ, "outputs/logs/09_nma.log"))
cat("NMA DONE:", length(NMA), "networks;", nrow(nma_df), "contrasts\n")
