# Render reference-style forest plots and an image-linked HTML report.

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(metafor)
  library(jsonlite)
})

PROJ <- normalizePath(".", winslash = "/", mustWork = TRUE)
source(file.path(PROJ, "R", "00_renderer_settings.R"))
SETTINGS <- load_renderer_settings(PROJ)
MANIFEST <- read.csv(file.path(PROJ, "data_derived/luna_manifest.csv"), check.names = FALSE,
                     na.strings = c("", "NA", "NR"))
RESULTS <- read.csv(file.path(PROJ, "outputs/tables/model_registry.csv"), check.names = FALSE,
                    na.strings = c("", "NA", "NR"))
FIGDIR <- SETTINGS$figdir
dir.create(FIGDIR, recursive = TRUE, showWarnings = FALSE)
LOGDIR <- file.path(PROJ, "outputs/logs")
dir.create(LOGDIR, recursive = TRUE, showWarnings = FALSE)
warning_lines <- character()
format_label <- paste(toupper(SETTINGS$formats), collapse = ", ")

fmt <- function(x, digits = 2) ifelse(is.na(x), "NA", formatC(x, format = "f", digits = digits))
fmt_effect_value <- function(x) {
  if (!is.finite(x)) return("NA")
  if (abs(x) >= 100 || (abs(x) > 0 && abs(x) < 0.001)) {
    return(formatC(x, format = "e", digits = 2))
  }
  fmt(x)
}
effect_text <- function(measure, est, lo, hi) {
  if (measure == "RR") return(sprintf("RR %s (%s, %s)", fmt_effect_value(est), fmt_effect_value(lo), fmt_effect_value(hi)))
  if (measure == "MD") return(sprintf("MD %s (%s, %s)", fmt_effect_value(est), fmt_effect_value(lo), fmt_effect_value(hi)))
  if (measure == "RD") return(sprintf("RD %s (%s, %s)", fmt_effect_value(est), fmt_effect_value(lo), fmt_effect_value(hi)))
  sprintf("Proportion %s (%s, %s)", fmt_effect_value(est), fmt_effect_value(lo), fmt_effect_value(hi))
}

study_effects <- function(dat) {
  measure <- dat$measure[1]
  if (measure == "proportion") {
    out <- lapply(seq_len(nrow(dat)), function(i) {
      bt <- binom.test(dat$ei[i], dat$ni[i])
      data.frame(est = dat$ei[i] / dat$ni[i], lo = bt$conf.int[1], hi = bt$conf.int[2], vi = NA_real_)
    })
    return(do.call(rbind, out))
  }
  if (measure == "MD") {
    es <- escalc(measure = "MD", m1i = dat$m1i, sd1i = dat$sd1i, n1i = dat$n1i,
                 m2i = dat$m2i, sd2i = dat$sd2i, n2i = dat$n2i)
    se <- sqrt(es$vi)
    return(data.frame(est = es$yi, lo = es$yi - 1.96 * se, hi = es$yi + 1.96 * se, vi = es$vi))
  }
  if (measure == "RR") {
    es <- escalc(measure = "RR", ai = dat$ei, bi = dat$ni - dat$ei,
                 ci = dat$ec, di = dat$nc - dat$ec, add = 0.5, to = "only0")
    se <- sqrt(es$vi)
    return(data.frame(est = exp(es$yi), lo = exp(es$yi - 1.96 * se), hi = exp(es$yi + 1.96 * se), vi = es$vi))
  }
  p1 <- dat$ei / dat$ni; p2 <- dat$ec / dat$nc
  se <- sqrt(p1 * (1 - p1) / dat$ni + p2 * (1 - p2) / dat$nc)
  data.frame(est = p1 - p2, lo = p1 - p2 - 1.96 * se, hi = p1 - p2 + 1.96 * se,
             vi = se^2)
}

pooled_effect <- function(dat) {
  measure <- dat$measure[1]
  if (measure == "proportion") {
    es <- escalc(measure = "PLO", xi = dat$ei, ni = dat$ni, add = 0.5, to = "only0")
    fit <- rma(yi, vi, data = es, method = "REML", test = "knha")
    pr <- predict(fit, transf = plogis)
    return(list(est = pr$pred, lo = pr$ci.lb, hi = pr$ci.ub, fit = fit))
  }
  if (measure == "MD") {
    es <- escalc(measure = "MD", m1i = dat$m1i, sd1i = dat$sd1i, n1i = dat$n1i,
                 m2i = dat$m2i, sd2i = dat$sd2i, n2i = dat$n2i)
    fit <- rma(yi, vi, data = es, method = "REML", test = "knha")
    pr <- predict(fit)
    return(list(est = pr$pred, lo = pr$ci.lb, hi = pr$ci.ub, fit = fit))
  }
  if (measure == "RR") {
    es <- escalc(measure = "RR", ai = dat$ei, bi = dat$ni - dat$ei,
                 ci = dat$ec, di = dat$nc - dat$ec, add = 0.5, to = "only0")
    fit <- rma(yi, vi, data = es, method = "REML", test = "knha")
    pr <- predict(fit, transf = exp)
    return(list(est = pr$pred, lo = pr$ci.lb, hi = pr$ci.ub, fit = fit))
  }
  fit <- suppressWarnings(rma.mh(ai = dat$ei, bi = dat$ni - dat$ei,
                                  ci = dat$ec, di = dat$nc - dat$ec,
                                  measure = "RD", add = 0, to = "none", correct = FALSE))
  pr <- predict(fit)
  list(est = pr$pred, lo = pr$ci.lb, hi = pr$ci.ub, fit = fit)
}

draw_one <- function(model_id, title, favour = NULL) {
  dat <- MANIFEST[MANIFEST$analysis_id == model_id, , drop = FALSE]
  if (!nrow(dat)) return(NULL)
  pooled <- pooled_effect(dat)
  se <- study_effects(dat)
  measure <- dat$measure[1]
  n_study <- nrow(dat)
  groups <- if (model_id == "ecl6_core") {
    list("1.1. Randomized studies" = which(dat$study_id == "Zhang 2026"),
         "1.2. Nonrandomized studies" = which(dat$study_id != "Zhang 2026"))
  } else if (measure == "proportion") {
    list("Supportive single-arm cohorts" = seq_len(nrow(dat)))
  } else {
    list("Comparative studies" = seq_len(nrow(dat)))
  }
  ultra_space <- SETTINGS$whitespace_mode == "ultra-compact"
  compact_space <- SETTINGS$whitespace_mode != "standard"
  row_step <- if (ultra_space) 0.40 else if (compact_space) 0.44 else 0.55
  group_gap <- if (ultra_space) 0.02 else if (compact_space) 0.05 else 0.10
  detail_step <- if (ultra_space) 0.39 else if (compact_space) 0.43 else 0.55
  model_step <- if (ultra_space) 0.35 else if (compact_space) 0.39 else 0.48
  study_start_y <- n_study + if (ultra_space) 5.87 else if (compact_space) 5.80 else 5.25
  header_shift <- if (compact_space) 0.35 else 0
  y_after_studies <- study_start_y
  for (idx in groups) {
    y_after_studies <- y_after_studies - row_step * (1 + length(idx))
    if (length(idx) > 1) y_after_studies <- y_after_studies - group_gap
  }
  total_y <- y_after_studies - 0.05
  summary_y <- total_y - detail_step - model_step
  axis_y <- summary_y - if (ultra_space) 0.20 else if (compact_space) 0.27 else 0.55
  has_direction_labels <- SETTINGS$show_direction_labels && !is.null(favour) && measure == "RR"
  favour_y <- axis_y - if (ultra_space) 1.10 else if (compact_space) 1.18 else 1.05
  bottom_y <- if (compact_space) {
    if (has_direction_labels) favour_y - if (ultra_space) 0.34 else 0.40 else axis_y - if (ultra_space) 0.72 else 0.78
  } else favour_y - 0.55
  top_y <- n_study + 7.9
  height_in <- if (SETTINGS$height_mode == "fixed") SETTINGS$height_in else if (ultra_space) {
    max(2.95, 2.15 + 0.22 * n_study + if (has_direction_labels) 0.33 else 0)
  } else if (compact_space) {
    max(3.05, 2.30 + 0.24 * n_study + if (has_direction_labels) 0.35 else 0)
  } else max(4.8, 3.8 + 0.38 * n_study)
  output_paths <- setNames(file.path(FIGDIR, paste0("forest_", model_id, ".", SETTINGS$formats)), SETTINGS$formats)

  draw_device <- function() {
    # Keep only a small device-edge inset so the frame stroke is never clipped.
    par(mai = rep(0.02, 4), xaxs = "i", yaxs = "i", xpd = NA)
    plot.new(); plot.window(xlim = c(0, 1), ylim = c(bottom_y, top_y), xaxs = "i", yaxs = "i")
    rect(0, bottom_y, 1, top_y, col = SETTINGS$background, border = SETTINGS$border, lwd = 1.2)
    fs <- SETTINGS$layout$font_scale
    effect_x <- SETTINGS$layout$effect_x
    forest_left <- SETTINGS$layout$forest_left
    forest_right <- SETTINGS$layout$forest_right
    text(0.5, n_study + 7.45, title, font = 2, cex = 1.15 * fs)
    text(0.02, n_study + 6.55 + header_shift, "Study or subgroup", adj = 0, font = 2, cex = 0.78 * fs)
    text(0.282, n_study + 6.55 + header_shift, ifelse(measure == "proportion", "Single arm", "Bio-PFO"), font = 2, cex = 0.64 * fs)
    text(0.282, n_study + 6.15 + header_shift, ifelse(measure == "MD", "Mean    n", "Events   Total"), font = 2, cex = 0.66 * fs)
    text(0.388, n_study + 6.55 + header_shift, ifelse(measure == "proportion", "", "Metallic"), font = 2, cex = 0.64 * fs)
    text(0.388, n_study + 6.15 + header_shift, ifelse(measure == "MD", "Mean    n", ifelse(measure == "proportion", "", "Events   Total")), font = 2, cex = 0.66 * fs)
    text(0.455, n_study + 6.35 + header_shift, ifelse(SETTINGS$show_weights, "Weight", ""), font = 2, cex = 0.68 * fs)
    text(effect_x, n_study + 6.55 + header_shift, ifelse(measure == "proportion", "Proportion", paste(measure, "[95% CI]")),
         adj = 0, font = 2, cex = 0.72 * fs)
    text(effect_x, n_study + 6.15 + header_shift, "Estimate [95% CI]", adj = 0, font = 2, cex = 0.68 * fs)
    segments(0.01, n_study + 5.85 + header_shift, 0.99, n_study + 5.85 + header_shift,
             col = SETTINGS$line, lwd = 1)

    xleft <- forest_left; xright <- forest_right
    if (measure == "RR") {
      focused_ecl6 <- model_id == "ecl6_core" && SETTINGS$ecl6_axis_mode == "focused"
      if (focused_ecl6) {
        xlo <- floor((min(se$lo, pooled$lo, 1, na.rm = TRUE) - 0.02) / 0.05) * 0.05
        xhi <- ceiling((max(se$hi, pooled$hi, 1, na.rm = TRUE) + 0.02) / 0.05) * 0.05
        ticks <- sort(unique(c(xlo, 0.9, 1, 1.1, xhi)))
        ticks <- ticks[ticks >= xlo - 1e-9 & ticks <= xhi + 1e-9]
      } else {
        xlo <- if (SETTINGS$extreme_ci_mode == "clip") 0.2 else min(0.2, min(c(se$lo[is.finite(se$lo)], pooled$lo), na.rm = TRUE) * 0.75)
        xhi <- if (SETTINGS$extreme_ci_mode == "clip") 5 else max(5, max(se$hi[is.finite(se$hi)], pooled$hi, na.rm = TRUE) * 1.15)
        ticks <- c(0.2, 0.5, 1, 2, 5); ticks <- ticks[ticks >= xlo & ticks <= xhi]
      }
      raw_map <- function(x) xleft + (log(pmax(x, .Machine$double.xmin)) - log(xlo)) / (log(xhi) - log(xlo)) * (xright - xleft)
      null <- NA
    } else if (measure == "MD") {
      xlo <- if (SETTINGS$extreme_ci_mode == "clip") -5 else min(-1, min(se$lo, pooled$lo, na.rm = TRUE) * 1.15)
      xhi <- if (SETTINGS$extreme_ci_mode == "clip") 5 else max(1, max(se$hi, pooled$hi, na.rm = TRUE) * 1.15)
      raw_map <- function(x) xleft + (x - xlo) / (xhi - xlo) * (xright - xleft)
      null <- NA; ticks <- pretty(c(xlo, xhi), n = 5)
    } else if (measure == "RD") {
      xlo <- if (SETTINGS$extreme_ci_mode == "clip") -0.05 else min(-0.05, min(se$lo, pooled$lo, na.rm = TRUE) * 1.15)
      xhi <- if (SETTINGS$extreme_ci_mode == "clip") 0.05 else max(0.05, max(se$hi, pooled$hi, na.rm = TRUE) * 1.15)
      raw_map <- function(x) xleft + (x - xlo) / (xhi - xlo) * (xright - xleft)
      null <- NA; ticks <- pretty(c(xlo, xhi), n = 5)
    } else {
      if (SETTINGS$proportion_axis_mode == "auto") {
        visible_lo <- min(c(se$lo, pooled$lo), na.rm = TRUE)
        visible_hi <- max(c(se$hi, pooled$hi), na.rm = TRUE)
        xlo <- max(0, floor((visible_lo - 0.02) / 0.05) * 0.05)
        xhi <- min(1, ceiling((visible_hi + 0.02) / 0.05) * 0.05)
        if (xhi - xlo < 0.15) xlo <- max(0, xhi - 0.15)
        ticks <- pretty(c(xlo, xhi), n = 5)
        ticks <- ticks[ticks >= xlo - 1e-9 & ticks <= xhi + 1e-9]
      } else {
        xlo <- 0; xhi <- 1; ticks <- c(0, .25, .5, .75, 1)
      }
      map <- function(x) xleft + (x - xlo) / (xhi - xlo) * (xright - xleft)
      null <- NA
    }
    if (measure != "proportion") {
      map <- if (SETTINGS$extreme_ci_mode == "clip") {
        function(x) pmin(xright, pmax(xleft, raw_map(x)))
      } else raw_map
      null <- if (measure == "RR") map(1) else map(0)
    }
    if (is.finite(null)) {
      segments(null, axis_y, null, n_study + 5.85 + header_shift, col = SETTINGS$line, lwd = 1)
    }

    y <- study_start_y
    weights <- ifelse(is.finite(se$vi) & se$vi > 0, 1 / se$vi, NA_real_)
    weights <- if (all(is.na(weights))) rep(100 / n_study, n_study) else 100 * weights / sum(weights, na.rm = TRUE)
    draw_effect <- function(label, y_value, font = 1, base_cex = 0.60 * fs) {
      available <- forest_left - effect_x - 0.012
      label_width <- strwidth(label, cex = base_cex, units = "user")
      label_cex <- if (is.finite(label_width) && label_width > available) base_cex * available / label_width else base_cex
      text(effect_x, y_value, label, adj = 0, cex = label_cex, font = font)
    }
    for (g in names(groups)) {
      text(0.02, y, g, adj = 0, font = 2, cex = 0.78 * fs); y <- y - row_step
      idx <- groups[[g]]
      for (i in idx) {
        text(0.02, y, dat$study_id[i], adj = 0, cex = 0.73 * fs)
        if (measure == "MD") {
          text(0.255, y, ifelse(is.na(dat$m1i[i]), "-", fmt(dat$m1i[i])), cex = 0.70 * fs)
          text(0.31, y, ifelse(is.na(dat$n1i[i]), "-", dat$n1i[i]), cex = 0.70 * fs)
          text(0.36, y, ifelse(is.na(dat$m2i[i]), "-", fmt(dat$m2i[i])), cex = 0.70 * fs)
          text(0.415, y, ifelse(is.na(dat$n2i[i]), "-", dat$n2i[i]), cex = 0.70 * fs)
        } else {
          text(0.255, y, ifelse(is.na(dat$ei[i]), "-", dat$ei[i]), cex = 0.73 * fs)
          text(0.31, y, ifelse(is.na(dat$ni[i]), "-", dat$ni[i]), cex = 0.73 * fs)
          text(0.36, y, ifelse(is.na(dat$ec[i]), "-", dat$ec[i]), cex = 0.73 * fs)
          text(0.415, y, ifelse(is.na(dat$nc[i]), "-", dat$nc[i]), cex = 0.73 * fs)
        }
        if (SETTINGS$show_weights) {
          text(0.465, y, ifelse(is.na(weights[i]), "-", sprintf("%.1f%%", weights[i])), cex = 0.68 * fs)
        }
        draw_effect(effect_text(measure, se$est[i], se$lo[i], se$hi[i]), y)
        if (is.finite(se$lo[i]) && is.finite(se$hi[i])) {
          segments(map(se$lo[i]), y, map(se$hi[i]), y, lwd = 1.1, col = SETTINGS$line)
          segments(map(se$lo[i]), y - 0.10, map(se$lo[i]), y + 0.10, col = SETTINGS$line)
          segments(map(se$hi[i]), y - 0.10, map(se$hi[i]), y + 0.10, col = SETTINGS$line)
          if (SETTINGS$extreme_ci_mode == "clip" && se$lo[i] < xlo) {
            polygon(c(xleft, xleft + 0.012, xleft + 0.012), c(y, y + 0.08, y - 0.08), col = SETTINGS$line, border = SETTINGS$line)
          }
          if (SETTINGS$extreme_ci_mode == "clip" && se$hi[i] > xhi) {
            polygon(c(xright, xright - 0.012, xright - 0.012), c(y, y + 0.08, y - 0.08), col = SETTINGS$line, border = SETTINGS$line)
          }
        }
        points(map(se$est[i]), y, pch = 15, col = SETTINGS$accent, cex = 0.65 + 0.04 * sqrt(max(weights[i], 1)))
        y <- y - row_step
      }
      if (length(idx) > 1) y <- y - group_gap
    }

    y <- y - 0.05
    separator_y <- y + if (ultra_space) 0.27 else if (compact_space) 0.29 else 0.35
    segments(0.01, separator_y, 0.99, separator_y,
             col = SETTINGS$line, lwd = 1)
    text(0.02, y, sprintf("Total (%s)", measure), adj = 0, font = 2, cex = 0.75 * fs)
    if (measure == "MD") {
      text(0.255, y, "-", cex = 0.72 * fs); text(0.31, y, "-", cex = 0.72 * fs)
      text(0.36, y, "-", cex = 0.72 * fs); text(0.415, y, "-", cex = 0.72 * fs)
    } else if (measure == "proportion") {
      text(0.255, y, sum(dat$ei, na.rm = TRUE), cex = 0.72 * fs)
      text(0.31, y, sum(dat$ni, na.rm = TRUE), cex = 0.72 * fs)
      text(0.36, y, "-", cex = 0.72 * fs); text(0.415, y, "-", cex = 0.72 * fs)
    } else {
      text(0.255, y, sum(dat$ei, na.rm = TRUE), cex = 0.72 * fs)
      text(0.31, y, sum(dat$ni, na.rm = TRUE), cex = 0.72 * fs)
      text(0.36, y, sum(dat$ec, na.rm = TRUE), cex = 0.72 * fs)
      text(0.415, y, sum(dat$nc, na.rm = TRUE), cex = 0.72 * fs)
    }
    draw_effect(effect_text(measure, pooled$est, pooled$lo, pooled$hi), y, font = 2, base_cex = 0.62 * fs)
    polygon(c(map(pooled$lo), map(pooled$est), map(pooled$hi), map(pooled$est)),
            c(y, y + 0.20, y, y - 0.20), col = SETTINGS$line, border = SETTINGS$line)
    if (SETTINGS$extreme_ci_mode == "clip" && pooled$lo < xlo) {
      polygon(c(xleft, xleft + 0.012, xleft + 0.012), c(y, y + 0.10, y - 0.10), col = SETTINGS$line, border = SETTINGS$line)
    }
    if (SETTINGS$extreme_ci_mode == "clip" && pooled$hi > xhi) {
      polygon(c(xright, xright - 0.012, xright - 0.012), c(y, y + 0.10, y - 0.10), col = SETTINGS$line, border = SETTINGS$line)
    }
    y <- y - detail_step
    if (SETTINGS$show_heterogeneity && !is.null(pooled$fit) && !is.null(pooled$fit$QE)) {
      text(0.02, y, sprintf("Heterogeneity: Q = %.2f, df = %d, p = %.3f, I2 = %s%%",
                            pooled$fit$QE, pooled$fit$k - 1, pooled$fit$QEp,
                            ifelse(is.na(pooled$fit$I2), "NA", fmt(pooled$fit$I2, 0))), adj = 0, cex = 0.68 * fs)
    }
    y <- y - model_step
    if (SETTINGS$show_model) {
      text(0.02, y, sprintf("Model: %s; k = %d; N = %d", RESULTS$method[RESULTS$model_id == model_id][1],
                            nrow(dat), if (measure == "proportion") sum(dat$ni, na.rm = TRUE) else if (measure == "MD") sum(dat$n1i, dat$n2i, na.rm = TRUE) else sum(dat$ni, dat$nc, na.rm = TRUE)),
           adj = 0, cex = 0.68 * fs)
    }
    tick_digits <- if (measure == "proportion") {
      if (length(ticks) > 1 && min(diff(ticks)) < 0.095) 2 else 1
    } else if (measure == "RR" && model_id == "ecl6_core" && SETTINGS$ecl6_axis_mode == "focused") 2 else if (max(abs(ticks)) < 1) 2 else 1
    axis(1, at = map(ticks), labels = if (compact_space) FALSE else fmt(ticks, tick_digits),
         pos = axis_y, cex.axis = 0.68 * fs, tck = -0.015)
    if (compact_space) {
      text(map(ticks), axis_y - if (ultra_space) 0.39 else 0.42,
           labels = fmt(ticks, tick_digits), cex = 0.68 * fs)
    }
    if (SETTINGS$show_direction_labels && !is.null(favour) && measure == "RR") {
      favour_label <- function(label) {
        if (grepl("biodegradable", label, ignore.case = TRUE)) return("Favors Bio-PFO")
        if (grepl("permanent", label, ignore.case = TRUE)) return("Favors metallic")
        label
      }
      text(xleft, favour_y, favour_label(favour[1]), adj = 0, cex = 0.60 * fs)
      text(xright, favour_y, favour_label(favour[2]), adj = 1, cex = 0.60 * fs)
    }
  }

  open_device <- function(format, path) {
    if (format == "png") return(png(path, width = round(SETTINGS$width_in * SETTINGS$dpi), height = round(height_in * SETTINGS$dpi), res = SETTINGS$dpi, bg = SETTINGS$background))
    if (format == "jpg") return(jpeg(path, width = round(SETTINGS$width_in * SETTINGS$dpi), height = round(height_in * SETTINGS$dpi), res = SETTINGS$dpi, quality = 95, bg = SETTINGS$background))
    if (format == "pdf") return(pdf(path, width = SETTINGS$width_in, height = height_in, onefile = FALSE, useDingbats = FALSE, bg = SETTINGS$background))
    if (format == "svg") return(svg(path, width = SETTINGS$width_in, height = height_in, bg = SETTINGS$background))
    if (format == "tiff") return(tiff(path, width = round(SETTINGS$width_in * SETTINGS$dpi), height = round(height_in * SETTINGS$dpi), res = SETTINGS$dpi, compression = "lzw", bg = SETTINGS$background))
    stop("Unsupported output format: ", format, call. = FALSE)
  }
  for (format in SETTINGS$formats) {
    open_device(format, output_paths[[format]])
    draw_device(); dev.off()
  }
  if (SETTINGS$extreme_ci_mode == "full" && ((measure == "RR" && pooled$hi > 100) || (measure == "proportion" && pooled$hi - pooled$lo > 0.5))) {
    warning_lines <<- c(warning_lines, sprintf("%s: extremely wide pooled confidence interval; consider clip mode.", model_id))
  }
  if (nrow(dat) < 3) warning_lines <<- c(warning_lines, sprintf("%s: k=%d; inference should be interpreted cautiously.", model_id, nrow(dat)))
  output_paths
}

plot_specs <- list(
  ecl6_core = c("Effective closure at approximately 6 months, general PFO", "Favors biodegradable", "Favors permanent metallic"),
  ecl12_gen_core = c("Effective closure at 12 months, general PFO", "Favors biodegradable", "Favors permanent metallic"),
  ecl12_mig = c("Effective closure at 12 months, migraine PFO subgroup", "Favors biodegradable", "Favors permanent metallic"),
  ecl24 = c("Long-term effective closure, 24 months or later", "Favors biodegradable", "Favors permanent metallic"),
  strict_early = c("Strict complete closure, early follow-up", "Favors biodegradable", "Favors permanent metallic"),
  cve_late_core = c("Recurrent stroke/TIA/CVE at longest follow-up", "Favors biodegradable", "Favors permanent metallic"),
  cve_late_with_luani = c("Recurrent stroke/TIA/CVE with conditional Luani sensitivity", "Favors biodegradable", "Favors permanent metallic"),
  midas_12m = c("MIDAS score at 12 months, migraine PFO", "Lower score", "Higher score"),
  af_late_core = c("New atrial fibrillation, AF-only family", "Favors biodegradable", "Favors permanent metallic"),
  af_late_with_luani = c("New atrial fibrillation with conditional Luani sensitivity", "Favors biodegradable", "Favors permanent metallic"),
  drt_any = c("Device-related thrombosis, corrected arm-level set", "Favors biodegradable", "Favors permanent metallic"),
  embolization = c("Device embolization/detachment, corrected arm-level set", "Favors biodegradable", "Favors permanent metallic"),
  pericardial = c("Pericardial tamponade, corrected arm-level set", "Favors biodegradable", "Favors permanent metallic"),
  sa_proc_success = c("Supportive single-arm procedural success", NULL, NULL),
  sa_closure_12m = c("Supportive single-arm closure at 12 months", NULL, NULL),
  sa_drt = c("Supportive single-arm device-related thrombus", NULL, NULL),
  close12_mig_strict = c("Strict complete closure at 12 months, migraine PFO subgroup", "Favors biodegradable", "Favors permanent metallic"),
  shunt_g6_general = c("Moderate-or-large residual shunt (grade 2-3) at 6 months, general PFO", "Favors biodegradable", "Favors permanent metallic"),
  af_late_full = c("New atrial fibrillation, AF-only family (full set incl. Zhang 2026)", "Favors biodegradable", "Favors permanent metallic"),
  af_late_full_with_luani = c("New atrial fibrillation, full set with conditional Luani sensitivity", "Favors biodegradable", "Favors permanent metallic"),
  af_any_6m = c("Any atrial arrhythmia at ~6 months", "Favors biodegradable", "Favors permanent metallic"),
  af_any_late = c("Any atrial arrhythmia, late follow-up", "Favors biodegradable", "Favors permanent metallic"),
  af_any_late_with_darsaklis = c("Any atrial arrhythmia, late follow-up with conditional Darsaklis sensitivity", "Favors biodegradable", "Favors permanent metallic"),
  ecl6_expanded_early = c("Expanded early effective closure (mixed ~4-6 months, sensitivity)", "Favors biodegradable", "Favors permanent metallic"),
  ecl12_gen_with_luani = c("Effective closure at 12 months with conditional Luani sensitivity", "Favors biodegradable", "Favors permanent metallic"),
  bleeding = c("Bleeding/vascular complications, corrected arm-level set", "Favors biodegradable", "Favors permanent metallic")
)

for (id in names(plot_specs)) draw_one(id, plot_specs[[id]][1], plot_specs[[id]][2:3])
if (!length(warning_lines)) warning_lines <- "No renderer warnings."
writeLines(unique(warning_lines), file.path(LOGDIR, "renderer_warnings.log"))

outputs_root <- normalizePath(file.path(PROJ, "outputs"), winslash = "/", mustWork = FALSE)
figdir_abs <- normalizePath(FIGDIR, winslash = "/", mustWork = FALSE)
figure_href_dir <- if (figdir_abs == outputs_root) "." else substring(figdir_abs, nchar(outputs_root) + 2L)
figure_href_dir <- gsub("\\\\", "/", figure_href_dir)

result_lines <- vapply(seq_len(nrow(RESULTS)), function(i) {
  sprintf("| %s | %s | %d | %s | %s | %s | %s |",
          RESULTS$model_id[i], RESULTS$measure[i], RESULTS$k[i], RESULTS$n_participants[i],
          fmt(RESULTS$estimate[i]), fmt(RESULTS$ci_low[i]), fmt(RESULTS$ci_high[i]))
}, character(1))
figure_lines <- paste0("![", names(plot_specs), "](", figure_href_dir, "/forest_", names(plot_specs), ".png)")
md <- c(
  "# PFO Meta-analysis v5 (Luna)", "",
  "**Analysis status:** Rebuilt in an isolated v5 package; v4 was not overwritten.", "",
  sprintf("**Design:** Publication-style grouped forest plots modeled on the supplied blue-gray reference figure. Export formats: %s; PNGs are linked below and embedded in the HTML report.", format_label), "",
  "## Methods", "",
  "Comparative models use risk ratios for binary efficacy/safety outcomes and mean differences for MIDAS, with REML random effects and Hartung-Knapp confidence intervals. Rare-event primary models use Mantel-Haenszel risk differences without a universal continuity correction. Supportive single-arm analyses use a logit-normal REML model and are not comparative treatment effects.", "",
  "Rows are selected through `data_derived/luna_manifest.csv`. Narrative-only rows, the countless Zhang duplicate, and Darsaklis narrative rare-event rows are excluded from quantitative models. Derived and conditional rows are labeled and isolated in the manifest.", "",
  "## Results", "",
  "| Model | Measure | k | N | Estimate | CI low | CI high |", "|---|---:|---:|---:|---:|---:|---:|",
  result_lines, "", "## Figures", "", figure_lines
)
writeLines(enc2utf8(md), file.path(PROJ, "outputs/PFO_Meta-analysis_v5_Luna_report.md"), useBytes = TRUE)

html_rows <- paste0("<tr><td>", RESULTS$model_id, "</td><td>", RESULTS$measure, "</td><td>", RESULTS$k,
                    "</td><td>", RESULTS$n_participants, "</td><td>", fmt(RESULTS$estimate),
                    "</td><td>", fmt(RESULTS$ci_low), "</td><td>", fmt(RESULTS$ci_high), "</td></tr>", collapse = "\n")
html_figs <- paste0("<h3>", names(plot_specs), "</h3><img src=\"", figure_href_dir, "/forest_", names(plot_specs),
                    ".png\" alt=\"", names(plot_specs), " forest plot\">", collapse = "\n")
html <- c(
  "<!doctype html>", "<html><head><meta charset=\"utf-8\"><title>PFO Meta-analysis v5 (Luna)</title>",
  "<style>body{font-family:Arial,sans-serif;background:#ffffff;color:#1f2937;max-width:1200px;margin:24px auto;padding:0 18px}h1,h2{color:#17365d}table{border-collapse:collapse;width:100%;background:white}th,td{border:1px solid #b9c4d5;padding:6px;text-align:left}th{background:#ffffff}img{display:block;max-width:100%;height:auto;background:#ffffff;margin:8px 0 26px;border:1px solid #000000}</style></head><body>",
  "<h1>PFO Meta-analysis v5 (Luna)</h1><p><b>Status:</b> Rebuilt in an isolated v5 package; v4 was not overwritten.</p>",
  "<h2>Methods</h2><p>Comparative models use RR for binary outcomes and MD for MIDAS, with REML random effects and Hartung-Knapp confidence intervals. Rare-event primary models use Mantel-Haenszel RD without a universal continuity correction. Single-arm models are supportive absolute proportions only.</p>",
  "<h2>Results</h2><table><thead><tr><th>Model</th><th>Measure</th><th>k</th><th>N</th><th>Estimate</th><th>CI low</th><th>CI high</th></tr></thead><tbody>", html_rows, "</tbody></table>",
  "<h2>Figures</h2>", html_figs, "</body></html>"
)
writeLines(enc2utf8(html), file.path(PROJ, "outputs/PFO_Meta-analysis_v5_Luna_report.html"), useBytes = TRUE)
cat("FIGURES AND REPORT WRITTEN\n")
