# Shared renderer settings loader and validator.

renderer_default_settings <- function() {
  list(
    formats = c("png", "pdf"),
    width_in = 8.5,
    height_mode = "auto",
    height_in = 6.0,
    dpi = 300,
    layout_preset = "reference",
    whitespace_mode = "standard",
    proportion_axis_mode = "full",
    ecl6_axis_mode = "standard",
    background = "#FFFFFF",
    border = "#000000",
    line = "#263238",
    accent = "#000000",
    show_weights = TRUE,
    show_heterogeneity = TRUE,
    show_model = TRUE,
    show_direction_labels = TRUE,
    extreme_ci_mode = "clip",
    output_dir = "outputs/figures",
    report = TRUE
  )
}

renderer_layout_presets <- function() {
  list(
    reference = list(effect_x = 0.51, forest_left = 0.68, forest_right = 0.96,
                     font_scale = 1.00, row_scale = 1.00),
    compact = list(effect_x = 0.50, forest_left = 0.66, forest_right = 0.96,
                   font_scale = 0.94, row_scale = 0.90),
    presentation = list(effect_x = 0.50, forest_left = 0.67, forest_right = 0.96,
                        font_scale = 1.12, row_scale = 1.08)
  )
}

renderer_merge_settings <- function(defaults, overrides) {
  if (!length(overrides)) return(defaults)
  for (nm in names(overrides)) {
    if (!is.null(overrides[[nm]])) defaults[[nm]] <- overrides[[nm]]
  }
  defaults
}

renderer_is_hex <- function(x) is.character(x) && length(x) == 1L && grepl("^#[0-9A-Fa-f]{6}$", x)

load_renderer_settings <- function(proj) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("The renderer requires the jsonlite package to read config/renderer_settings.json", call. = FALSE)
  }
  settings_path <- file.path(proj, "config", "renderer_settings.json")
  overrides <- if (file.exists(settings_path)) {
    jsonlite::fromJSON(settings_path, simplifyVector = TRUE)
  } else {
    list()
  }
  settings <- renderer_merge_settings(renderer_default_settings(), overrides)
  settings$formats <- tolower(unique(as.character(settings$formats)))
  allowed_formats <- c("png", "jpg", "pdf", "svg", "tiff")
  if (!length(settings$formats) || any(!settings$formats %in% allowed_formats)) {
    stop("formats must contain only png, jpg, pdf, svg, and tiff", call. = FALSE)
  }
  if (!is.numeric(settings$width_in) || length(settings$width_in) != 1L || settings$width_in <= 0) {
    stop("width_in must be a positive number", call. = FALSE)
  }
  if (!settings$height_mode %in% c("auto", "fixed")) stop("height_mode must be auto or fixed", call. = FALSE)
  if (!is.numeric(settings$height_in) || length(settings$height_in) != 1L || settings$height_in <= 0) {
    stop("height_in must be a positive number", call. = FALSE)
  }
  if (!is.numeric(settings$dpi) || length(settings$dpi) != 1L || settings$dpi < 72 || settings$dpi > 1200) {
    stop("dpi must be between 72 and 1200", call. = FALSE)
  }
  if (!settings$layout_preset %in% names(renderer_layout_presets())) {
    stop("layout_preset must be reference, compact, or presentation", call. = FALSE)
  }
  if (!settings$whitespace_mode %in% c("standard", "compact", "ultra-compact")) {
    stop("whitespace_mode must be standard, compact, or ultra-compact", call. = FALSE)
  }
  if (!settings$proportion_axis_mode %in% c("full", "auto")) {
    stop("proportion_axis_mode must be full or auto", call. = FALSE)
  }
  if (!settings$ecl6_axis_mode %in% c("standard", "focused")) {
    stop("ecl6_axis_mode must be standard or focused", call. = FALSE)
  }
  for (nm in c("background", "border", "line", "accent")) {
    if (!renderer_is_hex(settings[[nm]])) stop(nm, " must be a six-digit hex color", call. = FALSE)
  }
  for (nm in c("show_weights", "show_heterogeneity", "show_model", "show_direction_labels", "report")) {
    settings[[nm]] <- isTRUE(settings[[nm]])
  }
  if (!settings$extreme_ci_mode %in% c("full", "clip")) stop("extreme_ci_mode must be full or clip", call. = FALSE)
  if (!is.character(settings$output_dir) || length(settings$output_dir) != 1L || !nzchar(settings$output_dir)) {
    stop("output_dir must be a non-empty path", call. = FALSE)
  }
  if (settings$report && !("png" %in% settings$formats)) {
    stop("PNG must be selected when report is TRUE because the report embeds PNG figures", call. = FALSE)
  }
  settings$layout <- renderer_layout_presets()[[settings$layout_preset]]
  settings$width_in <- as.numeric(settings$width_in)
  settings$height_in <- as.numeric(settings$height_in)
  settings$dpi <- as.integer(settings$dpi)
  settings$figdir <- if (grepl("^[A-Za-z]:[/\\\\]", settings$output_dir)) {
    normalizePath(settings$output_dir, winslash = "/", mustWork = FALSE)
  } else {
    normalizePath(file.path(proj, settings$output_dir), winslash = "/", mustWork = FALSE)
  }
  outputs_root <- normalizePath(file.path(proj, "outputs"), winslash = "/", mustWork = FALSE)
  if (!(settings$figdir == outputs_root || startsWith(settings$figdir, paste0(outputs_root, "/")))) {
    stop("output_dir must be inside the project outputs folder so report links remain portable", call. = FALSE)
  }
  settings
}
