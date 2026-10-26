# Unified wildfire DNA-qSIP workflow.
# Validated by static R grammar and dependency checks; no full-data R execution was available.
# Edit INPUT_DIR, FILES, and temperature settings before running.
# Run: source("Wildfire_qSIP_Unified.R", encoding = "UTF-8")
# Input abundance: bacterial 16S gene copies per g soil in EACH density fraction.
# Design: UB/F3M/F3Y x two temperatures x label/unlabel x three independent plots.
# Label/unlabel and temperature contrasts retain the original independent bootstrap.
# Community contributions use detected-active union/zero masking; this is not evidence
# of zero biological growth in an untestable or non-active ASV at the other temperature.
# ASV growth is in gene-copy units, not carbon units. RGR is Fnew / incubation time.
# CUE.txt/deltaCUE.txt supply measured macro metrics/differences; CUE is not recalculated from qSIP.
# Original 90% CIs, Tukey/Duncan tests, and strict 3-replicate AFE filter are retained.
# Signed-log displays use sign(x) * log10(1 + abs(x)); raw values are also exported.

INPUT_DIR <- "C:/Users/79840/Desktop/qSIP数据及文献"
OUTPUT_DIR <- file.path(INPUT_DIR, "wildfire_qSIP_results")
FILES <- c(
  metadata = "input.txt",
  features = "16s-otu-table.tsv",
  taxonomy = "16s-taxonomy.tsv",
  macro_absolute = "CUE.txt",
  macro_delta = "deltaCUE.txt",
  tree = "rooted-tree.nwk",
  field_temperature = "daily average temperature.txt"
)

# Match temperature codes to input files and provide their actual temperatures.
TEMP_INPUT_CODES <- c("T21", "T28")
TEMP_C <- c(21, 28)
INCUBATION_DAYS <- 7
SOIL_WATER_18O_ATOM_FRACTION <- 0.50
PROP_O_FROM_WATER <- 0.60
MIN_FRACTIONS <- 3L
MIN_REPS <- 3L  # Complete 3+3 gradients required by this merged bootstrap model.
N_BOOT <- 1000L
APPLY_BOTTOM10_CORRECTION <- TRUE  # Retains the uploaded scripts; raw AFE is exported.
NON_GROWER_PROP <- 0.10
RANDOM_SEED <- 20260911L
K_NSIM_RGR <- 499L
RUN_SHARED_SENSITIVITY <- FALSE
RUN_ASINH_SENSITIVITY <- FALSE
CORRELOGRAM_N_POINTS <- 30L
CORRELOGRAM_CI_BS <- 99L
CORRELOGRAM_CI_CONF <- 0.90
CORRELOGRAM_MAX_ASVS <- 200L
INSTALL_MISSING_PACKAGES <- FALSE
STOP_ON_ERROR <- TRUE
# Twenty figures follow the retained plotting order; numbering is consecutive.
RUN_MODULES <- c("temperature", "density", "afe", "cue", "growth", "community",
                 "strategy", "phylo_g", "phylo_rgr", "contributors")
# Each run creates a new folder to prevent mixing old and current results.
RUN_ID <- format(Sys.time(), "run_%Y%m%d_%H%M%S")

main <- function() {
  if (!dir.exists(INPUT_DIR)) stop("INPUT_DIR does not exist: ", INPUT_DIR)
  if (length(TEMP_INPUT_CODES) != 2L || anyDuplicated(TEMP_INPUT_CODES) || length(TEMP_C) != 2L || any(!is.finite(TEMP_C)) || TEMP_C[2L] <= TEMP_C[1L]) {
    stop("Set two distinct temperature codes and increasing numeric temperatures.")
  }
  if (MIN_REPS != 3L) stop("This workflow requires MIN_REPS = 3; relaxing it needs a missingness-aware bootstrap.")
  if (N_BOOT < 100L || MIN_FRACTIONS < 1L || INCUBATION_DAYS <= 0 || SOIL_WATER_18O_ATOM_FRACTION <= 0.002000429 || SOIL_WATER_18O_ATOM_FRACTION > 1 ||
      PROP_O_FROM_WATER <= 0 || PROP_O_FROM_WATER > 1 || NON_GROWER_PROP <= 0 || NON_GROWER_PROP > 1) {
    stop("Invalid analysis settings.")
  }
  module_order <- c("temperature", "density", "afe", "cue", "growth", "community",
                    "strategy", "phylo_g", "phylo_rgr", "contributors")
  if (length(setdiff(RUN_MODULES, module_order))) stop("Unknown RUN_MODULES entry.")
  dependencies <- list(growth = "afe", community = c("afe", "growth"),
                       strategy = c("afe", "growth"), phylo_g = c("afe", "growth", "strategy"),
                       phylo_rgr = c("afe", "growth"), contributors = c("afe", "growth"))
  needed <- unique(c(RUN_MODULES, unlist(dependencies[RUN_MODULES], use.names = FALSE)))
  selected <- module_order[module_order %in% needed]
  module_inputs <- list(temperature = "field_temperature", density = "metadata",
                        afe = c("metadata", "features", "taxonomy"),
                        cue = c("macro_absolute", "macro_delta"), phylo_g = "tree", phylo_rgr = "tree")
  input_keys <- unique(unlist(module_inputs[selected], use.names = FALSE))
  INPUTS <- as.list(setNames(file.path(normalizePath(INPUT_DIR, winslash = "/"), FILES), names(FILES)))
  missing_files <- unlist(INPUTS[input_keys], use.names = FALSE)
  missing_files <- missing_files[!file.exists(missing_files)]
  if (length(missing_files)) stop("Missing input files:\n", paste(missing_files, collapse = "\n"))
  
  packages <- c("data.table", "Matrix", "matrixStats", "ggplot2", "patchwork", "scales")
  if (any(c("temperature", "cue") %in% selected)) packages <- c(packages, "dplyr", "tidyr")
  if ("cue" %in% selected) packages <- c(packages, "agricolae")
  if ("community" %in% selected) packages <- c(packages, "multcompView", "stringr", "RColorBrewer")
  if ("strategy" %in% selected) packages <- c(packages, "ggalluvial")
  if (any(c("phylo_g", "phylo_rgr") %in% selected)) packages <- c(packages, "ape")
  if ("phylo_rgr" %in% selected) packages <- c(packages, "phytools", "phylolm")
  if ("phylo_g" %in% selected) packages <- c(packages, "ggtree", "ggnewscale")
  if ("phylo_rgr" %in% selected) packages <- c(packages, "phylosignal", "phylobase")
  packages <- unique(packages)
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) && INSTALL_MISSING_PACKAGES) {
    cran <- setdiff(missing, "ggtree")
    if (length(cran)) install.packages(cran, repos = "https://cloud.r-project.org")
    if ("ggtree" %in% missing) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")
      BiocManager::install("ggtree", ask = FALSE, update = FALSE)
    }
    missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  }
  if (length(missing)) stop("Install these packages, or set INSTALL_MISSING_PACKAGES = TRUE: ", paste(missing, collapse = ", "))
  for (pkg in packages) suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  old_locale <- Sys.getlocale("LC_TIME")
  on.exit(Sys.setlocale("LC_TIME", old_locale), add = TRUE)
  Sys.setlocale("LC_TIME", "C")
  
  OUTPUT_ROOT <- file.path(OUTPUT_DIR, RUN_ID)
  if (dir.exists(OUTPUT_ROOT)) stop("Run folder already exists; choose a new RUN_ID: ", OUTPUT_ROOT)
  dir.create(OUTPUT_ROOT, recursive = TRUE)
  CACHE <- as.list(setNames(file.path(OUTPUT_ROOT, "Shared_data", module_order), module_order))
  for (p in CACHE[selected]) dir.create(p, recursive = TRUE, showWarnings = FALSE)
  INPUTS <- lapply(INPUTS, normalizePath, winslash = "/", mustWork = FALSE)
  data.table::fwrite(data.table::data.table(input = names(INPUTS), path = unlist(INPUTS)),
                     file.path(OUTPUT_ROOT, "Input_files.tsv"), sep = "\t")
  config_names <- c("TEMP_INPUT_CODES", "TEMP_C", "INCUBATION_DAYS", "SOIL_WATER_18O_ATOM_FRACTION",
                    "PROP_O_FROM_WATER", "MIN_FRACTIONS", "MIN_REPS", "N_BOOT", "APPLY_BOTTOM10_CORRECTION",
                    "NON_GROWER_PROP", "RANDOM_SEED", "K_NSIM_RGR", "RUN_SHARED_SENSITIVITY",
                    "RUN_ASINH_SENSITIVITY", "CORRELOGRAM_N_POINTS", "CORRELOGRAM_CI_BS", "CORRELOGRAM_MAX_ASVS")
  config_values <- vapply(config_names, function(n) paste(get(n, inherits = TRUE), collapse = ","), character(1))
  data.table::fwrite(data.table::data.table(parameter = config_names, value = config_values),
                     file.path(OUTPUT_ROOT, "Analysis_settings.tsv"), sep = "\t")
  writeLines(c(
    "Inputs: input.txt (fraction qPCR), 16s-otu-table.tsv, 16s-taxonomy.tsv, rooted-tree.nwk, SY2.txt, SY.txt, Daily average temperature.txt.",
    "Metadata requires: sample ID, fraction, density, abundance, iso_trt, fire_trt, tem, rep.",
    "ASV table: ASV IDs in the first column; count columns must match metadata fraction IDs exactly.",
    "SY2.txt requires fire_trt, Tem_trt, Respiration, Growth, CUE. SY.txt contains the three supplied differences, fire_trt, rep, in that order.",
    "Field temperature file: one header row and two blocks of Date, burned temperature, unburned temperature, difference (2022 and 2025).",
    "Temperature inputs are standardized internally to T_low/T_high; their actual degrees are in Analysis_settings.tsv.",
    "Original inferential definitions retained: 90% percentile CI; active if AFE lower CI > 0; strict 3+3-gradient eligibility.",
    "Bottom-10% AFE correction is optional and enabled by default to reproduce the uploaded workflow; it assumes the selected low-AFE taxa are non-growers.",
    "Detected-active community sums use the active union with other sides masked to zero; this is an operational summary, not proof of biological absence.",
    "Response/phylogeny analyses require testability at both temperatures and activity at either temperature; untestable and inactive pools are labeled Not testable.",
    "Community/taxon G is bacterial 16S gene-copy production per g soil per day; it is not an ASV-level carbon flux or CUE.",
    "RGR is the linear Fnew/time estimate retained from the uploaded code; exponential estimates are separate diagnostic columns.",
    "SY2 and SY metrics are supplied observations; this script does not infer pairing or recompute differences/CUE from these files.",
    "Signed-log displays use sign(x)*log10(1+abs(x)); contributor magnitudes use log10(1+magnitude). All raw values are retained.",
    "Top15 shares use true nonnegative contributions, with NA for zero denominators; no artificial one-copy floor is applied.",
    "Pairwise bootstrap tail probabilities and original Tukey/Duncan tests are retained without adding multiple-testing correction.",
    "Phylogenetic analyses use all eligible stage-specific taxa; only the correlogram may subsample to the configured maximum.",
    "Shared-ASV sensitivity is optional and cannot block stage-specific analyses when too few shared taxa remain.",
    "Field temperatures are observations from the supplied two periods, not an estimate of a three-year annual mean.",
    "One PDF per numbered folder, with source plot tables and associated statistics. Shared_data retains reusable computation tables/RDS objects."
  ), file.path(OUTPUT_ROOT, "Analysis_notes.txt"))
  
  FIRE_LEVELS <- c("UB", "F3M", "F3Y")
  FIRE_LABELS <- c(UB = "Unburned", F3M = "3 months", F3Y = "3 years")
  TEMP_LEVELS <- c("T_low", "T_high")
  TEMP_LABELS <- setNames(paste0(TEMP_C, "\u00B0C"), TEMP_LEVELS)
  TEMP_COLORS <- setNames(c("#4C78A8", "#E45756"), TEMP_LEVELS)
  TEMP_COLORS_DISPLAY <- setNames(unname(TEMP_COLORS), unname(TEMP_LABELS))
  FIRE_COLORS <- c(UB = "#4C78A8", F3M = "#72B7B2", F3Y = "#54A24B")
  FIRE_COLORS_DISPLAY <- setNames(unname(FIRE_COLORS), unname(FIRE_LABELS))
  REP_LEVELS <- c("1", "2", "3")
  ISO_LEVELS <- c("unlabel", "label")
  CI_PROBS <- c(0.05, 0.50, 0.95)
  N_BOOT_COMBINE <- N_BOOT
  NATURAL_18O <- 0.002000429
  SOIL_WATER_AFE <- SOIL_WATER_18O_ATOM_FRACTION - NATURAL_18O
  MAX_LABELABLE_DNA_AFE <- SOIL_WATER_AFE * PROP_O_FROM_WATER
  G_DENOMINATOR <- MAX_LABELABLE_DNA_AFE * INCUBATION_DAYS
  RGR_DENOMINATOR <- G_DENOMINATOR
  
  # Normalize temperature codes at the file boundary only.
  fread <- function(...) {
    x <- data.table::fread(...)
    if (data.table::is.data.table(x)) data.table::setnames(x, trimws(names(x))) else names(x) <- trimws(names(x))
    for (nm in intersect(c("tem", "Tem_trt"), names(x))) {
      vals <- trimws(as.character(x[[nm]]))
      if (!all(vals %in% TEMP_LEVELS)) {
        if (!all(vals %in% TEMP_INPUT_CODES)) stop("Unexpected temperature codes in ", nm, ": ", paste(unique(vals), collapse = ", "))
        vals <- TEMP_LEVELS[match(vals, TEMP_INPUT_CODES)]
      }
      if (data.table::is.data.table(x)) data.table::set(x, j = nm, value = vals) else x[[nm]] <- vals
    }
    x
  }
  unit <- grid::unit
  blank_plot <- function(label) ggplot2::ggplot() + ggplot2::annotate("text", x = 0, y = 0, label = label) + ggplot2::theme_void()
  figure_catalog <- data.table::data.table(
    module = c("temperature", "density", "afe", "cue", "cue", "community", "community", "community", "community", "strategy", "strategy", "strategy", "phylo_g", "phylo_g", "phylo_rgr", "phylo_rgr", "phylo_rgr", "contributors", "contributors", "contributors"),
    source_pdf = c("Field_temperature.pdf", "Fig_01_replicate_density_profiles_with_horizontal_WAD_90CI.pdf", "Fig_AFE_ASV_forest_6panel.pdf", "01_Fig_SY2_Absolute_Metrics_Tukey.pdf", "02_Fig_SY_Macro_Delta_Duncan.pdf", "Fig_AB_macro_thermal_response.pdf", "Fig_C_top10_phylum_DeltaG_polished.pdf", "Fig_ASV_Top15_heatmap_cumulative_DeltaG_percent_polished.pdf", "Fig_D_phylum_stacked_DeltaG.pdf", "Fig_A_volcano_3panel_signed_log10.pdf", "Fig_B_alluvial_pre_fire_all_thermal_strategies.pdf", "Fig_C_alluvial_pre_fire_Pos_ASVs.pdf", "Fig_A_circular_phylogenetic_tree_DeltaG.pdf", "Fig_C_Phylum_PieCharts.pdf", "Fig_B_BlombergK_DeltaRGR_scatter.pdf", "Fig_C_phylogenetic_correlogram_DeltaRGR.pdf", "Fig_S1_PagelLambda_DeltaRGR_scatter.pdf", "Fig_F2_Addon_Flux_Decomposition.pdf", "Fig_F3_Pos_Contributors_PointRange.pdf", "Fig_F3_Neg_Contributors_PointRange.pdf"),
    slug = c("Field_temperature", "Density_profiles", "AFE_all_ASVs", "Community_growth_respiration_CUE", "Community_carbon_temperature_response", "Detected_active_growth_response", "Top10_phylum_response", "Top15_positive_contribution_heatmap", "Phylum_growth_contributions", "ASV_response_classification", "Response_transitions_reference_pool", "Response_transitions_reference_positive", "Phylogeny_absolute_growth_response", "Phylum_response_composition", "Blomberg_K_relative_response", "Relative_response_phylogenetic_correlogram", "Pagel_lambda_relative_response", "Positive_negative_growth_decomposition", "Top15_positive_contributors", "Top15_negative_contributors"),
    table_pattern = c(".*", ".*", ".*", "^(01_|Raw_absolute|ANOVA_absolute|Tukey_)", "^(02_|Raw_delta|ANOVA_delta|Duncan_)", "^(01_|02_|Community_)", "^(03_|04_|Phylum_)", "^(05_|ASV_)", "^(04b_|Phylum_)", "^(01_|02_|03_|11_)", "^(04_|05_|06_|07_|08_|09_|10|11_)", "^(FigC_|11_)", "^(01_|02_|07|08|09_|stage_specific|.*tree.*nwk)", "^(07_|Phylum_response_)", "^(01_|02_|03_|06_|07_phylogenetic|10_|stage_specific)", "^(02_|07a|07b|10_)", "^(01_|02_|03_|06_|07_phylogenetic|10_|stage_specific)", "^Table_F2", "^Table_F3", "^Table_F3")
  )
  figure_catalog[, number := sprintf("Fig%02d", seq_len(.N))]
  figure_catalog[, folder := file.path(OUTPUT_ROOT, paste(number, slug, sep = "_"))]
  figure_catalog[, status := ifelse(module %in% selected, "not_generated", "not_requested")]
  figure_catalog[, pdf := file.path(folder, paste0(number, ".pdf"))]
  figure_catalog[, note := ""]
  for (d in figure_catalog[module %in% selected, folder]) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  flush_catalog <- function() {
    data.table::fwrite(figure_catalog, file.path(OUTPUT_ROOT, "Figure_index.tsv"), sep = "\t")
    for (i in which(figure_catalog$module %in% selected)) {
      writeLines(c(paste("Figure:", figure_catalog$number[i]), paste("Status:", figure_catalog$status[i]),
                   paste("Note:", figure_catalog$note[i])), file.path(figure_catalog$folder[i], "Status.txt"))
    }
  }
  flush_catalog()
  
  # Export actual plot-source tables, including patchwork panels and layer data.
  export_plot_data <- function(plot, folder, prefix = "plot") {
    write_source <- function(x, suffix) {
      if (!is.data.frame(x) || !ncol(x)) return(invisible(NULL))
      x <- data.table::copy(data.table::as.data.table(x))
      for (nm in names(x)) if (is.list(x[[nm]])) {
        data.table::set(x, j = nm, value = vapply(x[[nm]], function(v) paste(v, collapse = ";"), character(1)))
      }
      data.table::fwrite(x, file.path(folder, paste0(prefix, suffix, ".tsv")), sep = "\t")
    }
    if (inherits(plot, "patchwork")) {
      for (i in seq_len(length(plot))) export_plot_data(plot[[i]], folder, paste0(prefix, "_panel", sprintf("%02d", i)))
    } else if (inherits(plot, "ggplot")) {
      write_source(plot$data, "_source")
      for (i in seq_along(plot$layers)) write_source(plot$layers[[i]]$data, paste0("_layer", sprintf("%02d", i)))
    }
    invisible(NULL)
  }
  save_figure <- function(filename, plot, ..., device = NULL) {
    idx <- match(basename(filename), figure_catalog$source_pdf)
    if (is.na(idx)) stop("Unregistered figure: ", filename)
    folder <- figure_catalog$folder[idx]
    dir.create(folder, recursive = TRUE, showWarnings = FALSE)
    if (is.null(plot)) plot <- blank_plot("No eligible data for this figure")
    pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
    ggplot2::ggsave(filename = figure_catalog$pdf[idx], plot = plot, ..., device = pdf_device)
    export_plot_data(plot, folder)
    figure_catalog[idx, status := "generated"]
    flush_catalog()
    invisible(figure_catalog$pdf[idx])
  }
  copy_module_tables <- function(module_name) {
    source_dir <- CACHE[[module_name]]
    files <- list.files(source_dir, pattern = "\\.(tsv|txt|nwk)$", full.names = TRUE)
    rows <- which(figure_catalog$module == module_name)
    for (i in rows) {
      matched <- files[grepl(figure_catalog$table_pattern[i], basename(files))]
      if (!length(matched)) next
      dir.create(figure_catalog$folder[i], recursive = TRUE, showWarnings = FALSE)
      ok <- file.copy(matched, file.path(figure_catalog$folder[i], basename(matched)), overwrite = TRUE)
      if (!all(ok)) stop("Failed to copy associated figure tables for ", figure_catalog$number[i])
    }
  }
  
  safe_quantiles <- function(x, probs = CI_PROBS) {
    x <- x[is.finite(x)]
    
    if (length(x) == 0) {
      return(rep(NA_real_, length(probs)))
    }
    
    as.numeric( quantile( x, probs = probs, na.rm = TRUE, names = FALSE, type = 7 ) )
  }
  
  signed_log10_abs <- function(x) sign(x) * log10(1 + abs(x))
  
  safe_numeric <- function(x) {
    suppressWarnings(as.numeric(x))
  }
  
  safe_asinh_scale <- function(x) {
    x <- x[is.finite(x)]
    
    if (length(x) == 0) {
      return(1)
    }
    
    s <- median(abs(x), na.rm = TRUE)
    
    if (!is.finite(s) || s <= 0) {
      s <- stats::sd(x, na.rm = TRUE)
    }
    
    if (!is.finite(s) || s <= 0) {
      s <- 1
    }
    
    s
  }
  
  classify_delta_strategy <- function(ci_low, ci_high) {
    out <- rep(NA_character_, length(ci_low))
    
    out[ is.finite(ci_low) & ci_low > 0 ] <- "Pos"
    
    out[ is.finite(ci_high) & ci_high < 0 ] <- "Neg"
    
    out[ is.na(out) & is.finite(ci_low) & is.finite(ci_high) ] <- "Neutral"
    
    out
  }
  
  combine_independent_bootstrap_difference <- function(
    mat_hot,
    mat_cold,
    n_boot = N_BOOT_COMBINE
  ) {
    if (nrow(mat_hot) != nrow(mat_cold)) {
      stop("mat_hot and mat_cold must have the same row count.")
    }
    
    n_hot <- ncol(mat_hot)
    n_cold <- ncol(mat_cold)
    
    if (n_hot <= 0 || n_cold <= 0) {
      stop("Bootstrap matrices must have at least one column.")
    }
    
    hot_idx <- sample.int( n_hot, size = n_boot, replace = TRUE )
    
    cold_idx <- sample.int( n_cold, size = n_boot, replace = TRUE )
    
    mat_hot[, hot_idx, drop = FALSE] - mat_cold[, cold_idx, drop = FALSE]
  }
  
  format_p <- function(p) {
    if (!is.finite(p)) {
      return("P = NA")
    }
    
    if (p < 0.001) {
      return("P < 0.001")
    }
    
    sprintf("P = %.3f", p)
  }
  
  make_signal_axis_labels <- function(max_value) {
    if (!is.finite(max_value) || max_value <= 0) {
      return(scales::label_number(accuracy = 0.01))
    }
    
    if (max_value < 0.001) {
      return(scales::label_scientific(digits = 2))
    }
    
    if (max_value < 0.01) {
      return(scales::label_number(accuracy = 0.001))
    }
    
    if (max_value < 0.1) {
      return(scales::label_number(accuracy = 0.01))
    }
    
    scales::label_number(accuracy = 0.1)
  }
  
  # Temperature
  run_temperature <- function() {
    data_raw <- read.table(INPUTS$field_temperature, skip = 1, fill = TRUE, stringsAsFactors = FALSE)
    data_2022 <- data_raw[, 1:4]
    colnames(data_2022) <- c("Date", "T_B", "T_UB", "Delta_T")
    data_2022$Date <- as.Date(data_2022$Date)
    data_2022 <- data_2022[!is.na(data_2022$Date), ]
    if (!nrow(data_2022)) stop("No valid 2022 dates in field-temperature data.")
    avg_TB_22 <- mean(data_2022$T_B, na.rm = TRUE)
    avg_TUB_22 <- mean(data_2022$T_UB, na.rm = TRUE)
    avg_Delta_22 <- mean(data_2022$Delta_T, na.rm = TRUE)
    label_2022 <- sprintf("Mean T-B (Burned): %.1f °C\nMean T-UB (Mature forest): %.1f °C\nMean ∆T (Difference): %.1f °C",
                          avg_TB_22, avg_TUB_22, avg_Delta_22)
    data_2025 <- data_raw[, 5:8]
    colnames(data_2025) <- c("Date", "T_B", "T_UB", "Delta_T")
    data_2025$Date <- as.Date(data_2025$Date)
    data_2025 <- data_2025[!is.na(data_2025$Date), ]
    if (!nrow(data_2025)) stop("No valid 2025 dates in field-temperature data.")
    avg_TB_25 <- mean(data_2025$T_B, na.rm = TRUE)
    avg_TUB_25 <- mean(data_2025$T_UB, na.rm = TRUE)
    avg_Delta_25 <- mean(data_2025$Delta_T, na.rm = TRUE)
    label_2025 <- sprintf("Mean T-B (Burned): %.1f °C\nMean T-UB (Mature forest): %.1f °C\nMean ∆T (Difference): %.1f °C",
                          avg_TB_25, avg_TUB_25, avg_Delta_25)
    p1 <- ggplot(data_2022, aes(x = Date)) + geom_line(aes(y = T_B), color = "red", alpha = 0.3) +
      geom_line(aes(y = T_UB), color = "blue", alpha = 0.3) +
      geom_point(aes(y = T_B, color = "T-B (Burned)"), size = 2) +
      geom_point(aes(y = T_UB, color = "T-UB (Mature forest)"), size = 2) +
      scale_color_manual(name = "Site Type", values = c("T-B (Burned)" = "red", "T-UB (Mature forest)" = "blue")) +
      annotate("text", x = min(data_2022$Date), y = max(data_2022$T_B, na.rm = TRUE),
               label = label_2022, hjust = 0, vjust = 1, size = 4, fontface = "bold") +
      scale_x_date(date_labels = "%b %d", date_breaks = "10 days") +
      scale_y_continuous(expand = expansion(mult = c(0.1, 0.25))) +
      labs(title = "3 Months Post-fire (July 1 - August 31, 2022)", x = "Date", y = "0-10 cm Mineral Soil Temp (°C)") +
      theme_bw() +
      theme(legend.position = "bottom")
    p2 <- ggplot(data_2025, aes(x = Date)) + geom_line(aes(y = T_B), color = "red", alpha = 0.3) +
      geom_line(aes(y = T_UB), color = "blue", alpha = 0.3) +
      geom_point(aes(y = T_B, color = "T-B (Burned)"), size = 2) +
      geom_point(aes(y = T_UB, color = "T-UB (Mature forest)"), size = 2) +
      scale_color_manual(name = "Site Type", values = c("T-B (Burned)" = "red", "T-UB (Mature forest)" = "blue")) +
      annotate("text", x = min(data_2025$Date), y = max(data_2025$T_B, na.rm = TRUE),
               label = label_2025, hjust = 0, vjust = 1, size = 4, fontface = "bold") +
      scale_x_date(date_labels = "%b %d", date_breaks = "10 days") +
      scale_y_continuous(expand = expansion(mult = c(0.1, 0.25))) +
      labs(title = "3 Years Post-fire (July 1 - August 31, 2025)", x = "Date", y = "0-10 cm Mineral Soil Temp (°C)") +
      theme_bw() +
      theme(legend.position = "bottom")
    final_plot <- p1 / p2 + plot_layout(guides = 'collect') & theme(legend.position = 'bottom')
    save_figure("Field_temperature.pdf", final_plot, width = 12, height = 8, units = "in")
    fwrite(data_2022, file.path(CACHE$temperature, "Daily_temperature_2022.tsv"), sep = "\t")
    fwrite(data_2025, file.path(CACHE$temperature, "Daily_temperature_2025.tsv"), sep = "\t")
    fwrite(data.table(stage = c("F3M", "F3Y"), burned_mean = c(avg_TB_22, avg_TB_25),
                      unburned_mean = c(avg_TUB_22, avg_TUB_25), difference_mean = c(avg_Delta_22, avg_Delta_25)),
           file.path(CACHE$temperature, "Mean_temperature_summary.tsv"), sep = "\t")
    invisible(NULL)
  }
  
  # Density
  run_density <- function() {
    INPUT_FILE <- INPUTS$metadata
    output_dir <- CACHE$density
    if (!file.exists(INPUT_FILE)) {
      stop("Missing input file: ", INPUT_FILE)
    }
    
    ISO_LABELS  <- c(unlabel = "16O control", label = "18O label")
    ISO_COLORS  <- c(unlabel = "#4C78A8", label = "#E45756")
    ISO_LINES   <- c(unlabel = "dashed", label = "solid")
    
    bootstrap_independent_mean_difference <- function(label_values, unlabel_values, n_boot = N_BOOT) {
      label_values <- label_values[is.finite(label_values)]
      unlabel_values <- unlabel_values[is.finite(unlabel_values)]
      if (length(label_values) == 0 || length(unlabel_values) == 0) {
        return(rep(NA_real_, n_boot))
      }
      
      out <- numeric(n_boot)
      for (b in seq_len(n_boot)) {
        lab_b <- sample(label_values, size = length(label_values), replace = TRUE)
        unl_b <- sample(unlabel_values, size = length(unlabel_values), replace = TRUE)
        out[b] <- mean(lab_b) - mean(unl_b)
      }
      out
    }
    
    format_wad <- function(x) {
      if (!is.finite(x)) return("NA")
      sprintf("%.4f", x)
    }
    
    # Read and validate fraction-level qPCR metadata.
    meta <- fread(INPUT_FILE, sep = "\t", header = TRUE, check.names = FALSE)
    setnames(meta, trimws(names(meta)))
    if ("OTU ID" %in% names(meta)) {
      setnames(meta, "OTU ID", "fraction_sample_id")
    } else if ("sample_id" %in% names(meta)) {
      setnames(meta, "sample_id", "fraction_sample_id")
    } else if (!"fraction_sample_id" %in% names(meta)) {
      setnames(meta, names(meta)[1], "fraction_sample_id")
    }
    
    required_cols <- c(
      "fraction_sample_id", "fraction", "density", "abundance",
      "iso_trt", "fire_trt", "tem", "rep"
    )
    missing_cols <- setdiff(required_cols, names(meta))
    if (length(missing_cols) > 0) {
      stop("input(1).txt missing required columns: ", paste(missing_cols, collapse = ", "))
    }
    
    meta[, fraction_sample_id := trimws(as.character(fraction_sample_id))]
    meta[, iso_trt := trimws(as.character(iso_trt))]
    meta[, fire_trt := trimws(as.character(fire_trt))]
    meta[, tem := trimws(as.character(tem))]
    meta[, rep := trimws(as.character(rep))]
    meta[, fraction := as.integer(fraction)]
    meta[, density := as.numeric(density)]
    meta[, abundance := as.numeric(abundance)]
    if (anyDuplicated(meta$fraction_sample_id)) {
      stop("Duplicate fraction_sample_id values detected in input(1).txt.")
    }
    
    if (any(!meta$fire_trt %chin% FIRE_LEVELS)) {
      stop("Unexpected fire_trt values: ", paste(setdiff(unique(meta$fire_trt), FIRE_LEVELS), collapse = ", "))
    }
    
    if (any(!meta$tem %chin% TEMP_LEVELS)) {
      stop("Unexpected tem values: ", paste(setdiff(unique(meta$tem), TEMP_LEVELS), collapse = ", "))
    }
    
    if (any(!meta$iso_trt %chin% ISO_LEVELS)) {
      stop("Unexpected iso_trt values: ", paste(setdiff(unique(meta$iso_trt), ISO_LEVELS), collapse = ", "))
    }
    
    if (any(!meta$rep %chin% REP_LEVELS)) {
      stop("Unexpected rep values: ", paste(setdiff(unique(meta$rep), REP_LEVELS), collapse = ", "))
    }
    
    if (any(!is.finite(meta$density)) || any(!is.finite(meta$abundance)) || any(meta$abundance < 0)) {
      stop("Non-finite density/abundance or negative abundance values detected.")
    }
    
    meta[, tube_id := paste(fire_trt, tem, iso_trt, rep, sep = "_")]
    meta[, comparison_id := paste(fire_trt, tem, sep = "_")]
    expected_tubes <- CJ(
      fire_trt = FIRE_LEVELS,
      tem = TEMP_LEVELS,
      iso_trt = ISO_LEVELS,
      rep = REP_LEVELS,
      unique = TRUE
    )
    expected_tubes[, tube_id := paste(fire_trt, tem, iso_trt, rep, sep = "_")]
    observed_tubes <- unique(meta[, .(fire_trt, tem, iso_trt, rep, tube_id)])
    missing_tubes <- expected_tubes[!observed_tubes, on = .(fire_trt, tem, iso_trt, rep)]
    extra_tubes <- observed_tubes[!expected_tubes, on = .(fire_trt, tem, iso_trt, rep)]
    if (nrow(missing_tubes) > 0) {
      fwrite(missing_tubes, file.path(output_dir, "QC_missing_expected_tubes.tsv"), sep = "\t")
      stop("Complete factorial design check failed: missing expected tubes.")
    }
    
    if (nrow(extra_tubes) > 0) {
      fwrite(extra_tubes, file.path(output_dir, "QC_extra_unexpected_tubes.tsv"), sep = "\t")
      stop("Complete factorial design check failed: unexpected tubes detected.")
    }
    
    if (nrow(observed_tubes) != 36L) {
      stop("Expected exactly 36 gradient tubes; observed ", nrow(observed_tubes), ".")
    }
    
    comparison_manifest <- CJ( fire_trt = FIRE_LEVELS, tem = TEMP_LEVELS, unique = TRUE )
    comparison_manifest[, comparison_id := paste(fire_trt, tem, sep = "_")]
    comparison_manifest[, unlabel_group := paste(fire_trt, tem, "unlabel", sep = "_")]
    comparison_manifest[, label_group := paste(fire_trt, tem, "label", sep = "_")]
    comparison_manifest[, bootstrap_design := "independent_unlabel_and_label"]
    comparison_manifest[, replicate_pairing := "none"]
    fwrite(
      comparison_manifest,
      file.path(output_dir, "01_exact_label_unlabel_comparison_manifest.tsv"),
      sep = "\t"
    )
    tube_structure <- meta[
      order(fraction),
      .(
        n_fraction = .N,
        min_fraction = min(fraction),
        max_fraction = max(fraction),
        min_density = min(density),
        max_density = max(density),
        total_qPCR = sum(abundance),
        density_monotonic_decreasing = all(diff(density) <= 0),
        density_monotonic_increasing = all(diff(density) >= 0)
      ),
      by = .(fire_trt, tem, iso_trt, rep, tube_id, comparison_id)
    ]
    fwrite( tube_structure, file.path(output_dir, "02_gradient_tube_structure.tsv"), sep = "\t" )
    if (any(!(tube_structure$density_monotonic_decreasing | tube_structure$density_monotonic_increasing))) {
      warning("At least one gradient has non-monotonic density order. Check 02_gradient_tube_structure.tsv.")
    }
    
    if (any(!is.finite(tube_structure$total_qPCR) | tube_structure$total_qPCR <= 0)) {
      bad_tubes <- tube_structure[!is.finite(total_qPCR) | total_qPCR <= 0]
      fwrite(bad_tubes, file.path(output_dir, "QC_nonpositive_total_qPCR_tubes.tsv"), sep = "\t")
      stop("At least one gradient tube has total qPCR <= 0. See QC_nonpositive_total_qPCR_tubes.tsv.")
    }
    
    if (any(tube_structure$n_fraction < 2L)) {
      stop("At least one gradient tube has fewer than 2 fractions.")
    }
    
    meta[, rel_qPCR := abundance / sum(abundance), by = tube_id]
    meta[, rel_qPCR_percent := 100 * rel_qPCR]
    if (any(!is.finite(meta$rel_qPCR_percent))) {
      stop("Non-finite relative qPCR values detected after within-tube normalization.")
    }
    
    fwrite(
      meta[order(fire_trt, tem, iso_trt, rep, fraction)],
      file.path(output_dir, "03_fraction_qPCR_normalized.tsv"),
      sep = "\t"
    )
    profile_replicates <- meta[
      order(fire_trt, tem, iso_trt, rep, density),
      .(
        fraction_sample_id,
        fire_trt,
        tem,
        iso_trt,
        rep,
        tube_id,
        comparison_id,
        fraction,
        density,
        abundance,
        rel_qPCR,
        rel_qPCR_percent
      )
    ]
    fwrite( profile_replicates, file.path(output_dir, "04_density_profile_replicates_for_plot.tsv"), sep = "\t" )
    community_wad <- meta[
      ,
      .( community_WAD = sum(density * abundance) / sum(abundance), total_qPCR = sum(abundance), n_fraction = .N ),
      by = .(fire_trt, tem, iso_trt, rep, tube_id, comparison_id)
    ]
    if (any(!is.finite(community_wad$community_WAD))) {
      stop("Non-finite community WAD detected.")
    }
    
    fwrite( community_wad, file.path(output_dir, "05_community_WAD_per_tube.tsv"), sep = "\t" )
    set.seed(RANDOM_SEED + 201L)
    wad_group_summary_list <- list()
    wad_group_boot_list <- list()
    for (current_fire in FIRE_LEVELS) {
      for (current_tem in TEMP_LEVELS) {
        for (current_iso in ISO_LEVELS) {
          wad_values <- community_wad[
            fire_trt == current_fire & tem == current_tem & iso_trt == current_iso,
            community_WAD
          ]
          if (length(wad_values) != 3L) {
            stop(
              "Expected 3 independent WAD values for ",
              current_fire, " / ", current_tem, " / ", current_iso, "."
            )
          }
          
          wad_boot <- numeric(N_BOOT)
          for (b in seq_len(N_BOOT)) {
            wad_boot[b] <- mean(sample(wad_values, size = length(wad_values), replace = TRUE))
          }
          
          wad_q <- safe_quantiles(wad_boot)
          wad_group_summary_list[[length(wad_group_summary_list) + 1L]] <- data.table(
            fire_trt = current_fire,
            tem = current_tem,
            iso_trt = current_iso,
            comparison_id = paste(current_fire, current_tem, sep = "_"),
            n_rep = length(wad_values),
            mean_WAD = mean(wad_values),
            sd_WAD = sd(wad_values),
            WAD_CI_low_90 = wad_q[1],
            WAD_boot_median = wad_q[2],
            WAD_CI_high_90 = wad_q[3],
            n_boot = N_BOOT
          )
          wad_group_boot_list[[length(wad_group_boot_list) + 1L]] <- data.table(
            fire_trt = current_fire,
            tem = current_tem,
            iso_trt = current_iso,
            bootstrap_iteration = seq_len(N_BOOT),
            WAD_boot_mean = wad_boot
          )
        }
      }
    }
    
    wad_group_summary <- rbindlist(wad_group_summary_list, use.names = TRUE, fill = TRUE)
    wad_group_boot <- rbindlist(wad_group_boot_list, use.names = TRUE, fill = TRUE)
    fwrite( wad_group_summary, file.path(output_dir, "06_community_WAD_group_summary_with_90CI.tsv"), sep = "\t" )
    fwrite( wad_group_boot, file.path(output_dir, "06b_community_WAD_group_bootstrap.tsv"), sep = "\t" )
    set.seed(RANDOM_SEED)
    wad_delta_summary_list <- list()
    wad_delta_boot_list <- list()
    for (g in seq_len(nrow(comparison_manifest))) {
      current_fire <- comparison_manifest$fire_trt[g]
      current_tem <- comparison_manifest$tem[g]
      current_id <- comparison_manifest$comparison_id[g]
      label_values <- community_wad[
        fire_trt == current_fire & tem == current_tem & iso_trt == "label",
        community_WAD
      ]
      unlabel_values <- community_wad[
        fire_trt == current_fire & tem == current_tem & iso_trt == "unlabel",
        community_WAD
      ]
      if (length(label_values) != 3L || length(unlabel_values) != 3L) {
        stop("Expected 3 label and 3 unlabel WAD values for ", current_id, ".")
      }
      
      delta_boot <- bootstrap_independent_mean_difference( label_values, unlabel_values, N_BOOT )
      delta_q <- safe_quantiles(delta_boot)
      wad_delta_summary_list[[length(wad_delta_summary_list) + 1L]] <- data.table(
        fire_trt = current_fire,
        tem = current_tem,
        comparison_id = current_id,
        mean_WAD_unlabel = mean(unlabel_values),
        mean_WAD_label = mean(label_values),
        delta_WAD_observed = mean(label_values) - mean(unlabel_values),
        delta_WAD_CI_low_90 = delta_q[1],
        delta_WAD_boot_median = delta_q[2],
        delta_WAD_CI_high_90 = delta_q[3],
        DeltaWAD_positive_90CI = delta_q[1] > 0,
        bootstrap_design = "independent_unlabel_and_label",
        n_boot = N_BOOT
      )
      wad_delta_boot_list[[length(wad_delta_boot_list) + 1L]] <- data.table(
        fire_trt = current_fire,
        tem = current_tem,
        comparison_id = current_id,
        bootstrap_iteration = seq_len(N_BOOT),
        delta_WAD = delta_boot
      )
    }
    
    wad_delta_summary <- rbindlist(wad_delta_summary_list, use.names = TRUE, fill = TRUE)
    wad_delta_boot <- rbindlist(wad_delta_boot_list, use.names = TRUE, fill = TRUE)
    fwrite( wad_delta_summary, file.path(output_dir, "07_community_WAD_exact_comparison_summary.tsv"), sep = "\t" )
    fwrite( wad_delta_boot, file.path(output_dir, "08_community_WAD_independent_bootstrap.tsv"), sep = "\t" )
    plot_profile <- copy(profile_replicates)
    plot_wad <- copy(wad_group_summary)
    plot_delta <- copy(wad_delta_summary)
    panel_scale <- plot_profile[
      ,
      .(
        curve_max = max(rel_qPCR_percent, na.rm = TRUE),
        density_min = min(density, na.rm = TRUE),
        density_max = max(density, na.rm = TRUE)
      ),
      by = .(fire_trt, tem)
    ]
    plot_wad <- merge( plot_wad, panel_scale, by = c("fire_trt", "tem"), all.x = TRUE, sort = FALSE )
    plot_wad[, WAD_CI_y := fifelse( iso_trt == "label", curve_max * 1.31, curve_max * 1.18 )]
    plot_wad[, WAD_cap_halfheight := pmax(curve_max * 0.020, 0.03)]
    plot_delta <- merge( plot_delta, panel_scale, by = c("fire_trt", "tem"), all.x = TRUE, sort = FALSE )
    plot_delta[, annot_x := density_min + 0.02 * (density_max - density_min)]
    plot_delta[, annot_y := curve_max * 1.49]
    plot_delta[, wad_label := sprintf(
      "DeltaWAD = %s\n90%% CI [%s, %s]",
      vapply(delta_WAD_observed, format_wad, character(1)),
      vapply(delta_WAD_CI_low_90, format_wad, character(1)),
      vapply(delta_WAD_CI_high_90, format_wad, character(1))
    )]
    plot_profile[, fire_trt := factor(fire_trt, levels = FIRE_LEVELS, labels = FIRE_LABELS[FIRE_LEVELS])]
    plot_profile[, tem := factor(tem, levels = TEMP_LEVELS, labels = TEMP_LABELS[TEMP_LEVELS])]
    plot_profile[, iso_trt := factor(iso_trt, levels = ISO_LEVELS, labels = ISO_LABELS[ISO_LEVELS])]
    plot_profile[, rep := factor(rep, levels = REP_LEVELS)]
    plot_wad[, fire_trt := factor(fire_trt, levels = FIRE_LEVELS, labels = FIRE_LABELS[FIRE_LEVELS])]
    plot_wad[, tem := factor(tem, levels = TEMP_LEVELS, labels = TEMP_LABELS[TEMP_LEVELS])]
    plot_wad[, iso_trt := factor(iso_trt, levels = ISO_LEVELS, labels = ISO_LABELS[ISO_LEVELS])]
    plot_delta[, fire_trt := factor(fire_trt, levels = FIRE_LEVELS, labels = FIRE_LABELS[FIRE_LEVELS])]
    plot_delta[, tem := factor(tem, levels = TEMP_LEVELS, labels = TEMP_LABELS[TEMP_LEVELS])]
    plot_colors <- setNames(ISO_COLORS[ISO_LEVELS], ISO_LABELS[ISO_LEVELS])
    plot_lines <- setNames(ISO_LINES[ISO_LEVELS], ISO_LABELS[ISO_LEVELS])
    plot_shapes <- c( "16O control" = 1, "18O label" = 15 )
    p_density <- ggplot(
      plot_profile,
      aes(
        x = density,
        y = rel_qPCR_percent,
        colour = iso_trt,
        linetype = iso_trt,
        shape = iso_trt,
        group = interaction(iso_trt, rep, drop = TRUE)
      )
    ) +
      geom_line( linewidth = 0.62, alpha = 0.68 ) +
      geom_point( size = 1.55, stroke = 0.65, alpha = 0.82 ) +
      geom_segment(
        data = plot_wad,
        aes( x = WAD_CI_low_90, xend = WAD_CI_high_90, y = WAD_CI_y, yend = WAD_CI_y, colour = iso_trt ),
        inherit.aes = FALSE,
        linewidth = 1.00
      ) +
      geom_segment(
        data = plot_wad,
        aes(
          x = WAD_CI_low_90,
          xend = WAD_CI_low_90,
          y = WAD_CI_y - WAD_cap_halfheight,
          yend = WAD_CI_y + WAD_cap_halfheight,
          colour = iso_trt
        ),
        inherit.aes = FALSE,
        linewidth = 0.85
      ) +
      geom_segment(
        data = plot_wad,
        aes(
          x = WAD_CI_high_90,
          xend = WAD_CI_high_90,
          y = WAD_CI_y - WAD_cap_halfheight,
          yend = WAD_CI_y + WAD_cap_halfheight,
          colour = iso_trt
        ),
        inherit.aes = FALSE,
        linewidth = 0.85
      ) +
      geom_point(
        data = plot_wad,
        aes( x = mean_WAD, y = WAD_CI_y, colour = iso_trt, shape = iso_trt ),
        inherit.aes = FALSE,
        size = 3.15,
        stroke = 1.0
      ) +
      geom_text(
        data = plot_delta,
        aes(x = annot_x, y = annot_y, label = wad_label),
        inherit.aes = FALSE,
        hjust = 0,
        vjust = 1,
        size = 2.85,
        lineheight = 1.03,
        colour = "black"
      ) +
      facet_grid(rows = vars(tem), cols = vars(fire_trt)) +
      scale_colour_manual(values = plot_colors) +
      scale_shape_manual(values = plot_shapes) +
      scale_linetype_manual(values = plot_lines) +
      scale_x_continuous( labels = label_number(accuracy = 0.01), expand = expansion(mult = c(0.02, 0.02)) ) +
      scale_y_continuous( expand = expansion(mult = c(0.02, 0.03)) ) +
      labs(
        x = expression(paste("Buoyant density (g mL"^{-1}, ")")),
        y = "Relative qPCR abundance (%)",
        colour = NULL,
        shape = NULL,
        linetype = NULL
      ) +
      theme_classic(base_size = 12) +
      theme(
        strip.background = element_rect( fill = "grey95", colour = "black", linewidth = 0.6 ),
        strip.text = element_text(size = 11, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 9.5, colour = "black"),
        legend.position = "top",
        legend.box = "horizontal",
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
        panel.spacing = grid::unit(0.65, "lines"),
        plot.margin = margin(8, 10, 8, 8)
      )
    save_figure(
      filename = file.path( output_dir, "Fig_01_replicate_density_profiles_with_horizontal_WAD_90CI.pdf" ),
      plot = p_density,
      width = 12.5,
      height = 7.8,
      units = "in"
    )
    parameter_table <- data.table(
      Parameter = c(
        "Expected fire stages",
        "Expected temperatures",
        "Expected isotope treatments",
        "Expected biological replicates",
        "Expected gradient tubes",
        "Allowed fraction-count variation",
        "DeltaWAD bootstrap design",
        "Bootstrap iterations",
        "Confidence interval",
        "Within-tube normalization",
        "Density-profile display",
        "WAD CI display",
        "Mean profile across fractions"
      ),
      Value = c(
        paste(FIRE_LEVELS, collapse = ","),
        paste(TEMP_LEVELS, collapse = ","),
        paste(ISO_LEVELS, collapse = ","),
        paste(REP_LEVELS, collapse = ","),
        "36",
        "retained; actual fraction count is used for each tube",
        "independent_unlabel_and_label",
        as.character(N_BOOT),
        "90% bootstrap percentile CI (5th-95th)",
        "fraction qPCR / sum fraction qPCR within the same tube",
        "three independent replicate profiles per isotope group",
        "horizontal CI bar for isotope-group mean community WAD",
        "not calculated; no direct averaging by fraction number"
      )
    )
    fwrite( parameter_table, file.path(output_dir, "09_analysis_parameters.tsv"), sep = "\t" )
    invisible(NULL)
  }
  
  # Afe
  run_afe <- function() {
    META_FILE <- INPUTS$metadata
    FEATURE_FILE <- INPUTS$features
    TAXONOMY_FILE <- INPUTS$taxonomy
    output_dir <- CACHE$afe
    bootstrap_dir <- file.path(output_dir, "bootstrap_objects")
    if (!dir.exists(bootstrap_dir)) {
      dir.create(bootstrap_dir, recursive = TRUE)
    }
    
    required_input_files <- c(META_FILE, FEATURE_FILE, TAXONOMY_FILE)
    missing_input_files <- required_input_files[!file.exists(required_input_files)]
    if (length(missing_input_files) > 0) {
      stop("Missing required input files:\n", paste(missing_input_files, collapse = "\n"))
    }
    
    WAD_TO_GC_SLOPE <- 0.083506
    WAD_TO_GC_INTERCEPT <- 1.646057
    MW_LIGHT_AT <- 307.691
    MW_GC_COEF <- 0.496
    DELTA_MW_18O_MAX <- 12.07747
    comparison_table <- CJ(fire_trt = FIRE_LEVELS, tem = TEMP_LEVELS, unique = TRUE)
    comparison_table[, fire_order := match(fire_trt, FIRE_LEVELS)]
    comparison_table[, temp_order := match(tem, TEMP_LEVELS)]
    setorder(comparison_table, fire_order, temp_order)
    comparison_table[, comparison_id := paste(fire_trt, tem, sep = "_")]
    comparison_table[, c("fire_order", "temp_order") := NULL]
    comparison_table[, unlabel_group := paste(fire_trt, tem, "unlabel", sep = "_")]
    comparison_table[, label_group := paste(fire_trt, tem, "label", sep = "_")]
    comparison_table[, bootstrap_design := "independent_unlabel_and_label"]
    fwrite( comparison_table, file.path(output_dir, "01_exact_label_unlabel_comparison_manifest.tsv"), sep = "\t" )
    
    make_boot_weights <- function(n, n_boot) {
      out <- matrix(0, nrow = n, ncol = n_boot)
      for (b in seq_len(n_boot)) {
        sampled <- sample.int(n, size = n, replace = TRUE)
        out[, b] <- tabulate(sampled, nbins = n) / n
      }
      out
    }
    
    bottom_non_grower_shift <- function(x, prop = NON_GROWER_PROP) {
      if (!APPLY_BOTTOM10_CORRECTION) return(0)
      x <- x[is.finite(x)]
      if (length(x) == 0) {
        return(NA_real_)
      }
      n_bottom <- max(1L, floor(length(x) * prop))
      median(sort(x)[seq_len(n_bottom)])
    }
    
    calc_AFE_18O <- function(Wlight, Wlabel) {
      GC_estimated <- (Wlight - WAD_TO_GC_INTERCEPT) / WAD_TO_GC_SLOPE
      Mlight <- MW_GC_COEF * GC_estimated + MW_LIGHT_AT
      Mlabel <- (((Wlabel - Wlight) / Wlight) + 1) * Mlight
      ((Mlabel - Mlight) / DELTA_MW_18O_MAX) * (1 - NATURAL_18O)
    }
    
    clean_tax_value <- function(x) {
      x <- trimws(as.character(x))
      x <- sub("^[A-Za-z]__", "", x)
      x[x %chin% c("", "NA", "Na", "na", "Unassigned", "unassigned")] <- NA_character_
      x
    }
    
    extract_rank_from_taxonomy <- function(x, rank_prefix) {
      parts <- strsplit(as.character(x), ";", fixed = TRUE)
      vapply(
        parts,
        function(z) {
          z <- trimws(z)
          hit <- z[grepl(paste0("^", rank_prefix, "__"), z)]
          if (length(hit) == 0) {
            return(NA_character_)
          }
          clean_tax_value(hit[1])
        },
        FUN.VALUE = character(1)
      )
    }
    
    read_feature_table <- function(file) {
      header <- readLines(file, n = 50L, warn = FALSE)
      if (!length(header)) stop("Empty ASV table: ", file)
      meta_ids <- as.character(meta$fraction_sample_id)
      fields <- strsplit(sub("^\ufeff", "", header), "\t", fixed = TRUE)
      match_count <- vapply(fields, function(z) sum(trimws(z) %in% meta_ids), integer(1))
      if (!any(match_count > 0L)) stop("No metadata sample IDs found in the ASV header.")
      header_line <- which.max(match_count)
      data.table::fread(file, sep = "\t", header = TRUE, skip = header_line - 1L, check.names = FALSE)
    }
    
    meta <- fread(META_FILE, sep = "\t", header = TRUE, check.names = FALSE)
    setnames(meta, trimws(names(meta)))
    if ("OTU ID" %in% names(meta)) {
      setnames(meta, "OTU ID", "fraction_sample_id")
    } else if ("sample_id" %in% names(meta)) {
      setnames(meta, "sample_id", "fraction_sample_id")
    } else if (!"fraction_sample_id" %in% names(meta)) {
      setnames(meta, names(meta)[1], "fraction_sample_id")
    }
    
    required_meta_cols <- c(
      "fraction_sample_id", "fraction", "density", "abundance",
      "iso_trt", "fire_trt", "tem", "rep"
    )
    missing_meta_cols <- setdiff(required_meta_cols, names(meta))
    if (length(missing_meta_cols) > 0) {
      stop("Metadata missing columns: ", paste(missing_meta_cols, collapse = ", "))
    }
    
    meta[, fraction_sample_id := trimws(as.character(fraction_sample_id))]
    meta[, iso_trt := trimws(as.character(iso_trt))]
    meta[, fire_trt := trimws(as.character(fire_trt))]
    meta[, tem := trimws(as.character(tem))]
    meta[, rep := trimws(as.character(rep))]
    meta[, fraction := as.integer(fraction)]
    meta[, density := as.numeric(density)]
    meta[, abundance := as.numeric(abundance)]
    if (anyDuplicated(meta$fraction_sample_id)) {
      stop("Duplicate fraction_sample_id values detected.")
    }
    
    if (any(!meta$fire_trt %chin% FIRE_LEVELS)) {
      stop("Unexpected fire_trt values in metadata.")
    }
    
    if (any(!meta$tem %chin% TEMP_LEVELS)) {
      stop("Unexpected tem values in metadata.")
    }
    
    if (any(!meta$iso_trt %chin% ISO_LEVELS)) {
      stop("Unexpected iso_trt values in metadata.")
    }
    
    if (any(!meta$rep %chin% REP_LEVELS)) {
      stop("Unexpected rep values in metadata.")
    }
    
    if (any(!is.finite(meta$density)) || any(!is.finite(meta$abundance)) || any(meta$abundance < 0)) {
      stop("Invalid density or qPCR abundance values in metadata.")
    }
    
    meta[, tube_id := paste(fire_trt, tem, iso_trt, rep, sep = "_")]
    meta[, comparison_id := paste(fire_trt, tem, sep = "_")]
    expected_tubes <- CJ(
      fire_trt = FIRE_LEVELS,
      tem = TEMP_LEVELS,
      iso_trt = ISO_LEVELS,
      rep = REP_LEVELS,
      unique = TRUE
    )
    observed_tubes <- unique(meta[, .(fire_trt, tem, iso_trt, rep)])
    missing_tubes <- expected_tubes[!observed_tubes, on = .(fire_trt, tem, iso_trt, rep)]
    if (nrow(missing_tubes) > 0 || nrow(observed_tubes) != 36L) {
      if (nrow(missing_tubes) > 0) {
        fwrite(missing_tubes, file.path(output_dir, "QC_missing_expected_tubes.tsv"), sep = "\t")
      }
      stop("Complete 36-tube factorial design validation failed.")
    }
    
    taxonomy <- fread(TAXONOMY_FILE, sep = "\t", header = TRUE, check.names = FALSE)
    setnames(taxonomy, trimws(names(taxonomy)))
    candidate_tax_id_cols <- c("Feature ID", "OTU_ID", "OTU ID", "#OTU ID", "ASV", "Zotu")
    tax_id_col <- candidate_tax_id_cols[candidate_tax_id_cols %chin% names(taxonomy)][1]
    if (is.na(tax_id_col) || length(tax_id_col) == 0) {
      tax_id_col <- names(taxonomy)[1]
    }
    
    setnames(taxonomy, tax_id_col, "asv_id")
    taxonomy[, asv_id := trimws(as.character(asv_id))]
    candidate_tax_string_cols <- c("Taxon", "Taxonomy", "taxonomy", "Consensus.Lineage")
    tax_string_col <- candidate_tax_string_cols[candidate_tax_string_cols %chin% names(taxonomy)][1]
    if (!is.na(tax_string_col) && length(tax_string_col) > 0) {
      taxonomy[, taxonomy_string := trimws(as.character(get(tax_string_col)))]
    } else {
      rank_cols_available <- intersect(
        c("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
        names(taxonomy)
      )
      if (length(rank_cols_available) == 0) {
        stop("Taxonomy file has neither a taxonomy string nor recognizable rank columns.")
      }
      taxonomy[, taxonomy_string := do.call( paste, c(.SD, sep = ";") ), .SDcols = rank_cols_available]
    }
    
    domain_text <- rep(NA_character_, nrow(taxonomy))
    if ("Domain" %in% names(taxonomy)) {
      domain_text <- as.character(taxonomy$Domain)
    } else if ("Kingdom" %in% names(taxonomy)) {
      domain_text <- as.character(taxonomy$Kingdom)
    }
    
    is_bacteria <- grepl("(^|;|\\s)(d__|k__)?Bacteria($|;|\\s)", taxonomy$taxonomy_string, ignore.case = TRUE) |
      grepl("Bacteria", domain_text, ignore.case = TRUE)
    is_archaea <- grepl("Archaea", taxonomy$taxonomy_string, ignore.case = TRUE)
    is_eukaryota <- grepl("Eukary", taxonomy$taxonomy_string, ignore.case = TRUE)
    is_mitochondria <- grepl("mitochond", taxonomy$taxonomy_string, ignore.case = TRUE)
    is_chloroplast <- grepl("chloroplast|chlorophy", taxonomy$taxonomy_string, ignore.case = TRUE)
    is_bacteria[is.na(is_bacteria)] <- FALSE
    is_archaea[is.na(is_archaea)] <- FALSE
    is_eukaryota[is.na(is_eukaryota)] <- FALSE
    is_mitochondria[is.na(is_mitochondria)] <- FALSE
    is_chloroplast[is.na(is_chloroplast)] <- FALSE
    is_unassigned <- is.na(taxonomy$taxonomy_string) | trimws(taxonomy$taxonomy_string) == "" |
      grepl("^unassigned$", trimws(taxonomy$taxonomy_string), ignore.case = TRUE) |
      !is_bacteria
    rank_prefix_map <- c(
      Domain = "d",
      Phylum = "p",
      Class = "c",
      Order = "o",
      Family = "f",
      Genus = "g",
      Species = "s"
    )
    for (rank_name in names(rank_prefix_map)) {
      if (rank_name %in% names(taxonomy)) {
        taxonomy[, (rank_name) := clean_tax_value(get(rank_name))]
      } else {
        taxonomy[, (rank_name) := extract_rank_from_taxonomy(taxonomy_string, rank_prefix_map[[rank_name]])]
      }
    }
    
    taxonomy[, keep_target_bacteria := is_bacteria & !is_archaea & !is_eukaryota & !is_mitochondria & !is_chloroplast]
    taxonomy[, removal_reason := fifelse(
      is_archaea, "Archaea",
      fifelse(
        is_eukaryota, "Eukaryota",
        fifelse(
          is_mitochondria, "Mitochondria",
          fifelse(
            is_chloroplast, "Chloroplast",
            fifelse(!is_bacteria, "Unassigned_or_non_Bacteria", "Keep_target_Bacteria")
          )
        )
      )
    )]
    taxonomy_filter_summary <- taxonomy[, .N, by = .(removal_reason, keep_target_bacteria)]
    setorder(taxonomy_filter_summary, keep_target_bacteria, removal_reason)
    fwrite(taxonomy_filter_summary, file.path(output_dir, "02_taxonomy_filter_summary.tsv"), sep = "\t")
    fwrite(
      taxonomy[keep_target_bacteria == FALSE, .(asv_id, taxonomy_string, removal_reason)],
      file.path(output_dir, "03_taxonomy_removed_ASVs.tsv"),
      sep = "\t"
    )
    tax_target <- taxonomy[
      keep_target_bacteria == TRUE,
      .(asv_id, Domain, Phylum, Class, Order, Family, Genus, Species, taxonomy_string)
    ]
    tax_target[is.na(Phylum) | Phylum == "", Phylum := "Unclassified"]
    
    # Align bacterial counts with metadata and taxonomy.
    feature <- read_feature_table(FEATURE_FILE)
    setnames(feature, trimws(names(feature)))
    feature_id_col <- names(feature)[1]
    setnames(feature, feature_id_col, "asv_id")
    feature[, asv_id := trimws(as.character(asv_id))]
    sample_cols <- setdiff(names(feature), "asv_id")
    meta_samples <- meta$fraction_sample_id
    samples_missing_from_feature <- setdiff(meta_samples, sample_cols)
    feature_samples_not_in_meta <- setdiff(sample_cols, meta_samples)
    sample_matching_qc <- data.table(
      Metric = c(
        "metadata_fraction_samples",
        "feature_table_sample_columns",
        "metadata_samples_missing_from_feature",
        "feature_samples_not_in_metadata"
      ),
      Value = c(
        length(meta_samples),
        length(sample_cols),
        length(samples_missing_from_feature),
        length(feature_samples_not_in_meta)
      )
    )
    fwrite(sample_matching_qc, file.path(output_dir, "04_sample_matching_QC.tsv"), sep = "\t")
    if (length(samples_missing_from_feature) > 0 || length(feature_samples_not_in_meta) > 0) {
      if (length(samples_missing_from_feature) > 0) {
        fwrite(data.table(sample_id = samples_missing_from_feature), file.path(output_dir, "QC_metadata_samples_missing_from_feature.tsv"), sep = "\t")
      }
      if (length(feature_samples_not_in_meta) > 0) {
        fwrite(data.table(sample_id = feature_samples_not_in_meta), file.path(output_dir, "QC_feature_samples_not_in_metadata.tsv"), sep = "\t")
      }
      stop("Feature-table samples and metadata fraction samples do not match exactly.")
    }
    
    missing_taxonomy <- setdiff(feature$asv_id, taxonomy$asv_id)
    if (length(missing_taxonomy) > 0) {
      fwrite(data.table(asv_id = missing_taxonomy), file.path(output_dir, "QC_ASVs_missing_taxonomy.tsv"), sep = "\t")
      stop("Feature table contains ASVs without taxonomy. See QC_ASVs_missing_taxonomy.tsv.")
    }
    
    feature_target <- feature[asv_id %chin% tax_target$asv_id]
    if (nrow(feature_target) == 0) {
      stop("No target bacterial ASVs remained after taxonomy filtering.")
    }
    
    feature_target <- feature_target[, c("asv_id", meta_samples), with = FALSE]
    count_matrix_dense <- as.matrix(feature_target[, ..meta_samples])
    storage.mode(count_matrix_dense) <- "double"
    rownames(count_matrix_dense) <- feature_target$asv_id
    if (any(!is.finite(count_matrix_dense)) || any(count_matrix_dense < 0)) {
      stop("Feature table contains invalid or negative counts.")
    }
    
    bacterial_reads <- colSums(count_matrix_dense)
    if (any(bacterial_reads <= 0)) {
      bad_samples <- names(bacterial_reads)[bacterial_reads <= 0]
      fwrite(data.table(sample_id = bad_samples), file.path(output_dir, "QC_zero_target_bacterial_reads.tsv"), sep = "\t")
      stop("Some fraction samples have zero target-bacterial reads after filtering.")
    }
    
    fraction_depth <- data.table(
      fraction_sample_id = names(bacterial_reads),
      bacterial_reads = as.numeric(bacterial_reads)
    )
    fraction_depth <- merge(
      fraction_depth,
      meta[, .(fraction_sample_id, fire_trt, tem, iso_trt, rep, fraction, density, abundance, tube_id, comparison_id)],
      by = "fraction_sample_id",
      all.x = TRUE,
      sort = FALSE
    )
    fwrite(fraction_depth, file.path(output_dir, "05_fraction_bacterial_read_depth.tsv"), sep = "\t")
    count_matrix_sparse <- Matrix(count_matrix_dense, sparse = TRUE)
    rm(count_matrix_dense)
    gc()
    sparse_summary <- summary(count_matrix_sparse)
    long <- data.table(
      asv_id = rownames(count_matrix_sparse)[sparse_summary$i],
      fraction_sample_id = colnames(count_matrix_sparse)[sparse_summary$j],
      seq_abund = as.numeric(sparse_summary$x)
    )
    rm(count_matrix_sparse, sparse_summary)
    gc()
    long <- merge( long, fraction_depth, by = "fraction_sample_id", all.x = TRUE, sort = FALSE )
    if (anyNA(long$tube_id) || anyNA(long$bacterial_reads)) {
      stop("Failed to join sparse feature counts to metadata/read-depth table.")
    }
    
    long[, rel_abund := seq_abund / bacterial_reads]
    long[, taxon_16S_fraction := rel_abund * abundance]
    freq_by_tube <- long[
      seq_abund > 0,
      .( n_fraction_detected = .N, total_reads_ASV_tube = sum(seq_abund) ),
      by = .(asv_id, fire_trt, tem, iso_trt, rep, tube_id, comparison_id)
    ]
    freq_summary <- freq_by_tube[
      ,
      .(
        ASV_tube_combinations = .N,
        pass_min_fraction = sum(n_fraction_detected >= MIN_FRACTIONS),
        fail_min_fraction = sum(n_fraction_detected < MIN_FRACTIONS)
      ),
      by = .(fire_trt, tem, iso_trt)
    ]
    fwrite(freq_summary, file.path(output_dir, "06_fraction_frequency_filter_summary.tsv"), sep = "\t")
    long <- merge(
      long,
      freq_by_tube[, .(asv_id, tube_id, n_fraction_detected)],
      by = c("asv_id", "tube_id"),
      all.x = TRUE,
      sort = FALSE
    )
    long_pass <- long[n_fraction_detected >= MIN_FRACTIONS]
    wads <- long_pass[
      ,
      .(
        wad_numerator = sum(density * taxon_16S_fraction),
        asv_abundance = sum(taxon_16S_fraction),
        seq_reads_tube = sum(seq_abund),
        n_fraction_detected = unique(n_fraction_detected)
      ),
      by = .(asv_id, fire_trt, tem, iso_trt, rep, tube_id, comparison_id)
    ]
    wads <- wads[asv_abundance > 0]
    wads[, wad := wad_numerator / asv_abundance]
    wads[, wad_numerator := NULL]
    long_pass <- merge(
      long_pass,
      wads[, .(asv_id, tube_id, wad)],
      by = c("asv_id", "tube_id"),
      all = FALSE,
      sort = FALSE
    )
    wvd_table <- long_pass[
      ,
      .( wvd = sum(((density - unique(wad))^2) * taxon_16S_fraction) / sum(taxon_16S_fraction) ),
      by = .(asv_id, tube_id)
    ]
    wads <- merge(wads, wvd_table, by = c("asv_id", "tube_id"), all.x = TRUE, sort = FALSE)
    wads <- merge(wads, tax_target, by = "asv_id", all.x = TRUE, sort = FALSE)
    fwrite(
      wads[order(fire_trt, tem, iso_trt, rep, asv_id)],
      file.path(output_dir, "08_WAD_per_ASV_per_tube.tsv"),
      sep = "\t"
    )
    
    # Require valid WAD estimates in all three gradients per isotope.
    eligibility <- wads[
      ,
      .(
        n_unlabel_wad = uniqueN(rep[iso_trt == "unlabel" & is.finite(wad)]),
        n_label_wad = uniqueN(rep[iso_trt == "label" & is.finite(wad)]),
        unlabel_reps = paste(sort(unique(rep[iso_trt == "unlabel" & is.finite(wad)])), collapse = ","),
        label_reps = paste(sort(unique(rep[iso_trt == "label" & is.finite(wad)])), collapse = ",")
      ),
      by = .(asv_id, fire_trt, tem, comparison_id)
    ]
    eligibility[, AFE_testable :=
                  n_unlabel_wad == MIN_REPS & n_label_wad == MIN_REPS & unlabel_reps == paste(REP_LEVELS, collapse = ",") &
                  label_reps == paste(REP_LEVELS, collapse = ",")
    ]
    eligibility <- merge( eligibility, tax_target[, .(asv_id, Phylum)], by = "asv_id", all.x = TRUE, sort = FALSE )
    eligibility_summary <- eligibility[
      ,
      .( ASVs_with_any_WAD = .N, qSIP_testable_ASVs = sum(AFE_testable), not_testable_ASVs = sum(!AFE_testable) ),
      by = .(fire_trt, tem, comparison_id)
    ]
    fwrite(eligibility, file.path(output_dir, "09_ASV_AFE_eligibility.tsv"), sep = "\t")
    fwrite(eligibility_summary, file.path(output_dir, "07_ASV_AFE_eligibility_summary.tsv"), sep = "\t")
    
    # Independently bootstrap labeled and control gradients.
    afe_result_list <- list()
    correction_list <- list()
    bootstrap_manifest_list <- list()
    for (g in seq_len(nrow(comparison_table))) {
      current_fire <- comparison_table$fire_trt[g]
      current_tem <- comparison_table$tem[g]
      current_id <- comparison_table$comparison_id[g]
      message("AFE bootstrap: ", current_id)
      eligible_asv <- eligibility[ fire_trt == current_fire & tem == current_tem & AFE_testable == TRUE, asv_id ]
      if (length(eligible_asv) == 0) {
        stop("No qSIP-testable ASVs in ", current_id, ".")
      }
      
      sub_wad <- wads[
        fire_trt == current_fire & tem == current_tem & asv_id %chin% eligible_asv,
        .(asv_id, iso_trt, rep, wad)
      ]
      unlabel_wide <- dcast( sub_wad[iso_trt == "unlabel", .(asv_id, rep, wad)], asv_id ~ rep, value.var = "wad" )
      label_wide <- dcast( sub_wad[iso_trt == "label", .(asv_id, rep, wad)], asv_id ~ rep, value.var = "wad" )
      rep_cols_unlabel <- setdiff(names(unlabel_wide), "asv_id")
      rep_cols_label <- setdiff(names(label_wide), "asv_id")
      rep_cols_unlabel <- as.character(sort(as.integer(rep_cols_unlabel)))
      rep_cols_label <- as.character(sort(as.integer(rep_cols_label)))
      if (!identical(rep_cols_unlabel, REP_LEVELS) || !identical(rep_cols_label, REP_LEVELS)) {
        stop("Expected exact rep columns 1,2,3 for both isotope treatments in ", current_id, ".")
      }
      
      setkey(unlabel_wide, asv_id)
      setkey(label_wide, asv_id)
      ordered_asv <- sort(eligible_asv)
      unlabel_wide <- unlabel_wide[J(ordered_asv)]
      label_wide <- label_wide[J(ordered_asv)]
      if (anyNA(unlabel_wide$asv_id) || anyNA(label_wide$asv_id)) {
        stop("Failed to align eligible ASVs across isotope treatments in ", current_id, ".")
      }
      
      Wunlabel <- as.matrix(unlabel_wide[, ..REP_LEVELS])
      Wlabel <- as.matrix(label_wide[, ..REP_LEVELS])
      storage.mode(Wunlabel) <- "double"
      storage.mode(Wlabel) <- "double"
      if (anyNA(Wunlabel) || anyNA(Wlabel)) {
        stop("NA WAD values remain after strict eligibility filtering in ", current_id, ".")
      }
      
      Wunlabel_mean <- rowMeans(Wunlabel)
      Wlabel_mean <- rowMeans(Wlabel)
      AFE_raw_observed <- calc_AFE_18O(Wunlabel_mean, Wlabel_mean)
      observed_shift <- bottom_non_grower_shift(AFE_raw_observed, NON_GROWER_PROP)
      AFE_corrected_observed <- AFE_raw_observed - observed_shift
      set.seed(RANDOM_SEED + g)
      unlabel_weights <- make_boot_weights(n = length(REP_LEVELS), n_boot = N_BOOT)
      label_weights <- make_boot_weights(n = length(REP_LEVELS), n_boot = N_BOOT)
      Wunlabel_boot <- Wunlabel %*% unlabel_weights
      Wlabel_boot <- Wlabel %*% label_weights
      AFE_boot_raw <- calc_AFE_18O(Wunlabel_boot, Wlabel_boot)
      correction_shift_by_iteration <- apply( AFE_boot_raw, 2, bottom_non_grower_shift, prop = NON_GROWER_PROP )
      AFE_boot_corrected <- sweep( AFE_boot_raw, MARGIN = 2, STATS = correction_shift_by_iteration, FUN = "-" )
      raw_q <- matrixStats::rowQuantiles(AFE_boot_raw, probs = CI_PROBS, na.rm = TRUE, drop = FALSE)
      corrected_q <- matrixStats::rowQuantiles(AFE_boot_corrected, probs = CI_PROBS, na.rm = TRUE, drop = FALSE)
      current_result <- data.table(
        asv_id = ordered_asv,
        fire_trt = current_fire,
        tem = current_tem,
        comparison_id = current_id,
        WAD_unlabel_mean = Wunlabel_mean,
        WAD_label_mean = Wlabel_mean,
        delta_WAD = Wlabel_mean - Wunlabel_mean,
        AFE_raw_observed = AFE_raw_observed,
        AFE_corrected_observed = AFE_corrected_observed,
        AFE_raw_CI_low_90 = raw_q[, 1],
        AFE_raw_median = raw_q[, 2],
        AFE_raw_CI_high_90 = raw_q[, 3],
        AFE_CI_low_90 = corrected_q[, 1],
        AFE_ASV = corrected_q[, 2],
        AFE_CI_high_90 = corrected_q[, 3],
        bootstrap_prob_AFE_le_0 = rowMeans(AFE_boot_corrected <= 0, na.rm = TRUE),
        Active_ASV = corrected_q[, 1] > 0,
        Detected_negative_AFE = corrected_q[, 3] < 0,
        n_boot = N_BOOT,
        bootstrap_design = "independent_unlabel_and_label"
      )
      current_result <- merge( current_result, tax_target, by = "asv_id", all.x = TRUE, sort = FALSE )
      afe_result_list[[length(afe_result_list) + 1L]] <- current_result
      correction_list[[length(correction_list) + 1L]] <- data.table(
        fire_trt = current_fire,
        tem = current_tem,
        comparison_id = current_id,
        observed_bottom10_shift = observed_shift,
        bootstrap_shift_median = median(correction_shift_by_iteration),
        bootstrap_shift_CI_low_90 = quantile(correction_shift_by_iteration, 0.05, names = FALSE),
        bootstrap_shift_CI_high_90 = quantile(correction_shift_by_iteration, 0.95, names = FALSE),
        non_grower_prop = NON_GROWER_PROP,
        n_testable_ASV = length(ordered_asv)
      )
      bootstrap_object <- list(
        comparison_id = current_id,
        fire_trt = current_fire,
        tem = current_tem,
        bootstrap_design = "independent_unlabel_and_label",
        asv_id = ordered_asv,
        rep_levels = REP_LEVELS,
        n_boot = N_BOOT,
        random_seed = RANDOM_SEED + g,
        unlabel_weights = unlabel_weights,
        label_weights = label_weights,
        correction_shift_by_iteration = correction_shift_by_iteration,
        AFE_boot_corrected = AFE_boot_corrected
      )
      bootstrap_file <- file.path(bootstrap_dir, paste0("AFE_bootstrap_", current_id, ".rds"))
      saveRDS(bootstrap_object, bootstrap_file, compress = TRUE)
      bootstrap_manifest_list[[length(bootstrap_manifest_list) + 1L]] <- data.table(
        fire_trt = current_fire,
        tem = current_tem,
        comparison_id = current_id,
        n_testable_ASV = length(ordered_asv),
        n_boot = N_BOOT,
        bootstrap_design = "independent_unlabel_and_label",
        file = bootstrap_file
      )
      rm(
        Wunlabel, Wlabel, Wunlabel_boot, Wlabel_boot,
        AFE_boot_raw, AFE_boot_corrected, raw_q, corrected_q,
        unlabel_weights, label_weights
      )
      gc()
    }
    
    afe_all <- rbindlist(afe_result_list, use.names = TRUE, fill = TRUE)
    correction_summary <- rbindlist(correction_list, use.names = TRUE, fill = TRUE)
    bootstrap_manifest <- rbindlist(bootstrap_manifest_list, use.names = TRUE, fill = TRUE)
    TOP_N_PHYLA <- 10L
    phylum_counts <- afe_all[ , .(qSIP_testable_ASVs = .N), by = Phylum ][order(-qSIP_testable_ASVs, Phylum)]
    top_phyla <- phylum_counts[seq_len(min(TOP_N_PHYLA, .N)), Phylum]
    afe_all[, Phylum_group := fifelse(Phylum %chin% top_phyla, Phylum, "Others")]
    phylum_levels <- c(top_phyla, "Others")
    afe_all[, Phylum_group := factor(Phylum_group, levels = phylum_levels)]
    known_colors <- c(
      Acidobacteriota = "#BC75AF",
      Actinobacteriota = "#E96B73",
      Proteobacteria = "#6FA8DC",
      Alphaproteobacteria = "#C98E73",
      Gammaproteobacteria = "#F2A23A",
      Bacteroidota = "#73A9C2",
      Firmicutes = "#C8B991",
      Bacillota = "#C8B991",
      Gemmatimonadota = "#8E74C9",
      Myxococcota = "#8B4513",
      Planctomycetota = "#9B8DB6",
      Verrucomicrobiota = "#245316",
      Chloroflexi = "#7FA776",
      Patescibacteria = "#A8CB98",
      Bdellovibrionota = "#A7A437",
      Armatimonadota = "#C987C6",
      `WPS-2` = "#D65158",
      Unclassified = "#666666",
      Others = "#999999"
    )
    missing_palette <- setdiff(phylum_levels, names(known_colors))
    if (length(missing_palette) > 0) {
      extra_colors <- grDevices::hcl.colors(length(missing_palette), palette = "Dark 3")
      names(extra_colors) <- missing_palette
      known_colors <- c(known_colors, extra_colors)
    }
    
    phylum_colors <- known_colors[phylum_levels]
    palette_table <- data.table(
      Phylum_group = phylum_levels,
      color = unname(phylum_colors),
      display_order = seq_along(phylum_levels)
    )
    afe_group_summary <- afe_all[
      ,
      .(
        qSIP_testable_ASVs = .N,
        Active_ASVs = sum(Active_ASV),
        Active_percent = 100 * mean(Active_ASV),
        Detected_negative_ASVs = sum(Detected_negative_AFE),
        median_AFE_all = median(AFE_ASV),
        median_AFE_active = if (any(Active_ASV)) median(AFE_ASV[Active_ASV]) else NA_real_
      ),
      by = .(fire_trt, tem, comparison_id)
    ]
    fwrite(correction_summary, file.path(output_dir, "10_non_grower_correction_summary.tsv"), sep = "\t")
    fwrite(afe_group_summary, file.path(output_dir, "11_AFE_ASV_group_summary.tsv"), sep = "\t")
    fwrite(afe_all, file.path(output_dir, "12_AFE_ASV_all.tsv"), sep = "\t")
    fwrite(afe_all[Active_ASV == TRUE], file.path(output_dir, "13_AFE_ASV_active.tsv"), sep = "\t")
    fwrite(afe_all[Detected_negative_AFE == TRUE], file.path(output_dir, "14_QC_detected_negative_AFE.tsv"), sep = "\t")
    fwrite(palette_table, file.path(output_dir, "15_phylum_palette.tsv"), sep = "\t")
    fwrite(bootstrap_manifest, file.path(output_dir, "16_AFE_bootstrap_object_manifest.tsv"), sep = "\t")
    
    prepare_panel_data <- function(panel_dt, gap = 8L) {
      x <- copy(panel_dt)
      setorderv( x, cols = c("Phylum_group", "AFE_ASV", "asv_id"), order = c(1, 1, 1) )
      x[, within_block_rank := seq_len(.N), by = Phylum_group]
      block <- x[ , .( n_total = .N, n_active = sum(Active_ASV) ), by = Phylum_group ]
      setorder(block, Phylum_group)
      block[, block_offset := {
        if (.N == 1) {
          0L
        } else {
          c(0L, cumsum(head(n_total + gap, -1)))
        }
      }]
      x <- merge(x, block[, .(Phylum_group, block_offset)], by = "Phylum_group", all.x = TRUE, sort = FALSE)
      setorderv(x, c("Phylum_group", "AFE_ASV", "asv_id"), c(1, 1, 1))
      x[, within_block_rank := seq_len(.N), by = Phylum_group]
      x[, y_pos := block_offset + within_block_rank]
      block_summary <- x[
        ,
        .(
          y_min = min(y_pos),
          y_max = max(y_pos),
          y_mid = (min(y_pos) + max(y_pos)) / 2,
          n_total = .N,
          n_active = sum(Active_ASV)
        ),
        by = Phylum_group
      ]
      block_summary[, active_percent := 100 * n_active / n_total]
      block_summary[, label_text := sprintf("%d/%d (%.1f%%)", n_active, n_total, active_percent)]
      setorder(block_summary, y_min)
      separator <- copy(block_summary)
      if (nrow(separator) > 1) {
        separator <- separator[1:(.N - 1)]
        separator[, y_sep := y_max + gap / 2]
      } else {
        separator <- separator[0]
      }
      
      list(data = x, block = block_summary, separator = separator)
    }
    
    AFE_x_min_raw <- min(afe_all$AFE_CI_low_90, na.rm = TRUE)
    AFE_x_max_raw <- max(afe_all$AFE_CI_high_90, na.rm = TRUE)
    AFE_x_span <- AFE_x_max_raw - AFE_x_min_raw
    if (!is.finite(AFE_x_span) || AFE_x_span <= 0) {
      AFE_x_span <- 0.05
    }
    
    AFE_label_x <- AFE_x_max_raw + 0.04 * AFE_x_span
    AFE_x_limits <- c( min(0, AFE_x_min_raw) - 0.05 * AFE_x_span, AFE_x_max_raw + 0.34 * AFE_x_span )
    panel_plot_list <- list()
    panel_data_list <- list()
    panel_summary_list <- list()
    block_summary_list <- list()
    for (g in seq_len(nrow(comparison_table))) {
      current_fire <- comparison_table$fire_trt[g]
      current_tem <- comparison_table$tem[g]
      current_id <- comparison_table$comparison_id[g]
      panel_dt <- afe_all[fire_trt == current_fire & tem == current_tem]
      prep <- prepare_panel_data(panel_dt, gap = 8L)
      plot_dt <- prep$data
      block_dt <- prep$block
      separator_dt <- prep$separator
      n_total <- nrow(plot_dt)
      n_active <- sum(plot_dt$Active_ASV)
      active_percent <- 100 * n_active / n_total
      panel_summary_list[[length(panel_summary_list) + 1L]] <- data.table(
        fire_trt = current_fire,
        tem = current_tem,
        comparison_id = current_id,
        qSIP_testable_ASVs = n_total,
        Active_ASVs = n_active,
        Active_percent = active_percent
      )
      tmp_plot <- copy(plot_dt)
      tmp_plot[, comparison_id := current_id]
      panel_data_list[[length(panel_data_list) + 1L]] <- tmp_plot
      tmp_block <- copy(block_dt)
      tmp_block[, comparison_id := current_id]
      block_summary_list[[length(block_summary_list) + 1L]] <- tmp_block
      title_text <- sprintf(
        "%s | %s\nActive/testable ASVs: %d/%d (%.1f%%)",
        FIRE_LABELS[[current_fire]],
        TEMP_LABELS[[current_tem]],
        n_active,
        n_total,
        active_percent
      )
      p <- ggplot() + geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.55, colour = "grey40")
      if (nrow(separator_dt) > 0) {
        p <- p + geom_hline(
          data = separator_dt,
          aes(yintercept = y_sep),
          inherit.aes = FALSE,
          linewidth = 0.25,
          colour = "grey86"
        )
      }
      
      p <- p +
        geom_segment(
          data = plot_dt,
          aes( x = AFE_CI_low_90, xend = AFE_CI_high_90, y = y_pos, yend = y_pos, colour = Phylum_group ),
          linewidth = 0.20,
          alpha = 0.48
        ) +
        geom_point(
          data = plot_dt[Active_ASV == FALSE],
          aes(x = AFE_ASV, y = y_pos, colour = Phylum_group),
          shape = 1,
          size = 0.55,
          stroke = 0.35,
          alpha = 0.75
        ) +
        geom_point(
          data = plot_dt[Active_ASV == TRUE],
          aes(x = AFE_ASV, y = y_pos, colour = Phylum_group),
          shape = 16,
          size = 0.58,
          alpha = 0.95
        ) +
        geom_text(
          data = block_dt,
          aes(x = AFE_label_x, y = y_mid, label = label_text, colour = Phylum_group),
          hjust = 0,
          size = 2.35,
          show.legend = FALSE
        ) +
        scale_colour_manual(values = phylum_colors, drop = FALSE) +
        coord_cartesian(xlim = AFE_x_limits, clip = "off") +
        labs(
          title = title_text,
          x = expression(paste("Atom fraction excess (", {}^{18}, "O)")),
          y = "Ranked qSIP-testable ASVs",
          colour = "Phylum"
        ) +
        theme_classic(base_size = 11) +
        theme(
          plot.title = element_text(size = 10.5, face = "bold", hjust = 0.5),
          axis.title = element_text(size = 11.5, face = "bold"),
          axis.text.x = element_text(size = 8.5, colour = "black"),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.65),
          legend.position = "right",
          legend.title = element_text(size = 10, face = "bold"),
          legend.text = element_text(size = 8.5),
          plot.margin = margin(8, 64, 8, 8)
        )
      panel_plot_list[[current_id]] <- p
    }
    
    AFE_combined <- wrap_plots(panel_plot_list, ncol = 2, guides = "collect") & theme(legend.position = "right")
    save_figure(
      filename = file.path(output_dir, "Fig_AFE_ASV_forest_6panel.pdf"),
      plot = AFE_combined,
      width = 14,
      height = 18,
      units = "in",
      limitsize = FALSE
    )
    fwrite(rbindlist(panel_data_list, use.names = TRUE, fill = TRUE), file.path(output_dir, "17_AFE_forest_plot_data.tsv"), sep = "\t")
    fwrite(rbindlist(panel_summary_list, use.names = TRUE, fill = TRUE), file.path(output_dir, "18_AFE_forest_panel_summary.tsv"), sep = "\t")
    fwrite(rbindlist(block_summary_list, use.names = TRUE, fill = TRUE), file.path(output_dir, "19_AFE_forest_phylum_block_summary.tsv"), sep = "\t")
    parameter_table <- data.table(
      Parameter = c(
        "Fire stages",
        "Temperatures",
        "Isotope treatments",
        "Biological replicates",
        "Minimum detected fractions per ASV per tube",
        "Required valid unlabel replicates per comparison",
        "Required valid label replicates per comparison",
        "Bootstrap design",
        "Bootstrap iterations",
        "Bootstrap confidence interval",
        "Non-grower correction",
        "Primary active ASV criterion",
        "Natural abundance 18O",
        "Multiple-testing correction",
        "Forest plot phylum grouping"
      ),
      Value = c(
        paste(FIRE_LEVELS, collapse = ","),
        paste(TEMP_LEVELS, collapse = ","),
        paste(ISO_LEVELS, collapse = ","),
        paste(REP_LEVELS, collapse = ","),
        as.character(MIN_FRACTIONS),
        as.character(MIN_REPS),
        as.character(MIN_REPS),
        "independent_unlabel_and_label",
        as.character(N_BOOT),
        "90% (5th-95th percentile)",
        "median of bottom 10% AFE shifted to zero within each exact fire x temperature comparison",
        "AFE_CI_low_90 > 0",
        as.character(NATURAL_18O),
        "None",
        paste0("Global top ", TOP_N_PHYLA, " phyla among qSIP-testable ASVs + Others")
      )
    )
    fwrite(parameter_table, file.path(output_dir, "20_AFE_analysis_parameters.tsv"), sep = "\t")
    invisible(NULL)
  }
  
  # Cue
  run_cue <- function() {
    output_dir <- CACHE$cue
    df1 <- fread(INPUTS$macro_absolute, data.table = FALSE)
    fwrite(df1, file.path(output_dir, "Raw_absolute_observations.tsv"), sep = "\t")
    df_long1 <- df1 %>%
      pivot_longer(cols = c("Respiration", "Growth", "CUE"), names_to = "Metric", values_to = "Value")
    summary_df1 <- df_long1 %>% group_by(fire_trt, Tem_trt, Metric) %>%
      summarise( Mean = mean(Value), SE = sd(Value) / sqrt(n()), .groups = 'drop' )
    letter_list1 <- list()
    metrics1 <- unique(df_long1$Metric)
    temps1 <- unique(df_long1$Tem_trt)
    for (m in metrics1) {
      for (t in temps1) {
        sub_data <- df_long1 %>% filter(Metric == m, Tem_trt == t)
        model <- aov(Value ~ fire_trt, data = sub_data)
        tukey <- HSD.test(model, "fire_trt", group = TRUE)
        fwrite(as.data.table(summary(model)[[1]], keep.rownames = "term"),
               file.path(output_dir, paste0("ANOVA_absolute_", m, "_", t, ".tsv")), sep = "\t")
        fwrite(as.data.table(tukey$groups, keep.rownames = "fire_trt"),
               file.path(output_dir, paste0("Tukey_", m, "_", t, ".tsv")), sep = "\t")
        letters_df <- data.frame(
          fire_trt = rownames(tukey$groups),
          Tem_trt = t,
          Metric = m,
          Letter = tukey$groups$groups,
          stringsAsFactors = FALSE
        )
        if (t == "T_low") {
          letters_df$Letter <- tolower(letters_df$Letter)
        } else {
          letters_df$Letter <- toupper(letters_df$Letter)
        }
        letter_list1[[length(letter_list1) + 1]] <- letters_df
      }
    }
    
    all_letters1 <- bind_rows(letter_list1)
    plot_data1 <- left_join(summary_df1, all_letters1, by = c("fire_trt", "Tem_trt", "Metric"))
    plot_data1$fire_trt <- factor(plot_data1$fire_trt, levels = c("UB", "F3M", "F3Y"), labels = c("Unburned", "3 months", "3 years"))
    plot_data1$Tem_trt <- factor(plot_data1$Tem_trt, levels = c("T_low", "T_high"), labels = c(unname(TEMP_LABELS[1L]), unname(TEMP_LABELS[2L])))
    plot_data1$Metric <- factor(plot_data1$Metric, levels = c("Respiration", "Growth", "CUE"))
    plot_data1 <- plot_data1 %>% group_by(Metric) %>%
      mutate(y_max = max(Mean + SE, na.rm = TRUE), letter_y = Mean + SE + (y_max * 0.05)) %>%
      ungroup()
    p_abs <- ggplot(plot_data1, aes(x = fire_trt, y = Mean, fill = Tem_trt)) +
      geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black", linewidth = 0.6) +
      geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), position = position_dodge(width = 0.8), width = 0.15, linewidth = 0.6) +
      geom_text(aes(y = letter_y, label = Letter), position = position_dodge(width = 0.8), size = 5, vjust = 0, color = "black") +
      facet_wrap(~ Metric, scales = "free_y", ncol = 3) +
      scale_fill_manual(values = TEMP_COLORS_DISPLAY) +
      labs(x = NULL, y = "Values", fill = "Incubation temperature") +
      scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
      theme_classic(base_size = 14) +
      theme(
        legend.position = "top",
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 13),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        axis.line = element_blank(),
        axis.text.x = element_text(size = 13, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 14, face = "bold"),
        strip.background = element_rect(fill = "grey95", color = "black", linewidth = 1),
        strip.text = element_text(size = 14, face = "bold", margin = margin(t=8, b=8))
      )
    save_figure(file.path(output_dir, "01_Fig_SY2_Absolute_Metrics_Tukey.pdf"), plot = p_abs, width = 12, height = 4.5, units = "in", device = cairo_pdf)
    fwrite(plot_data1, file.path(output_dir, "01_Table_SY2_Absolute_Stats.tsv"), sep = "\t")
    dt2 <- fread(INPUTS$macro_delta)
    if (ncol(dt2) != 5L) stop("SY.txt must have exactly five columns: Delta_Respiration, Delta_Growth, Delta_CUE, fire_trt, rep.")
    setnames(dt2, 1:5, c("Delta_Respiration", "Delta_Growth", "Delta_CUE", "fire_trt", "rep"))
    fwrite(dt2, file.path(output_dir, "Raw_delta_observations.tsv"), sep = "\t")
    params2 <- c("Delta_Respiration", "Delta_Growth", "Delta_CUE")
    
    get_stats_duncan <- function(data, param) {
      sub_dt <- data[, .(fire_trt = as.character(fire_trt), value = get(param))]
      sub_dt <- sub_dt[!is.na(value)]
      summary_dt <- sub_dt[, .(Mean = mean(value), SE = sd(value) / sqrt(.N)), by = fire_trt]
      aov_res <- aov(value ~ fire_trt, data = sub_dt)
      duncan_res <- duncan.test(aov_res, "fire_trt")
      fwrite(as.data.table(summary(aov_res)[[1]], keep.rownames = "term"),
             file.path(output_dir, paste0("ANOVA_delta_", param, ".tsv")), sep = "\t")
      fwrite(as.data.table(duncan_res$groups, keep.rownames = "fire_trt"),
             file.path(output_dir, paste0("Duncan_", param, ".tsv")), sep = "\t")
      letter_dt <- as.data.table(duncan_res$groups, keep.rownames = "fire_trt")
      letter_dt <- letter_dt[, .(fire_trt, sig_letter = groups)]
      res <- merge(summary_dt, letter_dt, by = "fire_trt")
      res[, Parameter := param]
      return(res)
    }
    
    stat_list2 <- lapply(params2, function(p) get_stats_duncan(dt2, p))
    stat_dt2 <- rbindlist(stat_list2)
    stat_dt2[, fire_trt_label := factor(fire_trt, levels = FIRE_LEVELS, labels = FIRE_LABELS)]
    
    make_bar_plot <- function(data, param_name, y_label, plot_title) {
      plot_data <- data[Parameter == param_name]
      y_max <- max(plot_data$Mean + plot_data$SE, na.rm = TRUE) * 1.15
      ggplot(plot_data, aes(x = fire_trt_label, y = Mean, fill = fire_trt)) +
        geom_col(width = 0.6, color = "black", linewidth = 0.6) +
        geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.15, linewidth = 0.6) +
        geom_text(aes(y = Mean + SE + (y_max * 0.04), label = sig_letter),
                  size = 5.5, fontface = "bold", color = "black") +
        scale_fill_manual(values = FIRE_COLORS, guide = "none") +
        scale_y_continuous(expand = expansion(mult = c(0, 0))) +
        coord_cartesian(ylim = c(0, y_max)) +
        labs(x = NULL, y = y_label, title = plot_title) +
        theme_classic(base_size = 14) +
        theme(
          plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
          axis.title.y = element_text(face = "bold", size = 13),
          axis.text.x = element_text(color = "black", size = 12),
          axis.text.y = element_text(color = "black", size = 11),
          axis.line = element_line(linewidth = 0.6, color = "black"),
          axis.ticks = element_line(color = "black")
        )
    }
    
    p1 <- make_bar_plot(stat_dt2, "Delta_Respiration", expression(Delta * "Respiration"), expression(Delta * "Respiration"))
    p2 <- make_bar_plot(stat_dt2, "Delta_Growth", expression(Delta * "Growth"), expression(Delta * "Growth"))
    p3 <- make_bar_plot(stat_dt2, "Delta_CUE", expression(Delta * "CUE"), expression(Delta * "CUE"))
    combined_plot <- p1 + p2 + p3 + plot_layout(ncol = 3)
    save_figure(file.path(output_dir, "02_Fig_SY_Macro_Delta_Duncan.pdf"), plot = combined_plot, width = 12, height = 4.5, units = "in", device = cairo_pdf)
    fwrite(stat_dt2, file.path(output_dir, "02_Table_SY_Delta_Stats.tsv"), sep = "\t")
    invisible(NULL)
  }
  
  # Growth
  run_growth <- function() {
    AFE_FILE <- file.path(CACHE$afe, "12_AFE_ASV_all.tsv")
    WAD_FILE <- file.path(CACHE$afe, "08_WAD_per_ASV_per_tube.tsv")
    PALETTE_FILE <- file.path(CACHE$afe, "15_phylum_palette.tsv")
    AFE_BOOTSTRAP_DIR <- file.path(CACHE$afe, "bootstrap_objects")
    output_dir <- CACHE$growth
    growth_bootstrap_dir <- file.path(output_dir, "bootstrap_objects")
    if (!dir.exists(growth_bootstrap_dir)) {
      dir.create(growth_bootstrap_dir, recursive = TRUE)
    }
    
    required_files <- c(AFE_FILE, WAD_FILE, PALETTE_FILE)
    missing_files <- required_files[!file.exists(required_files)]
    if (length(missing_files) > 0) {
      stop("Missing required step-02 files:\n", paste(missing_files, collapse = "\n"))
    }
    
    comparison_table <- CJ(fire_trt = FIRE_LEVELS, tem = TEMP_LEVELS, unique = TRUE)
    comparison_table[, fire_order := match(fire_trt, FIRE_LEVELS)]
    comparison_table[, temp_order := match(tem, TEMP_LEVELS)]
    setorder(comparison_table, fire_order, temp_order)
    comparison_table[, comparison_id := paste(fire_trt, tem, sep = "_")]
    comparison_table[, c("fire_order", "temp_order") := NULL]
    comparison_table[, bootstrap_design := "independent_unlabel_and_label"]
    for (current_id in comparison_table$comparison_id) {
      bootstrap_file <- file.path(AFE_BOOTSTRAP_DIR, paste0("AFE_bootstrap_", current_id, ".rds"))
      if (!file.exists(bootstrap_file)) {
        stop("Missing step-02 bootstrap object: ", bootstrap_file)
      }
    }
    
    afe <- fread(AFE_FILE, sep = "\t", header = TRUE)
    wads <- fread(WAD_FILE, sep = "\t", header = TRUE)
    palette_table <- fread(PALETTE_FILE, sep = "\t", header = TRUE)
    required_afe_cols <- c(
      "asv_id", "fire_trt", "tem", "comparison_id",
      "AFE_ASV", "AFE_CI_low_90", "AFE_CI_high_90",
      "Active_ASV", "Phylum", "Phylum_group"
    )
    required_wad_cols <- c( "asv_id", "fire_trt", "tem", "iso_trt", "rep", "asv_abundance" )
    required_palette_cols <- c("Phylum_group", "color", "display_order")
    missing_afe_cols <- setdiff(required_afe_cols, names(afe))
    missing_wad_cols <- setdiff(required_wad_cols, names(wads))
    missing_palette_cols <- setdiff(required_palette_cols, names(palette_table))
    if (length(missing_afe_cols) > 0) {
      stop("AFE table missing columns: ", paste(missing_afe_cols, collapse = ", "))
    }
    
    if (length(missing_wad_cols) > 0) {
      stop("WAD table missing columns: ", paste(missing_wad_cols, collapse = ", "))
    }
    
    if (length(missing_palette_cols) > 0) {
      stop("Palette table missing columns: ", paste(missing_palette_cols, collapse = ", "))
    }
    
    afe[, asv_id := as.character(asv_id)]
    afe[, fire_trt := as.character(fire_trt)]
    afe[, tem := as.character(tem)]
    afe[, comparison_id := as.character(comparison_id)]
    afe[, AFE_ASV := as.numeric(AFE_ASV)]
    afe[, AFE_CI_low_90 := as.numeric(AFE_CI_low_90)]
    afe[, AFE_CI_high_90 := as.numeric(AFE_CI_high_90)]
    afe[, Active_ASV := as.logical(Active_ASV)]
    afe[, Phylum := as.character(Phylum)]
    afe[, Phylum_group := as.character(Phylum_group)]
    wads[, asv_id := as.character(asv_id)]
    wads[, fire_trt := as.character(fire_trt)]
    wads[, tem := as.character(tem)]
    wads[, iso_trt := as.character(iso_trt)]
    wads[, rep := as.character(rep)]
    wads[, asv_abundance := as.numeric(asv_abundance)]
    palette_table[, Phylum_group := as.character(Phylum_group)]
    palette_table[, color := as.character(color)]
    palette_table[, display_order := as.integer(display_order)]
    setorder(palette_table, display_order)
    phylum_levels <- palette_table$Phylum_group
    phylum_colors <- setNames(palette_table$color, palette_table$Phylum_group)
    afe[, Active_check := is.finite(AFE_CI_low_90) & AFE_CI_low_90 > 0]
    if (any(afe$Active_ASV != afe$Active_check, na.rm = TRUE)) {
      stop("Active_ASV does not agree with finalized rule AFE_CI_low_90 > 0.")
    }
    
    afe[, Active_check := NULL]
    afe[, F_new_DNA := AFE_ASV / MAX_LABELABLE_DNA_AFE]
    afe[, F_new_DNA_CI_low_90 := AFE_CI_low_90 / MAX_LABELABLE_DNA_AFE]
    afe[, F_new_DNA_CI_high_90 := AFE_CI_high_90 / MAX_LABELABLE_DNA_AFE]
    afe[, RGR_ASV := AFE_ASV / RGR_DENOMINATOR]
    afe[, RGR_CI_low_90 := AFE_CI_low_90 / RGR_DENOMINATOR]
    afe[, RGR_CI_high_90 := AFE_CI_high_90 / RGR_DENOMINATOR]
    afe[, g_exponential_day := fifelse(
      is.finite(F_new_DNA) & F_new_DNA >= 0 & F_new_DNA < 1,
      -log1p(-F_new_DNA) / INCUBATION_DAYS,
      NA_real_
    )]
    
    # Propagate AFE and abundance uncertainty into growth estimates.
    growth_result_list <- list()
    growth_qc_list <- list()
    growth_bootstrap_manifest_list <- list()
    for (g in seq_len(nrow(comparison_table))) {
      current_fire <- comparison_table$fire_trt[g]
      current_tem <- comparison_table$tem[g]
      current_id <- comparison_table$comparison_id[g]
      message("Growth calculation: ", current_id)
      afe_group <- afe[fire_trt == current_fire & tem == current_tem]
      if (nrow(afe_group) == 0) {
        stop("No finalized AFE results in ", current_id, ".")
      }
      
      bootstrap_file <- file.path(AFE_BOOTSTRAP_DIR, paste0("AFE_bootstrap_", current_id, ".rds"))
      boot <- readRDS(bootstrap_file)
      if (!identical(as.character(boot$comparison_id), current_id)) {
        stop("Bootstrap comparison_id mismatch for ", current_id, ".")
      }
      if (!identical(as.character(boot$bootstrap_design), "independent_unlabel_and_label")) {
        stop("Unexpected bootstrap design in ", current_id, ".")
      }
      
      current_asv <- as.character(boot$asv_id)
      AFE_boot <- boot$AFE_boot_corrected
      unlabel_weights <- boot$unlabel_weights
      if (nrow(AFE_boot) != length(current_asv)) {
        stop("AFE bootstrap row count mismatch in ", current_id, ".")
      }
      if (ncol(AFE_boot) != ncol(unlabel_weights)) {
        stop("Bootstrap iteration mismatch in ", current_id, ".")
      }
      
      unlabel_abund_wide <- dcast(
        wads[
          fire_trt == current_fire & tem == current_tem & iso_trt == "unlabel" & asv_id %chin% current_asv,
          .(asv_id, rep, asv_abundance)
        ],
        asv_id ~ rep,
        value.var = "asv_abundance"
      )
      rep_cols <- setdiff(names(unlabel_abund_wide), "asv_id")
      rep_cols <- as.character(sort(as.integer(rep_cols)))
      if (!identical(rep_cols, REP_LEVELS)) {
        stop("Expected exact unlabel replicate columns 1,2,3 in ", current_id, ".")
      }
      
      setkey(unlabel_abund_wide, asv_id)
      unlabel_abund_wide <- unlabel_abund_wide[J(current_asv)]
      if (anyNA(unlabel_abund_wide$asv_id)) {
        stop("Some qSIP-testable ASVs lack exact unlabel endpoint abundance in ", current_id, ".")
      }
      
      Ntotal_matrix <- as.matrix(unlabel_abund_wide[, ..REP_LEVELS])
      storage.mode(Ntotal_matrix) <- "double"
      if (anyNA(Ntotal_matrix) || any(!is.finite(Ntotal_matrix)) || any(Ntotal_matrix < 0)) {
        stop("Invalid unlabel endpoint abundance matrix in ", current_id, ".")
      }
      
      Ntotal_boot <- Ntotal_matrix %*% unlabel_weights
      RGR_boot <- AFE_boot / RGR_DENOMINATOR
      Fnew_boot <- AFE_boot / MAX_LABELABLE_DNA_AFE
      G_boot <- Ntotal_boot * RGR_boot
      RGR_boot_q <- matrixStats::rowQuantiles(RGR_boot, probs = CI_PROBS, na.rm = TRUE, drop = FALSE)
      G_boot_q <- matrixStats::rowQuantiles(G_boot, probs = CI_PROBS, na.rm = TRUE, drop = FALSE)
      Ntotal_boot_q <- matrixStats::rowQuantiles(Ntotal_boot, probs = CI_PROBS, na.rm = TRUE, drop = FALSE)
      setkey(afe_group, asv_id)
      afe_group_ordered <- afe_group[J(current_asv)]
      if (anyNA(afe_group_ordered$asv_id)) {
        stop("Failed to align finalized AFE summary to bootstrap ASVs in ", current_id, ".")
      }
      
      Ntotal_mean <- rowMeans(Ntotal_matrix)
      Ntotal_sd <- matrixStats::rowSds(Ntotal_matrix)
      current_result <- copy(afe_group_ordered)
      current_result[, N_TOTAL_unlabel_mean := Ntotal_mean]
      current_result[, N_TOTAL_unlabel_sd := Ntotal_sd]
      current_result[, N_TOTAL_boot_CI_low_90 := Ntotal_boot_q[, 1]]
      current_result[, N_TOTAL_boot_median := Ntotal_boot_q[, 2]]
      current_result[, N_TOTAL_boot_CI_high_90 := Ntotal_boot_q[, 3]]
      current_result[, RGR_ASV := AFE_ASV / RGR_DENOMINATOR]
      current_result[, RGR_CI_low_90 := AFE_CI_low_90 / RGR_DENOMINATOR]
      current_result[, RGR_CI_high_90 := AFE_CI_high_90 / RGR_DENOMINATOR]
      current_result[, RGR_boot_CI_low_90 := RGR_boot_q[, 1]]
      current_result[, RGR_boot_median := RGR_boot_q[, 2]]
      current_result[, RGR_boot_CI_high_90 := RGR_boot_q[, 3]]
      current_result[, G_CI_low_90 := G_boot_q[, 1]]
      current_result[, G_ASV := G_boot_q[, 2]]
      current_result[, G_CI_high_90 := G_boot_q[, 3]]
      current_result[, G_from_meanN_times_RGR := N_TOTAL_unlabel_mean * RGR_ASV]
      current_result[, G_boot_vs_simple_relative_difference :=
                       (G_ASV - G_from_meanN_times_RGR) / G_from_meanN_times_RGR]
      current_result[, F_new_DNA := AFE_ASV / MAX_LABELABLE_DNA_AFE]
      current_result[, F_new_DNA_CI_low_90 := AFE_CI_low_90 / MAX_LABELABLE_DNA_AFE]
      current_result[, F_new_DNA_CI_high_90 := AFE_CI_high_90 / MAX_LABELABLE_DNA_AFE]
      current_result[, Growth_model_valid_point := is.finite(F_new_DNA) & F_new_DNA >= 0 & F_new_DNA < 1]
      current_result[, F_new_exceeds_1 := is.finite(F_new_DNA) & F_new_DNA >= 1]
      current_result[, prop_boot_Fnew_lt_0 := rowMeans(Fnew_boot < 0, na.rm = TRUE)]
      current_result[, prop_boot_Fnew_ge_1 := rowMeans(Fnew_boot >= 1, na.rm = TRUE)]
      current_result[, Active_RGR := Active_ASV]
      current_result[, Active_G := Active_ASV]
      growth_result_list[[length(growth_result_list) + 1L]] <- current_result
      growth_qc_list[[length(growth_qc_list) + 1L]] <- data.table(
        fire_trt = current_fire,
        tem = current_tem,
        comparison_id = current_id,
        qSIP_testable_ASVs = nrow(current_result),
        Active_ASVs = sum(current_result$Active_ASV),
        point_Fnew_ge_1 = sum(current_result$F_new_exceeds_1, na.rm = TRUE),
        active_point_Fnew_ge_1 = sum(current_result$Active_ASV & current_result$F_new_exceeds_1, na.rm = TRUE),
        active_G_CI_low_le_0 = sum(current_result$Active_ASV & current_result$G_CI_low_90 <= 0, na.rm = TRUE),
        median_prop_boot_Fnew_lt_0_active = median(current_result$prop_boot_Fnew_lt_0[current_result$Active_ASV], na.rm = TRUE),
        median_prop_boot_Fnew_ge_1_active = median(current_result$prop_boot_Fnew_ge_1[current_result$Active_ASV], na.rm = TRUE)
      )
      growth_boot_object <- list(
        comparison_id = current_id,
        fire_trt = current_fire,
        tem = current_tem,
        bootstrap_design = "inherited_from_step02_independent_unlabel_and_label",
        asv_id = current_asv,
        rep_levels = REP_LEVELS,
        n_boot = ncol(AFE_boot),
        RGR_boot = RGR_boot,
        G_boot = G_boot,
        Ntotal_unlabel_boot = Ntotal_boot,
        Fnew_boot = Fnew_boot
      )
      growth_boot_file <- file.path(growth_bootstrap_dir, paste0("G_bootstrap_", current_id, ".rds"))
      saveRDS(growth_boot_object, growth_boot_file, compress = TRUE)
      growth_bootstrap_manifest_list[[length(growth_bootstrap_manifest_list) + 1L]] <- data.table(
        fire_trt = current_fire,
        tem = current_tem,
        comparison_id = current_id,
        qSIP_testable_ASVs = length(current_asv),
        active_ASVs = sum(current_result$Active_ASV),
        n_boot = ncol(AFE_boot),
        file = growth_boot_file
      )
      rm(
        boot, AFE_boot, unlabel_weights, Ntotal_matrix, Ntotal_boot,
        RGR_boot, Fnew_boot, G_boot, RGR_boot_q, G_boot_q, Ntotal_boot_q
      )
      gc()
    }
    
    growth <- rbindlist(growth_result_list, use.names = TRUE, fill = TRUE)
    growth_qc <- rbindlist(growth_qc_list, use.names = TRUE, fill = TRUE)
    growth_bootstrap_manifest <- rbindlist(growth_bootstrap_manifest_list, use.names = TRUE, fill = TRUE)
    growth_group_summary <- growth[
      ,
      .(
        qSIP_testable_ASVs = .N,
        Active_ASVs = sum(Active_ASV),
        Active_percent = 100 * mean(Active_ASV),
        median_RGR_active = if (any(Active_ASV)) median(RGR_ASV[Active_ASV], na.rm = TRUE) else NA_real_,
        median_G_active = if (any(Active_ASV)) median(G_ASV[Active_ASV], na.rm = TRUE) else NA_real_,
        sum_G_active_ASV_medians = if (any(Active_ASV)) sum(G_ASV[Active_ASV], na.rm = TRUE) else NA_real_
      ),
      by = .(fire_trt, tem, comparison_id)
    ]
    fwrite(growth, file.path(output_dir, "01_RGR_G_ASV_all.tsv"), sep = "\t")
    fwrite(growth[Active_ASV == TRUE], file.path(output_dir, "02_RGR_G_ASV_active.tsv"), sep = "\t")
    fwrite(growth_group_summary, file.path(output_dir, "03_RGR_G_group_summary.tsv"), sep = "\t")
    fwrite(growth[F_new_exceeds_1 == TRUE], file.path(output_dir, "04_QC_Fnew_ge_1.tsv"), sep = "\t")
    fwrite(growth_qc, file.path(output_dir, "05_growth_bootstrap_QC.tsv"), sep = "\t")
    fwrite(growth_bootstrap_manifest, file.path(output_dir, "06_growth_bootstrap_object_manifest.tsv"), sep = "\t")
    
    rgr_consistency <- growth[
      ,
      .(
        max_abs_RGR_difference = max(abs(RGR_ASV - RGR_boot_median), na.rm = TRUE),
        median_abs_RGR_difference = median(abs(RGR_ASV - RGR_boot_median), na.rm = TRUE)
      ),
      by = .(fire_trt, tem, comparison_id)
    ]
    exact_manifest <- copy(comparison_table)
    exact_manifest[, unlabel_source_for_Ntotal := paste(fire_trt, tem, "unlabel", sep = "_")]
    exact_manifest[, AFE_source := file.path(CACHE$afe, "bootstrap_objects", paste0("AFE_bootstrap_", comparison_id, ".rds"))]
    fwrite(rgr_consistency, file.path(output_dir, "11_RGR_bootstrap_consistency_QC.tsv"), sep = "\t")
    fwrite(exact_manifest, file.path(output_dir, "12_exact_growth_comparison_manifest.tsv"), sep = "\t")
    parameter_table <- data.table(
      Parameter = c(
        "Final soil-water 18O atom fraction",
        "Natural abundance 18O",
        "Soil-water 18O AFE",
        "Fraction DNA oxygen from water (propO)",
        "Incubation days",
        "Maximum labelable DNA AFE",
        "RGR formula",
        "G formula",
        "N_total source",
        "AFE bootstrap source",
        "Bootstrap design",
        "Primary active ASV criterion",
        "Absolute growth unit",
        "Multiple-testing correction"
      ),
      Value = c(
        as.character(SOIL_WATER_18O_ATOM_FRACTION),
        as.character(NATURAL_18O),
        as.character(SOIL_WATER_AFE),
        as.character(PROP_O_FROM_WATER),
        as.character(INCUBATION_DAYS),
        as.character(MAX_LABELABLE_DNA_AFE),
        "AFE_ASV / (soil_water_AFE * propO * incubation_days)",
        "median of N_total_unlabel_boot * RGR_boot",
        "exact corresponding unlabel fire x temperature endpoint ASV abundance",
        "exact corrected AFE bootstrap draws from step 02",
        "independent_unlabel_and_label inherited from step 02",
        "AFE_CI_low_90 > 0",
        "16S copies g^-1 dry soil day^-1",
        "None"
      )
    )
    fwrite(parameter_table, file.path(output_dir, "13_growth_analysis_parameters.tsv"), sep = "\t")
    invisible(NULL)
  }
  
  # Community
  run_community <- function() {
    output_dir <- CACHE$community
    comparison_table <- CJ(fire_trt = FIRE_LEVELS, tem = TEMP_LEVELS, unique = TRUE)
    comparison_table[, comparison_id := paste(fire_trt, tem, sep = "_")]
    
    safe_quantiles <- function(x, probs = CI_PROBS) {
      x <- x[is.finite(x)]
      if (length(x) == 0) return(rep(NA_real_, length(probs)))
      as.numeric(quantile(x, probs = probs, na.rm = TRUE, names = FALSE, type = 7))
    }
    
    safe_log10_positive <- function(x) {
      x <- as.numeric(x)
      out <- rep(NA_real_, length(x))
      ok <- is.finite(x) & x > 0
      out[ok] <- log10(x[ok])
      out
    }
    
    bootstrap_tail_probability <- function(x, null = 0) {
      x <- x[is.finite(x)]
      if (length(x) == 0) return(NA_real_)
      p_lower <- mean(x <= null)
      p_upper <- mean(x >= null)
      min(1, 2 * min(p_lower, p_upper))
    }
    
    normalize_asv_taxon_rank <- function(x) {
      x <- trimws(as.character(x))
      x <- sub("^[A-Za-z]__", "", x)
      low <- tolower(x)
      cond1 <- is.na(x) | x == ""
      cond2 <- low %in% c("na", "n/a", "unclassified", "unassigned", "unknown", "uncultured", "uncultured bacterium", "uncultured_bacterium")
      cond3 <- grepl("^unclassified($|[ _-])", low)
      cond4 <- grepl("^unassigned($|[ _-])", low)
      unresolved <- cond1 | cond2 | cond3 | cond4
      x[unresolved] <- NA_character_
      return(x)
    }
    
    first_resolved_asv_taxon <- function(x) {
      x <- normalize_asv_taxon_rank(x)
      x <- unique(x[!is.na(x)])
      if (length(x) == 0L) return("NA")
      sort(x)[1]
    }
    
    make_asv_taxonomy_display_label <- function(phylum, class, order, family, genus, asv_id) {
      ranks <- c(Phylum = as.character(phylum), Class = as.character(class), Order = as.character(order), Family = as.character(family), Genus = as.character(genus))
      ranks[is.na(ranks) | ranks == ""] <- "NA"
      resolved <- which(ranks != "NA")
      if (ranks[["Genus"]] != "NA") return(paste(ranks[["Genus"]], as.character(asv_id)))
      if (length(resolved) == 0L) return(paste(c(rep("NA", 5L), as.character(asv_id)), collapse = " "))
      deepest <- max(resolved)
      lower_na_n <- 5L - deepest
      paste(c(ranks[deepest], if (lower_na_n > 0L) rep("NA", lower_na_n) else character(0), as.character(asv_id)), collapse = " ")
    }
    
    make_compact_letters <- function(dt_one_taxon) {
      stage_split <- split(dt_one_taxon$DeltaG_group, dt_one_taxon$fire_trt)
      stage_vectors <- setNames(lapply(FIRE_LEVELS, function(st) {
        vals <- stage_split[[st]]
        if (is.null(vals) || length(vals) == 0) rep(0, max(dt_one_taxon$bootstrap_iteration, na.rm = TRUE)) else vals[is.finite(vals)]
      }), FIRE_LEVELS)
      pair_list <- combn(FIRE_LEVELS, 2, simplify = FALSE)
      pvals <- vapply(pair_list, function(pair) {
        x <- stage_vectors[[pair[1]]]
        y <- stage_vectors[[pair[2]]]
        n_use <- min(length(x), length(y))
        if (n_use == 0) return(1)
        bootstrap_tail_probability(x[seq_len(n_use)] - y[seq_len(n_use)], null = 0)
      }, FUN.VALUE = numeric(1))
      names(pvals) <- vapply(pair_list, function(pair) paste(pair, collapse = "-"), character(1))
      letters <- multcompView::multcompLetters(pvals, threshold = 0.05)$Letters
      data.table(fire_trt = names(letters), sig_letter = unname(letters))
    }
    
    growth <- fread(file.path(CACHE$growth, "01_RGR_G_ASV_all.tsv"), sep = "\t", header = TRUE)
    growth[, `:=`(asv_id = as.character(asv_id), fire_trt = as.character(fire_trt), tem = as.character(tem), Active_ASV = as.logical(Active_ASV), G_ASV = as.numeric(G_ASV))]
    asv_taxonomy_lookup <- growth[, .(
      Phylum = first_resolved_asv_taxon(Phylum),
      Class = first_resolved_asv_taxon(Class),
      Order = first_resolved_asv_taxon(Order),
      Family = first_resolved_asv_taxon(Family),
      Genus = first_resolved_asv_taxon(Genus)
    ), by = asv_id]
    asv_taxonomy_lookup[, ASV_display_label := mapply(make_asv_taxonomy_display_label, Phylum, Class, Order, Family, Genus, asv_id, USE.NAMES = FALSE)]
    
    # Sum detected-active gene-copy production using shared growth draws.
    community_boot_list <- list()
    community_summary_list <- list()
    tax_group_boot_list <- list()
    asv_active_G_boot_list <- list()
    for (g in seq_len(nrow(comparison_table))) {
      current_fire <- comparison_table$fire_trt[g]
      current_tem <- comparison_table$tem[g]
      current_id <- comparison_table$comparison_id[g]
      boot <- readRDS(file.path(CACHE$growth, "bootstrap_objects", paste0("G_bootstrap_", current_id, ".rds")))
      current_asv <- as.character(boot$asv_id)
      active_mask <- growth[fire_trt == current_fire & tem == current_tem][match(current_asv, asv_id)]$Active_ASV
      if (anyNA(active_mask)) stop("Growth/active-ASV alignment failed: ", current_id)
      G_boot_linear <- boot$G_boot
      active_G_boot <- G_boot_linear[active_mask, , drop = FALSE]
      rownames(active_G_boot) <- current_asv[active_mask]
      asv_active_G_boot_list[[current_id]] <- list(asv_id = current_asv[active_mask], G_boot = active_G_boot)
      community_G_boot <- colSums(active_G_boot)
      community_q <- safe_quantiles(community_G_boot, CI_PROBS)
      community_boot_list[[length(community_boot_list) + 1L]] <- data.table(fire_trt = current_fire, tem = current_tem, comparison_id = current_id, bootstrap_iteration = seq_along(community_G_boot), Total_active_G = community_G_boot)
      community_summary_list[[length(community_summary_list) + 1L]] <- data.table(fire_trt = current_fire, tem = current_tem, Total_active_G_CI_low_90 = community_q[1], Total_active_G = community_q[2], Total_active_G_CI_high_90 = community_q[3])
      tax_lookup <- unique(growth[fire_trt == current_fire & tem == current_tem & asv_id %chin% current_asv, .(asv_id, Phylum)])
      setkey(tax_lookup, asv_id)
      tax_lookup <- tax_lookup[J(current_asv)]
      phylum_labels <- trimws(as.character(tax_lookup$Phylum[active_mask]))
      phylum_labels[is.na(phylum_labels) | phylum_labels == "" | phylum_labels == "NA"] <- "Unclassified"
      group_mat <- rowsum(active_G_boot, group = phylum_labels, reorder = FALSE)
      group_dt <- as.data.table(group_mat, keep.rownames = "taxon")
      long_dt <- melt(group_dt, id.vars = "taxon", variable.name = "bootstrap_col", value.name = "Group_G")
      long_dt[, `:=`(bootstrap_iteration = as.integer(gsub("^V", "", bootstrap_col)), taxon_level = "Phylum", fire_trt = current_fire, tem = current_tem)]
      tax_group_boot_list[[length(tax_group_boot_list) + 1L]] <- long_dt[, .(taxon_level, taxon, fire_trt, tem, bootstrap_iteration, Group_G)]
    }
    
    community_boot <- rbindlist(community_boot_list, use.names = TRUE, fill = TRUE)
    community_summary <- rbindlist(community_summary_list, use.names = TRUE, fill = TRUE)
    tax_group_boot <- rbindlist(tax_group_boot_list, use.names = TRUE, fill = TRUE)
    delta_boot_list <- list()
    delta_summary_list <- list()
    for (current_fire in FIRE_LEVELS) {
      g21 <- community_boot[fire_trt == current_fire & tem == "T_low"][order(bootstrap_iteration)]
      g28 <- community_boot[fire_trt == current_fire & tem == "T_high"][order(bootstrap_iteration)]
      delta_G <- g28$Total_active_G - g21$Total_active_G
      delta_q <- safe_quantiles(delta_G, CI_PROBS)
      delta_summary_list[[length(delta_summary_list) + 1L]] <- data.table(fire_trt = current_fire, DeltaG_CI_low_90 = delta_q[1], DeltaG_community = delta_q[2], DeltaG_CI_high_90 = delta_q[3])
      delta_boot_list[[length(delta_boot_list) + 1L]] <- data.table(fire_trt = current_fire, bootstrap_iteration = seq_along(delta_G), DeltaG_community = delta_G)
    }
    
    delta_boot <- rbindlist(delta_boot_list, use.names = TRUE, fill = TRUE)
    delta_summary <- rbindlist(delta_summary_list, use.names = TRUE, fill = TRUE)
    stage_split_delta <- split(delta_boot$DeltaG_community, delta_boot$fire_trt)
    stage_vectors_delta <- setNames(lapply(FIRE_LEVELS, function(st) {
      vals <- stage_split_delta[[st]]
      if (is.null(vals) || length(vals) == 0) rep(0, max(delta_boot$bootstrap_iteration, na.rm = TRUE)) else vals[is.finite(vals)]
    }), FIRE_LEVELS)
    pair_list_delta <- combn(FIRE_LEVELS, 2, simplify = FALSE)
    pvals_delta <- vapply(pair_list_delta, function(pair) {
      x <- stage_vectors_delta[[pair[1]]]
      y <- stage_vectors_delta[[pair[2]]]
      n_use <- min(length(x), length(y))
      if (n_use == 0) return(1)
      bootstrap_tail_probability(x[seq_len(n_use)] - y[seq_len(n_use)], null = 0)
    }, FUN.VALUE = numeric(1))
    names(pvals_delta) <- vapply(pair_list_delta, function(pair) paste(pair, collapse = "-"), character(1))
    letters_res <- multcompView::multcompLetters(pvals_delta, threshold = 0.05)$Letters
    community_delta_letters <- data.table(fire_trt = names(letters_res), sig_letter = unname(letters_res))
    plot_A <- copy(community_summary)
    plot_A[, `:=`(log10_Total_active_G = safe_log10_positive(Total_active_G), log10_Total_active_G_CI_low_90 = safe_log10_positive(Total_active_G_CI_low_90), log10_Total_active_G_CI_high_90 = safe_log10_positive(Total_active_G_CI_high_90), fire_trt = factor(fire_trt, levels = FIRE_LEVELS, labels = unname(FIRE_LABELS[FIRE_LEVELS])), tem = factor(tem, levels = TEMP_LEVELS, labels = unname(TEMP_LABELS[TEMP_LEVELS])))]
    plot_B <- copy(delta_summary)
    plot_B <- merge(plot_B, community_delta_letters, by = "fire_trt", all.x = TRUE)
    plot_B[is.na(sig_letter), sig_letter := "a"]
    plot_B[, `:=`(signed_log10_DeltaG_community = signed_log10_abs(DeltaG_community), signed_log10_DeltaG_CI_low_90 = signed_log10_abs(DeltaG_CI_low_90), signed_log10_DeltaG_CI_high_90 = signed_log10_abs(DeltaG_CI_high_90))]
    plot_B[, fire_trt := factor(fire_trt, levels = FIRE_LEVELS, labels = unname(FIRE_LABELS[FIRE_LEVELS]))]
    plot_B[, letter_y_pos := signed_log10_DeltaG_CI_high_90 + 0.18]
    all_taxa <- sort(unique(tax_group_boot$taxon))
    n_boot_global <- tax_group_boot[, uniqueN(bootstrap_iteration)]
    delta_complete_list <- list()
    delta_summary_list <- list()
    for (current_fire in FIRE_LEVELS) {
      dt21 <- tax_group_boot[fire_trt == current_fire & tem == "T_low", .(taxon, bootstrap_iteration, G21 = Group_G)]
      dt28 <- tax_group_boot[fire_trt == current_fire & tem == "T_high", .(taxon, bootstrap_iteration, G28 = Group_G)]
      merged_delta <- merge(dt21, dt28, by = c("taxon", "bootstrap_iteration"), all = TRUE)
      merged_delta[is.na(G21), G21 := 0]; merged_delta[is.na(G28), G28 := 0]
      merged_delta[, DeltaG_group := G28 - G21]
      complete_idx <- CJ(taxon = all_taxa, bootstrap_iteration = seq_len(n_boot_global), unique = TRUE)
      merged_delta <- merge(complete_idx, merged_delta[, .(taxon, bootstrap_iteration, DeltaG_group)], by = c("taxon", "bootstrap_iteration"), all.x = TRUE)
      merged_delta[is.na(DeltaG_group), DeltaG_group := 0]
      merged_delta[, c("taxon_level", "fire_trt") := .("Phylum", current_fire)]
      delta_complete_list[[length(delta_complete_list) + 1L]] <- merged_delta
      delta_summary_list[[length(delta_summary_list) + 1L]] <- merged_delta[, {
        qv <- safe_quantiles(DeltaG_group, CI_PROBS)
        .(DeltaG_CI_low_90 = qv[1], DeltaG_group = qv[2], DeltaG_CI_high_90 = qv[3])
      }, by = .(taxon_level, fire_trt, taxon)]
    }
    
    tax_delta_boot <- rbindlist(delta_complete_list, use.names = TRUE, fill = TRUE)
    tax_delta_summary <- rbindlist(delta_summary_list, use.names = TRUE, fill = TRUE)
    top10_phylum <- tax_delta_summary[, .(rank_score = sum(abs(DeltaG_group), na.rm = TRUE)), by = taxon][order(-rank_score, taxon)][1:min(.N, 10L)]
    phylum_letters_list <- list()
    for (taxon_i in top10_phylum$taxon) {
      dt_i <- tax_delta_boot[taxon == taxon_i, .(fire_trt, bootstrap_iteration, DeltaG_group)]
      letters_i <- make_compact_letters(dt_i)
      letters_i[, `:=`(taxon = taxon_i, taxon_level = "Phylum")]
      phylum_letters_list[[length(phylum_letters_list) + 1L]] <- letters_i
    }
    
    phylum_letters <- rbindlist(phylum_letters_list, use.names = TRUE, fill = TRUE)
    top10_phylum[, taxon_plot_label := stringr::str_wrap(paste0(sprintf("%02d", seq_len(.N)), ". ", taxon), width = 22)]
    plot_C_data <- merge(tax_delta_summary[taxon %chin% top10_phylum$taxon], phylum_letters[, .(taxon, fire_trt, sig_letter)], by = c("taxon", "fire_trt"), all.x = TRUE)
    plot_C_data <- merge(plot_C_data, top10_phylum[, .(taxon, taxon_plot_label)], by = "taxon", all.x = TRUE)
    plot_C_data[, fire_trt_label := factor(fire_trt, levels = FIRE_LEVELS, labels = unname(FIRE_LABELS[FIRE_LEVELS]))]
    plot_C_data[, taxon_plot_label := factor(taxon_plot_label, levels = top10_phylum$taxon_plot_label)]
    panel_stats <- plot_C_data[, .(panel_ymin = min(c(0, DeltaG_CI_low_90), na.rm = TRUE), panel_ymax = max(c(0, DeltaG_CI_high_90), na.rm = TRUE)), by = .(taxon, taxon_plot_label)]
    panel_stats[, panel_span := pmax(panel_ymax - panel_ymin, abs(panel_ymax), abs(panel_ymin), 1)]
    plot_C_data <- merge(plot_C_data, panel_stats, by = c("taxon", "taxon_plot_label"), all.x = TRUE)
    plot_C_data[, label_y := ifelse(DeltaG_group >= 0, DeltaG_CI_high_90 + 0.10 * panel_span, 0 + 0.08 * panel_span)]
    plot_C_data[!is.finite(label_y), label_y := 0]
    plot_C_data[is.na(sig_letter), sig_letter := "a"]
    plot_C_data[, `:=`(signed_log10_DeltaG_group = signed_log10_abs(DeltaG_group), signed_log10_DeltaG_CI_low_90 = signed_log10_abs(DeltaG_CI_low_90), signed_log10_DeltaG_CI_high_90 = signed_log10_abs(DeltaG_CI_high_90), signed_log10_label_y = signed_log10_abs(label_y))]
    top_n <- 15L
    summary_list_15 <- list()
    for (current_fire in FIRE_LEVELS) {
      obj21 <- asv_active_G_boot_list[[paste(current_fire, "T_low", sep = "_")]]
      obj28 <- asv_active_G_boot_list[[paste(current_fire, "T_high", sep = "_")]]
      union_ids <- sort(unique(c(obj21$asv_id, obj28$asv_id)))
      if (length(union_ids) == 0L) next
      n_boot <- ncol(obj21$G_boot)
      G21_union <- matrix(0, nrow = length(union_ids), ncol = n_boot, dimnames = list(union_ids, NULL))
      G28_union <- matrix(0, nrow = length(union_ids), ncol = n_boot, dimnames = list(union_ids, NULL))
      G21_union[match(obj21$asv_id, union_ids), ] <- obj21$G_boot
      G28_union[match(obj28$asv_id, union_ids), ] <- obj28$G_boot
      delta_mat <- G28_union - G21_union
      delta_q <- t(apply(delta_mat, 1, function(x) safe_quantiles(x, CI_PROBS)))
      selected <- data.table(
        taxon = union_ids,
        DeltaG_low90 = as.numeric(delta_q[, 1]),
        DeltaG_group = as.numeric(delta_q[, 2])
      )[is.finite(DeltaG_low90) & DeltaG_low90 > 0][order(-DeltaG_group, taxon)]
      if (nrow(selected) == 0L) next
      selected[, rank := seq_len(.N)]
      selected <- selected[rank <= min(top_n, .N)]
      positive_all <- pmax(delta_mat, 0)
      denominator_boot <- colSums(positive_all)
      valid_denominator <- is.finite(denominator_boot) & denominator_boot > 0
      selected_idx <- match(selected$taxon, union_ids)
      selected_positive <- positive_all[selected_idx, , drop = FALSE]
      selected_cumulative <- apply(selected_positive, 2, cumsum)
      if (is.null(dim(selected_cumulative))) selected_cumulative <- matrix(selected_cumulative, nrow = nrow(selected_positive), ncol = n_boot)
      individual_percent <- matrix(NA_real_, nrow = nrow(selected_positive), ncol = n_boot)
      cumulative_percent <- matrix(NA_real_, nrow = nrow(selected_positive), ncol = n_boot)
      individual_percent[, valid_denominator] <- sweep(selected_positive[, valid_denominator, drop = FALSE], 2, denominator_boot[valid_denominator], "/") * 100
      cumulative_percent[, valid_denominator] <- sweep(selected_cumulative[, valid_denominator, drop = FALSE], 2, denominator_boot[valid_denominator], "/") * 100
      summary_current <- rbindlist(lapply(seq_len(nrow(selected)), function(rr) {
        q_ind <- safe_quantiles(individual_percent[rr, ], CI_PROBS)
        q_cum <- safe_quantiles(cumulative_percent[rr, ], CI_PROBS)
        data.table(fire_trt = current_fire, taxon = selected$taxon[rr], rank = selected$rank[rr], individual_percent_median = q_ind[2], cumulative_percent_CI_low_90 = q_cum[1], cumulative_percent_median = q_cum[2], cumulative_percent_CI_high_90 = q_cum[3])
      }), use.names = TRUE, fill = TRUE)
      summary_current <- merge(summary_current, asv_taxonomy_lookup[, .(asv_id, ASV_display_label)], by.x = "taxon", by.y = "asv_id", all.x = TRUE, sort = FALSE)
      summary_list_15[[length(summary_list_15) + 1L]] <- summary_current
    }
    
    plot_ASV_15 <- if (length(summary_list_15)) rbindlist(summary_list_15, use.names = TRUE, fill = TRUE) else data.table(
      fire_trt = character(), taxon = character(), rank = integer(), ASV_display_label = character(),
      individual_percent_median = numeric(), cumulative_percent_CI_low_90 = numeric(),
      cumulative_percent_median = numeric(), cumulative_percent_CI_high_90 = numeric())
    plot_ASV_15[, fire_trt_label := factor(fire_trt, levels = FIRE_LEVELS, labels = unname(FIRE_LABELS[FIRE_LEVELS]))]
    plot_ASV_15[, taxon_rank_label := stringr::str_wrap(paste0(sprintf("%02d", rank), ". ", ASV_display_label), width = 36)]
    p_A <- ggplot(plot_A, aes(x = fire_trt, y = log10_Total_active_G, fill = tem)) +
      geom_col(position = position_dodge(width = 0.72), width = 0.62, colour = "black", linewidth = 0.45) +
      geom_errorbar(aes(ymin = log10_Total_active_G_CI_low_90, ymax = log10_Total_active_G_CI_high_90), position = position_dodge(width = 0.72), width = 0.16, linewidth = 0.55) +
      scale_fill_manual(values = TEMP_COLORS_DISPLAY, name = "Incubation temperature") +
      scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
      labs(x = NULL, y = expression(log[10] * " total detected active growth (16S copies " * g^{-1} * " dry soil " * d^{-1} * ")")) +
      theme_classic(base_size = 12) + theme(axis.title.y = element_text(size = 12, face = "bold"), axis.text.x = element_text(size = 11, colour = "black"), axis.text.y = element_text(size = 10, colour = "black"), axis.line = element_line(linewidth = 0.6, colour = "black"), axis.ticks = element_line(colour = "black"), legend.position = "top", legend.title = element_text(size = 10, face = "bold"), legend.text = element_text(size = 10), plot.margin = margin(t = 8, r = 8, b = 8, l = 8))
    p_B <- ggplot(plot_B, aes(x = fire_trt, y = signed_log10_DeltaG_community, fill = fire_trt)) +
      geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.6, colour = "grey45") +
      geom_col(width = 0.62, colour = "black", linewidth = 0.45) +
      geom_errorbar(aes(ymin = signed_log10_DeltaG_CI_low_90, ymax = signed_log10_DeltaG_CI_high_90), width = 0.16, linewidth = 0.55) +
      geom_text(aes(y = letter_y_pos, label = sig_letter), size = 6.0, fontface = "bold", colour = "black", vjust = 0) +
      scale_fill_manual(values = FIRE_COLORS_DISPLAY, guide = "none") +
      scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
      labs(x = NULL, y = "Signed log10(1 + |community DeltaG|)") +
      theme_classic(base_size = 12) + theme(axis.title.y = element_text(size = 12, face = "bold"), axis.text.x = element_text(size = 11, colour = "black"), axis.text.y = element_text(size = 10, colour = "black"), axis.line = element_line(linewidth = 0.6, colour = "black"), axis.ticks = element_line(colour = "black"), plot.margin = margin(t = 8, r = 8, b = 8, l = 8))
    combined_plot <- p_A + p_B + plot_layout(ncol = 2, widths = c(1.15, 1)) + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 16, face = "bold"))
    p_C <- ggplot(plot_C_data, aes(x = fire_trt_label, y = signed_log10_DeltaG_group, fill = fire_trt_label)) +
      geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5, colour = "grey45") +
      geom_col(width = 0.68, colour = "black", linewidth = 0.4) +
      geom_errorbar(aes(ymin = signed_log10_DeltaG_CI_low_90, ymax = signed_log10_DeltaG_CI_high_90), width = 0.15, linewidth = 0.48) +
      geom_text(aes(y = signed_log10_label_y, label = sig_letter), size = 3.7, fontface = "bold", vjust = 0, colour = "black") +
      facet_wrap(~ taxon_plot_label, scales = "free_y", ncol = 4) +
      scale_fill_manual(values = FIRE_COLORS_DISPLAY, name = "Fire stage") +
      scale_y_continuous(expand = expansion(mult = c(0.12, 0.22))) +
      coord_cartesian(clip = "off") +
      labs(x = NULL, y = "Signed log10(1 + |phylum DeltaG|)", title = "Top 10 phyla contributing to warming-induced community growth change") +
      theme_classic(base_size = 11) + theme(strip.background = element_rect(fill = "grey95", colour = "black", linewidth = 0.5), strip.text = element_text(size = 9, face = "bold", colour = "black", lineheight = 0.95), axis.text.x = element_text(size = 8.5, angle = 25, hjust = 1, colour = "black"), axis.text.y = element_text(size = 8.5, colour = "black"), axis.title.y = element_text(size = 11.5, face = "bold"), axis.line = element_line(linewidth = 0.6, colour = "black"), axis.ticks = element_line(linewidth = 0.45, colour = "black"), legend.position = "top", legend.title = element_text(size = 10, face = "bold"), legend.text = element_text(size = 9.5), panel.spacing = unit(0.85, "lines"), plot.title = element_text(size = 12, face = "bold", hjust = 0.5), plot.margin = margin(t = 8, r = 16, b = 8, l = 8))
    
    make_stage_panel_15 <- function(stage_id) {
      stage_dt <- plot_ASV_15[fire_trt == stage_id][order(rank)]
      if (nrow(stage_dt) == 0L) return(NULL)
      stage_label <- unname(FIRE_LABELS[stage_id])
      stage_color <- unname(FIRE_COLORS[stage_id])
      stage_dt[, taxon_rank_label := factor(taxon_rank_label, levels = rev(taxon_rank_label))]
      heat_max <- max(plot_ASV_15$individual_percent_median, na.rm = TRUE)
      if (!is.finite(heat_max) || heat_max <= 0) heat_max <- 1
      p_heat <- ggplot(stage_dt, aes(x = 1, y = taxon_rank_label, fill = individual_percent_median)) +
        geom_tile(width = 0.92, height = 0.92, colour = "black", linewidth = 0.35) +
        scale_fill_gradient(low = "grey98", high = "red3", name = "Positive DeltaG\nshare (%)", limits = c(0, heat_max), oob = scales::squish, labels = function(x) sprintf("%.1f%%", x)) +
        scale_x_continuous(expand = c(0, 0), breaks = NULL) +
        labs(title = stage_label, subtitle = paste0("Top ", nrow(stage_dt), " positive ASVs (90% CI > 0)"), x = NULL, y = NULL) +
        theme_classic(base_size = 10.5) + theme(plot.title = element_text(size = 11.5, face = "bold", hjust = 0.5), plot.subtitle = element_text(size = 9.0, hjust = 0.5), axis.text.y = element_text(size = 7.3, colour = "black", hjust = 1, lineheight = 0.92), axis.ticks = element_blank(), axis.line = element_blank(), legend.position = "right", legend.title = element_text(size = 8.6, face = "bold"), legend.text = element_text(size = 8.1), plot.margin = margin(t = 5, r = 6, b = 2, l = 5))
      final_dt <- stage_dt[rank == max(rank)]
      p_cum <- ggplot(stage_dt, aes(x = rank, y = cumulative_percent_median)) +
        geom_ribbon(aes(ymin = cumulative_percent_CI_low_90, ymax = cumulative_percent_CI_high_90), fill = stage_color, alpha = 0.20, colour = NA) +
        geom_line(linewidth = 0.90, colour = stage_color) + geom_point(size = 1.9, colour = stage_color) +
        geom_text(data = final_dt, aes(label = sprintf("Top%d = %.1f%%\n90%% CI %.1f-%.1f%%", rank, cumulative_percent_median, cumulative_percent_CI_low_90, cumulative_percent_CI_high_90)), size = 2.8, fontface = "bold", hjust = 1, vjust = -0.25, colour = stage_color) +
        scale_x_continuous(breaks = c(1, seq(5, 15, by = 5)), limits = c(1, 15), minor_breaks = NULL, expand = expansion(mult = c(0.01, 0.03))) +
        scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 25), labels = function(x) paste0(x, "%"), expand = expansion(mult = c(0.02, 0.10))) +
        labs(x = "Rank within stage", y = "Cumulative % of total positive\nDeltaG (bootstrap median, 90% CI)") +
        theme_classic(base_size = 10) + theme(axis.title = element_text(size = 9.5, face = "bold"), axis.text = element_text(size = 8.3, colour = "black"), axis.line = element_line(linewidth = 0.55, colour = "black"), axis.ticks = element_line(linewidth = 0.4, colour = "black"), plot.margin = margin(t = 0, r = 6, b = 5, l = 5))
      (p_heat / p_cum) + patchwork::plot_layout(heights = c(4.4, 1.7))
    }
    
    panel_list <- Filter(Negate(is.null), lapply(FIRE_LEVELS, make_stage_panel_15))
    p_top15 <- if (!length(panel_list)) blank_plot("No positive ASVs pass the 90% CI criterion") else patchwork::wrap_plots(panel_list, nrow = 1, guides = "collect") +
      patchwork::plot_annotation(
        title = "Stage-specific concentration of warming-stimulated growth at the ASV level: Top 15",
        subtitle = "Upper: bootstrap median share of total positive DeltaG for strictly positive Top 15 ASVs. Lower: cumulative contribution with 90% percentile CI.",
        theme = theme(plot.title = element_text(size = 12.5, face = "bold", hjust = 0.5), plot.subtitle = element_text(size = 10, hjust = 0.5), legend.position = "right")
      )
    stacked_dt <- copy(tax_delta_summary)
    stacked_dt[, Phylum_plot := fifelse(taxon %in% top10_phylum$taxon, taxon, "Others")]
    stacked_summary <- stacked_dt[, .(DeltaG_sum = sum(DeltaG_group, na.rm = TRUE)), by = .(fire_trt, Phylum_plot)]
    stacked_summary[, fire_trt_label := factor(fire_trt, levels = FIRE_LEVELS, labels = unname(FIRE_LABELS[FIRE_LEVELS]))]
    phylum_levels_stacked <- c("Others", rev(top10_phylum$taxon))
    stacked_summary[, Phylum_plot := factor(Phylum_plot, levels = phylum_levels_stacked)]
    my_phylum_colors <- setNames(
      c("grey80", brewer.pal(min(10, length(top10_phylum$taxon)), "Paired")),
      phylum_levels_stacked
    )
    p_D <- ggplot(stacked_summary, aes(x = fire_trt_label, y = DeltaG_sum, fill = Phylum_plot)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.8) +
      geom_col(position = "stack", color = "black", linewidth = 0.3, width = 0.55) +
      scale_fill_manual(values = my_phylum_colors, name = "Phylum") +
      scale_y_continuous(labels = scales::scientific, expand = expansion(mult = c(0.1, 0.1))) +
      labs(x = NULL,
           y = expression(Delta * "G Flux (16S copies " * g^{-1} * " dry soil " * d^{-1} * ")"),
           title = "Phylum-level Composition of Net Thermal Response") +
      theme_classic(base_size = 12) +
      theme(
        axis.title.y = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 11, color = "black", face = "bold"),
        axis.text.y = element_text(size = 10, color = "black"),
        legend.position = "right",
        legend.title = element_text(face = "bold"),
        plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
      )
    save_figure(filename = file.path(output_dir, "Fig_AB_macro_thermal_response.pdf"), plot = combined_plot, width = 13.5, height = 5.8, units = "in", device = cairo_pdf)
    save_figure(filename = file.path(output_dir, "Fig_C_top10_phylum_DeltaG_polished.pdf"), plot = p_C, width = 10, height = 10, units = "in", device = cairo_pdf)
    save_figure(filename = file.path(output_dir, "Fig_ASV_Top15_heatmap_cumulative_DeltaG_percent_polished.pdf"), plot = p_top15, width = 17.0, height = 10.2, units = "in", device = cairo_pdf)
    save_figure(filename = file.path(output_dir, "Fig_D_phylum_stacked_DeltaG.pdf"), plot = p_D, width = 7.5, height = 5.5, units = "in", device = cairo_pdf)
    fwrite(plot_A, file.path(output_dir, "01_FigA_total_detected_active_growth_summary.tsv"), sep = "\t")
    fwrite(plot_B, file.path(output_dir, "02_FigB_community_DeltaG_summary.tsv"), sep = "\t")
    fwrite(top10_phylum, file.path(output_dir, "03_FigC_top10_phylum_selection.tsv"), sep = "\t")
    fwrite(plot_C_data, file.path(output_dir, "04_FigC_top10_phylum_plot_data.tsv"), sep = "\t")
    fwrite(plot_ASV_15, file.path(output_dir, "05_Fig_ASV_Top15_summary.tsv"), sep = "\t")
    fwrite(stacked_summary, file.path(output_dir, "04b_FigD_phylum_stacked_plot_data.tsv"), sep = "\t")
    fwrite(community_boot, file.path(output_dir, "Community_growth_bootstrap.tsv"), sep = "\t")
    fwrite(delta_boot, file.path(output_dir, "Community_delta_bootstrap.tsv"), sep = "\t")
    fwrite(data.table(comparison = names(pvals_delta), bootstrap_tail_probability = unname(pvals_delta)),
           file.path(output_dir, "Community_pairwise_tail_probabilities.tsv"), sep = "\t")
    fwrite(tax_delta_summary, file.path(output_dir, "Phylum_delta_summary.tsv"), sep = "\t")
    invisible(NULL)
  }
  
  # Strategy
  run_strategy <- function() {
    GROWTH_FILE <- file.path( CACHE$growth, "01_RGR_G_ASV_all.tsv" )
    GROWTH_BOOT_DIR <- file.path( CACHE$growth, "bootstrap_objects" )
    output_dir <- CACHE$strategy
    required_files <- c(GROWTH_FILE)
    missing_files <- required_files[!file.exists(required_files)]
    if (length(missing_files) > 0) {
      stop( "Missing required input files:\n", paste(missing_files, collapse = "\n") )
    }
    
    VOLCANO_COLORS <- c(
      Pos = "#D73027",
      Neg = "#4575B4",
      Neutral = "#8C8C8C"
    )
    STATUS_LEVELS <- c( "Pos", "Neg", "Neutral", "Not testable" )
    STATUS_COLORS <- c(
      Pos = "#D73027",
      Neg = "#4575B4",
      Neutral = "#8C8C8C",
      `Not testable` = "#EAFF9F"
    )
    
    safe_boot_z <- function(x) {
      x <- x[is.finite(x)]
      if (length(x) < 2) {
        return(NA_real_)
      }
      
      s <- sd(x)
      if (!is.finite(s) || s <= 0) {
        return(NA_real_)
      }
      
      abs(median(x) / s)
    }
    
    pct_text <- function(n, denom) {
      if (!is.finite(denom) || denom <= 0) {
        return("NA")
      }
      sprintf("%d (%.1f%%)", n, 100 * n / denom)
    }
    
    make_stage_label <- function(stage_code) {
      FIRE_LABELS[[stage_code]]
    }
    
    growth <- fread( GROWTH_FILE, sep = "\t", header = TRUE )
    required_growth_cols <- c(
      "asv_id",
      "fire_trt",
      "tem",
      "comparison_id",
      "Active_ASV",
      "G_ASV",
      "Phylum",
      "Phylum_group"
    )
    missing_growth_cols <- setdiff(required_growth_cols, names(growth))
    if (length(missing_growth_cols) > 0) {
      stop( "Growth table missing columns: ", paste(missing_growth_cols, collapse = ", ") )
    }
    
    growth[, asv_id := as.character(asv_id)]
    growth[, fire_trt := as.character(fire_trt)]
    growth[, tem := as.character(tem)]
    growth[, comparison_id := as.character(comparison_id)]
    growth[, Active_ASV := as.logical(Active_ASV)]
    growth[, G_ASV := as.numeric(G_ASV)]
    growth[, Phylum := as.character(Phylum)]
    growth[, Phylum_group := as.character(Phylum_group)]
    if (any(!growth$fire_trt %chin% FIRE_LEVELS)) {
      stop("Unexpected fire_trt values in growth table.")
    }
    
    if (any(!growth$tem %chin% TEMP_LEVELS)) {
      stop("Unexpected tem values in growth table.")
    }
    
    # Classify temperature responses only for comparable ASVs.
    stage_results_list <- list()
    volcano_plot_data_list <- list()
    stage_summary_list <- list()
    stage_qc_list <- list()
    for (i in seq_along(FIRE_LEVELS)) {
      current_fire <- FIRE_LEVELS[i]
      message("Calculating DeltaG strategies for ", current_fire)
      dt21 <- growth[
        fire_trt == current_fire & tem == "T_low",
        .(
          asv_id,
          Active_T_low = Active_ASV,
          G_T_low = G_ASV,
          Phylum_T_low = Phylum,
          Phylum_group_T_low = Phylum_group
        )
      ]
      dt28 <- growth[
        fire_trt == current_fire & tem == "T_high",
        .(
          asv_id,
          Active_T_high = Active_ASV,
          G_T_high = G_ASV,
          Phylum_T_high = Phylum,
          Phylum_group_T_high = Phylum_group
        )
      ]
      merged_stage <- merge( dt21, dt28, by = "asv_id", all = TRUE, sort = FALSE )
      merged_stage[, fire_trt := current_fire]
      merged_stage[, stage_label := make_stage_label(current_fire)]
      merged_stage[, testable_T_low := !is.na(G_T_low)]
      merged_stage[, testable_T_high := !is.na(G_T_high)]
      merged_stage[, testable_both := testable_T_low & testable_T_high]
      merged_stage[, Active_T_low := fifelse(is.na(Active_T_low), FALSE, Active_T_low)]
      merged_stage[, Active_T_high := fifelse(is.na(Active_T_high), FALSE, Active_T_high)]
      merged_stage[, active_any := Active_T_low | Active_T_high]
      merged_stage[, Phylum := fifelse(!is.na(Phylum_T_low), Phylum_T_low, Phylum_T_high)]
      merged_stage[, Phylum_group := fifelse(!is.na(Phylum_group_T_low), Phylum_group_T_low, Phylum_group_T_high)]
      boot_file_21 <- file.path( GROWTH_BOOT_DIR, paste0("G_bootstrap_", current_fire, "_T_low.rds") )
      boot_file_28 <- file.path( GROWTH_BOOT_DIR, paste0("G_bootstrap_", current_fire, "_T_high.rds") )
      if (!file.exists(boot_file_21) || !file.exists(boot_file_28)) {
        stop( "Missing G bootstrap object for ", current_fire, ". Expected: ", boot_file_21, " and ", boot_file_28 )
      }
      
      boot21 <- readRDS(boot_file_21)
      boot28 <- readRDS(boot_file_28)
      req_boot_items <- c("asv_id", "G_boot")
      if (!all(req_boot_items %in% names(boot21)) || !all(req_boot_items %in% names(boot28))) {
        stop("Bootstrap objects must contain asv_id and G_boot.")
      }
      
      common_asv <- intersect( as.character(boot21$asv_id), as.character(boot28$asv_id) )
      common_asv <- sort(common_asv)
      if (length(common_asv) == 0) {
        stop("No overlapping ASVs between T_low and T_high bootstrap objects for ", current_fire)
      }
      
      idx21 <- match(common_asv, as.character(boot21$asv_id))
      idx28 <- match(common_asv, as.character(boot28$asv_id))
      G21_boot <- boot21$G_boot[idx21, , drop = FALSE]
      G28_boot <- boot28$G_boot[idx28, , drop = FALSE]
      storage.mode(G21_boot) <- "double"
      storage.mode(G28_boot) <- "double"
      set.seed(RANDOM_SEED + i)
      delta_boot <- combine_independent_bootstrap_difference(
        mat_hot = G28_boot,
        mat_cold = G21_boot,
        n_boot = N_BOOT_COMBINE
      )
      delta_q <- t( apply( delta_boot, 1, safe_quantiles, probs = CI_PROBS ) )
      colnames(delta_q) <- c( "DeltaG_CI_low_90", "DeltaG_median", "DeltaG_CI_high_90" )
      delta_dt <- data.table(
        asv_id = common_asv,
        DeltaG_CI_low_90 = delta_q[, 1],
        DeltaG = delta_q[, 2],
        DeltaG_CI_high_90 = delta_q[, 3],
        DeltaG_boot_z = apply(delta_boot, 1, safe_boot_z),
        prop_boot_DeltaG_gt_0 = rowMeans(delta_boot > 0, na.rm = TRUE),
        prop_boot_DeltaG_lt_0 = rowMeans(delta_boot < 0, na.rm = TRUE),
        DeltaG_boot_sd = apply(delta_boot, 1, sd, na.rm = TRUE)
      )
      delta_dt[, Delta_strategy := classify_delta_strategy(DeltaG_CI_low_90, DeltaG_CI_high_90)]
      merged_stage <- merge( merged_stage, delta_dt, by = "asv_id", all.x = TRUE, sort = FALSE )
      merged_stage[ testable_both == FALSE | active_any == FALSE, Stage_status := "Not testable" ]
      merged_stage[ testable_both == TRUE & active_any == TRUE, Stage_status := Delta_strategy ]
      merged_stage[ , Stage_status := factor(Stage_status, levels = STATUS_LEVELS) ]
      comparable_dt <- merged_stage[ testable_both == TRUE & active_any == TRUE ]
      comparable_dt[, Volcano_group := as.character(Stage_status)]
      comparable_dt[, Volcano_group := factor(Volcano_group, levels = c("Pos", "Neg", "Neutral"))]
      counts_dt <- comparable_dt[ , .N, by = Volcano_group ]
      counts_all <- data.table(
        Volcano_group = factor(c("Pos", "Neg", "Neutral"), levels = c("Pos", "Neg", "Neutral"))
      )
      counts_dt <- merge( counts_all, counts_dt, by = "Volcano_group", all.x = TRUE, sort = FALSE )
      counts_dt[is.na(N), N := 0L]
      denom_n <- nrow(comparable_dt)
      ann_label <- paste0(
        "Pos: ", pct_text(counts_dt[Volcano_group == "Pos", N], denom_n), "\n",
        "Neg: ", pct_text(counts_dt[Volcano_group == "Neg", N], denom_n), "\n",
        "Neutral: ", pct_text(counts_dt[Volcano_group == "Neutral", N], denom_n), "\n",
        "N comparable = ", denom_n
      )
      volcano_plot_data_list[[length(volcano_plot_data_list) + 1L]] <- comparable_dt
      stage_summary_list[[length(stage_summary_list) + 1L]] <- data.table(
        fire_trt = current_fire,
        stage_label = make_stage_label(current_fire),
        comparable_ASVs = denom_n,
        Pos_ASVs = counts_dt[Volcano_group == "Pos", N],
        Neg_ASVs = counts_dt[Volcano_group == "Neg", N],
        Neutral_ASVs = counts_dt[Volcano_group == "Neutral", N],
        Pos_percent = if (denom_n > 0) 100 * counts_dt[Volcano_group == "Pos", N] / denom_n else NA_real_,
        Neg_percent = if (denom_n > 0) 100 * counts_dt[Volcano_group == "Neg", N] / denom_n else NA_real_,
        Neutral_percent = if (denom_n > 0) 100 * counts_dt[Volcano_group == "Neutral", N] / denom_n else NA_real_,
        facet_label = ann_label
      )
      stage_qc_list[[length(stage_qc_list) + 1L]] <- data.table(
        fire_trt = current_fire,
        testable_both_ASVs = sum(merged_stage$testable_both, na.rm = TRUE),
        active_any_among_testable_both = sum(merged_stage$testable_both & merged_stage$active_any, na.rm = TRUE),
        inactive_among_testable_both = sum(merged_stage$testable_both & !merged_stage$active_any, na.rm = TRUE),
        not_testable_ASVs = sum(!merged_stage$testable_both, na.rm = TRUE),
        shared_bootstrap_ASVs = length(common_asv),
        n_boot_DeltaG = N_BOOT_COMBINE
      )
      stage_results_list[[length(stage_results_list) + 1L]] <- merged_stage
      rm(boot21, boot28, G21_boot, G28_boot, delta_boot, delta_q, delta_dt)
      gc()
    }
    
    stage_results <- rbindlist(stage_results_list, use.names = TRUE, fill = TRUE)
    volcano_plot_data <- rbindlist(volcano_plot_data_list, use.names = TRUE, fill = TRUE)
    stage_summary <- rbindlist(stage_summary_list, use.names = TRUE, fill = TRUE)
    stage_qc <- rbindlist(stage_qc_list, use.names = TRUE, fill = TRUE)
    volcano_plot_data[, fire_trt := factor(fire_trt, levels = FIRE_LEVELS)]
    volcano_plot_data[, stage_label := factor(stage_label, levels = FIRE_LABELS[FIRE_LEVELS])]
    stage_summary[, fire_trt := factor(fire_trt, levels = FIRE_LEVELS)]
    stage_summary[, stage_label := factor(stage_label, levels = FIRE_LABELS[FIRE_LEVELS])]
    fwrite( stage_results, file.path(output_dir, "01_DeltaG_strategy_all.tsv"), sep = "\t" )
    fwrite( stage_summary, file.path(output_dir, "02_DeltaG_strategy_summary_by_stage.tsv"), sep = "\t" )
    fwrite( stage_qc, file.path(output_dir, "03_DeltaG_bootstrap_summary_QC.tsv"), sep = "\t" )
    volcano_plot_data[, DeltaG_signed_log10 := signed_log10_abs(DeltaG)]
    ann_dt <- copy(stage_summary)
    p_volcano <- ggplot(
      volcano_plot_data,
      aes( x = DeltaG_signed_log10, y = DeltaG_boot_z, color = Volcano_group )
    ) +
      geom_vline( xintercept = 0, linetype = "dashed", linewidth = 0.5, color = "grey45" ) +
      geom_point( alpha = 0.80, size = 1.4 ) +
      geom_text(
        data = ann_dt,
        aes( x = Inf, y = Inf, label = facet_label ),
        inherit.aes = FALSE,
        hjust = 1.05,
        vjust = 1.15,
        size = 3.1,
        lineheight = 1.05
      ) +
      facet_wrap( ~ stage_label, nrow = 1, scales = "free" ) +
      scale_color_manual( values = VOLCANO_COLORS, drop = FALSE ) +
      scale_x_continuous( breaks = scales::breaks_pretty(n = 6), expand = expansion(mult = c(0.1, 0.1)) ) +
      scale_y_continuous( expand = expansion(mult = c(0.02, 0.22)) ) +
      labs(
        x = expression("Signed log"[10] * "(1 + |" * Delta * "G|)"),
        y = expression("|Bootstrap z-score|"),
        color = NULL
      ) +
      theme_classic(base_size = 12) +
      theme(
        strip.background = element_rect( fill = "grey95", color = "black", linewidth = 0.6 ),
        strip.text = element_text( size = 11, face = "bold" ),
        axis.title = element_text( size = 12, face = "bold" ),
        axis.text = element_text( size = 10, color = "black" ),
        legend.position = "top",
        panel.border = element_rect( fill = NA, color = "black", linewidth = 0.6 ),
        plot.margin = margin(8, 12, 8, 8)
      )
    save_figure(
      filename = file.path(output_dir, "Fig_A_volcano_3panel_signed_log10.pdf"),
      plot = p_volcano,
      width = 12.5,
      height = 4.8,
      units = "in"
    )
    TRACKED_UB_LEVELS <- c("Pos", "Neg", "Neutral")
    LABEL_MIN_PERCENT <- 2.0
    stage_status_wide <- dcast(
      stage_results[ , .( asv_id, Phylum, Phylum_group, fire_trt, Stage_status ) ],
      asv_id + Phylum + Phylum_group ~ fire_trt,
      value.var = "Stage_status"
    )
    required_cast_stage_cols <- c("UB", "F3M", "F3Y")
    missing_cast_stage_cols <- setdiff(required_cast_stage_cols, names(stage_status_wide))
    if (length(missing_cast_stage_cols) > 0) {
      stop(
        "dcast() did not create all expected stage columns: ",
        paste(missing_cast_stage_cols, collapse = ", "),
        ". Available columns are: ",
        paste(names(stage_status_wide), collapse = ", ")
      )
    }
    
    setnames(
      stage_status_wide,
      old = required_cast_stage_cols,
      new = paste0("Stage_status_", required_cast_stage_cols)
    )
    required_status_cols <- paste0("Stage_status_", required_cast_stage_cols)
    if (!all(required_status_cols %in% names(stage_status_wide))) {
      stop("Failed to create Stage_status_UB / Stage_status_F3M / Stage_status_F3Y columns.")
    }
    
    pre_fire_classified_asv <- stage_results[
      fire_trt == "UB" & as.character(Stage_status) %chin% TRACKED_UB_LEVELS,
      unique(asv_id)
    ]
    if (length(pre_fire_classified_asv) == 0) {
      stop( "No pre-fire Pos/Neg/Neutral ASVs were identified; ", "cannot build alluvial plot." )
    }
    
    tracked_wide <- stage_status_wide[ asv_id %chin% pre_fire_classified_asv ]
    if (anyNA(tracked_wide$Stage_status_UB)) {
      stop("Some tracked ASVs unexpectedly lack a pre-fire Stage_status.")
    }
    
    if (any(!as.character(tracked_wide$Stage_status_UB) %chin% TRACKED_UB_LEVELS)) {
      stop("Tracked UB pool contains statuses outside Pos/Neg/Neutral.")
    }
    
    tracked_wide[ is.na(Stage_status_F3M), Stage_status_F3M := "Not testable" ]
    tracked_wide[ is.na(Stage_status_F3Y), Stage_status_F3Y := "Not testable" ]
    tracked_wide[, Stage_status_UB := factor( as.character(Stage_status_UB), levels = STATUS_LEVELS )]
    tracked_wide[, Stage_status_F3M := factor( as.character(Stage_status_F3M), levels = STATUS_LEVELS )]
    tracked_wide[, Stage_status_F3Y := factor( as.character(Stage_status_F3Y), levels = STATUS_LEVELS )]
    N_TRACKED <- nrow(tracked_wide)
    transition_UB_F3M <- tracked_wide[ , .N, by = .( UB = Stage_status_UB, F3M = Stage_status_F3M ) ][order(UB, F3M)]
    transition_UB_F3M[, percent_of_tracked_pool := 100 * N / N_TRACKED]
    transition_UB_F3M[ , percent_within_UB_origin := 100 * N / sum(N), by = UB ]
    transition_F3M_F3Y <- tracked_wide[ , .N, by = .( F3M = Stage_status_F3M, F3Y = Stage_status_F3Y ) ][order(F3M, F3Y)]
    transition_F3M_F3Y[, percent_of_tracked_pool := 100 * N / N_TRACKED]
    transition_F3M_F3Y[ , percent_within_F3M_status := 100 * N / sum(N), by = F3M ]
    transition_UB_F3Y <- tracked_wide[ , .N, by = .( UB = Stage_status_UB, F3Y = Stage_status_F3Y ) ][order(UB, F3Y)]
    transition_UB_F3Y[, percent_of_tracked_pool := 100 * N / N_TRACKED]
    transition_UB_F3Y[ , percent_within_UB_origin := 100 * N / sum(N), by = UB ]
    path_table <- tracked_wide[
      ,
      .N,
      by = .( UB = Stage_status_UB, F3M = Stage_status_F3M, F3Y = Stage_status_F3Y )
    ][order(UB, F3M, F3Y)]
    path_table[, percent_of_tracked_pool := 100 * N / N_TRACKED]
    path_table[ , percent_within_UB_origin := 100 * N / sum(N), by = UB ]
    path_table[, alluvium_id := seq_len(.N)]
    stage_count_long <- rbindlist(
      list(
        data.table(Stage = "UB",  Status = tracked_wide$Stage_status_UB),
        data.table(Stage = "F3M", Status = tracked_wide$Stage_status_F3M),
        data.table(Stage = "F3Y", Status = tracked_wide$Stage_status_F3Y)
      ),
      use.names = TRUE,
      fill = TRUE
    )
    stage_node_counts <- stage_count_long[ , .N, by = .(Stage, Status) ]
    stage_grid <- CJ( Stage = c("F3M", "UB", "F3Y"), Status = STATUS_LEVELS, unique = TRUE )
    stage_grid[, Stage := factor(Stage, levels = c("F3M", "UB", "F3Y"))]
    stage_grid[, Status := factor(Status, levels = STATUS_LEVELS)]
    stage_node_counts[, Stage := factor(Stage, levels = c("F3M", "UB", "F3Y"))]
    stage_node_counts[, Status := factor(as.character(Status), levels = STATUS_LEVELS)]
    stage_node_counts <- merge(
      stage_grid,
      stage_node_counts,
      by = c("Stage", "Status"),
      all.x = TRUE,
      sort = FALSE
    )
    stage_node_counts[is.na(N), N := 0L]
    stage_node_counts[, percent_of_tracked_pool := 100 * N / N_TRACKED]
    stage_node_counts[, show_label := N > 0 & percent_of_tracked_pool >= LABEL_MIN_PERCENT]
    stage_node_counts[, label := fifelse(
      show_label,
      sprintf( "%s\nn = %d\n%.1f%%", as.character(Status), N, percent_of_tracked_pool ),
      ""
    )]
    stage_node_counts[, label_colour := fifelse(
      as.character(Status) %chin% c("Pos", "Neg", "Neutral"),
      "white",
      "black"
    )]
    setorder(stage_node_counts, Stage, Status)
    stage_node_counts[ , ymin := c(0, head(cumsum(N), -1)), by = Stage ]
    stage_node_counts[, ymax := cumsum(N), by = Stage]
    stage_node_counts[, ymid := (ymin + ymax) / 2]
    origin_fate_summary <- tracked_wide[
      ,
      {
        origin_status <- as.character(Stage_status_UB[1])
        .(
          origin_n = .N,
          same_strategy_F3M_n = sum(as.character(Stage_status_F3M) == origin_status),
          same_strategy_F3Y_n = sum(as.character(Stage_status_F3Y) == origin_status),
          Neutral_F3M_n = sum(as.character(Stage_status_F3M) == "Neutral"),
          Neutral_F3Y_n = sum(as.character(Stage_status_F3Y) == "Neutral"),
          Not_testable_F3M_n = sum(as.character(Stage_status_F3M) == "Not testable"),
          Not_testable_F3Y_n = sum(as.character(Stage_status_F3Y) == "Not testable")
        )
      },
      by = .(UB = Stage_status_UB)
    ]
    origin_fate_summary[, same_strategy_F3M_percent := 100 * same_strategy_F3M_n / origin_n]
    origin_fate_summary[, same_strategy_F3Y_percent := 100 * same_strategy_F3Y_n / origin_n]
    origin_fate_summary[, Neutral_F3M_percent := 100 * Neutral_F3M_n / origin_n]
    origin_fate_summary[, Neutral_F3Y_percent := 100 * Neutral_F3Y_n / origin_n]
    origin_fate_summary[, Not_testable_F3M_percent := 100 * Not_testable_F3M_n / origin_n]
    origin_fate_summary[, Not_testable_F3Y_percent := 100 * Not_testable_F3Y_n / origin_n]
    alluvial_conservation_qc <- data.table(
      check = c( "Tracked pool", "Sum path_table N", "UB node total", "F3M node total", "F3Y node total" ),
      N = c(
        N_TRACKED,
        sum(path_table$N),
        stage_node_counts[Stage == "UB", sum(N)],
        stage_node_counts[Stage == "F3M", sum(N)],
        stage_node_counts[Stage == "F3Y", sum(N)]
      )
    )
    if (any(alluvial_conservation_qc$N != N_TRACKED)) {
      fwrite( alluvial_conservation_qc, file.path(output_dir, "QC_alluvial_conservation_FAILED.tsv"), sep = "\t" )
      stop("Alluvial conservation check failed: ASV counts are not conserved across stages.")
    }
    
    fwrite( tracked_wide, file.path(output_dir, "04_tracked_pre_fire_classified_ASV_status_table.tsv"), sep = "\t" )
    fwrite(
      transition_UB_F3M,
      file.path(output_dir, "05_transition_matrix_UB_to_F3M_allStrategies.tsv"),
      sep = "\t"
    )
    fwrite(
      transition_F3M_F3Y,
      file.path(output_dir, "06_transition_matrix_F3M_to_F3Y_allStrategies.tsv"),
      sep = "\t"
    )
    fwrite(
      transition_UB_F3Y,
      file.path(output_dir, "07_transition_matrix_UB_to_F3Y_allStrategies.tsv"),
      sep = "\t"
    )
    fwrite( path_table, file.path(output_dir, "08_transition_paths_UB_F3M_F3Y_allStrategies.tsv"), sep = "\t" )
    fwrite( origin_fate_summary, file.path(output_dir, "09_origin_strategy_fate_summary.tsv"), sep = "\t" )
    fwrite( stage_node_counts, file.path(output_dir, "10_alluvial_node_counts.tsv"), sep = "\t" )
    fwrite( alluvial_conservation_qc, file.path(output_dir, "10b_alluvial_conservation_QC.tsv"), sep = "\t" )
    PLOT_STAGE_LEVELS <- c("F3M", "UB", "F3Y")
    left_flow_table <- copy(transition_UB_F3M)
    left_flow_table <- left_flow_table[N > 0]
    left_flow_table[, z_order := fcase(
      as.character(F3M) == "Not testable", 1L,
      as.character(F3M) == "Neutral",      2L,
      as.character(F3M) == "Neg",          3L,
      as.character(F3M) == "Pos",          4L,
      default = 0L
    )]
    setorder(left_flow_table, z_order)
    left_flow_table[, flow_id := factor(paste0("L", seq_len(.N)), levels = paste0("L", seq_len(.N)))]
    left_flow_table[, UB := factor(as.character(UB), levels = STATUS_LEVELS)]
    left_flow_table[, F3M := factor(as.character(F3M), levels = STATUS_LEVELS)]
    left_lodes <- rbindlist(
      list(
        left_flow_table[, .(Stage = "F3M", Status = as.character(F3M), alluvium_id = flow_id, N = N, Flow_status = as.character(F3M))],
        left_flow_table[, .(Stage = "UB",  Status = as.character(UB),  alluvium_id = flow_id, N = N, Flow_status = as.character(F3M))]
      ),
      use.names = TRUE, fill = TRUE
    )
    right_flow_table <- copy(transition_UB_F3Y)
    right_flow_table <- right_flow_table[N > 0]
    right_flow_table[, z_order := fcase(
      as.character(F3Y) == "Not testable", 1L,
      as.character(F3Y) == "Neutral",      2L,
      as.character(F3Y) == "Neg",          3L,
      as.character(F3Y) == "Pos",          4L,
      default = 0L
    )]
    setorder(right_flow_table, z_order)
    right_flow_table[, flow_id := factor(paste0("R", seq_len(.N)), levels = paste0("R", seq_len(.N)))]
    right_flow_table[, UB := factor(as.character(UB), levels = STATUS_LEVELS)]
    right_flow_table[, F3Y := factor(as.character(F3Y), levels = STATUS_LEVELS)]
    right_lodes <- rbindlist(
      list(
        right_flow_table[, .(Stage = "UB",  Status = as.character(UB),  alluvium_id = flow_id, N = N, Flow_status = as.character(F3Y))],
        right_flow_table[, .(Stage = "F3Y", Status = as.character(F3Y), alluvium_id = flow_id, N = N, Flow_status = as.character(F3Y))]
      ),
      use.names = TRUE, fill = TRUE
    )
    strata_input <- copy(path_table)
    strata_input[, UB := factor(as.character(UB), levels = STATUS_LEVELS)]
    strata_input[, F3M := factor(as.character(F3M), levels = STATUS_LEVELS)]
    strata_input[, F3Y := factor(as.character(F3Y), levels = STATUS_LEVELS)]
    strata_lodes <- data.table::as.data.table(
      ggalluvial::to_lodes_form(
        strata_input,
        axes = c("F3M", "UB", "F3Y"),
        key = "Stage",
        value = "Status",
        id = "alluvium_id"
      )
    )
    for (nm in c("left_lodes", "right_lodes", "strata_lodes")) {
      x <- get(nm)
      x[, Stage := factor(as.character(Stage), levels = PLOT_STAGE_LEVELS)]
      x[, Status := factor(as.character(Status), levels = STATUS_LEVELS)]
      if ("Flow_status" %in% names(x)) {
        x[, Flow_status := factor(as.character(Flow_status), levels = STATUS_LEVELS)]
      }
      assign(nm, x)
    }
    
    node_label_dt <- stage_node_counts[show_label == TRUE]
    node_label_dt[, Stage := factor(as.character(Stage), levels = PLOT_STAGE_LEVELS)]
    subtitle_text <- sprintf(
      paste0(
        "Tracking pre-fire classified thermal strategies (n = %d); ",
        "Pre-fire centered; ribbons show direct fates at 3 months and 3 years; ",
        "labels shown for >= %.1f%% of tracked pool"
      ),
      N_TRACKED,
      LABEL_MIN_PERCENT
    )
    p_alluvial <- ggplot() +
      geom_flow(
        data = left_lodes,
        aes(x = Stage, stratum = Status, alluvium = alluvium_id, y = N, fill = Flow_status),
        alpha = 0.70,
        color = "transparent",
        linewidth = 0,
        width = 0.35,
        curve_type = "sine",
        reverse = FALSE,
        show.legend = FALSE
      ) +
      geom_flow(
        data = right_lodes,
        aes(x = Stage, stratum = Status, alluvium = alluvium_id, y = N, fill = Flow_status),
        alpha = 0.70,
        color = "transparent",
        linewidth = 0,
        width = 0.35,
        curve_type = "sine",
        reverse = FALSE,
        show.legend = FALSE
      ) +
      geom_stratum(
        data = strata_lodes,
        aes(x = Stage, stratum = Status, y = N, fill = Status),
        color = "black",
        linewidth = 0.40,
        width = 0.35,
        reverse = FALSE
      ) +
      geom_text(
        data = node_label_dt,
        aes(x = Stage, y = ymid, label = label, color = label_colour),
        inherit.aes = FALSE,
        size = 2.9,
        lineheight = 0.95,
        fontface = "plain"
      ) +
      scale_fill_manual( values = STATUS_COLORS, breaks = STATUS_LEVELS, drop = FALSE ) +
      scale_color_identity() +
      scale_x_discrete(
        limits = PLOT_STAGE_LEVELS,
        labels = c(F3M = "3 months", UB = "Unburned", F3Y = "3 years"),
        expand = expansion(add = c(0.22, 0.22))
      ) +
      scale_y_continuous( labels = label_number(big.mark = ","), expand = expansion(mult = c(0.01, 0.04)) ) +
      labs( x = NULL, y = "Number of ASVs", fill = NULL, subtitle = subtitle_text ) +
      theme_classic(base_size = 12) +
      theme(
        axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10, color = "black"),
        plot.subtitle = element_text(size = 9.8, hjust = 0.5, margin = margin(b = 7)),
        legend.position = "top",
        legend.key.width = grid::unit(0.9, "cm"),
        legend.key.height = grid::unit(0.45, "cm"),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 0.6),
        plot.margin = margin(8, 12, 8, 8)
      ) +
      guides(fill = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(alpha = 1)))
    save_figure(
      filename = file.path(output_dir, "Fig_B_alluvial_pre_fire_all_thermal_strategies.pdf"),
      plot = p_alluvial,
      width = 11.5,
      height = 7.0,
      units = "in"
    )
    tracked_pos_wide <- tracked_wide[as.character(Stage_status_UB) == "Pos"]
    N_TRACKED_POS <- nrow(tracked_pos_wide)
    if (N_TRACKED_POS > 0L) {
      pos_transition_UB_F3M <- tracked_pos_wide[, .N, by = .(UB = Stage_status_UB, F3M = Stage_status_F3M)][order(UB, F3M)]
      pos_transition_UB_F3M[, percent_of_pre_fire_Pos_pool := 100 * N / N_TRACKED_POS]
      pos_transition_UB_F3Y <- tracked_pos_wide[, .N, by = .(UB = Stage_status_UB, F3Y = Stage_status_F3Y)][order(UB, F3Y)]
      pos_transition_UB_F3Y[, percent_of_pre_fire_Pos_pool := 100 * N / N_TRACKED_POS]
      pos_path_table <- tracked_pos_wide[, .N, by = .(UB = Stage_status_UB, F3M = Stage_status_F3M, F3Y = Stage_status_F3Y)][order(UB, F3M, F3Y)]
      pos_path_table[, percent_of_pre_fire_Pos_pool := 100 * N / N_TRACKED_POS]
      pos_path_table[, alluvium_id := seq_len(.N)]
      pos_stage_count_long <- rbindlist(list(
        data.table(Stage = "UB", Status = tracked_pos_wide$Stage_status_UB),
        data.table(Stage = "F3M", Status = tracked_pos_wide$Stage_status_F3M),
        data.table(Stage = "F3Y", Status = tracked_pos_wide$Stage_status_F3Y)
      ), use.names = TRUE, fill = TRUE)
      pos_stage_node_counts <- pos_stage_count_long[, .N, by = .(Stage, Status)]
      pos_stage_grid <- CJ(Stage = PLOT_STAGE_LEVELS, Status = STATUS_LEVELS, unique = TRUE)
      pos_stage_grid[, Stage := factor(Stage, levels = PLOT_STAGE_LEVELS)]
      pos_stage_grid[, Status := factor(Status, levels = STATUS_LEVELS)]
      pos_stage_node_counts[, Stage := factor(Stage, levels = PLOT_STAGE_LEVELS)]
      pos_stage_node_counts[, Status := factor(as.character(Status), levels = STATUS_LEVELS)]
      pos_stage_node_counts <- merge(pos_stage_grid, pos_stage_node_counts, by = c("Stage", "Status"), all.x = TRUE, sort = FALSE)
      pos_stage_node_counts[is.na(N), N := 0L]
      pos_stage_node_counts[, percent_of_pre_fire_Pos_pool := 100 * N / N_TRACKED_POS]
      pos_stage_node_counts[, show_label := N > 0 & percent_of_pre_fire_Pos_pool >= LABEL_MIN_PERCENT]
      pos_stage_node_counts[, label := fifelse(show_label, sprintf("%s\nn = %d\n%.1f%%", as.character(Status), N, percent_of_pre_fire_Pos_pool), "")]
      pos_stage_node_counts[, label_colour := fifelse(as.character(Status) %chin% c("Pos", "Neg", "Neutral"), "white", "black")]
      setorder(pos_stage_node_counts, Stage, Status)
      pos_stage_node_counts[, ymin := c(0, head(cumsum(N), -1)), by = Stage]
      pos_stage_node_counts[, ymax := cumsum(N), by = Stage]
      pos_stage_node_counts[, ymid := (ymin + ymax) / 2]
      pos_alluvial_conservation_qc <- data.table(
        check = c("Pre-fire Pos tracked pool", "Sum pos_path_table N", "UB node total", "F3M node total", "F3Y node total"),
        N = c(N_TRACKED_POS, sum(pos_path_table$N), pos_stage_node_counts[Stage == "UB", sum(N)], pos_stage_node_counts[Stage == "F3M", sum(N)], pos_stage_node_counts[Stage == "F3Y", sum(N)])
      )
      if (any(pos_alluvial_conservation_qc$N != N_TRACKED_POS)) {
        fwrite(pos_alluvial_conservation_qc, file.path(output_dir, "QC_FigC_pre_fire_Pos_conservation_FAILED.tsv"), sep = "\t")
        stop("Fig. C conservation check failed.")
      }
      
      fwrite(tracked_pos_wide, file.path(output_dir, "FigC_01_tracked_pre_fire_Pos_ASVs.tsv"), sep = "\t")
      fwrite(pos_transition_UB_F3M, file.path(output_dir, "FigC_02_transition_UB_to_F3M_pre_fire_Pos.tsv"), sep = "\t")
      fwrite(pos_transition_UB_F3Y, file.path(output_dir, "FigC_03_transition_UB_to_F3Y_pre_fire_Pos.tsv"), sep = "\t")
      fwrite(pos_path_table, file.path(output_dir, "FigC_04_three_stage_paths_pre_fire_Pos.tsv"), sep = "\t")
      fwrite(pos_stage_node_counts, file.path(output_dir, "FigC_05_node_counts_pre_fire_Pos.tsv"), sep = "\t")
      fwrite(pos_alluvial_conservation_qc, file.path(output_dir, "FigC_06_conservation_QC.tsv"), sep = "\t")
      pos_left_flow_table <- copy(pos_transition_UB_F3M)
      pos_left_flow_table <- pos_left_flow_table[N > 0]
      pos_left_flow_table[, z_order := fcase(
        as.character(F3M) == "Not testable", 1L,
        as.character(F3M) == "Neutral",      2L,
        as.character(F3M) == "Neg",          3L,
        as.character(F3M) == "Pos",          4L,
        default = 0L
      )]
      setorder(pos_left_flow_table, z_order)
      pos_left_flow_table[, flow_id := factor(paste0("PL", seq_len(.N)), levels = paste0("PL", seq_len(.N)))]
      pos_left_flow_table[, UB := factor(as.character(UB), levels = STATUS_LEVELS)]
      pos_left_flow_table[, F3M := factor(as.character(F3M), levels = STATUS_LEVELS)]
      pos_left_lodes <- rbindlist(
        list(
          pos_left_flow_table[, .(Stage = "F3M", Status = as.character(F3M), alluvium_id = flow_id, N = N, Flow_status = as.character(F3M))],
          pos_left_flow_table[, .(Stage = "UB",  Status = as.character(UB),  alluvium_id = flow_id, N = N, Flow_status = as.character(F3M))]
        ),
        use.names = TRUE, fill = TRUE
      )
      pos_right_flow_table <- copy(pos_transition_UB_F3Y)
      pos_right_flow_table <- pos_right_flow_table[N > 0]
      pos_right_flow_table[, z_order := fcase(
        as.character(F3Y) == "Not testable", 1L,
        as.character(F3Y) == "Neutral",      2L,
        as.character(F3Y) == "Neg",          3L,
        as.character(F3Y) == "Pos",          4L,
        default = 0L
      )]
      setorder(pos_right_flow_table, z_order)
      pos_right_flow_table[, flow_id := factor(paste0("PR", seq_len(.N)), levels = paste0("PR", seq_len(.N)))]
      pos_right_flow_table[, UB := factor(as.character(UB), levels = STATUS_LEVELS)]
      pos_right_flow_table[, F3Y := factor(as.character(F3Y), levels = STATUS_LEVELS)]
      pos_right_lodes <- rbindlist(
        list(
          pos_right_flow_table[, .(Stage = "UB",  Status = as.character(UB),  alluvium_id = flow_id, N = N, Flow_status = as.character(F3Y))],
          pos_right_flow_table[, .(Stage = "F3Y", Status = as.character(F3Y), alluvium_id = flow_id, N = N, Flow_status = as.character(F3Y))]
        ),
        use.names = TRUE, fill = TRUE
      )
      pos_strata_input <- copy(pos_path_table)
      pos_strata_input[, UB := factor(as.character(UB), levels = STATUS_LEVELS)]
      pos_strata_input[, F3M := factor(as.character(F3M), levels = STATUS_LEVELS)]
      pos_strata_input[, F3Y := factor(as.character(F3Y), levels = STATUS_LEVELS)]
      pos_strata_lodes <- data.table::as.data.table(
        ggalluvial::to_lodes_form(
          pos_strata_input,
          axes = c("F3M", "UB", "F3Y"),
          key = "Stage",
          value = "Status",
          id = "alluvium_id"
        )
      )
      for (nm in c("pos_left_lodes", "pos_right_lodes", "pos_strata_lodes")) {
        x <- get(nm)
        x[, Stage := factor(as.character(Stage), levels = PLOT_STAGE_LEVELS)]
        x[, Status := factor(as.character(Status), levels = STATUS_LEVELS)]
        if ("Flow_status" %in% names(x)) {
          x[, Flow_status := factor(as.character(Flow_status), levels = STATUS_LEVELS)]
        }
        assign(nm, x)
      }
      
      pos_node_label_dt <- pos_stage_node_counts[show_label == TRUE]
      pos_node_label_dt[, Stage := factor(as.character(Stage), levels = PLOT_STAGE_LEVELS)]
      pos_subtitle_text <- sprintf(
        paste0(
          "Tracking pre-fire Pos-ASVs only (n = %d); ",
          "Pre-fire centered; ribbons show direct fates at 3 months and 3 years; ",
          "labels shown for >= %.1f%% of the pre-fire Pos pool"
        ),
        N_TRACKED_POS,
        LABEL_MIN_PERCENT
      )
      p_alluvial_pos <- ggplot() +
        geom_flow(
          data = pos_left_lodes,
          aes(x = Stage, stratum = Status, alluvium = alluvium_id, y = N, fill = Flow_status),
          alpha = 0.70,
          color = "transparent",
          linewidth = 0,
          width = 0.35,
          curve_type = "sine",
          reverse = FALSE,
          show.legend = FALSE
        ) +
        geom_flow(
          data = pos_right_lodes,
          aes(x = Stage, stratum = Status, alluvium = alluvium_id, y = N, fill = Flow_status),
          alpha = 0.70,
          color = "transparent",
          linewidth = 0,
          width = 0.35,
          curve_type = "sine",
          reverse = FALSE,
          show.legend = FALSE
        ) +
        geom_stratum(
          data = pos_strata_lodes,
          aes(x = Stage, stratum = Status, y = N, fill = Status),
          color = "black",
          linewidth = 0.40,
          width = 0.35,
          reverse = FALSE
        ) +
        geom_text(
          data = pos_node_label_dt,
          aes(x = Stage, y = ymid, label = label, color = label_colour),
          inherit.aes = FALSE,
          size = 2.9,
          lineheight = 0.95,
          fontface = "plain"
        ) +
        scale_fill_manual( values = STATUS_COLORS, breaks = STATUS_LEVELS, drop = FALSE ) +
        scale_color_identity() +
        scale_x_discrete(
          limits = PLOT_STAGE_LEVELS,
          labels = c(F3M = "3 months", UB = "Unburned", F3Y = "3 years"),
          expand = expansion(add = c(0.22, 0.22))
        ) +
        scale_y_continuous( labels = label_number(big.mark = ","), expand = expansion(mult = c(0.01, 0.04)) ) +
        labs( x = NULL, y = "Number of ASVs", fill = NULL, subtitle = pos_subtitle_text ) +
        theme_classic(base_size = 12) +
        theme(
          axis.title = element_text(size = 12, face = "bold"),
          axis.text = element_text(size = 10, color = "black"),
          plot.subtitle = element_text(size = 9.8, hjust = 0.5, margin = margin(b = 7)),
          legend.position = "top",
          legend.key.width = grid::unit(0.9, "cm"),
          legend.key.height = grid::unit(0.45, "cm"),
          panel.border = element_rect(fill = NA, color = "black", linewidth = 0.6),
          plot.margin = margin(8, 12, 8, 8)
        ) +
        guides(fill = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(alpha = 1)))
      save_figure(
        filename = file.path(output_dir, "Fig_C_alluvial_pre_fire_Pos_ASVs.pdf"),
        plot = p_alluvial_pos,
        width = 11.5,
        height = 7.0,
        units = "in"
      )
    } else {
      message("No unburned positive responders; the positive-only transition figure is skipped.")
    }
    parameter_table <- data.table(
      Parameter = c(
        "Input growth table",
        "Input bootstrap directory",
        "Volcano x-axis",
        "Volcano y-axis",
        "Volcano denominator",
        "DeltaG formula",
        "DeltaG bootstrap construction",
        "Bootstrap iterations for DeltaG",
        "DeltaG CI",
        "Pos-ASV definition",
        "Neg-ASV definition",
        "Neutral-ASV definition",
        "Alluvial tracked pool",
        "Alluvial aggregation",
        "Alluvial flow color",
        "Alluvial statuses",
        "Alluvial node percentages denominator",
        "Alluvial in-node label threshold",
        "Fig C tracked pool",
        "Fig C layout",
        "Fig C node percentages denominator"
      ),
      Value = c(
        GROWTH_FILE,
        GROWTH_BOOT_DIR,
        "DeltaG = G_ASV(T_high) - G_ASV(T_low)",
        "|bootstrap z-score| = |median(DeltaG_boot) / sd(DeltaG_boot)|",
        "testable in both T_low and T_high AND active in at least one temperature",
        "G_ASV(T_high) - G_ASV(T_low)",
        "independent resampling of T_high and T_low bootstrap draws",
        as.character(N_BOOT_COMBINE),
        "90% (5th-95th percentile)",
        "DeltaG_CI_low_90 > 0",
        "DeltaG_CI_high_90 < 0",
        "CI overlaps 0",
        "all pre-fire Pos/Neg/Neutral ASVs",
        "Pre-fire-centered direct pairwise flows: F3M <- UB -> F3Y; ribbon width = ASV count",
        "outer-stage fate status: F3M on left and F3Y on right",
        paste(STATUS_LEVELS, collapse = ", "),
        "all tracked pre-fire classified ASVs",
        paste0(">= ", LABEL_MIN_PERCENT, "% of tracked pool"),
        "pre-fire Pos-ASVs only",
        "Pre-fire-centered direct pairwise flows: F3M <- UB -> F3Y",
        "pre-fire Pos-ASV pool"
      )
    )
    fwrite( parameter_table, file.path(output_dir, "11_analysis_parameters.tsv"), sep = "\t" )
    invisible(NULL)
  }
  
  # Phylo G
  run_phylo_g <- function() {
    TREE_FILE <- INPUTS$tree
    delta_g_candidates <- c( file.path( CACHE$strategy, "01_DeltaG_strategy_all.tsv" ) )
    delta_g_existing <- delta_g_candidates[file.exists(delta_g_candidates)]
    if (length(delta_g_existing) == 0L) {
      stop(
        "Cannot find 01_DeltaG_strategy_all.tsv in any supported Step-06 output folder. " ,
        "Run Step 06 first or update delta_g_candidates."
      )
    }
    
    DELTA_G_FILE <- delta_g_existing[1]
    GROWTH_FILE <- file.path( CACHE$growth, "01_RGR_G_ASV_all.tsv" )
    PALETTE_FILE <- file.path( CACHE$afe, "15_phylum_palette.tsv" )
    output_dir <- CACHE$phylo_g
    required_files <- c( TREE_FILE, DELTA_G_FILE, GROWTH_FILE )
    missing_files <- required_files[ !file.exists(required_files) ]
    if (length(missing_files) > 0) {
      stop( "Missing required files:\n", paste(missing_files, collapse = "\n") )
    }
    
    TREE_HEATMAP_ABS_QUANTILE <- 0.98
    ZERO_EDGE_MULTIPLIER <- 0.5
    set.seed(RANDOM_SEED)
    
    tree_raw <- ape::read.tree(TREE_FILE)
    if (is.null(tree_raw$tip.label) || length(tree_raw$tip.label) == 0) {
      stop("No tip labels found in rooted-tree.nwk.")
    }
    
    if (anyDuplicated(tree_raw$tip.label)) {
      dup <- unique(tree_raw$tip.label[duplicated(tree_raw$tip.label)])
      fwrite( data.table(asv_id = dup), file.path(output_dir, "QC_duplicate_tree_tip_labels.tsv"), sep = "\t" )
      stop("Duplicate tip labels detected in rooted-tree.nwk.")
    }
    
    if (is.null(tree_raw$edge.length)) {
      stop("The supplied tree has no branch lengths; the circular tree requires branch lengths.")
    }
    
    if (any(!is.finite(tree_raw$edge.length))) {
      stop("Non-finite branch lengths detected in rooted-tree.nwk.")
    }
    
    rooted_status <- ape::is.rooted(tree_raw)
    if (!rooted_status) {
      stop(
        paste0(
          "ape::is.rooted() reports that rooted-tree.nwk is not rooted. ",
          "Do not automatically root the tree here; verify the supplied tree."
        )
      )
    }
    
    nonpositive_edges <- which(tree_raw$edge.length <= 0)
    edge_replacement <- NA_real_
    tree <- tree_raw
    if (length(nonpositive_edges) > 0) {
      positive_edges <- tree$edge.length[ is.finite(tree$edge.length) & tree$edge.length > 0 ]
      if (length(positive_edges) == 0) {
        stop("Tree contains no positive branch lengths.")
      }
      
      edge_replacement <- min(positive_edges) * ZERO_EDGE_MULTIPLIER
      if (!is.finite(edge_replacement) || edge_replacement <= 0) {
        edge_replacement <- 1e-10
      }
      
      tree$edge.length[nonpositive_edges] <- edge_replacement
      warning(
        length(nonpositive_edges),
        " non-positive branch lengths were replaced by ",
        format(edge_replacement, scientific = TRUE),
        " for the tree display."
      )
    }
    
    stage_dt <- fread( DELTA_G_FILE, sep = "\t", header = TRUE )
    growth <- fread( GROWTH_FILE, sep = "\t", header = TRUE )
    required_stage_cols <- c(
      "asv_id",
      "fire_trt",
      "testable_both",
      "active_any",
      "DeltaG",
      "Stage_status",
      "Phylum",
      "Phylum_group"
    )
    missing_stage_cols <- setdiff( required_stage_cols, names(stage_dt) )
    if (length(missing_stage_cols) > 0) {
      stop( "06 DeltaG table missing columns: ", paste(missing_stage_cols, collapse = ", ") )
    }
    
    required_growth_cols <- c(
      "asv_id",
      "fire_trt",
      "tem",
      "Active_ASV",
      "N_TOTAL_unlabel_mean",
      "Phylum",
      "Phylum_group"
    )
    missing_growth_cols <- setdiff( required_growth_cols, names(growth) )
    if (length(missing_growth_cols) > 0) {
      stop( "04 growth table missing columns: ", paste(missing_growth_cols, collapse = ", ") )
    }
    
    stage_dt[, asv_id := as.character(asv_id)]
    stage_dt[, fire_trt := as.character(fire_trt)]
    stage_dt[, testable_both := as.logical(testable_both)]
    stage_dt[, active_any := as.logical(active_any)]
    stage_dt[, DeltaG := safe_numeric(DeltaG)]
    stage_dt[, Stage_status := as.character(Stage_status)]
    stage_dt[, Phylum := as.character(Phylum)]
    stage_dt[, Phylum_group := as.character(Phylum_group)]
    growth[, asv_id := as.character(asv_id)]
    growth[, fire_trt := as.character(fire_trt)]
    growth[, tem := as.character(tem)]
    growth[, Active_ASV := as.logical(Active_ASV)]
    growth[, N_TOTAL_unlabel_mean := safe_numeric(N_TOTAL_unlabel_mean)]
    growth[, Phylum := as.character(Phylum)]
    growth[, Phylum_group := as.character(Phylum_group)]
    stage_dt[ , in_tree := asv_id %chin% tree$tip.label ]
    stage_dt[ , phylo_eligible := testable_both == TRUE & active_any == TRUE & is.finite(DeltaG) & in_tree == TRUE ]
    stage_eligible <- stage_dt[ phylo_eligible == TRUE ]
    eligible_summary <- stage_dt[
      ,
      .(
        total_ASVs_in_stage_table = .N,
        testable_both = sum(testable_both, na.rm = TRUE),
        active_any = sum(testable_both & active_any, na.rm = TRUE),
        finite_DeltaG = sum(testable_both & active_any & is.finite(DeltaG), na.rm = TRUE),
        present_in_tree = sum(testable_both & active_any & is.finite(DeltaG) & in_tree, na.rm = TRUE),
        phylo_eligible_ASVs = sum(phylo_eligible, na.rm = TRUE)
      ),
      by = fire_trt
    ]
    eligible_summary[ , stage_label := FIRE_LABELS[fire_trt] ]
    fwrite( eligible_summary, file.path(output_dir, "01_tree_matching_and_eligibility_QC.tsv"), sep = "\t" )
    fwrite( stage_eligible, file.path(output_dir, "02_stage_specific_phylo_eligible_ASVs.tsv"), sep = "\t" )
    all_active_visual_ids <- sort(unique(stage_eligible$asv_id))
    visual_traits_long <- stage_eligible[ , .( asv_id, fire_trt, DeltaG, Stage_status, Phylum, Phylum_group ) ]
    visual_taxonomy <- unique( visual_traits_long[, .(asv_id, Phylum, Phylum_group)], by = "asv_id" )
    visual_traits_wide <- dcast( visual_traits_long, asv_id ~ fire_trt, value.var = c("DeltaG", "Stage_status") )
    visual_traits_wide <- merge( visual_traits_wide, visual_taxonomy, by = "asv_id", all.x = TRUE, sort = FALSE )
    for (st in FIRE_LEVELS) {
      raw_col <- paste0("DeltaG_", st)
      status_col <- paste0("Stage_status_", st)
      trans_col <- paste0("DeltaG_asinh_", st)
      if (!raw_col %in% names(visual_traits_wide)) {
        visual_traits_wide[, (raw_col) := NA_real_]
      }
      if (!status_col %in% names(visual_traits_wide)) {
        visual_traits_wide[, (status_col) := NA_character_]
      }
    }
    
    global_DeltaG_scale <- safe_asinh_scale(stage_eligible$DeltaG)
    for (st in FIRE_LEVELS) {
      raw_col <- paste0("DeltaG_", st)
      trans_col <- paste0("DeltaG_asinh_", st)
      visual_traits_wide[, (trans_col) := asinh(get(raw_col) / global_DeltaG_scale)]
    }
    
    shared_traits_wide <- visual_traits_wide
    abundance_rank <- growth[
      asv_id %chin% all_active_visual_ids,
      .( mean_endpoint_abundance = mean( N_TOTAL_unlabel_mean, na.rm = TRUE ) ),
      by = asv_id
    ]
    abundance_rank[ !is.finite(mean_endpoint_abundance), mean_endpoint_abundance := -Inf ]
    setorder( abundance_rank, -mean_endpoint_abundance, asv_id )
    plot_ids <- all_active_visual_ids
    n_plot <- length(plot_ids)
    plot_tree <- ape::drop.tip( tree, setdiff(tree$tip.label, plot_ids) )
    plot_meta <- shared_traits_wide[ asv_id %chin% plot_ids ]
    plot_meta <- merge( plot_meta, abundance_rank, by = "asv_id", all.x = TRUE, sort = FALSE )
    plot_meta[ is.na(Phylum_group) | Phylum_group == "", Phylum_group := "Unclassified" ]
    tree_signal_match_qc <- rbindlist(lapply(FIRE_LEVELS, function(st) {
      raw_col <- paste0("DeltaG_", st)
      ring_ids <- plot_meta[is.finite(get(raw_col)), unique(asv_id)]
      signal_ids <- stage_eligible[fire_trt == st, unique(asv_id)]
      data.table(
        fire_trt = st,
        eligible_ASVs = length(signal_ids),
        tree_ring_nonNA_ASVs = length(ring_ids),
        exact_same_ASV_set = setequal(ring_ids, signal_ids)
      )
    }))
    union_signal_ids <- sort(unique(stage_eligible$asv_id))
    tree_union_qc <- data.table(
      fire_trt = "UNION_TREE",
      eligible_ASVs = length(union_signal_ids),
      tree_ring_nonNA_ASVs = length(plot_tree$tip.label),
      exact_same_ASV_set = setequal(plot_tree$tip.label, union_signal_ids)
    )
    tree_signal_match_qc <- rbind(tree_signal_match_qc, tree_union_qc, fill = TRUE)
    if (!all(tree_signal_match_qc$exact_same_ASV_set)) {
      fwrite(tree_signal_match_qc, file.path(output_dir, "QC_tree_vs_eligible_ASV_sets_FAILED.tsv"), sep = "\t")
      stop("Circular-tree taxon set does not exactly match the stage-specific eligible sets.")
    }
    
    fwrite( tree_signal_match_qc, file.path(output_dir, "07b_tree_vs_eligible_ASV_set_QC.tsv"), sep = "\t" )
    for (st in FIRE_LEVELS) {
      raw_col <- paste0("DeltaG_", st)
      status_col <- paste0("Stage_status_", st)
      display_col <- paste0("DeltaG_display_", st)
      sig_col <- paste0("DeltaG_significant_", st)
      plot_meta[ , (display_col) := asinh(get(raw_col) / global_DeltaG_scale) ]
      plot_meta[ , (sig_col) := get(status_col) %chin% c("Pos", "Neg") ]
    }
    
    fwrite( plot_meta, file.path(output_dir, "07_circular_tree_display_ASVs.tsv"), sep = "\t" )
    plot_phyla <- sort(unique(plot_meta$Phylum_group))
    if (file.exists(PALETTE_FILE)) {
      palette_dt <- fread( PALETTE_FILE, sep = "\t", header = TRUE )
      if (all(c("Phylum_group", "color") %in% names(palette_dt))) {
        palette_dt[, Phylum_group := as.character(Phylum_group)]
        palette_dt[, color := as.character(color)]
      } else {
        palette_dt <- data.table( Phylum_group = character(), color = character() )
      }
    } else {
      palette_dt <- data.table( Phylum_group = character(), color = character() )
    }
    
    existing_colors <- setNames( palette_dt$color, palette_dt$Phylum_group )
    missing_phyla <- setdiff( plot_phyla, names(existing_colors) )
    if (length(missing_phyla) > 0) {
      extra_cols <- grDevices::hcl.colors( length(missing_phyla), palette = "Dark 3" )
      names(extra_cols) <- missing_phyla
      existing_colors <- c(existing_colors, extra_cols)
    }
    
    phylum_colors <- existing_colors[plot_phyla]
    fwrite(
      data.table( Phylum_group = plot_phyla, color = unname(phylum_colors) ),
      file.path(output_dir, "08_circular_tree_phylum_palette.tsv"),
      sep = "\t"
    )
    
    build_branch_phylum_groups <- function(tree_obj, tip_group_map) {
      n_tip <- length(tree_obj$tip.label)
      n_total <- n_tip + tree_obj$Nnode
      node_group <- rep(NA_character_, n_total)
      tip_group <- unname(tip_group_map[tree_obj$tip.label])
      tip_group[is.na(tip_group) | tip_group == ""] <- "Unclassified"
      node_group[seq_len(n_tip)] <- tip_group
      node_depth_all <- ape::node.depth.edgelength(tree_obj)
      internal_nodes <- seq.int(n_tip + 1L, n_total)
      internal_nodes <- internal_nodes[ order(node_depth_all[internal_nodes], decreasing = TRUE) ]
      for (node_i in internal_nodes) {
        child_nodes <- tree_obj$edge[ tree_obj$edge[, 1] == node_i, 2 ]
        if (length(child_nodes) == 0L) {
          node_group[node_i] <- "Mixed"
          next
        }
        
        child_groups <- node_group[child_nodes]
        if (
          length(child_groups) == length(child_nodes) && all(!is.na(child_groups)) &&
          length(unique(child_groups)) == 1L &&
          child_groups[1] != "Mixed"
        ) {
          node_group[node_i] <- child_groups[1]
        } else {
          node_group[node_i] <- "Mixed"
        }
      }
      
      out <- data.table(
        node = seq_len(n_total),
        Branch_group = node_group,
        node_type = ifelse( seq_len(n_total) <= n_tip, "tip", "internal" ),
        tip_label = NA_character_
      )
      out[ seq_len(n_tip), tip_label := tree_obj$tip.label ]
      out
    }
    
    tip_phylum_map <- setNames( as.character(plot_meta$Phylum_group), plot_meta$asv_id )
    missing_tip_phylum <- setdiff( plot_tree$tip.label, names(tip_phylum_map) )
    if (length(missing_tip_phylum) > 0L) {
      stop( "Missing Phylum_group for ", length(missing_tip_phylum), " Figure-A tree tips." )
    }
    
    branch_phylum_dt <- build_branch_phylum_groups( tree_obj = plot_tree, tip_group_map = tip_phylum_map )
    branch_colors <- c(
      phylum_colors,
      Mixed = "#BDBDBD"
    )
    if (any(!branch_phylum_dt$Branch_group %chin% names(branch_colors))) {
      bad_groups <- unique( branch_phylum_dt[ !Branch_group %chin% names(branch_colors), Branch_group ] )
      stop( "Tree branch groups missing from branch color palette: ", paste(bad_groups, collapse = ", ") )
    }
    
    fwrite( branch_phylum_dt, file.path(output_dir, "08b_circular_tree_branch_phylum_groups.tsv"), sep = "\t" )
    phylum_df <- data.frame(
      Phylum = factor(plot_meta$Phylum_group, levels = plot_phyla),
      row.names = plot_meta$asv_id,
      check.names = FALSE
    )
    
    get_display_values <- function(stage) {
      display_col <- paste0("DeltaG_display_", stage)
      as.numeric(plot_meta[[display_col]])
    }
    
    delta_display_df <- data.frame(
      `Unburned` = get_display_values("UB"),
      `3 months` = get_display_values("F3M"),
      `3 years` = get_display_values("F3Y"),
      row.names = plot_meta$asv_id,
      check.names = FALSE
    )
    all_shared_sig_display <- numeric(0)
    for (st in FIRE_LEVELS) {
      trans_col <- paste0("DeltaG_asinh_", st)
      x <- as.numeric(shared_traits_wide[[trans_col]])
      all_shared_sig_display <- c(all_shared_sig_display, x[is.finite(x)])
    }
    
    if (length(all_shared_sig_display) == 0) {
      all_shared_sig_display <- unlist(shared_traits_wide[, .(DeltaG_asinh_UB, DeltaG_asinh_F3M, DeltaG_asinh_F3Y)], use.names = FALSE)
      all_shared_sig_display <- all_shared_sig_display[is.finite(all_shared_sig_display)]
    }
    
    heat_abs_limit <- as.numeric(quantile(abs(all_shared_sig_display), probs = TREE_HEATMAP_ABS_QUANTILE, na.rm = TRUE, names = FALSE, type = 7))
    if (!is.finite(heat_abs_limit) || heat_abs_limit <= 0) heat_abs_limit <- max(abs(all_shared_sig_display), na.rm = TRUE)
    if (!is.finite(heat_abs_limit) || heat_abs_limit <= 0) heat_abs_limit <- 1
    tip_depth <- ape::node.depth.edgelength(plot_tree)[seq_along(plot_tree$tip.label)]
    tree_depth <- max(tip_depth, na.rm = TRUE)
    if (!is.finite(tree_depth) || tree_depth <= 0) tree_depth <- 1
    p_tree <- ggtree::ggtree(
      plot_tree,
      layout = "fan",
      open.angle = 25,
      linewidth = 0.08,
      color = "#BDBDBD"
    ) %>% ggtree::rotate_tree(90)
    p_tree$data$Branch_group <- branch_phylum_dt$Branch_group[match(p_tree$data$node, branch_phylum_dt$node)]
    p_tree <- p_tree + ggtree::geom_tree(aes(color = Branch_group), linewidth = 0.08) +
      scale_color_manual(values = branch_colors, na.value = "#BDBDBD", guide = "none")
    p_tree <- ggtree::gheatmap(
      p_tree, phylum_df,
      offset = 0.015 * tree_depth, width = 0.035,
      colnames = FALSE, color = NA
    ) +
      scale_fill_manual(values = phylum_colors, name = "Phylum", drop = FALSE)
    p_tree <- p_tree + ggnewscale::new_scale_fill()
    p_tree <- ggtree::gheatmap(
      p_tree, delta_display_df,
      offset = 0.055 * tree_depth, width = 0.18,
      colnames = FALSE,
      color = NA
    ) +
      scale_fill_gradient2(
        low = "#0000A7", mid = "white", high = "#B3150E",
        midpoint = 0, limits = c(-heat_abs_limit, heat_abs_limit),
        oob = scales::squish, na.value = "#F2F2F2", name = "Scaled DeltaG\n(asinh)"
      ) +
      theme(legend.position = "right", plot.margin = margin(12, 40, 12, 12))
    ring_start_x <- tree_depth + 0.055 * tree_depth
    ring_step_x  <- (0.18 * tree_depth) / 3
    x_UB  <- ring_start_x + ring_step_x * 0.5
    x_F3M <- ring_start_x + ring_step_x * 1.5
    x_F3Y <- ring_start_x + ring_step_x * 2.5
    gap_y <- length(plot_tree$tip.label) + 2
    p_tree <- p_tree + annotate("text", x = x_UB, y = gap_y, label = "UB", size = 3.5, fontface = "bold") +
      annotate("text", x = x_F3M, y = gap_y, label = "F3M", size = 3.5, fontface = "bold") +
      annotate("text", x = x_F3Y, y = gap_y, label = "F3Y", size = 3.5, fontface = "bold")
    save_figure(
      filename = file.path(output_dir, "Fig_A_circular_phylogenetic_tree_DeltaG.pdf"),
      plot = p_tree, width = 14, height = 14, units = "in", limitsize = FALSE
    )
    
    # Record tree-display settings; absolute-response K/lambda analyses were removed.
    tree_qc <- data.table(
      Metric = c("Tree file", "Original tree tips", "Tree is rooted",
                 "Original non-positive edge lengths", "Replacement edge length",
                 "Eligible ASVs: UB", "Eligible ASVs: F3M", "Eligible ASVs: F3Y",
                 "Circular-tree tips", "Heatmap transformation", "Heatmap color quantile"),
      Value = c(TREE_FILE, as.character(length(tree_raw$tip.label)), as.character(rooted_status),
                as.character(length(nonpositive_edges)), as.character(edge_replacement),
                as.character(stage_eligible[fire_trt == "UB", uniqueN(asv_id)]),
                as.character(stage_eligible[fire_trt == "F3M", uniqueN(asv_id)]),
                as.character(stage_eligible[fire_trt == "F3Y", uniqueN(asv_id)]),
                as.character(length(plot_ids)), "asinh(DeltaG / median(abs(DeltaG)))",
                as.character(TREE_HEATMAP_ABS_QUANTILE)))
    fwrite( tree_qc, file.path(output_dir, "09_analysis_parameters_and_tree_QC.tsv"), sep = "\t" )
    message(">>> Generating Figure C: Phylum-level ASV response pie charts...")
    pie_list <- list()
    for (st in FIRE_LEVELS) {
      status_col <- paste0("Stage_status_", st)
      dt <- plot_meta[, .(asv_id, Phylum_group, Stage_status = get(status_col))]
      dt[, Status := "Not testable"]
      dt[Stage_status == "Neutral", Status := "Neutral"]
      dt[Stage_status == "Pos", Status := "Positive"]
      dt[Stage_status == "Neg", Status := "Negative"]
      dt_count <- dt[, .N, by = .(Phylum_group, Status)]
      dt_count[, Treatment := st]
      pie_list[[length(pie_list) + 1]] <- dt_count
    }
    
    pie_data <- rbindlist(pie_list)
    pie_data[, Treatment := factor(Treatment, levels = FIRE_LEVELS)]
    pie_data <- pie_data[!is.na(Phylum_group) & Phylum_group != ""]
    plot_height <- max(6, length(unique(pie_data$Phylum_group)) * 0.9)
    p_pie <- ggplot(pie_data, aes(x = "", y = N, fill = Status)) +
      geom_bar(stat = "identity", width = 1, color = "black", linewidth = 0.3) +
      coord_polar(theta = "y") +
      facet_grid(Phylum_group ~ Treatment, switch = "y") +
      scale_fill_manual(
        values = c(
          "Positive" = "#B3150E",
          "Negative" = "#0000A7",
          "Neutral" = "#E0E0E0", "Not testable" = "#F4E6A2"
        )
      ) +
      labs( title = "Number of Responding ASVs by Phylum", fill = "Response Type" ) +
      theme_void(base_size = 12) +
      theme(
        strip.text.y.left = element_text(angle = 0, hjust = 1, size = 10, face = "bold"),
        strip.text.x      = element_text(size = 12, face = "bold", margin = margin(b = 8)),
        legend.position   = "bottom",
        legend.title      = element_text(face = "bold"),
        plot.title        = element_text(hjust = 0.5, face = "bold", size = 14, margin = margin(b = 15)),
        plot.margin       = margin(20, 20, 20, 20),
        panel.border      = element_rect(color = "black", fill = NA, linewidth = 0.5)
      )
    save_figure(
      filename = file.path(output_dir, "Fig_C_Phylum_PieCharts.pdf"),
      plot = p_pie,
      width = 9,
      height = plot_height,
      units = "in"
    )
    fwrite(pie_data, file.path(output_dir, "Phylum_response_counts.tsv"), sep = "\t")
    invisible(NULL)
  }
  
  # Phylo Rgr
  run_phylo_rgr <- function() {
    K_NSIM <- K_NSIM_RGR
    TREE_FILE <- INPUTS$tree
    GROWTH_FILE <- file.path( CACHE$growth, "01_RGR_G_ASV_all.tsv" )
    GROWTH_BOOT_DIR <- file.path( CACHE$growth, "bootstrap_objects" )
    PALETTE_FILE <- file.path( CACHE$afe, "15_phylum_palette.tsv" )
    output_dir <- CACHE$phylo_rgr
    required_files <- c( TREE_FILE, GROWTH_FILE )
    missing_files <- required_files[ !file.exists(required_files) ]
    if (length(missing_files) > 0) {
      stop( "Missing required files:\n", paste(missing_files, collapse = "\n") )
    }
    
    if (!dir.exists(GROWTH_BOOT_DIR)) {
      stop("Missing bootstrap directory: ", GROWTH_BOOT_DIR)
    }
    
    RUN_PHYLO_CORRELOGRAM <- TRUE
    CORRELOGRAM_STRATIFY_BY_PHYLUM <- TRUE
    CORRELOGRAM_DISTANCE <- "patristic"
    VALIDATE_CUSTOM_K <- TRUE
    MIN_PHYLO_TAXA <- 10L
    ZERO_EDGE_MULTIPLIER <- 0.5
    CORRELOGRAM_SAMPLE_SEED <- RANDOM_SEED + 9000L
    set.seed(RANDOM_SEED)
    
    safe_boot_z <- function(x) {
      x <- x[is.finite(x)]
      if (length(x) < 2) {
        return(NA_real_)
      }
      
      s <- stats::sd(x, na.rm = TRUE)
      if (!is.finite(s) || s <= 0) {
        return(NA_real_)
      }
      
      abs(stats::median(x, na.rm = TRUE) / s)
    }
    
    prepare_tree_trait <- function(
    tree,
    trait,
    min_taxa = MIN_PHYLO_TAXA
    ) {
      trait <- trait[is.finite(trait)]
      trait <- trait[!duplicated(names(trait))]
      keep <- intersect( tree$tip.label, names(trait) )
      if (length(keep) < min_taxa) {
        return(NULL)
      }
      
      tr <- ape::drop.tip( tree, setdiff(tree$tip.label, keep) )
      x <- trait[tr$tip.label]
      if (
        length(x) < min_taxa || !all(is.finite(x)) || !is.finite(stats::sd(x)) || stats::sd(x) <= 0
      ) {
        return(NULL)
      }
      
      list( tree = tr, trait = x )
    }
    
    build_k_cache <- function(tree) {
      if (is.null(tree$edge.length)) {
        stop("Tree must contain branch lengths for Blomberg K.")
      }
      
      tr <- ape::reorder.phylo(tree, order = "cladewise")
      edge <- tr$edge
      edge_len <- tr$edge.length
      n_tip <- length(tr$tip.label)
      n_all <- n_tip + tr$Nnode
      parent_nodes <- unique(edge[, 1])
      child_nodes <- unique(edge[, 2])
      root_node <- setdiff(parent_nodes, child_nodes)
      if (length(root_node) != 1L) {
        stop("Could not identify exactly one root node.")
      }
      
      node_depth <- ape::node.depth.edgelength(tr)
      internal_order <- parent_nodes[order(node_depth[parent_nodes], decreasing = TRUE)]
      child_list <- vector("list", length(internal_order))
      branch_list <- vector("list", length(internal_order))
      node_var <- numeric(n_all)
      for (ii in seq_along(internal_order)) {
        parent <- internal_order[ii]
        rows <- which(edge[, 1] == parent)
        ch <- edge[rows, 2]
        bl <- edge_len[rows]
        child_list[[ii]] <- ch
        branch_list[[ii]] <- bl
        eff_var <- node_var[ch] + bl
        if (length(eff_var) == 1L) {
          node_var[parent] <- eff_var[1]
        } else {
          v <- eff_var[1]
          for (jj in 2:length(eff_var)) {
            v2 <- eff_var[jj]
            v <- v * v2 / (v + v2)
          }
          node_var[parent] <- v
        }
      }
      
      root_var <- node_var[root_node]
      trace_C <- sum(node_depth[seq_len(n_tip)])
      expected_ratio <- (trace_C - n_tip * root_var) / (n_tip - 1)
      if (!is.finite(root_var) || root_var <= 0) {
        stop("Invalid root GLS variance in K cache.")
      }
      
      if (!is.finite(expected_ratio) || expected_ratio <= 0) {
        stop("Invalid Brownian expectation term in K cache.")
      }
      
      list(
        tree = tr,
        n_tip = n_tip,
        n_all = n_all,
        root_node = root_node,
        internal_order = internal_order,
        child_list = child_list,
        branch_list = branch_list,
        node_var = node_var,
        expected_ratio = expected_ratio
      )
    }
    
    k_value_ordered <- function(x_ordered, cache) {
      if (length(x_ordered) != cache$n_tip) {
        stop("Trait length does not match tree tips.")
      }
      
      if (any(!is.finite(x_ordered))) {
        stop("Non-finite trait values passed to K.")
      }
      
      mu <- numeric(cache$n_all)
      mu[seq_len(cache$n_tip)] <- x_ordered
      sse_phylo <- 0
      for (ii in seq_along(cache$internal_order)) {
        parent <- cache$internal_order[ii]
        ch <- cache$child_list[[ii]]
        bl <- cache$branch_list[[ii]]
        m <- mu[ch[1]]
        v <- cache$node_var[ch[1]] + bl[1]
        if (length(ch) > 1L) {
          for (jj in 2:length(ch)) {
            m2 <- mu[ch[jj]]
            v2 <- cache$node_var[ch[jj]] + bl[jj]
            denom <- v + v2
            sse_phylo <- sse_phylo + (m - m2)^2 / denom
            m <- (m * v2 + m2 * v) / denom
            v <- v * v2 / denom
          }
        }
        
        mu[parent] <- m
      }
      
      a_hat <- mu[cache$root_node]
      sse_raw <- sum((x_ordered - a_hat)^2)
      if (!is.finite(sse_phylo) || sse_phylo <= 0) {
        return(NA_real_)
      }
      
      as.numeric((sse_raw / sse_phylo) / cache$expected_ratio)
    }
    
    custom_blomberg_k <- function(tree, trait, nsim = K_NSIM, seed = RANDOM_SEED) {
      cache <- build_k_cache(tree)
      x_ordered <- as.numeric(trait[cache$tree$tip.label])
      if (any(!is.finite(x_ordered))) {
        stop("Trait/tree alignment produced non-finite values in K.")
      }
      
      K_obs <- k_value_ordered(x_ordered, cache)
      if (!is.finite(K_obs)) {
        return(list(K = NA_real_, P = NA_real_))
      }
      
      if (nsim <= 1L) {
        return(list(K = K_obs, P = NA_real_))
      }
      
      set.seed(seed)
      sim_K <- numeric(nsim)
      sim_K[1] <- K_obs
      if (nsim > 1L) {
        for (bb in 2:nsim) {
          sim_K[bb] <- k_value_ordered(sample(x_ordered, replace = FALSE), cache)
        }
      }
      
      P_raw <- mean(sim_K >= K_obs, na.rm = TRUE)
      list( K = K_obs, P = P_raw )
    }
    
    validate_custom_k <- function(tree, trait, seed = RANDOM_SEED, n_validate = 80L) {
      if (!VALIDATE_CUSTOM_K) {
        return(invisible(TRUE))
      }
      
      ids <- intersect(tree$tip.label, names(trait))
      ids <- ids[is.finite(trait[ids])]
      if (length(ids) < 10L) {
        return(invisible(TRUE))
      }
      
      set.seed(seed)
      ids <- sample(ids, size = min(n_validate, length(ids)), replace = FALSE)
      tr_small <- ape::drop.tip(tree, setdiff(tree$tip.label, ids))
      x_small <- trait[tr_small$tip.label]
      K_custom <- custom_blomberg_k( tree = tr_small, trait = x_small, nsim = 1L, seed = seed )$K
      K_ref <- suppressWarnings(
        as.numeric( phytools::phylosig( tree = tr_small, x = x_small, method = "K", test = FALSE ) )
      )
      rel_err <- abs(K_custom - K_ref) / max(abs(K_ref), 1e-12)
      if (!is.finite(rel_err) || rel_err > 1e-5) {
        stop(
          "Custom-K validation failed: K_custom = ", K_custom,
          ", K_phytools = ", K_ref,
          ", relative error = ", rel_err
        )
      }
      
      message(
        ">>> Custom-K validation passed on ", length(ids),
        " tips; relative error = ", format(rel_err, scientific = TRUE)
      )
      invisible(TRUE)
    }
    
    custom_pagel_lambda <- function(tree, trait) {
      x_ordered <- as.numeric(trait[tree$tip.label])
      if (any(!is.finite(x_ordered))) {
        stop("Trait/tree alignment produced non-finite values in lambda fit.")
      }
      
      dat_lambda <- data.frame(trait_value = x_ordered)
      rownames(dat_lambda) <- tree$tip.label
      fit_lambda <- suppressWarnings(
        phylolm::phylolm(
          trait_value ~ 1,
          data = dat_lambda,
          phy = tree,
          model = "lambda",
          lower.bound = 1e-8,
          upper.bound = 1,
          starting.value = 0.5,
          boot = 0,
          full.matrix = FALSE
        )
      )
      fit_null <- suppressWarnings(
        phylolm::phylolm(
          trait_value ~ 1,
          data = dat_lambda,
          phy = tree,
          model = "lambda",
          lower.bound = 1e-10,
          upper.bound = 1e-8,
          starting.value = 1e-9,
          boot = 0,
          full.matrix = FALSE
        )
      )
      lambda_hat <- suppressWarnings(as.numeric(fit_lambda$optpar[1]))
      logL1 <- suppressWarnings(as.numeric(fit_lambda$logLik[1]))
      logL0 <- suppressWarnings(as.numeric(fit_null$logLik[1]))
      if (!is.finite(lambda_hat)) {
        lambda_hat <- NA_real_
      }
      
      if (!is.finite(logL1) || !is.finite(logL0)) {
        return( list( lambda = lambda_hat, P = NA_real_, logL = logL1, logL0 = logL0 ) )
      }
      
      LR <- max(0, 2 * (logL1 - logL0))
      P_raw <- stats::pchisq(LR, df = 1, lower.tail = FALSE)
      list( lambda = lambda_hat, P = P_raw, logL = logL1, logL0 = logL0 )
    }
    
    calc_phylo_signal <- function(
    tree,
    trait,
    stage,
    stage_label,
    trait_version,
    analysis_set,
    k_nsim = K_NSIM
    ) {
      prep <- prepare_tree_trait( tree = tree, trait = trait, min_taxa = MIN_PHYLO_TAXA )
      if (is.null(prep)) {
        return(
          data.table(
            analysis_set = analysis_set,
            fire_trt = stage,
            stage_label = stage_label,
            trait = "DeltaRGR",
            trait_version = trait_version,
            N_ASV = sum(is.finite(trait)),
            Blomberg_K = NA_real_,
            K_P_raw = NA_real_,
            Pagel_lambda = NA_real_,
            lambda_P_raw = NA_real_,
            K_nsim = k_nsim
          )
        )
      }
      
      tr <- prep$tree
      x <- prep$trait
      message(
        ">>> Phylogenetic signal: ",
        analysis_set,
        " | ", stage,
        " | ", "DeltaRGR",
        " | ", trait_version,
        " | N = ", length(x),
        " | K permutations = ", k_nsim
      )
      validate_custom_k( tree = tr, trait = x, seed = RANDOM_SEED + match(stage, FIRE_LEVELS) )
      kres <- tryCatch(
        custom_blomberg_k(
          tree = tr,
          trait = x,
          nsim = k_nsim,
          seed = RANDOM_SEED + 100L * match(stage, FIRE_LEVELS)
        ),
        error = function(e) {
          warning(
            "Blomberg K failed for ",
            analysis_set, " / ", stage, " / ", "DeltaRGR", " / ", trait_version,
            ": ", conditionMessage(e)
          )
          NULL
        }
      )
      lres <- tryCatch(
        custom_pagel_lambda( tree = tr, trait = x ),
        error = function(e) {
          warning(
            "Pagel lambda failed for ",
            analysis_set, " / ", stage, " / ", "DeltaRGR", " / ", trait_version,
            ": ", conditionMessage(e)
          )
          NULL
        }
      )
      K_value <- if (!is.null(kres)) kres$K else NA_real_
      K_p <- if (!is.null(kres)) kres$P else NA_real_
      lambda_value <- if (!is.null(lres)) lres$lambda else NA_real_
      lambda_p <- if (!is.null(lres)) lres$P else NA_real_
      data.table(
        analysis_set = analysis_set,
        fire_trt = stage,
        stage_label = stage_label,
        trait = "DeltaRGR",
        trait_version = trait_version,
        N_ASV = length(x),
        Blomberg_K = K_value,
        K_P_raw = K_p,
        Pagel_lambda = lambda_value,
        lambda_P_raw = lambda_p,
        K_nsim = k_nsim
      )
    }
    
    format_signal_value <- function(x) {
      if (!is.finite(x)) {
        return("NA")
      }
      
      ax <- abs(x)
      if (ax == 0) {
        return("0")
      }
      
      if (ax < 0.001) {
        return(formatC(x, format = "e", digits = 2))
      }
      
      if (ax < 0.01) {
        return(formatC(x, format = "f", digits = 4))
      }
      
      formatC(x, format = "f", digits = 2)
    }
    
    select_correlogram_subset <- function(
    dt,
    max_n = CORRELOGRAM_MAX_ASVS,
    stratify = CORRELOGRAM_STRATIFY_BY_PHYLUM,
    seed = CORRELOGRAM_SAMPLE_SEED
    ) {
      x <- copy(dt)
      x <- unique(x, by = "asv_id")
      if (!is.finite(max_n) || nrow(x) <= max_n) {
        return(x)
      }
      
      max_n <- as.integer(max_n)
      if (max_n < MIN_PHYLO_TAXA) {
        stop("CORRELOGRAM_MAX_ASVS must be >= MIN_PHYLO_TAXA.")
      }
      
      set.seed(seed)
      if (!stratify || !"Phylum_group" %in% names(x)) {
        keep_ids <- sample(x$asv_id, size = max_n, replace = FALSE)
        return(x[asv_id %chin% keep_ids])
      }
      
      x[, cor_stratum := as.character(Phylum_group)]
      x[is.na(cor_stratum) | cor_stratum == "", cor_stratum := "Unclassified"]
      strata <- sort(unique(x$cor_stratum))
      if (length(strata) <= max_n) {
        mandatory_ids <- unlist(
          lapply(strata, function(gr) {
            ids <- x[cor_stratum == gr, asv_id]
            sample(ids, size = 1L)
          }),
          use.names = FALSE
        )
      } else {
        chosen_strata <- sample(strata, size = max_n, replace = FALSE)
        mandatory_ids <- unlist(
          lapply(chosen_strata, function(gr) {
            ids <- x[cor_stratum == gr, asv_id]
            sample(ids, size = 1L)
          }),
          use.names = FALSE
        )
      }
      
      remaining_n <- max_n - length(mandatory_ids)
      remaining_ids <- setdiff(x$asv_id, mandatory_ids)
      if (remaining_n > 0L) {
        extra_ids <- sample( remaining_ids, size = min(remaining_n, length(remaining_ids)), replace = FALSE )
        keep_ids <- c(mandatory_ids, extra_ids)
      } else {
        keep_ids <- mandatory_ids
      }
      
      out <- x[asv_id %chin% keep_ids]
      out[, cor_stratum := NULL]
      setorder(out, asv_id)
      out
    }
    
    calc_phylo_correlogram <- function(
    tree,
    trait,
    stage,
    stage_label,
    n_points = CORRELOGRAM_N_POINTS,
    ci_bs = CORRELOGRAM_CI_BS,
    ci_conf = CORRELOGRAM_CI_CONF,
    dist_method = CORRELOGRAM_DISTANCE
    ) {
      prep <- prepare_tree_trait( tree = tree, trait = trait, min_taxa = MIN_PHYLO_TAXA )
      if (is.null(prep)) {
        warning("Skipping phylogenetic correlogram for ", stage, ": insufficient valid taxa.")
        return(NULL)
      }
      
      tr <- prep$tree
      x <- prep$trait
      trait_df <- data.frame( DeltaRGR = as.numeric(x), row.names = names(x), check.names = FALSE )
      p4d <- phylobase::phylo4d( tr, tip.data = trait_df, missing.data = "fail", extra.data = "fail" )
      message(
        ">>> Phylogenetic correlogram: ", stage,
        " | N = ", length(x),
        " | distance = ", dist_method,
        " | points = ", n_points,
        " | CI bootstrap = ", ci_bs,
        " | CI = ", 100 * ci_conf, "%"
      )
      pc <- phylosignal::phyloCorrelogram(
        p4d = p4d,
        trait = "DeltaRGR",
        dist.phylo = dist_method,
        sigma = NULL,
        n.points = n_points,
        ci.bs = ci_bs,
        ci.conf = ci_conf
      )
      res <- as.data.table(pc$res)
      if (ncol(res) != 4L) {
        stop("Unexpected phylosignal::phyloCorrelogram result structure for ", stage, ".")
      }
      
      setnames( res, c( "phylogenetic_distance", "CI_low", "CI_high", "Moran_I" ) )
      h0 <- -1 / (pc$n - 1)
      res[
        ,
        `:=`(
          fire_trt = stage,
          stage_label = stage_label,
          N_ASV = as.integer(pc$n),
          H0 = h0,
          sigma = as.numeric(pc$sigma),
          ci_boot = as.integer(pc$ci.bs),
          ci_conf = as.numeric(pc$ci.conf)
        )
      ]
      res[
        ,
        correlogram_class := fifelse(
          is.finite(CI_low) & is.finite(CI_high) & CI_low > H0 & CI_high > H0,
          "Positive autocorrelation",
          fifelse(
            is.finite(CI_low) & is.finite(CI_high) & CI_low < H0 & CI_high < H0,
            "Negative autocorrelation",
            "Not significant"
          )
        )
      ]
      list( object = pc, table = res, tree = tr, trait = x )
    }
    
    message(">>> Reading rooted phylogenetic tree...")
    tree_raw <- ape::read.tree(TREE_FILE)
    if (
      is.null(tree_raw$tip.label) || length(tree_raw$tip.label) == 0
    ) {
      stop("No tip labels found in rooted-tree.nwk.")
    }
    
    if (anyDuplicated(tree_raw$tip.label)) {
      dup <- unique( tree_raw$tip.label[ duplicated(tree_raw$tip.label) ] )
      fwrite( data.table(asv_id = dup), file.path( output_dir, "QC_duplicate_tree_tip_labels.tsv" ), sep = "\t" )
      stop("Duplicate tip labels detected in rooted-tree.nwk.")
    }
    
    if (is.null(tree_raw$edge.length)) {
      stop( "The supplied tree has no branch lengths; K/lambda require branch lengths." )
    }
    
    if (any(!is.finite(tree_raw$edge.length))) {
      stop("Non-finite branch lengths detected in rooted-tree.nwk.")
    }
    
    if (!ape::is.rooted(tree_raw)) {
      stop(
        paste0(
          "ape::is.rooted() reports that rooted-tree.nwk is not rooted. ",
          "Verify the supplied tree rather than automatically rooting it here."
        )
      )
    }
    
    tree <- tree_raw
    nonpositive_edges <- which(tree$edge.length <= 0)
    edge_replacement <- NA_real_
    if (length(nonpositive_edges) > 0) {
      positive_edges <- tree$edge.length[ is.finite(tree$edge.length) & tree$edge.length > 0 ]
      if (length(positive_edges) == 0) {
        stop("Tree contains no positive branch lengths.")
      }
      
      edge_replacement <- min(positive_edges) * ZERO_EDGE_MULTIPLIER
      if (!is.finite(edge_replacement) || edge_replacement <= 0) {
        edge_replacement <- 1e-10
      }
      
      tree$edge.length[nonpositive_edges] <- edge_replacement
      warning(
        length(nonpositive_edges),
        " non-positive branch lengths were replaced by ",
        format(edge_replacement, scientific = TRUE),
        "."
      )
    }
    
    message(">>> Reading Step-04 growth table...")
    growth <- fread( GROWTH_FILE, sep = "\t", header = TRUE )
    required_growth_cols <- c(
      "asv_id",
      "fire_trt",
      "tem",
      "RGR_ASV",
      "Active_ASV",
      "N_TOTAL_unlabel_mean",
      "Phylum",
      "Phylum_group"
    )
    missing_growth_cols <- setdiff( required_growth_cols, names(growth) )
    if (length(missing_growth_cols) > 0) {
      stop( "04 growth table missing columns: ", paste(missing_growth_cols, collapse = ", ") )
    }
    
    growth[, asv_id := as.character(asv_id)]
    growth[, fire_trt := as.character(fire_trt)]
    growth[, tem := as.character(tem)]
    growth[, RGR_ASV := safe_numeric(RGR_ASV)]
    growth[, Active_ASV := as.logical(Active_ASV)]
    growth[, N_TOTAL_unlabel_mean := safe_numeric(N_TOTAL_unlabel_mean)]
    growth[, Phylum := as.character(Phylum)]
    growth[, Phylum_group := as.character(Phylum_group)]
    if (any(!growth$fire_trt %chin% FIRE_LEVELS)) {
      stop("Unexpected fire_trt values in growth table.")
    }
    
    if (any(!growth$tem %chin% TEMP_LEVELS)) {
      stop("Unexpected tem values in growth table.")
    }
    
    key_count <- growth[ , .N, by = .( asv_id, fire_trt, tem ) ]
    if (any(key_count$N > 1)) {
      fwrite( key_count[N > 1], file.path( output_dir, "QC_duplicate_ASV_fire_temp_rows.tsv" ), sep = "\t" )
      stop( "Duplicate asv_id x fire_trt x tem rows detected in Step-04 growth table." )
    }
    
    message(">>> Reconstructing DeltaRGR bootstrap distributions...")
    
    # Classify temperature responses only for comparable ASVs.
    stage_results_list <- list()
    for (i in seq_along(FIRE_LEVELS)) {
      current_fire <- FIRE_LEVELS[[i]]
      message( ">>> DeltaRGR stage: ", current_fire, " (", FIRE_LABELS[[current_fire]], ")" )
      dt21 <- growth[
        fire_trt == current_fire & tem == "T_low",
        .(
          asv_id,
          RGR_T_low = RGR_ASV,
          Active_T_low = Active_ASV,
          N_TOTAL_T_low = N_TOTAL_unlabel_mean,
          Phylum_T_low = Phylum,
          Phylum_group_T_low = Phylum_group
        )
      ]
      dt28 <- growth[
        fire_trt == current_fire & tem == "T_high",
        .(
          asv_id,
          RGR_T_high = RGR_ASV,
          Active_T_high = Active_ASV,
          N_TOTAL_T_high = N_TOTAL_unlabel_mean,
          Phylum_T_high = Phylum,
          Phylum_group_T_high = Phylum_group
        )
      ]
      merged_stage <- merge( dt21, dt28, by = "asv_id", all = TRUE, sort = FALSE )
      merged_stage[, fire_trt := current_fire]
      merged_stage[, stage_label := FIRE_LABELS[[current_fire]]]
      merged_stage[, testable_T_low := is.finite(RGR_T_low)]
      merged_stage[, testable_T_high := is.finite(RGR_T_high)]
      merged_stage[, testable_both := testable_T_low & testable_T_high]
      merged_stage[ is.na(Active_T_low), Active_T_low := FALSE ]
      merged_stage[ is.na(Active_T_high), Active_T_high := FALSE ]
      merged_stage[, active_any := Active_T_low | Active_T_high]
      merged_stage[ , Phylum := fifelse( !is.na(Phylum_T_low) & Phylum_T_low != "", Phylum_T_low, Phylum_T_high ) ]
      merged_stage[
        ,
        Phylum_group := fifelse(
          !is.na(Phylum_group_T_low) & Phylum_group_T_low != "",
          Phylum_group_T_low,
          Phylum_group_T_high
        )
      ]
      merged_stage[ , mean_endpoint_abundance := rowMeans( cbind( N_TOTAL_T_low, N_TOTAL_T_high ), na.rm = TRUE ) ]
      merged_stage[ !is.finite(mean_endpoint_abundance), mean_endpoint_abundance := NA_real_ ]
      merged_stage[ testable_both == TRUE, DeltaRGR_point := RGR_T_high - RGR_T_low ]
      boot_file_21 <- file.path( GROWTH_BOOT_DIR, paste0( "G_bootstrap_", current_fire, "_T_low.rds" ) )
      boot_file_28 <- file.path( GROWTH_BOOT_DIR, paste0( "G_bootstrap_", current_fire, "_T_high.rds" ) )
      if (
        !file.exists(boot_file_21) || !file.exists(boot_file_28)
      ) {
        stop(
          "Missing Step-04 bootstrap objects for ",
          current_fire,
          ". Expected:\n",
          boot_file_21,
          "\n",
          boot_file_28
        )
      }
      
      boot21 <- readRDS(boot_file_21)
      boot28 <- readRDS(boot_file_28)
      required_boot_items <- c( "asv_id", "RGR_boot" )
      if (
        !all(required_boot_items %in% names(boot21)) || !all(required_boot_items %in% names(boot28))
      ) {
        stop( "Step-04 bootstrap objects must contain asv_id and RGR_boot." )
      }
      
      common_asv <- intersect( as.character(boot21$asv_id), as.character(boot28$asv_id) )
      common_asv <- sort(unique(common_asv))
      if (length(common_asv) == 0) {
        stop( "No overlapping ASVs between T_low and T_high RGR bootstrap objects for ", current_fire, "." )
      }
      
      idx21 <- match( common_asv, as.character(boot21$asv_id) )
      idx28 <- match( common_asv, as.character(boot28$asv_id) )
      RGR21_boot <- boot21$RGR_boot[ idx21, , drop = FALSE ]
      RGR28_boot <- boot28$RGR_boot[ idx28, , drop = FALSE ]
      storage.mode(RGR21_boot) <- "double"
      storage.mode(RGR28_boot) <- "double"
      set.seed(RANDOM_SEED + i)
      delta_boot <- combine_independent_bootstrap_difference(
        mat_hot = RGR28_boot,
        mat_cold = RGR21_boot,
        n_boot = N_BOOT_COMBINE
      )
      delta_q <- t( apply( delta_boot, 1, safe_quantiles, probs = CI_PROBS ) )
      colnames(delta_q) <- c( "DeltaRGR_CI_low_90", "DeltaRGR", "DeltaRGR_CI_high_90" )
      delta_dt <- data.table(
        asv_id = common_asv,
        DeltaRGR_CI_low_90 = delta_q[, 1],
        DeltaRGR = delta_q[, 2],
        DeltaRGR_CI_high_90 = delta_q[, 3],
        DeltaRGR_boot_z = apply( delta_boot, 1, safe_boot_z ),
        prop_boot_DeltaRGR_gt_0 = rowMeans( delta_boot > 0, na.rm = TRUE ),
        prop_boot_DeltaRGR_lt_0 = rowMeans( delta_boot < 0, na.rm = TRUE ),
        DeltaRGR_boot_sd = apply( delta_boot, 1, stats::sd, na.rm = TRUE )
      )
      delta_dt[ , DeltaRGR_strategy := classify_delta_strategy( DeltaRGR_CI_low_90, DeltaRGR_CI_high_90 ) ]
      merged_stage <- merge( merged_stage, delta_dt, by = "asv_id", all.x = TRUE, sort = FALSE )
      merged_stage[ testable_both == FALSE, Stage_status := "Not testable" ]
      merged_stage[ testable_both == TRUE & active_any == FALSE, Stage_status := "Inactive" ]
      merged_stage[ testable_both == TRUE & active_any == TRUE & DeltaRGR_strategy == "Pos", Stage_status := "Pos" ]
      merged_stage[ testable_both == TRUE & active_any == TRUE & DeltaRGR_strategy == "Neg", Stage_status := "Neg" ]
      merged_stage[
        testable_both == TRUE & active_any == TRUE & DeltaRGR_strategy == "Neutral",
        Stage_status := "Neutral"
      ]
      stage_results_list[[length(stage_results_list) + 1L]] <- merged_stage
    }
    
    stage_dt <- rbindlist( stage_results_list, use.names = TRUE, fill = TRUE )
    stage_dt[ , fire_trt := factor( fire_trt, levels = FIRE_LEVELS ) ]
    setorder( stage_dt, fire_trt, asv_id )
    fwrite( stage_dt, file.path( output_dir, "01_DeltaRGR_strategy_all.tsv" ), sep = "\t" )
    stage_dt[, fire_trt_chr := as.character(fire_trt)]
    stage_dt[, in_tree := asv_id %chin% tree$tip.label]
    stage_dt[
      ,
      phylo_eligible := testable_both == TRUE & active_any == TRUE & is.finite(DeltaRGR) & in_tree == TRUE
    ]
    stage_eligible <- stage_dt[ phylo_eligible == TRUE ]
    eligible_summary <- stage_dt[
      ,
      .(
        total_ASVs_in_stage_table = .N,
        testable_both = sum(testable_both, na.rm = TRUE),
        active_any = sum(testable_both & active_any, na.rm = TRUE),
        finite_DeltaRGR = sum( testable_both & active_any & is.finite(DeltaRGR), na.rm = TRUE ),
        present_in_tree = sum( testable_both & active_any & is.finite(DeltaRGR) & in_tree, na.rm = TRUE ),
        phylo_eligible_ASVs = sum( phylo_eligible, na.rm = TRUE ),
        Pos = sum(Stage_status == "Pos", na.rm = TRUE),
        Neg = sum(Stage_status == "Neg", na.rm = TRUE),
        Neutral = sum(Stage_status == "Neutral", na.rm = TRUE)
      ),
      by = fire_trt_chr
    ]
    eligible_summary[ , stage_label := FIRE_LABELS[fire_trt_chr] ]
    setnames( eligible_summary, "fire_trt_chr", "fire_trt" )
    fwrite( eligible_summary, file.path( output_dir, "02_tree_matching_and_eligibility_QC.tsv" ), sep = "\t" )
    fwrite( stage_eligible, file.path( output_dir, "03_stage_specific_phylo_eligible_ASVs.tsv" ), sep = "\t" )
    shared_ids <- Reduce(
      intersect,
      lapply(
        FIRE_LEVELS,
        function(st) {
          unique( stage_eligible[ fire_trt_chr == st, asv_id ] )
        }
      )
    )
    shared_ids <- sort(unique(shared_ids))
    shared_tree <- NULL
    shared_traits_wide <- data.table()
    global_DeltaRGR_scale <- safe_asinh_scale(stage_eligible$DeltaRGR)
    if (length(shared_ids) >= MIN_PHYLO_TAXA) {
      shared_stage <- stage_eligible[ asv_id %chin% shared_ids ]
      shared_counts <- shared_stage[ , .N, by = fire_trt_chr ]
      if (
        nrow(shared_counts) != 3L || any(shared_counts$N != length(shared_ids))
      ) {
        stop( "Shared-ASV construction failed: stage counts are inconsistent." )
      }
      
      fwrite(
        data.table(asv_id = shared_ids),
        file.path( output_dir, "04_shared_DeltaRGR_ASV_list.tsv" ),
        sep = "\t"
      )
      shared_traits_long <- shared_stage[
        ,
        .(
          asv_id,
          fire_trt = fire_trt_chr,
          DeltaRGR,
          DeltaRGR_point,
          DeltaRGR_CI_low_90,
          DeltaRGR_CI_high_90,
          Stage_status,
          mean_endpoint_abundance,
          Phylum,
          Phylum_group
        )
      ]
      taxonomy_dt <- unique( shared_traits_long[ , .( asv_id, Phylum, Phylum_group ) ], by = "asv_id" )
      shared_traits_wide <- dcast(
        shared_traits_long,
        asv_id ~ fire_trt,
        value.var = c(
          "DeltaRGR",
          "DeltaRGR_point",
          "DeltaRGR_CI_low_90",
          "DeltaRGR_CI_high_90",
          "Stage_status",
          "mean_endpoint_abundance"
        )
      )
      shared_traits_wide <- merge( shared_traits_wide, taxonomy_dt, by = "asv_id", all.x = TRUE, sort = FALSE )
      expected_shared_cols <- c( "DeltaRGR_UB", "DeltaRGR_F3M", "DeltaRGR_F3Y" )
      missing_shared_cols <- setdiff( expected_shared_cols, names(shared_traits_wide) )
      if (length(missing_shared_cols) > 0) {
        stop( "Shared DeltaRGR table missing columns: ", paste(missing_shared_cols, collapse = ", ") )
      }
      
      global_DeltaRGR_scale <- safe_asinh_scale(
        unlist( shared_traits_wide[ , .( DeltaRGR_UB, DeltaRGR_F3M, DeltaRGR_F3Y ) ], use.names = FALSE )
      )
      for (st in FIRE_LEVELS) {
        raw_col <- paste0("DeltaRGR_", st)
        trans_col <- paste0("DeltaRGR_asinh_", st)
        shared_traits_wide[ , (trans_col) := asinh( get(raw_col) / global_DeltaRGR_scale ) ]
      }
      
      fwrite( shared_traits_wide, file.path( output_dir, "05_shared_DeltaRGR_trait_table.tsv" ), sep = "\t" )
      shared_tree <- ape::drop.tip( tree, setdiff( tree$tip.label, shared_ids ) )
      ape::write.tree( shared_tree, file = file.path( output_dir, "shared_DeltaRGR_tree.nwk" ) )
      message(">>> Main phylogenetic-signal analysis: ALL stage-specific active/comparable ASVs...")
    } else {
      warning("Shared-ASV sensitivity unavailable: fewer than MIN_PHYLO_TAXA shared taxa.")
      RUN_SHARED_SENSITIVITY <- FALSE
    }
    # Estimate signal from each stage-specific eligible ASV set.
    main_signal_list <- list()
    for (st in FIRE_LEVELS) {
      current <- stage_eligible[fire_trt_chr == st]
      if (nrow(current) < MIN_PHYLO_TAXA) {
        stop(
          "Only ", nrow(current),
          " phylogenetically eligible active ASVs at ", st,
          "; fewer than MIN_PHYLO_TAXA = ", MIN_PHYLO_TAXA, "."
        )
      }
      
      message( ">>> Main DeltaRGR phylogeny: ", st, " | all eligible active ASVs = ", nrow(current) )
      stage_tree <- ape::drop.tip( tree, setdiff(tree$tip.label, current$asv_id) )
      ape::write.tree(
        stage_tree,
        file = file.path(output_dir, paste0("stage_specific_all_active_tree_", st, ".nwk"))
      )
      x_raw <- current$DeltaRGR
      names(x_raw) <- current$asv_id
      main_signal_list[[length(main_signal_list) + 1L]] <-
        calc_phylo_signal(
          tree = stage_tree,
          trait = x_raw,
          stage = st,
          stage_label = FIRE_LABELS[[st]],
          trait_version = "raw",
          analysis_set = "stage_specific_all_active"
        )
      if (RUN_ASINH_SENSITIVITY) {
        stage_scale <- safe_asinh_scale(x_raw)
        x_asinh <- asinh(x_raw / stage_scale)
        names(x_asinh) <- current$asv_id
        main_signal_list[[length(main_signal_list) + 1L]] <-
          calc_phylo_signal(
            tree = stage_tree,
            trait = x_asinh,
            stage = st,
            stage_label = FIRE_LABELS[[st]],
            trait_version = "asinh_sensitivity",
            analysis_set = "stage_specific_all_active"
          )
      }
    }
    
    main_signal <- rbindlist(main_signal_list, use.names = TRUE, fill = TRUE)
    main_signal[, fire_trt := factor(fire_trt, levels = FIRE_LEVELS)]
    setorder(main_signal, trait_version, fire_trt)
    fwrite( main_signal, file.path(output_dir, "06_phylogenetic_signal_all_active_ASVs.tsv"), sep = "\t" )
    
    # Estimate response autocorrelation across phylogenetic distances.
    correlogram_table <- data.table()
    correlogram_sampling_qc_list <- list()
    correlogram_sample_ids_list <- list()
    correlogram_object_dir <- file.path(output_dir, "correlogram_objects")
    if (RUN_PHYLO_CORRELOGRAM) {
      if (!dir.exists(correlogram_object_dir)) {
        dir.create(correlogram_object_dir, recursive = TRUE)
      }
      
      correlogram_list <- list()
      for (ii in seq_along(FIRE_LEVELS)) {
        st <- FIRE_LEVELS[[ii]]
        current_full <- stage_eligible[fire_trt_chr == st]
        current <- select_correlogram_subset(
          dt = current_full,
          max_n = CORRELOGRAM_MAX_ASVS,
          stratify = CORRELOGRAM_STRATIFY_BY_PHYLUM,
          seed = CORRELOGRAM_SAMPLE_SEED + ii
        )
        correlogram_sampling_qc_list[[length(correlogram_sampling_qc_list) + 1L]] <- data.table(
          fire_trt = st,
          stage_label = FIRE_LABELS[[st]],
          N_all_Klambda = nrow(current_full),
          N_correlogram = nrow(current),
          sample_fraction = nrow(current) / nrow(current_full),
          max_ASVs_setting = as.character(CORRELOGRAM_MAX_ASVS),
          stratified_by_phylum = CORRELOGRAM_STRATIFY_BY_PHYLUM,
          n_phyla_all = uniqueN(current_full$Phylum_group),
          n_phyla_correlogram = uniqueN(current$Phylum_group)
        )
        correlogram_sample_ids_list[[length(correlogram_sample_ids_list) + 1L]] <- current[
          ,
          .( fire_trt = st, stage_label = FIRE_LABELS[[st]], asv_id, Phylum, Phylum_group, DeltaRGR )
        ]
        message(
          ">>> Correlogram subset: ", st,
          " | all K/lambda N = ", nrow(current_full),
          " | correlogram N = ", nrow(current)
        )
        stage_tree <- ape::drop.tip( tree, setdiff(tree$tip.label, current$asv_id) )
        x_cor <- current$DeltaRGR
        names(x_cor) <- current$asv_id
        set.seed(RANDOM_SEED + 7000L + ii)
        cor_res <- tryCatch(
          calc_phylo_correlogram( tree = stage_tree, trait = x_cor, stage = st, stage_label = FIRE_LABELS[[st]] ),
          error = function(e) {
            warning( "Phylogenetic correlogram failed for ", st, ": ", conditionMessage(e) )
            NULL
          }
        )
        if (!is.null(cor_res)) {
          cor_res$table[, N_ASV_full_Klambda := nrow(current_full)]
          cor_res$table[, sample_fraction := nrow(current) / nrow(current_full)]
          correlogram_list[[length(correlogram_list) + 1L]] <- cor_res$table
          saveRDS(
            cor_res$object,
            file = file.path( correlogram_object_dir, paste0("phylogenetic_correlogram_DeltaRGR_", st, ".rds") )
          )
        }
      }
      
      if (length(correlogram_sampling_qc_list) > 0L) {
        correlogram_sampling_qc <- rbindlist( correlogram_sampling_qc_list, use.names = TRUE, fill = TRUE )
        fwrite(
          correlogram_sampling_qc,
          file.path(output_dir, "07a_phylogenetic_correlogram_sampling_QC.tsv"),
          sep = "\t"
        )
      }
      
      if (length(correlogram_sample_ids_list) > 0L) {
        correlogram_sample_ids <- rbindlist( correlogram_sample_ids_list, use.names = TRUE, fill = TRUE )
        fwrite(
          correlogram_sample_ids,
          file.path(output_dir, "07a2_phylogenetic_correlogram_sampled_ASVs.tsv"),
          sep = "\t"
        )
      }
      
      if (length(correlogram_list) > 0L) {
        correlogram_table <- rbindlist( correlogram_list, use.names = TRUE, fill = TRUE )
        correlogram_table[ , stage_label := factor( stage_label, levels = FIRE_LABELS[FIRE_LEVELS] ) ]
        setorder(correlogram_table, fire_trt, phylogenetic_distance)
        fwrite( correlogram_table, file.path(output_dir, "07b_phylogenetic_correlogram_DeltaRGR.tsv"), sep = "\t" )
      }
    } else {
      message(">>> Phylogenetic correlogram skipped: RUN_PHYLO_CORRELOGRAM = FALSE")
    }
    
    if (RUN_SHARED_SENSITIVITY) {
      message(">>> Running optional shared-ASV DeltaRGR sensitivity...")
      shared_signal_list <- list()
      for (st in FIRE_LEVELS) {
        x_raw <- shared_traits_wide[[paste0("DeltaRGR_", st)]]
        names(x_raw) <- shared_traits_wide$asv_id
        shared_signal_list[[length(shared_signal_list) + 1L]] <-
          calc_phylo_signal(
            tree = shared_tree,
            trait = x_raw,
            stage = st,
            stage_label = FIRE_LABELS[[st]],
            trait_version = "raw",
            analysis_set = "shared_ASV_sensitivity"
          )
      }
      
      shared_signal <- rbindlist(shared_signal_list, use.names = TRUE, fill = TRUE)
      shared_signal[, fire_trt := factor(fire_trt, levels = FIRE_LEVELS)]
      setorder(shared_signal, fire_trt)
      fwrite(
        shared_signal,
        file.path(output_dir, "07_phylogenetic_signal_shared_ASV_sensitivity.tsv"),
        sep = "\t"
      )
    } else {
      message(">>> Shared-ASV sensitivity skipped: RUN_SHARED_SENSITIVITY = FALSE")
    }
    
    plot_K <- main_signal[ trait_version == "raw" & analysis_set == "stage_specific_all_active" ]
    plot_K[ , stage_label := factor( stage_label, levels = FIRE_LABELS[FIRE_LEVELS] ) ]
    plot_K[ , p_label := vapply( K_P_raw, format_p, FUN.VALUE = character(1) ) ]
    plot_K[
      ,
      sig_class := factor(
        ifelse( is.finite(K_P_raw) & K_P_raw < 0.05, "P < 0.05", "P >= 0.05" ),
        levels = c( "P >= 0.05", "P < 0.05" )
      )
    ]
    plot_K[
      ,
      point_label := ifelse(
        is.finite(Blomberg_K),
        paste0(
          "K = ",
          vapply( Blomberg_K, format_signal_value, FUN.VALUE = character(1) ),
          "\n",
          p_label,
          "\nN = ", N_ASV
        ),
        "K = NA"
      )
    ]
    k_values <- plot_K$Blomberg_K[ is.finite(plot_K$Blomberg_K) ]
    k_data_max <- if (length(k_values) > 0) {
      max(k_values, na.rm = TRUE)
    } else {
      NA_real_
    }
    
    if (!is.finite(k_data_max) || k_data_max <= 0) {
      k_upper <- 0.01
    } else {
      k_upper <- k_data_max * 1.85
    }
    
    k_axis_labeller <- make_signal_axis_labels( k_data_max )
    k_label_nudge <- k_upper * 0.055
    p_K <- ggplot( plot_K, aes( x = stage_label, y = Blomberg_K ) ) +
      geom_hline( yintercept = 0, linewidth = 0.55, colour = "grey35" ) +
      geom_point( aes(fill = sig_class), shape = 21, size = 5.8, stroke = 1.0, colour = "black" ) +
      geom_text(
        aes( y = Blomberg_K + k_label_nudge, label = point_label ),
        vjust = 0,
        size = 3.9,
        lineheight = 1.0
      ) +
      scale_fill_manual( values = c( "P >= 0.05" = "white", "P < 0.05" = "grey25" ), drop = FALSE, guide = "none" ) +
      scale_y_continuous(
        limits = c( 0, k_upper ),
        labels = k_axis_labeller,
        breaks = scales::breaks_pretty(n = 5),
        expand = expansion(mult = c(0, 0))
      ) +
      labs(
        x = NULL,
        y = expression("Blomberg's K for " * Delta * RGR),
        caption = paste0(
          "K = 1 is the Brownian-motion expectation and lies outside ",
          "the displayed y-range when not shown. Filled points indicate raw P < 0.05."
        )
      ) +
      coord_cartesian(clip = "off") +
      theme_classic(base_size = 12) +
      theme(
        axis.title = element_text( size = 12, face = "bold" ),
        axis.text.x = element_text( size = 10.5, colour = "black", face = "bold" ),
        axis.text.y = element_text( size = 10, colour = "black" ),
        panel.grid.major.y = element_line( colour = "grey90", linewidth = 0.38 ),
        panel.grid.minor = element_blank(),
        panel.border = element_rect( fill = NA, colour = "black", linewidth = 0.7 ),
        plot.caption = element_text( size = 8.8, colour = "grey35", hjust = 0 ),
        plot.margin = margin( 10, 12, 10, 10 )
      )
    save_figure(
      filename = file.path( output_dir, "Fig_B_BlombergK_DeltaRGR_scatter.pdf" ),
      plot = p_K,
      width = 6.4,
      height = 5.2,
      units = "in"
    )
    if (RUN_PHYLO_CORRELOGRAM && nrow(correlogram_table) > 0L) {
      cor_y_values <- c(
        correlogram_table$CI_low,
        correlogram_table$CI_high,
        correlogram_table$Moran_I,
        correlogram_table$H0
      )
      cor_y_values <- cor_y_values[is.finite(cor_y_values)]
      if (length(cor_y_values) == 0L) {
        cor_y_min <- -0.1
        cor_y_max <- 0.1
      } else {
        cor_y_min_data <- min(cor_y_values)
        cor_y_max_data <- max(cor_y_values)
        cor_span <- cor_y_max_data - cor_y_min_data
        if (!is.finite(cor_span) || cor_span <= 0) {
          cor_span <- max(abs(cor_y_values), 0.1)
        }
        cor_y_min <- max(-1, cor_y_min_data - 0.08 * cor_span)
        cor_y_max <- min(1, cor_y_max_data + 0.08 * cor_span)
      }
      
      correlogram_table[
        ,
        stage_label_plot := paste0( as.character(stage_label), "\nN = ", N_ASV, " / ", N_ASV_full_Klambda )
      ]
      stage_plot_levels <- vapply(
        FIRE_LEVELS,
        function(st) {
          row <- correlogram_table[fire_trt == st][1]
          if (nrow(row) == 0L) return(FIRE_LABELS[[st]])
          paste0(FIRE_LABELS[[st]], "\nN = ", row$N_ASV, " / ", row$N_ASV_full_Klambda)
        },
        FUN.VALUE = character(1)
      )
      correlogram_table[, stage_label_plot := factor(stage_label_plot, levels = stage_plot_levels)]
      cor_h0_dt <- unique( correlogram_table[ , .( fire_trt, stage_label_plot, H0, N_ASV, N_ASV_full_Klambda ) ] )
      p_cor <- ggplot( correlogram_table, aes( x = phylogenetic_distance, y = Moran_I ) ) +
        geom_ribbon( aes( ymin = CI_low, ymax = CI_high ), fill = "grey82", alpha = 0.70, colour = NA ) +
        geom_hline(
          data = cor_h0_dt,
          aes(yintercept = H0),
          linetype = "dashed",
          linewidth = 0.55,
          colour = "grey35"
        ) +
        geom_line( linewidth = 0.85, colour = "black" ) +
        geom_point( aes(colour = correlogram_class), size = 1.65 ) +
        facet_wrap( ~ stage_label_plot, nrow = 1, scales = "fixed" ) +
        scale_colour_manual(
          values = c(
            "Positive autocorrelation" = "#D73027",
            "Not significant" = "grey45",
            "Negative autocorrelation" = "#4575B4"
          ),
          breaks = c( "Positive autocorrelation", "Not significant", "Negative autocorrelation" ),
          drop = FALSE,
          name = NULL
        ) +
        scale_x_continuous(
          breaks = scales::breaks_pretty(n = 5),
          labels = scales::label_number(accuracy = 0.01),
          expand = expansion(mult = c(0.01, 0.02))
        ) +
        scale_y_continuous(
          breaks = scales::breaks_pretty(n = 6),
          labels = scales::label_number(accuracy = 0.01),
          expand = expansion(mult = c(0, 0))
        ) +
        coord_cartesian( ylim = c(cor_y_min, cor_y_max), clip = "off" ) +
        labs(
          x = "Patristic phylogenetic distance",
          y = expression("Moran's I for " * Delta * RGR),
          caption = paste0(
            "Grey ribbon = ",
            round(100 * CORRELOGRAM_CI_CONF),
            "% bootstrap confidence envelope; dashed line = H0 = -1/(N-1). ",
            "Red/blue points indicate confidence intervals entirely above/below H0. ",
            "Facet N is correlogram N / full K-lambda N. ",
            "When subsampling is required, ASVs are selected reproducibly with Phylum stratification. ",
            "Adjacent phylogenetic-distance points are not independent tests."
          )
        ) +
        theme_classic(base_size = 12) +
        theme(
          strip.background = element_rect( fill = "grey95", colour = "black", linewidth = 0.65 ),
          strip.text = element_text( size = 11, face = "bold", colour = "black" ),
          axis.title = element_text( size = 12, face = "bold" ),
          axis.text = element_text( size = 9.5, colour = "black" ),
          panel.grid.major.y = element_line( colour = "grey92", linewidth = 0.35 ),
          panel.grid.minor = element_blank(),
          panel.border = element_rect( fill = NA, colour = "black", linewidth = 0.7 ),
          legend.position = "top",
          legend.text = element_text(size = 9.5),
          plot.caption = element_text( size = 8.5, colour = "grey35", hjust = 0 ),
          plot.margin = margin(10, 12, 10, 10)
        )
      save_figure(
        filename = file.path( output_dir, "Fig_C_phylogenetic_correlogram_DeltaRGR.pdf" ),
        plot = p_cor,
        width = 10.8,
        height = 4.6,
        units = "in"
      )
    }
    
    plot_lambda <- copy(plot_K)
    plot_lambda[ , p_label := vapply( lambda_P_raw, format_p, FUN.VALUE = character(1) ) ]
    plot_lambda[
      ,
      sig_class := factor(
        ifelse( is.finite(lambda_P_raw) & lambda_P_raw < 0.05, "P < 0.05", "P >= 0.05" ),
        levels = c( "P >= 0.05", "P < 0.05" )
      )
    ]
    plot_lambda[
      ,
      point_label := ifelse(
        is.finite(Pagel_lambda),
        paste0(
          "lambda = ",
          vapply( Pagel_lambda, format_signal_value, FUN.VALUE = character(1) ),
          "\n",
          p_label,
          "\nN = ", N_ASV
        ),
        "lambda = NA"
      )
    ]
    lambda_values <- plot_lambda$Pagel_lambda[ is.finite(plot_lambda$Pagel_lambda) ]
    lambda_data_max <- if (length(lambda_values) > 0) {
      max(lambda_values, na.rm = TRUE)
    } else {
      NA_real_
    }
    
    if (!is.finite(lambda_data_max) || lambda_data_max <= 0) {
      lambda_upper <- 0.10
    } else {
      lambda_upper <- min( 1, max( lambda_data_max * 1.65, lambda_data_max + 0.06 ) )
    }
    
    lambda_label_nudge <- lambda_upper * 0.055
    p_lambda <- ggplot( plot_lambda, aes( x = stage_label, y = Pagel_lambda ) ) +
      geom_hline( yintercept = 0, linewidth = 0.55, colour = "grey35" ) +
      geom_point( aes(fill = sig_class), shape = 21, size = 5.8, stroke = 1.0, colour = "black" ) +
      geom_text(
        aes( y = Pagel_lambda + lambda_label_nudge, label = point_label ),
        vjust = 0,
        size = 3.9,
        lineheight = 1.0
      ) +
      scale_fill_manual( values = c( "P >= 0.05" = "white", "P < 0.05" = "grey25" ), drop = FALSE, guide = "none" ) +
      scale_y_continuous(
        limits = c( 0, lambda_upper ),
        breaks = scales::breaks_pretty(n = 5),
        labels = scales::label_number(accuracy = 0.01),
        expand = expansion(mult = c(0, 0))
      ) +
      labs(
        x = NULL,
        y = expression("Pagel's " * lambda * " for " * Delta * RGR),
        caption = paste0(
          "lambda = 1 is the Brownian-motion reference and lies outside ",
          "the displayed y-range when not shown. Filled points indicate raw P < 0.05."
        )
      ) +
      coord_cartesian(clip = "off") +
      theme_classic(base_size = 12) +
      theme(
        axis.title = element_text( size = 12, face = "bold" ),
        axis.text.x = element_text( size = 10.5, colour = "black", face = "bold" ),
        axis.text.y = element_text( size = 10, colour = "black" ),
        panel.grid.major.y = element_line( colour = "grey90", linewidth = 0.38 ),
        panel.grid.minor = element_blank(),
        panel.border = element_rect( fill = NA, colour = "black", linewidth = 0.7 ),
        plot.caption = element_text( size = 8.8, colour = "grey35", hjust = 0 ),
        plot.margin = margin( 10, 12, 10, 10 )
      )
    save_figure(
      filename = file.path( output_dir, "Fig_S1_PagelLambda_DeltaRGR_scatter.pdf" ),
      plot = p_lambda,
      width = 6.4,
      height = 5.2,
      units = "in"
    )
    parameter_dt <- data.table(
      parameter = c(
        "primary_trait",
        "DeltaRGR_definition",
        "DeltaRGR_estimator_for_phylogeny",
        "bootstrap_combination",
        "N_BOOT_COMBINE",
        "DeltaRGR_CI",
        "active_filter",
        "main_inference_taxon_set",
        "stage_specific_N_UB",
        "stage_specific_N_F3M",
        "stage_specific_N_F3Y",
        "shared_ASV_count_optional_sensitivity",
        "Blomberg_K_nsim",
        "Pagel_lambda_test",
        "multiple_testing_correction",
        "RUN_ASINH_SENSITIVITY",
        "RUN_SHARED_SENSITIVITY",
        "tree_tip_count_raw",
        "tree_tip_count_shared_optional_sensitivity",
        "nonpositive_edges_replaced",
        "edge_replacement_value",
        "RUN_PHYLO_CORRELOGRAM",
        "correlogram_trait",
        "correlogram_distance",
        "correlogram_max_ASVs",
        "correlogram_stratify_by_phylum",
        "correlogram_n_points",
        "correlogram_CI_bootstrap",
        "correlogram_CI_confidence",
        "correlogram_null_expectation",
        "correlogram_interpretation",
        "random_seed"
      ),
      value = c(
        "DeltaRGR",
        "RGR_T_high - RGR_T_low",
        "median of independently combined T_high-minus-T_low RGR bootstrap distribution",
        "independent resampling of T_high and T_low bootstrap columns",
        as.character(N_BOOT_COMBINE),
        "90% CI = 5th to 95th percentile",
        "Active_T_low OR Active_T_high",
        "all eligible active/comparable ASVs within each stage",
        as.character(stage_eligible[fire_trt_chr == "UB", uniqueN(asv_id)]),
        as.character(stage_eligible[fire_trt_chr == "F3M", uniqueN(asv_id)]),
        as.character(stage_eligible[fire_trt_chr == "F3Y", uniqueN(asv_id)]),
        as.character(length(shared_ids)),
        as.character(K_NSIM),
        "phylolm::phylolm(model='lambda'); LRT versus near-zero lambda null",
        "none; raw P values",
        as.character(RUN_ASINH_SENSITIVITY),
        as.character(RUN_SHARED_SENSITIVITY),
        as.character(length(tree_raw$tip.label)),
        as.character(length(shared_tree$tip.label)),
        as.character(length(nonpositive_edges)),
        ifelse( is.finite(edge_replacement), format(edge_replacement, scientific = TRUE), "NA" ),
        as.character(RUN_PHYLO_CORRELOGRAM),
        "raw bootstrap-median DeltaRGR; reproducible Phylum-stratified capped subset when stage N exceeds correlogram_max_ASVs",
        CORRELOGRAM_DISTANCE,
        as.character(CORRELOGRAM_MAX_ASVS),
        as.character(CORRELOGRAM_STRATIFY_BY_PHYLUM),
        as.character(CORRELOGRAM_N_POINTS),
        as.character(CORRELOGRAM_CI_BS),
        as.character(CORRELOGRAM_CI_CONF),
        "Moran I H0 = -1/(N-1)",
        "CI entirely above H0 = positive autocorrelation; below H0 = negative; adjacent lags are not independent tests",
        as.character(RANDOM_SEED)
      )
    )
    fwrite( parameter_dt, file.path( output_dir, "10_analysis_parameters_and_tree_QC.tsv" ), sep = "\t" )
    invisible(NULL)
  }
  
  # Contributors
  run_contributors <- function() {
    output_dir <- CACHE$contributors
    
    safe_quantiles <- function(x, probs = CI_PROBS) {
      x <- x[is.finite(x)]
      if(length(x) == 0) return(rep(NA_real_, length(probs)))
      as.numeric(quantile(x, probs = probs, type = 7))
    }
    
    growth <- fread(file.path(CACHE$growth, "01_RGR_G_ASV_all.tsv"))
    boot_dir <- file.path(CACHE$growth, "bootstrap_objects")
    delta_boot_list <- list()
    for (st in FIRE_LEVELS) {
      boot21 <- readRDS(file.path(boot_dir, paste0("G_bootstrap_", st, "_T_low.rds")))
      boot28 <- readRDS(file.path(boot_dir, paste0("G_bootstrap_", st, "_T_high.rds")))
      act21 <- growth[fire_trt == st & tem == "T_low" & Active_ASV == TRUE, asv_id]
      act28 <- growth[fire_trt == st & tem == "T_high" & Active_ASV == TRUE, asv_id]
      union_asvs <- sort(unique(c(act21, act28)))
      n_boot <- ncol(boot21$G_boot)
      G21_union <- matrix(0, nrow = length(union_asvs), ncol = n_boot, dimnames = list(union_asvs, NULL))
      G28_union <- matrix(0, nrow = length(union_asvs), ncol = n_boot, dimnames = list(union_asvs, NULL))
      idx21 <- match(act21, boot21$asv_id)
      G21_union[match(act21, union_asvs), ] <- boot21$G_boot[idx21, , drop = FALSE]
      idx28 <- match(act28, boot28$asv_id)
      G28_union[match(act28, union_asvs), ] <- boot28$G_boot[idx28, , drop = FALSE]
      delta_mat <- G28_union - G21_union
      delta_boot_list[[st]] <- delta_mat
    }
    
    # Decompose detected-active growth changes into positive and negative parts.
    decomp_list <- list()
    for (st in FIRE_LEVELS) {
      mat <- delta_boot_list[[st]]
      P_boot <- colSums(pmax(mat, 0), na.rm = TRUE)
      L_boot <- colSums(pmax(-mat, 0), na.rm = TRUE)
      Net_boot <- P_boot - L_boot
      P_q <- safe_quantiles(P_boot)
      L_q <- safe_quantiles(L_boot)
      Net_q <- safe_quantiles(Net_boot)
      decomp_list[[length(decomp_list) + 1]] <- data.table(
        fire_trt = st,
        Component = c("Positive Flux", "Negative Flux", "Net Flux"),
        Median = c(P_q[2], -L_q[2], Net_q[2]),
        Low90  = c(P_q[1], -L_q[3], Net_q[1]),
        High90 = c(P_q[3], -L_q[1], Net_q[3])
      )
    }
    
    decomp_dt <- rbindlist(decomp_list)
    decomp_dt[, fire_trt := factor(fire_trt, levels = FIRE_LEVELS, labels = unname(FIRE_LABELS[FIRE_LEVELS]))]
    p_decomp <- ggplot(decomp_dt, aes(x = fire_trt)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
      geom_col(data = decomp_dt[Component %chin% c("Positive Flux", "Negative Flux")],
               aes(y = Median, fill = Component), width = 0.5, alpha = 0.8, color = "black", linewidth = 0.4) +
      geom_errorbar(data = decomp_dt[Component %chin% c("Positive Flux", "Negative Flux")],
                    aes(ymin = Low90, ymax = High90), width = 0.15, linewidth = 0.6) +
      geom_pointrange(data = decomp_dt[Component == "Net Flux"],
                      aes(y = Median, ymin = Low90, ymax = High90, shape = Component),
                      size = 0.8, color = "black", fill = "white", stroke = 1) +
      scale_fill_manual(values = c("Positive Flux" = "#D73027", "Negative Flux" = "#4575B4")) +
      scale_shape_manual(values = c("Net Flux" = 21)) +
      scale_y_continuous(labels = scales::scientific, expand = expansion(mult = c(0.1, 0.1))) +
      labs(x = NULL,
           y = bquote(Delta * "G (16S copies "*g^{-1}*" "*d^{-1}*")"),
           title = "Detected active-growth change",
           subtitle = bquote(Delta*G[Net] == sum(max(Delta*G[i],0)) - sum(max(-Delta*G[i],0)))) +
      theme_classic(base_size = 12) +
      theme(legend.position = "right",
            legend.title = element_blank(),
            axis.text.x = element_text(size = 11, color="black", face="bold"))
    save_figure(file.path(output_dir, "Fig_F2_Addon_Flux_Decomposition.pdf"), p_decomp, width = 5.7, height = 4)
    fwrite(decomp_dt, file.path(output_dir, "Table_F2_Flux_Decomposition.tsv"), sep = "\t")
    top_n <- 15L
    
    # Select credible contributors and quantify rank stability.
    contributor_list <- list()
    for (st in FIRE_LEVELS) {
      mat <- delta_boot_list[[st]]
      if (!nrow(mat)) next
      n_boot <- ncol(mat)
      dG_median <- apply(mat, 1, median, na.rm = TRUE)
      dG_low90  <- apply(mat, 1, function(x) safe_quantiles(x)[1])
      dG_high90 <- apply(mat, 1, function(x) safe_quantiles(x)[3])
      P_boot <- colSums(pmax(mat, 0), na.rm = TRUE)
      L_boot <- colSums(pmax(-mat, 0), na.rm = TRUE)
      valid_pos_asvs <- names(dG_median)[dG_low90 > 0]
      pos_candidates <- dG_median[valid_pos_asvs]
      pos_asvs <- head(names(sort(pos_candidates, decreasing = TRUE)), top_n)
      ranks_pos <- apply(mat, 2, function(x) rank(-x, ties.method = "min"))
      if (is.null(dim(ranks_pos))) ranks_pos <- matrix(ranks_pos, nrow = nrow(mat), dimnames = list(rownames(mat), NULL))
      for (i in seq_along(pos_asvs)) {
        asv <- pos_asvs[i]
        flux_b <- pmax(mat[asv, ], 0)
        pct_b <- ifelse(P_boot > 0, flux_b / P_boot * 100, NA_real_)
        freq_in_top15 <- mean(ranks_pos[asv, ] <= top_n & mat[asv, ] > 0) * 100
        q_abs <- safe_quantiles(flux_b); q_pct <- safe_quantiles(pct_b)
        contributor_list[[length(contributor_list) + 1]] <- data.table(
          fire_trt = st, Direction = "Positive", ASV_Rank = i, asv_id = asv,
          Abs_Median = q_abs[2], Abs_Low90 = q_abs[1], Abs_High90 = q_abs[3],
          Log10_Abs_Median = log10(1 + q_abs[2]), Log10_Abs_Low = log10(1 + q_abs[1]), Log10_Abs_High = log10(1 + q_abs[3]),
          Pct_Median = q_pct[2], Pct_Low = q_pct[1], Pct_High = q_pct[3], Stability_Freq_Top15_Pct = freq_in_top15
        )
      }
      
      valid_neg_asvs <- names(dG_median)[dG_high90 < 0]
      neg_candidates <- dG_median[valid_neg_asvs]
      neg_asvs <- head(names(sort(neg_candidates, decreasing = FALSE)), top_n)
      ranks_neg <- apply(mat, 2, function(x) rank(x, ties.method = "min"))
      if (is.null(dim(ranks_neg))) ranks_neg <- matrix(ranks_neg, nrow = nrow(mat), dimnames = list(rownames(mat), NULL))
      for (i in seq_along(neg_asvs)) {
        asv <- neg_asvs[i]
        flux_b <- pmax(-mat[asv, ], 0)
        pct_b <- ifelse(L_boot > 0, flux_b / L_boot * 100, NA_real_)
        freq_in_top15 <- mean(ranks_neg[asv, ] <= top_n & mat[asv, ] < 0) * 100
        q_abs <- safe_quantiles(flux_b); q_pct <- safe_quantiles(pct_b)
        contributor_list[[length(contributor_list) + 1]] <- data.table(
          fire_trt = st, Direction = "Negative", ASV_Rank = i, asv_id = asv,
          Abs_Median = q_abs[2], Abs_Low90 = q_abs[1], Abs_High90 = q_abs[3],
          Log10_Abs_Median = log10(1 + q_abs[2]), Log10_Abs_Low = log10(1 + q_abs[1]), Log10_Abs_High = log10(1 + q_abs[3]),
          Pct_Median = q_pct[2], Pct_Low = q_pct[1], Pct_High = q_pct[3], Stability_Freq_Top15_Pct = freq_in_top15
        )
      }
    }
    
    contributors_dt <- if (length(contributor_list)) rbindlist(contributor_list) else data.table(
      fire_trt = character(), Direction = character(), ASV_Rank = integer(), asv_id = character(),
      Abs_Median = numeric(), Abs_Low90 = numeric(), Abs_High90 = numeric(),
      Log10_Abs_Median = numeric(), Log10_Abs_Low = numeric(), Log10_Abs_High = numeric(),
      Pct_Median = numeric(), Pct_Low = numeric(), Pct_High = numeric(), Stability_Freq_Top15_Pct = numeric())
    asv_tax <- unique(growth[, .(asv_id, Phylum, Family, Genus)])
    asv_tax[, Genus := fifelse(is.na(Genus) | Genus == "" | Genus == "Unclassified", NA_character_, Genus)]
    asv_tax[, Family := fifelse(is.na(Family) | Family == "" | Family == "Unclassified", NA_character_, Family)]
    asv_tax[, Display_Label := fcase(
      !is.na(Genus), paste(Genus, asv_id),
      !is.na(Family), paste(Family, "NA", asv_id),
      default = paste(Phylum, "NA", asv_id)
    )]
    contributors_dt <- merge(contributors_dt, asv_tax[, .(asv_id, Display_Label)], by = "asv_id", all.x = TRUE)
    contributors_dt[, Full_X_Label := sprintf("%02d. %s", ASV_Rank, Display_Label)]
    contributors_dt[, Full_X_Label := factor(Full_X_Label, levels = unique(Full_X_Label[order(ASV_Rank)]))]
    contributors_dt[, fire_trt := factor(fire_trt, levels = FIRE_LEVELS, labels = unname(FIRE_LABELS[FIRE_LEVELS]))]
    p_pos_contributors <- ggplot(contributors_dt[Direction == "Positive"], aes(x = Full_X_Label, y = Log10_Abs_Median)) +
      geom_pointrange(aes(ymin = Log10_Abs_Low, ymax = Log10_Abs_High), color = "#D73027", size = 0.5) +
      geom_text(aes(label = sprintf("%.1f%%", Pct_Median), y = Log10_Abs_High + 0.2), size = 2.7, angle = 90, hjust = 0, color = "black") +
      facet_wrap(~ fire_trt, scales = "free_x") +
      scale_y_continuous(expand = expansion(mult = c(0.1, 0.35))) +
      labs(x = "ASV Rank & Identity", y = bquote("Positive "*Delta*"G: log"[10]*"(1 + magnitude)"), title = "Top 15 Positive Contributors") +
      theme_bw(base_size = 12) +
      theme(
        strip.background = element_rect(fill = "grey95"),
        axis.text.x = element_text(size = 9.5, angle = 50, hjust = 1, color = "black")
      )
    save_figure(file.path(output_dir, "Fig_F3_Pos_Contributors_PointRange.pdf"), if (any(contributors_dt$Direction == "Positive")) p_pos_contributors else blank_plot("No positive contributors pass the 90% CI criterion"), width = 12, height = 8)
    p_neg_contributors <- ggplot(contributors_dt[Direction == "Negative"], aes(x = Full_X_Label, y = Log10_Abs_Median)) +
      geom_pointrange(aes(ymin = Log10_Abs_Low, ymax = Log10_Abs_High), color = "#4575B4", size = 0.5) +
      geom_text(aes(label = sprintf("%.1f%%", Pct_Median), y = Log10_Abs_High + 0.2), size = 2.7, angle = 90, hjust = 0, color = "black") +
      facet_wrap(~ fire_trt, scales = "free_x") +
      scale_y_continuous(expand = expansion(mult = c(0.1, 0.35))) +
      labs(x = "ASV Rank & Identity", y = bquote("Negative "*Delta*"G: log"[10]*"(1 + magnitude)"), title = "Top 15 Negative Contributors") +
      theme_bw(base_size = 12) +
      theme(
        strip.background = element_rect(fill = "grey95"),
        axis.text.x = element_text(size = 9.5, angle = 50, hjust = 1, color = "black")
      )
    save_figure(file.path(output_dir, "Fig_F3_Neg_Contributors_PointRange.pdf"), if (any(contributors_dt$Direction == "Negative")) p_neg_contributors else blank_plot("No negative contributors pass the 90% CI criterion"), width = 12, height = 8)
    fwrite(contributors_dt, file.path(output_dir, "Table_F3_All_Top15_Contributors_and_Stability.tsv"), sep = "\t")
    invisible(NULL)
  }
  
  runners <- list(temperature = run_temperature, density = run_density, afe = run_afe,
                  cue = run_cue, growth = run_growth, community = run_community,
                  strategy = run_strategy, phylo_g = run_phylo_g,
                  phylo_rgr = run_phylo_rgr, contributors = run_contributors)
  module_log <- data.table::data.table(module = character(), status = character(),
                                       elapsed_seconds = numeric(), message = character())
  on.exit({
    flush_catalog()
    data.table::fwrite(module_log, file.path(OUTPUT_ROOT, "Module_status.tsv"), sep = "\t")
    capture.output(sessionInfo(), file = file.path(OUTPUT_ROOT, "Session_info.txt"))
  }, add = TRUE)
  for (module_name in selected) {
    message("Running ", module_name, " ...")
    start_time <- proc.time()[["elapsed"]]
    error_text <- NULL
    warnings <- character()
    tryCatch(withCallingHandlers({
      deps <- dependencies[[module_name]]
      failed_deps <- module_log[module %in% deps & status != "completed", module]
      if (length(failed_deps)) stop("Required module failed: ", paste(failed_deps, collapse = ", "))
      runners[[module_name]]()
    }, warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
    }), error = function(e) error_text <<- conditionMessage(e))
    copy_module_tables(module_name)
    if (length(warnings)) writeLines(unique(warnings), file.path(CACHE[[module_name]], "Warnings.txt"))
    module_log <- rbind(module_log, data.table::data.table(
      module = module_name, status = if (is.null(error_text)) "completed" else "failed",
      elapsed_seconds = proc.time()[["elapsed"]] - start_time,
      message = if (is.null(error_text)) "" else error_text))
    if (!is.null(error_text)) {
      figure_catalog[module == module_name & status == "not_generated", note := error_text]
      if (STOP_ON_ERROR) stop("Module ", module_name, " failed: ", error_text, call. = FALSE)
      warning("Module ", module_name, " failed: ", error_text)
    }
    data.table::fwrite(module_log, file.path(OUTPUT_ROOT, "Module_status.tsv"), sep = "\t")
    flush_catalog()
    gc(verbose = FALSE)
  }
  message("Finished. Results: ", normalizePath(OUTPUT_ROOT, winslash = "/"))
  invisible(OUTPUT_ROOT)
}

main()
