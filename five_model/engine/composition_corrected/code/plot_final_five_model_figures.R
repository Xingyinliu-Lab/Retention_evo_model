suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(readr)})
options(warn = -1)

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_file <- normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = TRUE)
package_root <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/", mustWork = TRUE)
model_file <- file.path(package_root, "output", "five_model_fits.tsv")
p_file <- file.path(package_root, "output", "paired_tests_vs_sequential_retention.tsv")
fig_dir <- file.path(package_root, "figures")
src_dir <- file.path(package_root, "source_data")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(src_dir, recursive = TRUE, showWarnings = FALSE)

tree_labels <- c(
  gtdb_r226_expanded = "GTDB r226",
  asgard971_addcore_gtdb_rooted = "Asgard971 add-core",
  asgard971_astral89_gtdb_rooted = "Asgard971 ASTRAL89",
  archaea1338_addcore_gtdb_rooted = "Archaea1338 add-core",
  archaea1338_denovo_gtdb_rooted = "Archaea1338 de novo",
  archaea1338_astral64_gtdb_rooted = "Archaea1338 ASTRAL64"
)
tree_cols <- c(
  "GTDB r226" = "#4D4D4D", "Asgard971 add-core" = "#3B73A3",
  "Asgard971 ASTRAL89" = "#56A6B8", "Archaea1338 add-core" = "#D08A3C",
  "Archaea1338 de novo" = "#C85A63", "Archaea1338 ASTRAL64" = "#8B6DAE"
)
tree_shapes <- c(
  "GTDB r226" = 16, "Asgard971 add-core" = 17, "Asgard971 ASTRAL89" = 15,
  "Archaea1338 add-core" = 18, "Archaea1338 de novo" = 8,
  "Archaea1338 ASTRAL64" = 7
)
model_labels <- c(
  adjacency_coexist_middle = "Sequential retention\n(coexist middle)",
  adjacency_dALP_middle = "dALP-middle\nadjacency",
  adjacency_cALP_middle = "cALP-middle\nadjacency",
  ancestral_coexist_divergence = "Ancestral coexist\ndivergence",
  ARD_unrestricted = "ARD unrestricted"
)

models <- read_tsv(model_file, col_types = cols()) %>%
  mutate(
    tree_label = factor(unname(tree_labels[tree]), levels = unname(tree_labels)),
    model_label = unname(model_labels[model]),
    y_plot = log10(delta_AIC + 1)
  )
p_values <- read_tsv(p_file, col_types = cols()) %>%
  transmute(
    branch_scale,
    model_label = recode(
      comparison_model,
      "ARD unrestricted" = "ARD unrestricted",
      "dALP-middle adjacency" = "dALP-middle\nadjacency",
      "cALP-middle adjacency" = "cALP-middle\nadjacency",
      "Ancestral coexist divergence" = "Ancestral coexist\ndivergence"
    ),
    exact_p_two_sided,
    p_label = paste0("P = ", formatC(exact_p_two_sided, format = "f", digits = 3))
  )

order_models <- function(x) {
  fixed <- c("ARD unrestricted", "Sequential retention\n(coexist middle)")
  remaining <- setdiff(unique(x$model_label), fixed)
  medians <- tapply(x$delta_AIC, x$model_label, median, na.rm = TRUE)[remaining]
  c(fixed, names(sort(medians)))
}

plot_scale <- function(scale_name) {
  x <- models %>% filter(branch_scale == scale_name)
  levels_x <- order_models(x)
  x <- x %>% mutate(model_label = factor(model_label, levels = levels_x))
  annotation_y <- max(x$y_plot, na.rm = TRUE) + 0.18
  annotations <- tibble(
    model_label = factor(levels_x, levels = levels_x),
    p_label = ifelse(levels_x == "Sequential retention\n(coexist middle)", "Reference", NA_character_)
  ) %>%
    left_join(
      p_values %>% filter(branch_scale == scale_name) %>%
        select(model_label, paired_p_label = p_label),
      by = "model_label"
    ) %>%
    mutate(p_label = coalesce(p_label, paired_p_label), x = seq_along(levels_x), y = annotation_y)

  write_tsv(x %>% mutate(model_label = as.character(model_label)),
            file.path(src_dir, paste0("five_model_", scale_name, "_source.tsv")))
  write_tsv(annotations %>% select(model_label, p_label, x, y),
            file.path(src_dir, paste0("five_model_", scale_name, "_p_annotations.tsv")))

  breaks_raw <- c(0, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000)
  title_text <- ifelse(scale_name == "unit_topology", "Unit topology", "Source branch lengths")
  ggplot(x, aes(as.numeric(model_label), y_plot, group = model_label)) +
    geom_hline(yintercept = 0, size = 0.30, colour = "#666666", linetype = "22") +
    geom_boxplot(width = 0.48, size = 0.30, outlier.shape = NA, fill = NA, colour = "#8F8F8F") +
    geom_point(
      aes(colour = tree_label, shape = tree_label),
      position = position_jitter(width = 0.075, height = 0, seed = 20260820), size = 1.85
    ) +
    geom_text(
      data = annotations, aes(x = x, y = y, label = p_label), inherit.aes = FALSE,
      family = "Arial", size = 2.25, colour = "#252525", vjust = 0
    ) +
    scale_x_continuous(breaks = seq_along(levels_x), labels = levels_x,
                       limits = c(0.65, length(levels_x) + 0.35)) +
    scale_y_continuous(breaks = log10(breaks_raw + 1), labels = breaks_raw,
                       minor_breaks = NULL, expand = expansion(mult = c(0.01, 0.13))) +
    scale_colour_manual(values = tree_cols, drop = FALSE) +
    scale_shape_manual(values = tree_shapes, drop = FALSE) +
    labs(
      x = NULL, y = expression(Delta * "AIC"), colour = NULL, shape = NULL,
      title = title_text,
      caption = paste(
        "P values: unadjusted two-sided exact paired Wilcoxon signed-rank tests",
        "across six matched species-tree reconstructions; each comparison is versus sequential retention.",
        "Trees are sensitivity reconstructions, not independent biological replicates."
      )
    ) +
    theme_classic(base_size = 7.2, base_family = "Arial") +
    theme(
      axis.line = element_line(size = 0.30), axis.ticks = element_line(size = 0.28),
      axis.text.x = element_text(size = 6.3, hjust = 0.5),
      plot.title = element_text(size = 8.2, face = "bold"),
      plot.caption = element_text(size = 5.2, colour = "#555555", hjust = 0),
      legend.position = "top", legend.text = element_text(size = 5.5),
      plot.margin = margin(3, 6, 2, 5)
    ) +
    guides(colour = guide_legend(nrow = 1), shape = guide_legend(nrow = 1))
}

outputs <- c(
  unit_topology = "model_AIC_root_corrected_five_models_unit_topology_log_spacing.pdf",
  source_branch_lengths = "model_AIC_root_corrected_five_models_source_branches_log_spacing.pdf"
)
for (scale_name in names(outputs)) {
  ggsave(file.path(fig_dir, unname(outputs[scale_name])), plot_scale(scale_name),
         width = 183, height = 96, units = "mm", device = cairo_pdf)
}
