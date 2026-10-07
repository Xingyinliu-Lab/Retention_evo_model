suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/")
root <- dirname(dirname(script_path))
input_file <- if (length(args) >= 1) args[[1]] else file.path(
  root, "source_data", "ARD_transition_contribution_summary_composition_corrected.tsv"
)
figure_stem <- if (length(args) >= 2) args[[2]] else
  "fig_ARD_transition_contribution_composition_corrected_20260901"

figure_dir <- file.path(root, "figures")
source_dir <- file.path(root, "source_data")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

dat <- read.delim(input_file, check.names = FALSE, stringsAsFactors = FALSE) %>%
  mutate(
    tree_label = recode(
      tree,
      gtdb_r226_expanded = "GTDB r226",
      asgard971_addcore_gtdb_rooted = "Asgard971 add-core",
      asgard971_astral89_gtdb_rooted = "Asgard971 ASTRAL89",
      archaea1338_addcore_gtdb_rooted = "Archaea1338 add-core",
      archaea1338_denovo_gtdb_rooted = "Archaea1338 de novo",
      archaea1338_astral64_gtdb_rooted = "Archaea1338 ASTRAL64"
    ),
    branch_label = recode(
      branch_scale,
      unit_topology = "unit topology",
      source_branch_lengths = "source branches"
    ),
    direct_percent = 100 * direct_fraction,
    coexist_percent = 100 - direct_percent,
    stable = history_identifiability == "stable"
  )

tree_ids <- c(
  "GTDB r226", "Asgard971 add-core", "Asgard971 ASTRAL89",
  "Archaea1338 add-core", "Archaea1338 de novo", "Archaea1338 ASTRAL64"
)
tree_short <- c(
  "GTDB r226", "A971 add-core", "A971 ASTRAL89",
  "A1338 add-core", "A1338 de novo", "A1338 ASTRAL64"
)
tree_cols <- setNames(
  c("#565656", "#3F7FB3", "#58AABD", "#D48A35", "#CB5B63", "#8D6BB3"),
  tree_ids
)
dat$tree_label <- factor(dat$tree_label, levels = tree_ids)
dat$branch_label <- factor(dat$branch_label, levels = c("unit topology", "source branches"))

long_dat <- bind_rows(
  dat %>% transmute(
    tree, tree_label, branch_scale, branch_label, stable, history_identifiability,
    route = "Via coexistence", contribution_percent = coexist_percent
  ),
  dat %>% transmute(
    tree, tree_label, branch_scale, branch_label, stable, history_identifiability,
    route = "Direct transition", contribution_percent = direct_percent
  )
) %>%
  mutate(
    route = factor(route, levels = c("Via coexistence", "Direct transition")),
    route_x = as.integer(route),
    history_key = interaction(tree_label, branch_label, drop = TRUE, lex.order = TRUE),
    jitter_offset = seq(-0.21, 0.21, length.out = 12)[as.integer(history_key)],
    x_plot = route_x + jitter_offset
  )

stable_dat <- filter(long_dat, stable)

if (nrow(dat) != 12L || sum(dat$stable) != 9L || sum(!dat$stable) != 3L) {
  stop("Expected 12 computed histories: 9 stable and 3 flagged")
}

write.table(
  long_dat,
  file.path(source_dir, paste0(figure_stem, "_source_data.tsv")),
  sep = "\t", quote = FALSE, row.names = FALSE, fileEncoding = "UTF-8"
)
write.table(
  filter(dat, !stable),
  file.path(source_dir, paste0(figure_stem, "_excluded_histories.tsv")),
  sep = "\t", quote = FALSE, row.names = FALSE, fileEncoding = "UTF-8"
)

p <- ggplot(stable_dat, aes(route_x, contribution_percent, group = route)) +
  geom_boxplot(
    width = 0.36, linewidth = 0.34, outlier.shape = NA,
    fill = NA, colour = "#8A8A8A"
  ) +
  geom_point(
    data = long_dat,
    aes(x = x_plot, y = contribution_percent, colour = tree_label, shape = branch_label),
    inherit.aes = FALSE,
    size = 1.85, alpha = 0.98, stroke = 0.52
  ) +
  annotate(
    "text", x = 1.50, y = 109,
    label = "12/12 computed   |   9 included in boxplots",
    size = 2.25, family = "Arial", fontface = "bold", colour = "#555555"
  ) +
  scale_x_continuous(
    breaks = 1:2,
    labels = c(
      "Via coexistence\n(two adjacent transitions)",
      "Direct dALP-only / cALP-only\ntransition"
    ),
    limits = c(0.72, 2.28)
  ) +
  scale_y_continuous(
    breaks = seq(0, 100, 20), limits = c(0, 112),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_colour_manual(
    values = tree_cols, breaks = tree_ids, labels = tree_short,
    limits = tree_ids, drop = FALSE, name = NULL
  ) +
  scale_shape_manual(
    values = c("unit topology" = 16, "source branches" = 1),
    name = NULL
  ) +
  labs(
    title = "ARD transition-route decomposition",
    subtitle = "Final composition-corrected unrestricted ARD fits",
    x = NULL,
    y = "Share of posterior expected transitions (%)",
    caption = paste0(
      "Points show all 12 tree/branch-scale histories; boxplots summarize the nine numerically stable histories.\n",
      "The three high-turnover or rate-boundary-saturated source-branch histories remain visible as ordinary open circles,\n",
      "but are not included in the boxplots or proportional interpretation."
    )
  ) +
  theme_classic(base_size = 6.5, base_family = "Arial") +
  theme(
    axis.line = element_line(linewidth = 0.32, colour = "#333333"),
    axis.ticks = element_line(linewidth = 0.30, colour = "#333333"),
    axis.text.x = element_text(size = 6.0, lineheight = 0.94, colour = "#222222"),
    axis.text.y = element_text(size = 5.8, colour = "#333333"),
    axis.title.y = element_text(size = 6.3),
    plot.title = element_text(size = 7.7, face = "bold", margin = margin(b = 1.5)),
    plot.subtitle = element_text(size = 6.1, colour = "#555555", margin = margin(b = 3)),
    plot.caption = element_text(size = 5.2, colour = "#555555", hjust = 0, lineheight = 1.10, margin = margin(t = 4)),
    legend.position = "top",
    legend.box = "vertical",
    legend.text = element_text(size = 5.4),
    legend.key.width = grid::unit(3.5, "mm"),
    legend.key.height = grid::unit(2.5, "mm"),
    legend.spacing.x = grid::unit(0.8, "mm"),
    plot.margin = margin(4, 6, 4, 5)
  ) +
  guides(
    colour = guide_legend(nrow = 2, byrow = TRUE, order = 1, override.aes = list(size = 2.0, shape = 16)),
    shape = guide_legend(nrow = 1, order = 2, override.aes = list(colour = "#444444", size = 2.0))
  )

pdf_file <- file.path(figure_dir, paste0(figure_stem, ".pdf"))
png_file <- file.path(figure_dir, paste0(figure_stem, ".png"))
ggsave(pdf_file, p, width = 122, height = 102, units = "mm", device = cairo_pdf)
ggsave(png_file, p, width = 122, height = 102, units = "mm", dpi = 600, bg = "white")

summary_out <- data.frame(
  metric = c(
    "computed_tree_scale_histories",
    "stable_tree_scale_histories",
    "flagged_source_branch_histories",
    "coexist_mediated_majority_stable_histories",
    "median_direct_endpoint_percent_stable"
  ),
  value = c(
    nrow(dat), sum(dat$stable), sum(!dat$stable),
    sum(dat$stable & dat$coexist_percent > 50),
    median(dat$direct_percent[dat$stable])
  )
)
write.table(
  summary_out,
  file.path(source_dir, paste0(figure_stem, "_summary.tsv")),
  sep = "\t", quote = FALSE, row.names = FALSE, fileEncoding = "UTF-8"
)

message(pdf_file)
message(png_file)
