suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
})

args <- commandArgs(FALSE)
script <- sub("^--file=", "", args[grep("^--file=", args)][1])
root <- normalizePath(file.path(dirname(script), ".."), mustWork = TRUE)
figure_dir <- file.path(root, "figures", "source_branch_lengths_20260821")
source_dir <- file.path(root, "source_data", "source_branch_lengths_20260821")

unit <- read_tsv(file.path(root, "results",
  "MreB_ALP_composition_corrected_five_models_six_trees.tsv"), show_col_types = FALSE) %>%
  filter(tree == "GTDB_r226") %>% mutate(branch_scale = "Unit topology")
source <- read_tsv(file.path(root, "results", "source_branch_lengths_20260820",
  "MreB_ALP_five_models_source_branch_lengths.tsv"), show_col_types = FALSE) %>%
  filter(tree == "GTDB_r226") %>% mutate(branch_scale = "Source branch lengths")

model_levels <- c("ARD unrestricted", "Sequential\nretention", "ALP-middle\nadjacency",
                  "M-middle\nadjacency", "Ancestral coexist\n(strict)")
dat <- bind_rows(unit, source) %>% mutate(
  model_label = recode(model,
    ARD_unrestricted = "ARD unrestricted",
    adjacency_coexist_middle = "Sequential\nretention",
    adjacency_ALP_middle = "ALP-middle\nadjacency",
    adjacency_M_middle = "M-middle\nadjacency",
    ancestral_coexistence_strict = "Ancestral coexist\n(strict)"),
  model_label = factor(model_label, levels = model_levels),
  branch_scale = factor(branch_scale, levels = c("Unit topology", "Source branch lengths")),
  label = sprintf("%.1f", delta_AIC)
)

unit_max <- max(dat$delta_AIC[dat$branch_scale == "Unit topology"])
source_max <- max(dat$delta_AIC[dat$branch_scale == "Source branch lengths"])
scale_factor <- source_max / unit_max
dat <- dat %>% mutate(
  y_plot = if_else(branch_scale == "Source branch lengths", delta_AIC / scale_factor, delta_AIC)
)
stopifnot(nrow(dat) == 10, is.finite(scale_factor), scale_factor > 0)
write_tsv(dat, file.path(source_dir,
  "MreB_ALP_five_model_AIC_GTDB_unit_vs_source_branch_lengths_dual_axis.tsv"))

cols <- c("Unit topology" = "#777777", "Source branch lengths" = "#4F857D")
p <- ggplot(dat, aes(model_label, y_plot, fill = branch_scale)) +
  geom_col(position = position_dodge(width = 0.74), width = 0.64, colour = NA) +
  geom_text(aes(label = label), position = position_dodge(width = 0.74),
            vjust = -0.45, size = 1.95, family = "Arial", colour = "#303030") +
  scale_fill_manual(values = cols, name = NULL) +
  scale_y_continuous(
    name = expression("Unit-topology "*Delta*"AIC"),
    limits = c(0, unit_max * 1.13),
    breaks = scales::pretty_breaks(6),
    expand = expansion(mult = c(0, 0)),
    sec.axis = sec_axis(~ . * scale_factor,
                        name = expression("Source branch-length "*Delta*"AIC"),
                        breaks = scales::pretty_breaks(6)(c(0, source_max)))
  ) +
  labs(
    title = "GTDB r226 model comparison", x = NULL,
    caption = "Composition-corrected root priors. Bars use separate linear axes; labels show raw delta-AIC values."
  ) +
  theme_classic(base_size = 7.2, base_family = "Arial") +
  theme(
    axis.line = element_line(linewidth = 0.30, colour = "black"),
    axis.ticks = element_line(linewidth = 0.28, colour = "black"),
    axis.text.x = element_text(size = 5.9, lineheight = 0.92),
    axis.text.y.left = element_text(size = 6.2, colour = cols[["Unit topology"]]),
    axis.text.y.right = element_text(size = 6.2, colour = cols[["Source branch lengths"]]),
    axis.title.y.left = element_text(size = 7.0, colour = cols[["Unit topology"]]),
    axis.title.y.right = element_text(size = 7.0, colour = cols[["Source branch lengths"]]),
    plot.title = element_text(size = 8.1, face = "bold"),
    plot.caption = element_text(size = 5.3, colour = "#555555", hjust = 0,
                                margin = margin(t = 4)),
    legend.position = "top", legend.text = element_text(size = 6.0),
    legend.key.width = grid::unit(4.5, "mm"),
    plot.margin = margin(4, 7, 4, 7)
  )

ggsave(file.path(figure_dir,
  "fig_mreb_alp_five_model_AIC_GTDB_unit_vs_source_branch_lengths_dual_axis.pdf"),
  p, width = 148, height = 98, units = "mm", device = cairo_pdf)
ggsave(file.path(figure_dir,
  "fig_mreb_alp_five_model_AIC_GTDB_unit_vs_source_branch_lengths_dual_axis.png"),
  p, width = 148, height = 98, units = "mm", dpi = 600, bg = "white")

writeLines(c("status=complete", "backend=R", "tree=GTDB_r226",
             paste0("source_to_unit_scale_factor=", signif(scale_factor, 12)), "figures=1"),
  file.path(figure_dir, "GTDB_FIVE_MODEL_DUAL_AXIS_COMPLETE.marker"))
