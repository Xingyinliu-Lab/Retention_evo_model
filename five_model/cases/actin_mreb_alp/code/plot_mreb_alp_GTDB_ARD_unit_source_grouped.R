suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
})

args <- commandArgs(FALSE)
script <- sub("^--file=", "", args[grep("^--file=", args)][1])
root <- normalizePath(file.path(dirname(script), ".."), mustWork = TRUE)
figure_dir <- file.path(root, "figures", "source_branch_lengths_20260821")
source_dir <- file.path(root, "source_data", "source_branch_lengths_20260821")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

unit <- read_tsv(
  file.path(root, "results", "MreB_ALP_ARD_transition_contribution_summary_six_trees.tsv"),
  show_col_types = FALSE
) %>% filter(tree == "GTDB_r226") %>%
  transmute(branch_scale = "Unit topology",
            `Coexist-mediated` = 100 * coexist_mediated_fraction,
            `Direct endpoint transition` = 100 * direct_endpoint_fraction)

source <- read_tsv(
  file.path(root, "results", "source_branch_lengths_20260820",
            "MreB_ALP_ARD_transition_contribution_source_branch_lengths.tsv"),
  show_col_types = FALSE
) %>% filter(tree == "GTDB_r226") %>%
  transmute(branch_scale = "Source branch lengths",
            `Coexist-mediated` = 100 * coexist_mediated_fraction,
            `Direct endpoint transition` = 100 * direct_endpoint_fraction)

dat <- bind_rows(unit, source) %>%
  pivot_longer(-branch_scale, names_to = "route", values_to = "percent") %>%
  mutate(
    route = factor(route, levels = c("Coexist-mediated", "Direct endpoint transition")),
    branch_scale = factor(branch_scale, levels = c("Unit topology", "Source branch lengths")),
    label = sprintf("%.1f%%", percent)
  )
stopifnot(nrow(dat) == 4, all(is.finite(dat$percent)),
          all(abs(dat %>% group_by(branch_scale) %>% summarise(x=sum(percent)) %>% pull(x) - 100) < 1e-7))
write_tsv(dat, file.path(source_dir,
  "MreB_ALP_ARD_transition_contribution_GTDB_unit_vs_source_branch_lengths.tsv"))

cols <- c("Unit topology" = "#777777", "Source branch lengths" = "#4F857D")
p <- ggplot(dat, aes(route, percent, fill = branch_scale)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.62, colour = NA) +
  geom_text(aes(label = label), position = position_dodge(width = 0.72),
            vjust = -0.45, size = 2.25, family = "Arial", colour = "#303030") +
  scale_fill_manual(values = cols, name = NULL) +
  scale_x_discrete(labels = c("Coexist-mediated", "Direct endpoint\ntransition")) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20),
                     expand = expansion(mult = c(0, 0))) +
  labs(
    title = "GTDB r226 ARD transition routes",
    x = NULL, y = "Share of inferred ARD transitions (%)",
    caption = "Fixed-Q stochastic-history posterior means; unit topology and source branch-length sensitivity."
  ) +
  theme_classic(base_size = 7.2, base_family = "Arial") +
  theme(
    axis.line = element_line(linewidth = 0.30, colour = "black"),
    axis.ticks = element_line(linewidth = 0.28, colour = "black"),
    axis.text.x = element_text(size = 6.2, lineheight = 0.92),
    axis.text.y = element_text(size = 6.2),
    axis.title.y = element_text(size = 7.0),
    plot.title = element_text(size = 8.1, face = "bold"),
    plot.caption = element_text(size = 5.3, colour = "#555555", hjust = 0,
                                margin = margin(t = 4)),
    legend.position = "top",
    legend.text = element_text(size = 6.0),
    legend.key.width = grid::unit(4.5, "mm"),
    plot.margin = margin(4, 6, 4, 5)
  )

ggsave(file.path(figure_dir,
  "fig_mreb_alp_ARD_transition_contribution_GTDB_unit_vs_source_branch_lengths.pdf"),
  p, width = 108, height = 92, units = "mm", device = cairo_pdf)
ggsave(file.path(figure_dir,
  "fig_mreb_alp_ARD_transition_contribution_GTDB_unit_vs_source_branch_lengths.png"),
  p, width = 108, height = 92, units = "mm", dpi = 600, bg = "white")

writeLines(c("status=complete", "backend=R", "tree=GTDB_r226",
             "branch_scales=unit_topology,source_branch_lengths", "figures=1"),
  file.path(figure_dir, "GTDB_ARD_UNIT_SOURCE_GROUPED_COMPLETE.marker"))
