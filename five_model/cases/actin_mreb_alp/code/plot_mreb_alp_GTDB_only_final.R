suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = TRUE)
root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
figure_dir <- file.path(root, "figures")
source_dir <- file.path(root, "source_data")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

theme_bar <- function() {
  theme_classic(base_size = 7.2, base_family = "Arial") +
    theme(
      axis.line = element_line(size = 0.30, colour = "black"),
      axis.ticks = element_line(size = 0.28, colour = "black"),
      axis.text.x = element_text(size = 6.0, lineheight = 0.92),
      axis.text.y = element_text(size = 6.2),
      axis.title.y = element_text(size = 7.0),
      plot.title = element_text(size = 8.1, face = "bold"),
      plot.caption = element_text(size = 5.3, colour = "#555555", hjust = 0,
                                  margin = margin(t = 4)),
      plot.margin = margin(4, 6, 4, 5),
      legend.position = "none"
    )
}

# Final composition-corrected five-model result, GTDB r226 only.
model_levels <- c(
  "ARD unrestricted",
  "Sequential\nretention",
  "ALP-middle\nadjacency",
  "M-middle\nadjacency",
  "Ancestral coexist\n(strict)"
)
model_dat <- read.delim(
  file.path(root, "results", "MreB_ALP_composition_corrected_five_models_six_trees.tsv"),
  check.names = FALSE, stringsAsFactors = FALSE
) %>%
  filter(tree == "GTDB_r226") %>%
  mutate(
    model_label = recode(
      model,
      ARD_unrestricted = "ARD unrestricted",
      adjacency_coexist_middle = "Sequential\nretention",
      adjacency_ALP_middle = "ALP-middle\nadjacency",
      adjacency_M_middle = "M-middle\nadjacency",
      ancestral_coexistence_strict = "Ancestral coexist\n(strict)"
    ),
    model_label = factor(model_label, levels = model_levels),
    value_label = sprintf("%.1f", delta_AIC)
  ) %>%
  arrange(model_label)

stopifnot(nrow(model_dat) == 5, all(model_dat$branch_scale == "unit_topology"))
write.table(
  model_dat,
  file.path(source_dir, "MreB_ALP_five_model_AIC_GTDB_only.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE, fileEncoding = "UTF-8"
)

model_cols <- c(
  "ARD unrestricted" = "#4D4D4D",
  "Sequential\nretention" = "#4F857D",
  "ALP-middle\nadjacency" = "#9A83B5",
  "M-middle\nadjacency" = "#D79A45",
  "Ancestral coexist\n(strict)" = "#B86168"
)
ymax <- max(model_dat$delta_AIC) * 1.12
p_aic <- ggplot(model_dat, aes(model_label, delta_AIC, fill = model_label)) +
  geom_col(width = 0.62, colour = NA) +
  geom_text(aes(label = value_label), vjust = -0.45, size = 2.20,
            family = "Arial", colour = "#303030") +
  scale_fill_manual(values = model_cols) +
  scale_y_continuous(
    limits = c(0, ymax), breaks = seq(0, 500, 100),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = "GTDB r226 model comparison",
    x = NULL,
    y = expression(Delta*"AIC"),
    caption = "Composition-corrected root priors; GTDB r226 expanded tree, unit topology."
  ) +
  theme_bar()

ggsave(file.path(figure_dir, "fig_mreb_alp_five_model_AIC_GTDB_only.pdf"),
       p_aic, width = 125, height = 88, units = "mm", device = cairo_pdf)
ggsave(file.path(figure_dir, "fig_mreb_alp_five_model_AIC_GTDB_only.png"),
       p_aic, width = 125, height = 88, units = "mm", dpi = 600, bg = "white")

# Final fixed-Q ARD route contribution, GTDB r226 only.
route_dat <- read.delim(
  file.path(root, "results", "MreB_ALP_ARD_transition_contribution_summary_six_trees.tsv"),
  check.names = FALSE, stringsAsFactors = FALSE
) %>%
  filter(tree == "GTDB_r226") %>%
  transmute(
    `Coexist-mediated` = 100 * coexist_mediated_fraction,
    `Direct endpoint\ntransition` = 100 * direct_endpoint_fraction
  ) %>%
  pivot_longer(cols = everything(), names_to = "route_label", values_to = "contribution_percent") %>%
  mutate(
    route_label = factor(route_label,
                         levels = c("Coexist-mediated", "Direct endpoint\ntransition")),
    value_label = sprintf("%.1f%%", contribution_percent)
  )

stopifnot(nrow(route_dat) == 2, abs(sum(route_dat$contribution_percent) - 100) < 1e-8)
write.table(
  route_dat,
  file.path(source_dir, "MreB_ALP_ARD_transition_contribution_GTDB_only.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE, fileEncoding = "UTF-8"
)

route_cols <- c(
  "Coexist-mediated" = "#4F857D",
  "Direct endpoint\ntransition" = "#B86168"
)
p_route <- ggplot(route_dat, aes(route_label, contribution_percent, fill = route_label)) +
  geom_col(width = 0.58, colour = NA) +
  geom_text(aes(label = value_label), vjust = -0.45, size = 2.35,
            family = "Arial", colour = "#303030") +
  scale_fill_manual(values = route_cols) +
  scale_y_continuous(
    limits = c(0, 100), breaks = seq(0, 100, 20),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = "GTDB r226 ARD transition routes",
    x = NULL,
    y = "Share of inferred ARD transitions (%)",
    caption = "Fixed-Q stochastic-history posterior means; unit topology."
  ) +
  theme_bar()

ggsave(file.path(figure_dir, "fig_mreb_alp_ARD_transition_contribution_GTDB_only.pdf"),
       p_route, width = 96, height = 88, units = "mm", device = cairo_pdf)
ggsave(file.path(figure_dir, "fig_mreb_alp_ARD_transition_contribution_GTDB_only.png"),
       p_route, width = 96, height = 88, units = "mm", dpi = 600, bg = "white")
