suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(stringr)
  library(patchwork)
})
args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
code_dir <- dirname(normalizePath(script_arg, winslash = "/"))
root <- dirname(code_dir)
data_dir <- file.path(root, "data")
fig_dir <- file.path(root, "fig")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
read_src <- function(name) read_tsv(file.path(data_dir, name), show_col_types = FALSE)
save_pair <- function(plot, stem, width_mm, height_mm) {
  ggsave(file.path(fig_dir, paste0(stem, ".pdf")), plot, width = width_mm, height = height_mm, units = "mm", device = cairo_pdf, bg = "white")
  ggsave(file.path(fig_dir, paste0(stem, ".png")), plot, width = width_mm, height = height_mm, units = "mm", dpi = 600, bg = "white")
}
palette_gate <- c(strong = "#D45A56", supported = "#2F958B", compatible = "#D6A235", not_supported = "#B9BEC8")
palette_state <- c("A-only" = "#5F83A7", "A+B coexist" = "#70A397", "B-only" = "#C57D5D")
theme_pub <- theme_classic(base_family = "Arial", base_size = 8) + theme(
  axis.line = element_line(linewidth = 0.35, colour = "#202733"),
  axis.ticks = element_line(linewidth = 0.35, colour = "#202733"),
  axis.text = element_text(colour = "#28364A"),
  axis.title = element_text(colour = "#1C293C", face = "bold"),
  panel.grid = element_blank(), plot.title = element_text(face = "bold", size = 11, colour = "#111820"),
  plot.subtitle = element_text(size = 8.5, colour = "#4A4F55"), legend.title = element_text(face = "bold"))

cases <- read_src("formal29_case_summary.tsv") %>% arrange(factor(evidence_class, levels = c("strong", "supported")), desc(n_AB_coexist_MAG), og_id)
strong_order_y <- rev((cases %>% filter(evidence_class == "strong"))$og_id)
fits <- read_src("formal11_strong_five_model_fits.tsv") %>% mutate(
  og_id = factor(og_id, levels = strong_order_y),
  branch_scale = factor(branch_scale, levels = c("unit_topology", "source_branch_lengths"), labels = c("Unit topology", "Source branch lengths")),
  model = factor(model, levels = c("ARD", "Retention", "A-middle", "B-middle", "Ancestral")),
  delta_cap = pmin(delta_AIC, 20), best = delta_AIC < 1e-8)
ard <- read_src("formal11_strong_ard_contribution.tsv") %>% mutate(
  og_id = factor(og_id, levels = strong_order_y),
  branch_scale = factor(branch_scale, levels = c("unit_topology", "source_branch_lengths"), labels = c("Unit topology", "Source branch lengths")))
plot_fit_block <- function(scale_value, title_text, show_y = FALSE) {
  d <- fits %>% filter(branch_scale == scale_value)
  ggplot(d, aes(model, og_id, fill = delta_cap)) +
    geom_tile(colour = "white", linewidth = 0.65) +
    geom_point(data = d %>% filter(best), size = 1.25, colour = "#24323C") +
    scale_fill_gradientn(colours = c("#F5F4EF", "#C8DAD7", "#82AFB1", "#2A7680"), values = rescale(c(0, 2, 8, 20)), limits = c(0, 20), name = "Delta AIC from best") +
    labs(x = NULL, y = NULL, title = title_text) + theme_pub +
    theme(axis.text.x = element_text(angle = 48, hjust = 1), axis.text.y = if (show_y) element_text(size = 6) else element_blank(),
          axis.ticks.y = if (show_y) element_line(linewidth = 0.35) else element_blank(), axis.line.y = if (show_y) element_line(linewidth = 0.35) else element_blank(),
          plot.title = element_text(size = 9, hjust = 0.5), legend.position = "bottom")
}
plot_ard_block <- function(scale_value, title_text) {
  d <- ard %>% filter(branch_scale == scale_value)
  ggplot(d, aes("ARD path", og_id, fill = coexist_mediated_fraction)) +
    geom_tile(colour = "#D1C39E", linewidth = 0.4) +
    geom_text(aes(label = percent(coexist_mediated_fraction, accuracy = 1)), size = 2.1, fontface = "bold", colour = "#5E5744") +
    scale_fill_gradient(low = "#F5F2E9", high = "#D4C08A", guide = "none") + labs(x = NULL, y = NULL, title = title_text) + theme_pub +
    theme(axis.text.x = element_text(angle = 48, hjust = 1), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(), plot.title = element_text(size = 8, hjust = 0.5))
}
p <- plot_fit_block("Unit topology", "Unit topology", TRUE) + plot_ard_block("Unit topology", "Unit ARD path %") +
  plot_fit_block("Source branch lengths", "Source branch lengths") + plot_ard_block("Source branch lengths", "Source ARD path %") +
  plot_layout(widths = c(5, 1.15, 5, 1.15), guides = "collect") & theme(legend.position = "bottom")
save_pair(p, "bacterial_strong11_five_model_ARD_matrix", 300, 123)

