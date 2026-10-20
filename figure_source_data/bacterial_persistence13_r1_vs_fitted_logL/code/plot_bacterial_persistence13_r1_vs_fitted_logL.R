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

burdens <- read_src("formal13_persistence_supported_burdens.tsv") %>% mutate(
  treatment = factor(treatment, levels = c("r=1", "fitted r")),
  support_group = if_else(persistence_status == "COEXIST_PERSISTENCE_SUPPORTED", "Persistence supported", "Null boundary"))
support_cols <- c("Persistence supported" = "#5D9CA7", "Null boundary" = "#DF945B")
paired_plot <- function(value_col, ylab, stem, title = NULL) {
  p <- ggplot(burdens, aes(x = treatment, y = .data[[value_col]], group = og_id, colour = support_group)) +
    geom_line(alpha = 0.58, linewidth = 0.4) +
    geom_boxplot(inherit.aes = FALSE, data = burdens, mapping = aes(x = treatment, y = .data[[value_col]], fill = treatment), width = 0.46, outlier.shape = NA, colour = "#2D2D2D", linewidth = 0.45) +
    geom_point(size = 1.85, shape = 21, fill = "white", stroke = 0.65, position = position_jitter(width = 0.035, height = 0, seed = 20261019)) +
    scale_colour_manual(values = support_cols, name = NULL) + scale_fill_manual(values = c("r=1" = "#C8D5E5", "fitted r" = "#EAC0A4"), guide = "none") +
    labs(title = title, x = NULL, y = ylab) + theme_pub + theme(legend.position = "top")
  save_pair(p, stem, 132, 126)
}

paired_plot("log_likelihood", "Raw log likelihood", "bacterial_persistence13_r1_vs_fitted_logL", "logL: r=1 versus fitted-r model")


