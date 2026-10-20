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
strong_cases <- cases %>% filter(evidence_class == "strong")
strong_order_x <- strong_cases$og_id
axis_labels <- setNames(paste0(strong_cases$og_id, "\n", str_wrap(strong_cases$functional_label, width = 24)), strong_cases$og_id)
state_counts <- read_src("formal11_strong_host_state_counts.tsv") %>% mutate(og_id = factor(og_id, levels = strong_order_x), state = factor(state, levels = c("A-only", "A+B coexist", "B-only")))
p <- ggplot(state_counts, aes(og_id, n_hosts, fill = state)) + geom_col(width = 0.70, colour = "white", linewidth = 0.18) +
  scale_fill_manual(values = palette_state, name = NULL, drop = FALSE) + scale_x_discrete(labels = axis_labels, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.035))) + labs(x = NULL, y = "Genome count") + theme_pub +
  theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 6), legend.position = "top", plot.margin = margin(5, 6, 28, 7))
save_pair(p, "bacterial_strong11_host_state_counts", 300, 125)

