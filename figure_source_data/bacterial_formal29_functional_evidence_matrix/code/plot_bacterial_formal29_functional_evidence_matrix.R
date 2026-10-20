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

cases <- read_src("formal29_case_summary.tsv") %>% arrange(factor(evidence_class, levels = c("strong", "supported")), desc(n_AB_coexist_MAG), og_id) %>% mutate(row = rev(row_number()), label = functional_label)
category_palette <- c("Metabolism" = "#75A082", "Energy/redox" = "#C4884B", "Membrane/transport" = "#628EAE", "Cell envelope" = "#6D9B9B", "Defense/regulation" = "#9B7DAE", "Stress/regulation" = "#AA8C70", "Genome maintenance" = "#8573A4", "Protein processing" = "#8094A1", "Signal regulation" = "#8D85B5", "Uncharacterized" = "#A4A4A4")
p <- ggplot(cases) + geom_rect(data = cases %>% filter(row %% 2 == 1), aes(xmin = 0, xmax = 100, ymin = row - 0.45, ymax = row + 0.45), fill = "#ECECEC", colour = NA, alpha = 0.55) +
  geom_text(aes(x = 1, y = row, label = og_id), hjust = 0, size = 2.45, family = "Arial", colour = "#3C4248") +
  geom_text(aes(x = 17, y = row, label = label), hjust = 0, size = 2.35, family = "Arial", colour = "#3C4248") +
  geom_rect(aes(xmin = 59, xmax = 77, ymin = row - 0.30, ymax = row + 0.30, fill = broad_category), colour = NA) +
  geom_text(aes(x = 68, y = row, label = broad_category), size = 2.1, family = "Arial", fontface = "bold", colour = "white") +
  geom_point(aes(x = 83, y = row, size = n_coexist_phyla), shape = 21, fill = "#3B86A6", colour = "white", stroke = 0.3, alpha = 0.9) +
  geom_text(aes(x = 83, y = row, label = n_coexist_phyla), size = 2.0, family = "Arial", fontface = "bold", colour = "white") +
  geom_rect(aes(xmin = 88, xmax = 99, ymin = row - 0.30, ymax = row + 0.30, fill = evidence_class), colour = NA) +
  geom_text(aes(x = 93.5, y = row, label = str_to_title(evidence_class)), size = 2.1, family = "Arial", fontface = "bold", colour = "white") +
  annotate("text", x = c(1, 17, 68, 83, 93.5), y = max(cases$row) + 1, label = c("OG", "Representative annotation", "Broad category", "Coexist\nphyla", "Gate class"), hjust = c(0, 0, 0.5, 0.5, 0.5), family = "Arial", fontface = "bold", size = 2.7) +
  scale_fill_manual(values = c(category_palette, palette_gate), guide = "none") + scale_size_continuous(range = c(3.0, 6.2), guide = "none") +
  coord_cartesian(xlim = c(0, 100), ylim = c(0.4, max(cases$row) + 1.7), clip = "off") + theme_void(base_family = "Arial") + theme(plot.margin = margin(6, 6, 6, 6))
save_pair(p, "bacterial_formal29_functional_evidence_matrix", 250, 205)

