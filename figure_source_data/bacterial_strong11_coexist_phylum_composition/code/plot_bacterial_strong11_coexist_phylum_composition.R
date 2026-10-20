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
phylum_counts <- read_src("formal11_strong_coexist_phylum_counts.tsv")
phylum_totals <- phylum_counts %>% group_by(phylum) %>% summarise(total = sum(n_hosts), .groups = "drop") %>% arrange(desc(total), phylum) %>%
  mutate(abundance_rank = row_number(), phylum_plot = case_when(abundance_rank <= 10 ~ phylum, abundance_rank <= 20 ~ "Other phyla (ranks 11-20)", abundance_rank <= 40 ~ "Other phyla (ranks 21-40)", TRUE ~ "Other phyla (rank >40)"))
kept_phyla <- head(phylum_totals$phylum, 10)
tail_groups <- c("Other phyla (ranks 11-20)", "Other phyla (ranks 21-40)", "Other phyla (rank >40)")
phylum_plot <- phylum_counts %>% left_join(phylum_totals %>% select(phylum, phylum_plot), by = "phylum") %>% group_by(og_id, phylum_plot) %>% summarise(n_hosts = sum(n_hosts), .groups = "drop") %>%
  mutate(og_id = factor(og_id, levels = strong_order_x), phylum_plot = factor(phylum_plot, levels = c(kept_phyla, tail_groups)))
phyla <- c(kept_phyla, tail_groups)
phy_cols <- c(setNames(hcl.colors(length(kept_phyla), "Set 3"), kept_phyla), "Other phyla (ranks 11-20)" = "#B8A58A", "Other phyla (ranks 21-40)" = "#AEB8C2", "Other phyla (rank >40)" = "#D4D7DA")
p <- ggplot(phylum_plot, aes(og_id, n_hosts, fill = phylum_plot)) + geom_col(width = 0.70, colour = "white", linewidth = 0.12) +
  scale_fill_manual(values = phy_cols, breaks = phyla, name = "Phylum", guide = guide_legend(nrow = 2, byrow = TRUE)) +
  scale_x_discrete(labels = axis_labels, drop = FALSE) + scale_y_continuous(expand = expansion(mult = c(0, 0.04))) + labs(x = NULL, y = "Coexist genomes") + theme_pub +
  theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 6), legend.position = "bottom", legend.text = element_text(size = 5.7), legend.key.size = grid::unit(3.8, "mm"), plot.margin = margin(5, 6, 28, 7))
save_pair(p, "bacterial_strong11_coexist_phylum_composition", 300, 150)

