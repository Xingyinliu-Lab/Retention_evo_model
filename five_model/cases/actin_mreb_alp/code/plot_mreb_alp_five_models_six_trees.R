suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
})

all_args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", all_args[grep("^--file=", all_args)][1])
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
input <- file.path(root, "results", "MreB_ALP_composition_corrected_five_models_six_trees.tsv")
p_file <- file.path(root, "results", "paired_tests_vs_sequential_retention.tsv")
figdir <- file.path(root, "figures")
srcdir <- file.path(root, "source_data")
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
dir.create(srcdir, recursive = TRUE, showWarnings = FALSE)

model_levels <- c("ARD unrestricted", "Sequential retention", "ALP-middle adjacency", "M-middle adjacency", "Ancestral coexist (strict)")
tree_levels <- c("GTDB r226", "A1338 add-core", "A1338 de novo", "A1338 ASTRAL64", "Asgard971 add-core", "Asgard971 ASTRAL89")

dat <- suppressMessages(read_tsv(input)) %>%
  mutate(
    model_label = recode(model,
      ARD_unrestricted = "ARD unrestricted",
      adjacency_coexist_middle = "Sequential retention",
      adjacency_ALP_middle = "ALP-middle adjacency",
      adjacency_M_middle = "M-middle adjacency",
      ancestral_coexistence_strict = "Ancestral coexist (strict)"),
    tree_label = recode(tree,
      GTDB_r226 = "GTDB r226", A1338_addcore = "A1338 add-core",
      A1338_denovo = "A1338 de novo", A1338_ASTRAL64 = "A1338 ASTRAL64",
      Asgard971_addcore = "Asgard971 add-core", Asgard971_ASTRAL89 = "Asgard971 ASTRAL89"),
    model_label = factor(model_label, levels = model_levels),
    tree_label = factor(tree_label, levels = tree_levels),
    y_plot = log10(delta_AIC + 1)
  )

stopifnot(nrow(dat) == 30, all(table(dat$model_label) == 6), all(table(dat$tree_label) == 5), all(dat$success))
write_tsv(dat, file.path(srcdir, "MreB_ALP_composition_corrected_five_models_six_trees_figure_data.tsv"))

pdat <- suppressMessages(read_tsv(p_file)) %>%
  mutate(model_label = recode(comparison_model,
    ARD_unrestricted = "ARD unrestricted",
    adjacency_ALP_middle = "ALP-middle adjacency",
    adjacency_M_middle = "M-middle adjacency",
    ancestral_coexistence_strict = "Ancestral coexist (strict)"),
    model_label = factor(model_label, levels=model_levels),
    label=paste0("P = ", format(exact_p_two_sided, digits=3)),
    y=max(dat$y_plot)*1.045)

tree_cols <- c("GTDB r226"="#333333", "A1338 add-core"="#D28E35", "A1338 de novo"="#B95B66",
               "A1338 ASTRAL64"="#8066A6", "Asgard971 add-core"="#2D8C8C", "Asgard971 ASTRAL89"="#4E78B4")
tree_shapes <- c("GTDB r226"=16, "A1338 add-core"=17, "A1338 de novo"=15,
                 "A1338 ASTRAL64"=18, "Asgard971 add-core"=8, "Asgard971 ASTRAL89"=3)
breaks_raw <- c(0, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000)

p <- ggplot(dat, aes(model_label, y_plot, group=model_label)) +
  geom_hline(yintercept=0, size=0.28, colour="#777777", linetype="22") +
  geom_boxplot(width=0.43, size=0.34, outlier.shape=NA, fill=NA, colour="#8B8B8B") +
  geom_jitter(aes(colour=tree_label, shape=tree_label), width=0.10, height=0,
              size=2.15, stroke=0.45, alpha=0.98) +
  geom_text(data=pdat, aes(x=model_label, y=y, label=label), inherit.aes=FALSE,
            size=2.05, family="Arial", colour="#333333") +
  scale_y_continuous(breaks=log10(breaks_raw+1), labels=breaks_raw, minor_breaks=NULL,
                     expand=expansion(mult=c(0.02,0.13))) +
  scale_colour_manual(values=tree_cols, drop=FALSE, name=NULL) +
  scale_shape_manual(values=tree_shapes, drop=FALSE, name=NULL) +
  labs(x=NULL, y=expression(Delta*"AIC"),
       title="MreB–ALP aggregate state models",
       caption="Composition-corrected root priors; unit topology. Y-axis spacing is logarithmic; labels show raw values.") +
  theme_classic(base_size=7.2, base_family="Arial") +
  theme(axis.line=element_line(size=0.30), axis.ticks=element_line(size=0.28),
        axis.text.x=element_text(size=6.1, lineheight=0.92), axis.text.y=element_text(size=6.3),
        plot.title=element_text(size=8.2, face="bold"), plot.caption=element_text(size=5.3, colour="#555555"),
        legend.position="top", legend.text=element_text(size=5.8), legend.key.width=grid::unit(3.8,"mm"),
        plot.margin=margin(3,6,2,5))

ggsave(file.path(figdir, "fig_mreb_alp_five_model_AIC_six_trees.pdf"), p,
       width=178, height=106, units="mm", device=cairo_pdf)
ggsave(file.path(figdir, "fig_mreb_alp_five_model_AIC_six_trees.png"), p,
       width=178, height=106, units="mm", dpi=600, bg="white")
