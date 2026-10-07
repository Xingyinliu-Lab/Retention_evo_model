suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
})

args <- commandArgs(FALSE)
script <- sub("^--file=", "", args[grep("^--file=", args)][1])
root <- normalizePath(file.path(dirname(script), ".."), mustWork=TRUE)
result_dir <- file.path(root, "results", "source_branch_lengths_20260820")
figure_dir <- file.path(root, "figures", "source_branch_lengths_20260821")
source_dir <- file.path(root, "source_data", "source_branch_lengths_20260821")
dir.create(figure_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(source_dir, recursive=TRUE, showWarnings=FALSE)

tree_levels <- c("GTDB r226", "A1338 add-core", "A1338 de novo", "A1338 ASTRAL64",
                 "Asgard971 add-core", "Asgard971 ASTRAL89")
tree_cols <- c("GTDB r226"="#333333", "A1338 add-core"="#D28E35",
               "A1338 de novo"="#B95B66", "A1338 ASTRAL64"="#8066A6",
               "Asgard971 add-core"="#2D8C8C", "Asgard971 ASTRAL89"="#4E78B4")
tree_shapes <- c("GTDB r226"=16, "A1338 add-core"=17, "A1338 de novo"=15,
                 "A1338 ASTRAL64"=18, "Asgard971 add-core"=8, "Asgard971 ASTRAL89"=3)
tree_labels <- c(GTDB_r226="GTDB r226", A1338_addcore="A1338 add-core",
                 A1338_denovo="A1338 de novo", A1338_ASTRAL64="A1338 ASTRAL64",
                 Asgard971_addcore="Asgard971 add-core", Asgard971_ASTRAL89="Asgard971 ASTRAL89")

theme_pub <- function() {
  theme_classic(base_size=7.2, base_family="Arial") +
    theme(axis.line=element_line(linewidth=0.30), axis.ticks=element_line(linewidth=0.28),
          axis.text.x=element_text(size=6.1, lineheight=0.92), axis.text.y=element_text(size=6.3),
          plot.title=element_text(size=8.2, face="bold"),
          plot.caption=element_text(size=5.3, colour="#555555"),
          legend.position="top", legend.text=element_text(size=5.8),
          legend.key.width=grid::unit(3.8, "mm"), plot.margin=margin(3,6,3,5))
}

# Five-model comparison: source branch lengths
fit <- suppressMessages(read_tsv(file.path(result_dir, "MreB_ALP_five_models_source_branch_lengths.tsv")))
ptest <- suppressMessages(read_tsv(file.path(result_dir, "paired_tests_vs_sequential_retention_source_branch_lengths.tsv")))
model_levels <- c("ARD unrestricted", "Sequential retention", "ALP-middle adjacency",
                  "M-middle adjacency", "Ancestral coexist (strict)")
fit_plot <- fit %>% mutate(
  model_label=recode(model,
    ARD_unrestricted="ARD unrestricted", adjacency_coexist_middle="Sequential retention",
    adjacency_ALP_middle="ALP-middle adjacency", adjacency_M_middle="M-middle adjacency",
    ancestral_coexistence_strict="Ancestral coexist (strict)"),
  tree_label=recode(tree, !!!tree_labels),
  model_label=factor(model_label, levels=model_levels),
  tree_label=factor(tree_label, levels=tree_levels),
  y_plot=log10(delta_AIC + 1)
)
stopifnot(nrow(fit_plot)==30, all(fit_plot$success), all(fit_plot$branch_scale=="source_branch_lengths"))
write_tsv(fit_plot, file.path(source_dir, "MreB_ALP_five_model_AIC_source_branch_lengths_figure_data.tsv"))

p_label <- ptest %>% mutate(
  model_label=recode(comparison_model,
    ARD_unrestricted="ARD unrestricted", adjacency_ALP_middle="ALP-middle adjacency",
    adjacency_M_middle="M-middle adjacency", ancestral_coexistence_strict="Ancestral coexist (strict)"),
  model_label=factor(model_label, levels=model_levels),
  label=paste0("P = ", format(exact_p_two_sided, digits=3)),
  y=max(fit_plot$y_plot)*1.045
)
raw_breaks <- c(0,1,2,5,10,20,50,100,200,500,1000)
p_five <- ggplot(fit_plot, aes(model_label, y_plot, group=model_label)) +
  geom_hline(yintercept=0, linewidth=0.28, colour="#777777", linetype="22") +
  geom_boxplot(width=0.43, linewidth=0.34, outlier.shape=NA, fill=NA, colour="#8B8B8B") +
  geom_jitter(aes(colour=tree_label, shape=tree_label), width=0.10, height=0,
              size=2.15, stroke=0.45, alpha=0.98) +
  geom_text(data=p_label, aes(x=model_label, y=y, label=label), inherit.aes=FALSE,
            size=2.05, family="Arial", colour="#333333") +
  scale_y_continuous(breaks=log10(raw_breaks+1), labels=raw_breaks, minor_breaks=NULL,
                     expand=expansion(mult=c(0.02,0.13))) +
  scale_colour_manual(values=tree_cols, drop=FALSE, name=NULL) +
  scale_shape_manual(values=tree_shapes, drop=FALSE, name=NULL) +
  labs(x=NULL, y=expression(Delta*"AIC"), title="MreB–ALP aggregate state models",
       caption="Composition-corrected root priors; source species-tree branch lengths. Y-axis spacing is logarithmic; labels show raw values.") +
  theme_pub()

ggsave(file.path(figure_dir, "fig_mreb_alp_five_model_AIC_six_trees_source_branch_lengths.pdf"),
       p_five, width=178, height=106, units="mm", device=cairo_pdf)
ggsave(file.path(figure_dir, "fig_mreb_alp_five_model_AIC_six_trees_source_branch_lengths.png"),
       p_five, width=178, height=106, units="mm", dpi=600, bg="white")

# ARD route contributions: source branch lengths
ard <- suppressMessages(read_tsv(file.path(result_dir, "MreB_ALP_ARD_transition_contribution_source_branch_lengths.tsv")))
ard_plot <- bind_rows(
  ard %>% transmute(tree, route="Coexist-mediated", percent=100*coexist_mediated_fraction),
  ard %>% transmute(tree, route="Direct endpoint transition", percent=100*direct_endpoint_fraction)
) %>% mutate(
  tree_label=factor(recode(tree, !!!tree_labels), levels=tree_levels),
  route=factor(route, levels=c("Coexist-mediated", "Direct endpoint transition"))
)
stopifnot(nrow(ard_plot)==12, all(is.finite(ard_plot$percent)))
write_tsv(ard_plot, file.path(source_dir, "MreB_ALP_ARD_transition_contribution_source_branch_lengths_figure_data.tsv"))

p_ard <- ggplot(ard_plot, aes(route, percent, group=route)) +
  geom_boxplot(width=0.42, linewidth=0.34, outlier.shape=NA, fill=NA, colour="#8B8B8B") +
  geom_jitter(aes(colour=tree_label, shape=tree_label), width=0.08, height=0,
              size=2.15, stroke=0.45) +
  scale_colour_manual(values=tree_cols, name=NULL) +
  scale_shape_manual(values=tree_shapes, name=NULL) +
  scale_y_continuous(limits=c(0,100), breaks=seq(0,100,20), expand=expansion(mult=c(0.01,0.03))) +
  scale_x_discrete(labels=c("Coexist-mediated\ntransitions", "Direct endpoint\ntransitions")) +
  labs(x=NULL, y="Share of inferred ARD transitions (%)",
       title="MreB–ALP ARD transition-route contributions",
       caption="Fixed-Q stochastic-history posterior means; six species-tree reconstructions; source branch lengths.") +
  theme_pub()

ggsave(file.path(figure_dir, "fig_mreb_alp_ARD_transition_contribution_six_trees_source_branch_lengths.pdf"),
       p_ard, width=142, height=96, units="mm", device=cairo_pdf)
ggsave(file.path(figure_dir, "fig_mreb_alp_ARD_transition_contribution_six_trees_source_branch_lengths.png"),
       p_ard, width=142, height=96, units="mm", dpi=600, bg="white")

writeLines(c("status=complete", "backend=R", "branch_scale=source_branch_lengths", "figures=2"),
           file.path(figure_dir, "FIGURES_COMPLETE.marker"))
