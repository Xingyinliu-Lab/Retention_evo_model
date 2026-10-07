suppressPackageStartupMessages({library(dplyr); library(readr)})

args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_file <- normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = TRUE)
package_root <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/", mustWork = TRUE)
input <- if (length(args) >= 1) args[[1]] else file.path(package_root, "output", "five_model_fits.tsv")
output <- if (length(args) >= 2) args[[2]] else file.path(package_root, "output", "paired_tests_vs_sequential_retention.tsv")

benchmark <- "adjacency_coexist_middle"
labels <- c(
  adjacency_coexist_middle = "Sequential retention (coexist middle)",
  adjacency_dALP_middle = "dALP-middle adjacency",
  adjacency_cALP_middle = "cALP-middle adjacency",
  ancestral_coexist_divergence = "Ancestral coexist divergence",
  ARD_unrestricted = "ARD unrestricted"
)

exact_signed_rank <- function(diff, alternative = c("two.sided", "greater")) {
  alternative <- match.arg(alternative)
  diff <- diff[is.finite(diff) & abs(diff) > 1e-12]
  n <- length(diff)
  if (n == 0) return(1)
  ranks <- rank(abs(diff), ties.method = "average")
  observed <- sum(sign(diff) * ranks)
  signs <- as.matrix(expand.grid(rep(list(c(-1, 1)), n)))
  null <- as.numeric(signs %*% ranks)
  if (alternative == "two.sided") mean(abs(null) >= abs(observed) - 1e-12)
  else mean(null >= observed - 1e-12)
}

d <- read_tsv(input, col_types = cols())
rows <- list(); ii <- 1
for (scale in unique(d$branch_scale)) {
  b <- d %>% filter(branch_scale == scale, model == benchmark) %>%
    select(tree, benchmark_AIC = AIC)
  for (model_id in setdiff(unique(d$model), benchmark)) {
    z <- d %>% filter(branch_scale == scale, model == model_id) %>%
      select(tree, comparison_AIC = AIC) %>% inner_join(b, by = "tree") %>%
      mutate(delta = comparison_AIC - benchmark_AIC)
    rows[[ii]] <- tibble(
      branch_scale = scale,
      benchmark = unname(labels[benchmark]),
      comparison_model = unname(labels[model_id]),
      n_paired_trees = nrow(z),
      benchmark_better_trees = sum(z$delta > 1e-12),
      comparison_better_trees = sum(z$delta < -1e-12),
      ties = sum(abs(z$delta) <= 1e-12),
      median_AIC_difference_other_minus_benchmark = median(z$delta),
      exact_p_two_sided = exact_signed_rank(z$delta, "two.sided"),
      exact_p_benchmark_better_one_sided = exact_signed_rank(z$delta, "greater")
    )
    ii <- ii + 1
  }
}
out <- bind_rows(rows) %>% arrange(branch_scale, exact_p_two_sided, comparison_model)
write_tsv(out, output)
print(out, n = Inf)
