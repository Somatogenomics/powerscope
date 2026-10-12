#' RNA-seq differential expression power for one gene (negative binomial)
#'
#' Counts in group g follow NB with mean `mu_g` and dispersion `phi`
#' (variance `mu + phi * mu^2`). The log fold change is estimated from the
#' group means; its variance is
#' `(1/mu1 + phi)/n1 + (1/mu2 + phi)/n2` on the natural-log scale
#' (Hart et al. 2013), reported here on the log2 scale. The test is a Wald
#' z-test with known dispersion, so with few replicates (where dispersion is
#' estimated and shrunk, as in DESeq2/edgeR) realised power is somewhat lower.
#'
#' @param log2fc True log2 fold change (group 2 / group 1).
#' @param mean_count Expected normalised count in group 1 at the current
#'   sequencing depth. DESeq2 `baseMean` (average over all samples) is a close
#'   substitute for modest fold changes.
#' @param dispersion NB dispersion `phi` (squared biological coefficient of
#'   variation). Estimate from your data, e.g. `mcols(dds)$dispersion`.
#' @param n1 Replicates in group 1.
#' @param n2 Replicates in group 2; defaults to `n1`.
#' @param depth_factor Multiplier on sequencing depth relative to the depth
#'   at which `mean_count` was measured.
#' @param alpha Per-gene two-sided significance level. For FDR control across
#'   genes, obtain a per-gene level with [ps_rnaseq_genes()] or
#'   [ps_fdr_alpha()].
#' @return A `ps_result` data.table; estimate scale is log2 fold change.
#' @references Hart SN, Therneau TM, Zhang Y, Poland GA, Kocher JP (2013).
#'   Calculating sample size estimates for RNA sequencing data. *Journal of
#'   Computational Biology* 20:970-978. \doi{10.1089/cmb.2012.0283}
#' @examples
#' ps_rnaseq(log2fc = 1, mean_count = 50, dispersion = 0.1, n1 = 3:8)
#' @export
ps_rnaseq <- function(log2fc, mean_count, dispersion, n1, n2 = NULL,
                      depth_factor = 1, alpha = 0.05) {
  if (any(log2fc == 0)) stop("`log2fc` must be nonzero")
  .check_pos(mean_count, "mean_count"); .check_pos(dispersion, "dispersion", TRUE)
  .check_pos(n1, "n1"); .check_pos(depth_factor, "depth_factor")
  g <- .ps_grid(list(log2fc = log2fc, mean_count = mean_count, dispersion = dispersion,
                     n1 = n1, n2 = n2, depth_factor = depth_factor, alpha = alpha))
  if (is.null(n2)) g[, n2 := n1]
  mu1 <- g$mean_count * g$depth_factor
  mu2 <- mu1 * 2^g$log2fc
  se_ln <- sqrt((1 / mu1 + g$dispersion) / g$n1 + (1 / mu2 + g$dispersion) / g$n2)
  out <- .ps_finish(g, "rnaseq_nb", "log2 fold change", g$log2fc, se_ln / log(2), Inf, g$alpha)
  out[, method := "NB Wald, known dispersion (Hart et al. 2013)"]
  out[]
}

#' RNA-seq power across genes, optionally under FDR control
#'
#' Applies [ps_rnaseq()] to every gene in a table of means and dispersions
#' measured from your own data, so that power reflects the observed spread of
#' expression and dispersion. If `q` and `pi1` are given, the per-gene level
#' is the Jung (2005) solution for Benjamini-Hochberg FDR `q`, assuming the
#' truly DE genes are a random draw from `genes`.
#'
#' @param genes A data.frame with columns `mean_count` and `dispersion`
#'   (one row per gene); optionally `gene_id`.
#' @param log2fc Assumed true |log2 fold change| for DE genes (scalar).
#' @param n1,n2,depth_factor See [ps_rnaseq()]. Scalars.
#' @param alpha Per-gene level, used when `q` is `NULL`.
#' @param q Target FDR, or `NULL`.
#' @param pi1 Fraction of genes truly DE (required with `q`; must be
#'   estimated from data or prior studies).
#' @return A list: `per_gene` (data.table of per-gene power, Type S,
#'   exaggeration), `summary` (one-row data.table: genes, per-gene alpha,
#'   mean and median power, expected true and false discoveries when `q` is
#'   used, and power-weighted exaggeration ratio and Type S rate, i.e. their
#'   expected values among DE genes that reach significance).
#' @examples
#' set.seed(1)
#' genes <- data.frame(mean_count = rlnorm(2000, log(100), 1.5),
#'                     dispersion = rlnorm(2000, log(0.1), 0.6))
#' ps_rnaseq_genes(genes, log2fc = 1, n1 = 5, q = 0.05, pi1 = 0.1)$summary
#' @export
ps_rnaseq_genes <- function(genes, log2fc, n1, n2 = n1, depth_factor = 1,
                            alpha = 0.05, q = NULL, pi1 = NULL) {
  genes <- data.table::as.data.table(genes)
  need <- c("mean_count", "dispersion")
  if (!all(need %in% names(genes))) stop("`genes` needs columns: ", paste(need, collapse = ", "))
  genes <- genes[is.finite(mean_count) & mean_count > 0 & is.finite(dispersion) & dispersion >= 0]
  if (!nrow(genes)) stop("no genes with finite positive mean_count and dispersion")
  if (length(log2fc) != 1 || length(n1) != 1 || length(n2) != 1) stop("log2fc, n1, n2 must be scalars")

  mu1 <- genes$mean_count * depth_factor
  mu2 <- mu1 * 2^log2fc
  se <- sqrt((1 / mu1 + genes$dispersion) / n1 + (1 / mu2 + genes$dispersion) / n2) / log(2)
  avg_power <- function(a) {
    z <- stats::qnorm(1 - a / 2); lam <- log2fc / se
    mean(stats::pnorm(z - lam, lower.tail = FALSE) + stats::pnorm(-z - lam))
  }
  exp_true <- exp_false <- NA_real_
  if (!is.null(q)) {
    if (is.null(pi1)) stop("`pi1` is required with `q`")
    a <- ps_fdr_alpha(avg_power, q = q, m = nrow(genes), pi1 = pi1)
    if (is.na(a)) stop("no per-gene alpha achieves FDR ", q, " at this design")
    exp_true <- attr(a, "expected_true"); exp_false <- attr(a, "expected_false")
    alpha <- as.numeric(a)
  }
  eng <- ps_retrodesign(effect = log2fc, se = se, alpha = alpha)
  per_gene <- cbind(genes, eng[, .(se, power, type_s, exaggeration, min_sig_ratio)])
  summary <- data.table::data.table(
    genes = nrow(genes), log2fc = log2fc, n1 = n1, n2 = n2, depth_factor = depth_factor,
    alpha_per_gene = alpha, fdr = if (is.null(q)) NA_real_ else q,
    pi1 = if (is.null(pi1)) NA_real_ else pi1,
    mean_power = mean(per_gene$power), median_power = stats::median(per_gene$power),
    expected_true_discoveries = exp_true, expected_false_discoveries = exp_false,
    # expected |estimate|/|true| among DE genes that reach significance:
    # power-weighted because each gene contributes to the discovery set in
    # proportion to its power
    exaggeration_among_discoveries = stats::weighted.mean(per_gene$exaggeration, per_gene$power, na.rm = TRUE),
    type_s_among_discoveries = stats::weighted.mean(per_gene$type_s, per_gene$power, na.rm = TRUE)
  )
  list(per_gene = per_gene[], summary = summary[])
}
