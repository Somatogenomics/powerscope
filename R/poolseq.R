#' Pooled-sequencing allele-frequency difference between two groups
#'
#' Each pool contains `pool_size` individuals of the given `ploidy`, so
#' `k = ploidy * pool_size` chromosomes, and is sequenced to read depth
#' `depth` at the site. Two-stage binomial sampling (chromosomes into the
#' pool, then reads from the pool) gives
#' `Var(p_hat) = p(1-p) * (1/k + 1/C - 1/(k*C))` for depth `C`
#' (e.g. Gautier et al. 2013). With `n_pools` independent replicate pools per
#' group, frequencies are averaged across pools. `overdispersion` (>= 1)
#' multiplies the total variance to account for unequal DNA contribution of
#' individuals, library and mapping effects; estimate it from the variance
#' between replicate pools relative to this binomial expectation.
#'
#' `target = "population"` tests a difference between the populations the
#' pools were drawn from. `target = "pool"` treats the pools' own chromosome
#' frequencies as the quantity of interest (read sampling only), e.g. when the
#' pools are the complete set of individuals of interest.
#'
#' @param p1,p2 True allele frequencies in groups 1 and 2.
#' @param pool_size1 Individuals per pool in group 1.
#' @param pool_size2 Individuals per pool in group 2; defaults to `pool_size1`.
#' @param depth1 Read depth per pool at the site, group 1.
#' @param depth2 Read depth per pool, group 2; defaults to `depth1`.
#' @param n_pools1,n_pools2 Replicate pools per group.
#' @param ploidy Chromosome sets per individual (2 for diploid; 1 for haploid
#'   or hemizygous loci).
#' @param overdispersion Variance inflation factor (>= 1).
#' @param alpha Per-site two-sided significance level.
#' @param target `"population"` or `"pool"`.
#' @return A `ps_result` data.table; estimate scale is the allele-frequency
#'   difference `p2 - p1`.
#' @references Gautier M, Foucaud J, Gharbi K, et al. (2013). Estimation of
#'   population allele frequencies from next-generation sequencing data:
#'   pool-versus individual-based genotyping. *Molecular Ecology*
#'   22:3766-3779. \doi{10.1111/mec.12360}
#' @examples
#' ps_poolseq(p1 = 0.5, p2 = 0.6, pool_size1 = 50, depth1 = c(30, 60, 100))
#' @export
ps_poolseq <- function(p1, p2, pool_size1, pool_size2 = NULL, depth1, depth2 = NULL,
                       n_pools1 = 1, n_pools2 = NULL, ploidy = 2, overdispersion = 1,
                       alpha = 0.05, target = c("population", "pool")) {
  target <- match.arg(target)
  .check_prob(p1, "p1"); .check_prob(p2, "p2")
  .check_pos(pool_size1, "pool_size1"); .check_pos(depth1, "depth1")
  .check_pos(n_pools1, "n_pools1"); .check_pos(ploidy, "ploidy")
  if (any(overdispersion < 1)) stop("`overdispersion` must be >= 1")
  g <- .ps_grid(list(p1 = p1, p2 = p2, pool_size1 = pool_size1, pool_size2 = pool_size2,
                     depth1 = depth1, depth2 = depth2, n_pools1 = n_pools1, n_pools2 = n_pools2,
                     ploidy = ploidy, overdispersion = overdispersion, alpha = alpha))
  if (is.null(pool_size2)) g[, pool_size2 := pool_size1]
  if (is.null(depth2)) g[, depth2 := depth1]
  if (is.null(n_pools2)) g[, n_pools2 := n_pools1]
  if (any(g$p1 == g$p2)) stop("`p1` and `p2` must differ")
  v <- function(p, N, C, k) {
    kN <- k * N
    if (target == "population") p * (1 - p) * (1 / kN + 1 / C - 1 / (kN * C)) else p * (1 - p) / C
  }
  var_tot <- g$overdispersion * (v(g$p1, g$pool_size1, g$depth1, g$ploidy) / g$n_pools1 +
                                 v(g$p2, g$pool_size2, g$depth2, g$ploidy) / g$n_pools2)
  out <- .ps_finish(g, paste0("poolseq_", target), "allele-frequency difference",
                    g$p2 - g$p1, sqrt(var_tot), Inf, g$alpha)
  out[, method := "two-stage binomial, Wald z"]
  out[]
}
