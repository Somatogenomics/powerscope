#' powerscope: power, sample size, and Type S/M errors for genomics designs
#'
#' Every design function reduces a study to an effect and its standard error
#' and passes them to [ps_retrodesign()], which returns power, the sign error
#' rate (Type S) and the exaggeration ratio (Type M). [ps_sample_size()] and
#' [ps_mde()] solve for sample size and minimum detectable effect;
#' [ps_fdr_alpha()] converts an FDR target into a per-test level.
#'
#' @import data.table
#' @keywords internal
"_PACKAGE"

utils::globalVariables(c(
  ".", "alpha", "beta_sd", "case_freq", "delta", "depth1", "depth_factor", "dispersion",
  "exaggeration", "maf", "mean_count", "method", "min_sig_ratio", "n1", "n_pools1",
  "pool_size1", "power", "r2", "se", "sig_cutoff", "type_s", "n_cases",
  "grp", "value", "variable", "y", ".data"
))
