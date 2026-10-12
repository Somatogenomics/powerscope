#' Power, Type S and Type M (exaggeration) for an effect and its standard error
#'
#' The shared engine behind every design in powerscope. The estimate is
#' modelled as `effect + se * T`, where `T` is standard normal (`df = Inf`)
#' or Student t with `df` degrees of freedom. A result is called significant
#' when `|estimate| >= q * se`, with `q` the two-sided critical value at
#' `alpha`.
#'
#' * `power`: P(significant).
#' * `type_s`: P(estimate has the opposite sign to `effect` | significant).
#' * `exaggeration`: E(|estimate| | significant) / |effect| (Type M).
#' * `min_sig_ratio`: `q * se / |effect|`, the smallest overestimation any
#'   significant result can have.
#'
#' For `df = Inf` all quantities are closed form. For finite `df` the
#' exaggeration ratio is obtained by numerical integration.
#'
#' @param effect Assumed true effect (nonzero), on the scale of the estimate.
#' @param se Standard error of the estimate (> 0).
#' @param alpha Two-sided significance level.
#' @param df Degrees of freedom of the reference distribution; `Inf` for normal.
#' @return A [data.table::data.table()] with one row per input combination
#'   (inputs are recycled to a common length).
#' @references Gelman A, Carlin J (2014). Beyond power calculations: assessing
#'   Type S (sign) and Type M (magnitude) errors. *Perspectives on
#'   Psychological Science* 9:641-651. \doi{10.1177/1745691614551642}
#' @examples
#' # Gelman, Hill & Vehtari (2020) Figure 16.1
#' ps_retrodesign(effect = 2, se = 8.1)
#' @export
ps_retrodesign <- function(effect, se, alpha = 0.05, df = Inf) {
  n <- max(length(effect), length(se), length(alpha), length(df))
  effect <- rep_len(effect, n); se <- rep_len(se, n)
  alpha <- rep_len(alpha, n); df <- rep_len(df, n)
  if (any(!is.finite(effect)) || any(effect == 0)) stop("`effect` must be finite and nonzero")
  if (any(!is.finite(se)) || any(se <= 0)) stop("`se` must be finite and > 0")
  if (any(alpha <= 0 | alpha >= 1)) stop("`alpha` must be in (0, 1)")

  res <- lapply(seq_len(n), function(i) .retro_one(effect[i], se[i], alpha[i], df[i]))
  out <- data.table::rbindlist(res)
  out[]
}

.retro_one <- function(D, s, alpha, df) {
  q <- if (is.infinite(df)) stats::qnorm(1 - alpha / 2) else stats::qt(1 - alpha / 2, df)
  lam <- D / s
  if (is.infinite(df)) {
    p_hi <- stats::pnorm(q - lam, lower.tail = FALSE)
    p_lo <- stats::pnorm(-q - lam)
  } else {
    p_hi <- stats::pt(q - lam, df, lower.tail = FALSE)
    p_lo <- stats::pt(-q - lam, df)
  }
  power <- p_hi + p_lo
  type_s <- if (D > 0) p_lo / power else p_hi / power

  cut <- q * s
  if (is.infinite(df)) {
    a <- (cut - D) / s
    b <- (-cut - D) / s
    # E[X; X > cut] and E[-X; X < -cut] for X ~ N(D, s^2)
    e_hi <- D * stats::pnorm(a, lower.tail = FALSE) + s * stats::dnorm(a)
    e_lo <- -D * stats::pnorm(b) + s * stats::dnorm(b)
  } else {
    dens <- function(x) stats::dt((x - D) / s, df) / s
    e_hi <- stats::integrate(function(x) x * dens(x), cut, Inf, rel.tol = 1e-8)$value
    e_lo <- stats::integrate(function(x) -x * dens(x), -Inf, -cut, rel.tol = 1e-8)$value
  }
  exaggeration <- if (power > 0) ((e_hi + e_lo) / power) / abs(D) else NA_real_

  data.table::data.table(
    effect = D, se = s, alpha = alpha, df = df, snr = abs(D) / s,
    crit = q, sig_cutoff = cut, power = power, type_s = type_s,
    exaggeration = exaggeration, min_sig_ratio = cut / abs(D)
  )
}

#' Simulate the sampling distribution of an estimate and summarise significance
#'
#' Monte Carlo counterpart of [ps_retrodesign()], used to check the analytic
#' values.
#'
#' @inheritParams ps_retrodesign
#' @param n_sims Number of draws.
#' @param seed Optional seed; if supplied the RNG state is restored afterwards.
#' @return A one-row data.table with `sim_power`, `sim_type_s`,
#'   `sim_exaggeration` and `sim_n_sig`.
#' @examples
#' ps_retrodesign_sim(2, 8.1, n_sims = 1e5, seed = 1)
#' @export
ps_retrodesign_sim <- function(effect, se, alpha = 0.05, df = Inf, n_sims = 1e6, seed = NULL) {
  if (!is.null(seed)) {
    old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
    on.exit(if (!is.null(old)) assign(".Random.seed", old, envir = globalenv()), add = TRUE)
    set.seed(seed)
  }
  q <- if (is.infinite(df)) stats::qnorm(1 - alpha / 2) else stats::qt(1 - alpha / 2, df)
  est <- effect + se * (if (is.infinite(df)) stats::rnorm(n_sims) else stats::rt(n_sims, df))
  sig <- abs(est) >= q * se
  data.table::data.table(
    sim_power = mean(sig),
    sim_type_s = if (any(sig)) mean(sign(est[sig]) != sign(effect)) else NA_real_,
    sim_exaggeration = if (any(sig)) mean(abs(est[sig])) / abs(effect) else NA_real_,
    sim_n_sig = sum(sig)
  )
}
