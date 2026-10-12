#!/usr/bin/env Rscript
# Analytic vs simulated power, Type S and exaggeration for every design.
# Simulation applies the test a user would run (estimated nuisance parameters),
# so differences show where the analytic approximation is accurate and where
# it is not (mainly very small samples).
#
# Usage (from the package root): Rscript validation/run_validation.R [n_sims]
suppressPackageStartupMessages({ library(data.table); devtools::load_all(quiet = TRUE) })
args <- commandArgs(TRUE)
n_sims <- if (length(args)) as.numeric(args[1]) else 1e5
SEED <- 20261012
out_dir <- "validation/results"; dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cases <- list(
  list("means", ps_means, list(delta = 0.5, sd = 1, n1 = 5,  n2 = 5)),
  list("means", ps_means, list(delta = 0.5, sd = 1, n1 = 20, n2 = 20)),
  list("means", ps_means, list(delta = 0.5, sd = 1, n1 = 64, n2 = 64)),
  list("means", ps_means, list(delta = 0.2, sd = 1, n1 = 10, n2 = 10)),
  list("props", ps_props, list(p1 = 0.5, p2 = 0.75, n1 = 30, n2 = 30)),
  list("props", ps_props, list(p1 = 0.1, p2 = 0.2,  n1 = 100, n2 = 150)),
  list("rnaseq", ps_rnaseq, list(log2fc = 1,   mean_count = 50,  dispersion = 0.1,  n1 = 3,  n2 = 3)),
  list("rnaseq", ps_rnaseq, list(log2fc = 1,   mean_count = 50,  dispersion = 0.1,  n1 = 10, n2 = 10)),
  list("rnaseq", ps_rnaseq, list(log2fc = 0.5, mean_count = 200, dispersion = 0.05, n1 = 6,  n2 = 6)),
  list("rnaseq", ps_rnaseq, list(log2fc = 0.5, mean_count = 10,  dispersion = 0.2,  n1 = 4,  n2 = 4)),
  list("poolseq", ps_poolseq, list(p1 = 0.5, p2 = 0.6, pool_size1 = 50,  depth1 = 60, n_pools1 = 1)),
  list("poolseq", ps_poolseq, list(p1 = 0.2, p2 = 0.4, pool_size1 = 30,  depth1 = 40, n_pools1 = 3)),
  list("poolseq", ps_poolseq, list(p1 = 0.5, p2 = 0.55, pool_size1 = 100, depth1 = 200, n_pools1 = 2)),
  list("casecontrol", ps_casecontrol, list(maf = 0.2, odds_ratio = 1.3, n_cases = 1000, n_controls = 1000, alpha = 0.05)),
  list("casecontrol", ps_casecontrol, list(maf = 0.05, odds_ratio = 1.5, n_cases = 500, n_controls = 1500, alpha = 0.05)),
  list("casecontrol", ps_casecontrol, list(maf = 0.3, odds_ratio = 1.2, n_cases = 3000, n_controls = 3000, alpha = 5e-8)),
  list("qtl", ps_qtl, list(maf = 0.3, beta_sd = 0.15, n = 300, alpha = 0.05)),
  list("qtl", ps_qtl, list(maf = 0.1, beta_sd = 0.3,  n = 150, alpha = 0.05))
)

res <- rbindlist(lapply(seq_along(cases), function(i) {
  cs <- cases[[i]]
  an <- do.call(cs[[2]], cs[[3]])
  ns <- if (cs[[1]] == "qtl") min(n_sims, 1e4) else n_sims
  sm <- ps_simulate(cs[[1]], cs[[3]], n_sims = ns, seed = SEED + i)
  data.table(case = i, design = cs[[1]],
             params = paste(names(cs[[3]]), unlist(cs[[3]]), sep = "=", collapse = ", "),
             power = an$power, sim_power = sm$sim_power, mc_se_power = sm$mc_se_power,
             power_z = (sm$sim_power - an$power) / sm$mc_se_power,
             type_s = an$type_s, sim_type_s = sm$sim_type_s,
             exaggeration = an$exaggeration, sim_exaggeration = sm$sim_exaggeration,
             n_sims = ns)
}))
res[, power_abs_diff := abs(sim_power - power)]
res[, exaggeration_rel_diff := sim_exaggeration / exaggeration - 1]
fwrite(res, file.path(out_dir, "validation_results.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
print(res[, .(case, design, params, power = round(power, 4), sim_power = round(sim_power, 4),
              power_z = round(power_z, 1), type_s = signif(type_s, 2), sim_type_s = signif(sim_type_s, 2),
              exag = round(exaggeration, 3), sim_exag = round(sim_exaggeration, 3))], nrows = 100)
