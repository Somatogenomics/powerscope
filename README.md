# powerscope

Power, sample size, minimum detectable effect, sign error rate (Type S) and
exaggeration ratio (Type M) for genomics study designs. Available as an R
package and as a browser-only app.

- App: <https://somatogenomics.github.io/powerscope/app/>
- Documentation: <https://somatogenomics.github.io/powerscope/>

## Why Type S and Type M

Power alone does not describe what a significant result from a weak design
means. When the true effect is small relative to the standard error,
significant estimates overstate the effect and can have the wrong sign
(Gelman & Carlin 2014). Every function here reports both alongside power.
For the Gelman, Hill & Vehtari (2020) Figure 16.1 design (true effect 2,
SE 8.1):

```r
library(powerscope)
ps_retrodesign(effect = 2, se = 8.1)
#> power 0.0570, type_s 0.240, exaggeration 9.54, min_sig_ratio 7.94
```

## Designs

| Function | Design | Basis |
|---|---|---|
| `ps_retrodesign()` | any effect + SE | Gelman & Carlin (2014); closed form (normal) or integration (t) |
| `ps_means()` | two-group or paired t-test | exact (noncentral t; SD integrated out for Type M) |
| `ps_props()` | two proportions, pooled or Wald z | exact enumeration; normal approximation for very large n |
| `ps_rnaseq()` | one gene, NB counts | Hart et al. (2013) variance, Wald test, known dispersion |
| `ps_rnaseq_genes()` | all genes from your own mean/dispersion table, BH FDR | per-gene alpha from Jung (2005) |
| `ps_poolseq()` | allele-frequency difference between pools | two-stage binomial (Gautier et al. 2013) + overdispersion |
| `ps_casecontrol()` | allelic case-control association | log-OR Wald, HWE, multiplicative model |
| `ps_qtl()` | additive SNP on a quantitative trait | linear regression |

Solvers and helpers: `ps_sample_size()`, `ps_mde()`, `ps_fdr_alpha()`,
`ps_simulate()` (simulates the test actually run), `ps_plot()`.

All vector arguments are fully crossed, so a single call gives a design grid:

```r
ps_rnaseq(log2fc = 1, mean_count = c(10, 50, 500), dispersion = 0.1, n1 = 2:12)
ps_sample_size(ps_means, "n1", 0.8, delta = 0.5, sd = 1, link = list(n2 = 1))   # n1 = n2 = 64
ps_mde(ps_casecontrol, "odds_ratio", 0.8, maf = 0.3, n_cases = 2000, interval = c(1.0001, 5))
```

RNA-seq across genes, using means and dispersions from your own data:

```r
genes <- data.frame(mean_count = mcols(dds)$baseMean, dispersion = mcols(dds)$dispersion)
ps_rnaseq_genes(genes, log2fc = 1, n1 = 4, q = 0.05, pi1 = 0.1)$summary
```

## Inputs must come from data

Effect sizes, SDs, dispersions, allele frequencies, overdispersion and the
fraction of true effects are model inputs. The package has no defaults for
them. Estimate them from pilot data, previous experiments or the literature,
and report results over a range of plausible values. The app's starting
values are examples only.

## Validation

`validation/run_validation.R` compares every analytic result with a
simulation of the test a user would actually run, with estimated nuisance
parameters (1e5 simulated studies per case; QTL 1e4). Results are in
`validation/results/validation_results.tsv`:

- t-test, proportions, RNA-seq NB and QTL: analytic power within Monte
  Carlo error of simulation (|z| ≤ 1.2); exaggeration within 3%.
- Pool-seq: the plug-in Wald test rejects slightly more often than the
  analytic value at p ≈ 0.5 (+0.9 points of power), because estimated
  `p(1-p)` is at most its true value there.
- Case-control with a rare allele (MAF 0.05, OR 1.5, unequal groups):
  analytic power 1.6 points above the simulated Wald test.

Where accuracy at a specific design matters, run `ps_simulate()` or the
app's "Simulation check" tab.

The unit tests (`tests/testthat`) check exact agreement with
`stats::power.t.test(strict = TRUE)` and `stats::power.prop.test(strict = TRUE)`,
the textbook n = 64 per group for d = 0.5, the Figure 16.1 values, FDR identities,
and simulation agreement for each design.

## Assumptions

- RNA-seq: Wald test with known dispersion. DESeq2/edgeR estimate and shrink
  dispersions, so realised power with 2-3 replicates is lower.
- Pool-seq: equal DNA contribution per individual unless `overdispersion` is
  set; estimate it from the variance between replicate pools.
- Association: Hardy-Weinberg equilibrium, no population structure;
  multiplicative (case-control) or additive (QTL) model. The default
  `alpha = 5e-8` is the human genome-wide convention; derive the threshold for
  your species and marker set from the effective number of independent tests.
- Multiple testing other than BH across genes is handled by setting `alpha`.

## Install and run

```r
# install.packages("remotes")
remotes::install_github("Somatogenomics/powerscope", build_vignettes = TRUE)
```

Run the app locally from a source checkout:

```sh
Rscript tools/build_app.R      # copies R/ into app/powerscope_core.R
Rscript tools/test_app.R       # server-side tests of every module
R -e 'shiny::runApp("app")'
```

Tested with R 4.3.3, data.table 1.14.10, ggplot2 3.4.4, shiny 1.8.0,
bslib 0.6.1, testthat 3.2.1.

## Deployment

`.github/workflows/pages.yaml` builds the pkgdown site into `docs/`, exports
the app with `shinylive::export()` into `docs/app/`, and deploys to GitHub
Pages. In the repository settings, set Pages > Source to "GitHub Actions".
`.github/workflows/R-CMD-check.yaml` runs `R CMD check` on Linux, macOS and
Windows.

## References

- Gelman A, Carlin J (2014). Beyond power calculations: assessing Type S (sign) and Type M (magnitude) errors. *Perspect Psychol Sci* 9:641-651. doi:10.1177/1745691614551642
- Gelman A, Hill J, Vehtari A (2020). *Regression and Other Stories*, ch. 16. Cambridge University Press.
- Hart SN, Therneau TM, Zhang Y, Poland GA, Kocher JP (2013). Calculating sample size estimates for RNA sequencing data. *J Comput Biol* 20:970-978. doi:10.1089/cmb.2012.0283
- Jung SH (2005). Sample size for FDR-control in microarray data analysis. *Bioinformatics* 21:3097-3104. doi:10.1093/bioinformatics/bti456
- Gautier M et al. (2013). Estimation of population allele frequencies from next-generation sequencing data: pool-versus individual-based genotyping. *Mol Ecol* 22:3766-3779. doi:10.1111/mec.12360
- Pe'er I, Yelensky R, Altshuler D, Daly MJ (2008). Estimation of the multiple testing burden for genomewide association studies of nearly all common variants. *Genet Epidemiol* 32:381-385. doi:10.1002/gepi.20303

## License

MIT © Abdulhadi Kobiowu
