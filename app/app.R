# powerscope Shiny app. Runs in the browser via Shinylive (webR) or locally:
#   Rscript tools/build_app.R && R -e 'shiny::runApp("app")'
library(shiny)
library(bslib)
library(data.table)
library(ggplot2)
if (file.exists("powerscope_core.R")) source("powerscope_core.R") else library(powerscope)

# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
fmt <- function(x, d = 3) {
  if (is.null(x) || length(x) == 0 || is.na(x)) return("NA")
  if (abs(x) < 1e-4 && x != 0) formatC(x, digits = 2, format = "e") else formatC(x, digits = d, format = "fg", flag = "#")
}
num <- function(id, label, value, min = NA, max = NA, step = NA) list(type = "num", id = id, label = label, value = value, min = min, max = max, step = step)
chk <- function(id, label, value) list(type = "chk", id = id, label = label, value = value)
sel <- function(id, label, choices, selected) list(type = "sel", id = id, label = label, choices = choices, selected = selected)

param_input <- function(ns, p) {
  switch(p$type,
    num = numericInput(ns(p$id), p$label, p$value, min = p$min, max = p$max, step = p$step),
    chk = checkboxInput(ns(p$id), p$label, p$value),
    sel = selectInput(ns(p$id), p$label, p$choices, p$selected))
}

# integer bisection: smallest n with power >= target (power assumed monotone in n)
solve_n <- function(f, target, lo = 2, hi = 1e6) {
  if (f(hi) < target) return(NA_integer_)
  if (f(lo) >= target) return(as.integer(lo))
  while (hi - lo > 1) { m <- (lo + hi) %/% 2; if (f(m) >= target) hi <- m else lo <- m }
  as.integer(hi)
}

r_call <- function(fun, args) {
  a <- vapply(names(args), function(n) paste0(n, " = ", deparse(args[[n]], width.cutoff = 500L)), "")
  paste0(fun, "(", paste(a, collapse = ", "), ")")
}

# ---------------------------------------------------------------------------
# design specifications: inputs -> arguments of the package function.
# `args(p, n, nvar)` builds the call; when `n` is given it replaces the
# sample-size variable `nvar` and any group-2 counterpart linked to it.
# Starting values are examples only; replace them with estimates from your data.
# ---------------------------------------------------------------------------
specs <- list(
  means = list(
    title = "Means (t-test)", fun = "ps_means", sim = "means",
    note = "Exact for the two-sample or paired t-test on normal data, including Type S and the exaggeration ratio.",
    params = list(num("delta", "True difference in means", 0.5, step = 0.05),
                  num("sd", "SD within group (or of paired differences)", 1, min = 1e-9, step = 0.05),
                  num("n1", "n in group 1 (or number of pairs)", 20, min = 2, step = 1),
                  num("ratio", "Allocation n2 / n1", 1, min = 0.05, step = 0.1),
                  chk("paired", "Paired design", FALSE)),
    nvars = c("Group 1 n (or pairs)" = "n1"), effect = "delta", default_alpha = 0.05,
    args = function(p, n = NULL, nvar = "n1") {
      n1 <- if (is.null(n)) p$n1 else n
      a <- list(delta = p$delta, sd = p$sd, n1 = n1, alpha = p$alpha, paired = isTRUE(p$paired))
      if (!isTRUE(p$paired)) a$n2 <- ceiling(n1 * p$ratio)
      a
    },
    interval = function(p) sign(p$delta) * c(1e-6 * p$sd, 50 * p$sd),
    curve = c(2, 100)),
  props = list(
    title = "Proportions", fun = "ps_props", sim = "props",
    note = "Exact enumeration of binomial outcomes for the z-test when n1 x n2 <= 4e6; normal approximation otherwise.",
    params = list(num("p1", "True proportion, group 1", 0.3, min = 1e-6, max = 1 - 1e-6, step = 0.01),
                  num("p2", "True proportion, group 2", 0.45, min = 1e-6, max = 1 - 1e-6, step = 0.01),
                  num("n1", "n in group 1", 80, min = 2, step = 1),
                  num("ratio", "Allocation n2 / n1", 1, min = 0.05, step = 0.1),
                  sel("method", "Test", c("Pooled z-test" = "pooled", "Wald (unpooled)" = "wald"), "pooled")),
    nvars = c("Group 1 n" = "n1"), effect = "p2", default_alpha = 0.05,
    args = function(p, n = NULL, nvar = "n1") {
      n1 <- if (is.null(n)) p$n1 else n
      list(p1 = p$p1, p2 = p$p2, n1 = n1, n2 = ceiling(n1 * p$ratio), alpha = p$alpha, method = p$method)
    },
    interval = function(p) if (p$p2 > p$p1) c(p$p1 + 1e-6, 1 - 1e-9) else c(1e-9, p$p1 - 1e-6),
    curve = c(10, 500)),
  rnaseq = list(
    title = "RNA-seq (one gene)", fun = "ps_rnaseq", sim = "rnaseq",
    note = "Negative binomial Wald test with known dispersion (Hart et al. 2013). With 2-3 replicates, estimated and shrunk dispersions make realised power lower.",
    params = list(num("log2fc", "True log2 fold change", 1, step = 0.1),
                  num("mean_count", "Mean normalised count, group 1", 50, min = 1e-6, step = 5),
                  num("dispersion", "NB dispersion (phi)", 0.1, min = 0, step = 0.01),
                  num("n1", "Replicates, group 1", 4, min = 2, step = 1),
                  num("ratio", "Allocation n2 / n1", 1, min = 0.05, step = 0.1),
                  num("depth_factor", "Depth multiplier vs. current", 1, min = 1e-3, step = 0.25)),
    nvars = c("Replicates per group" = "n1"), effect = "log2fc", default_alpha = 0.05,
    args = function(p, n = NULL, nvar = "n1") {
      n1 <- if (is.null(n)) p$n1 else n
      list(log2fc = p$log2fc, mean_count = p$mean_count, dispersion = p$dispersion,
           n1 = n1, n2 = ceiling(n1 * p$ratio), depth_factor = p$depth_factor, alpha = p$alpha)
    },
    interval = function(p) sign(p$log2fc) * c(1e-4, 15),
    curve = c(2, 30)),
  poolseq = list(
    title = "Pool-seq", fun = "ps_poolseq", sim = "poolseq",
    note = "Two-stage binomial sampling (chromosomes into the pool, then reads). Set overdispersion > 1 from replicate-pool variance. The plug-in Wald test can differ by about 1 point of power; use the simulation check.",
    params = list(num("p1", "Allele frequency, group 1", 0.5, min = 1e-6, max = 1 - 1e-6, step = 0.01),
                  num("p2", "Allele frequency, group 2", 0.6, min = 1e-6, max = 1 - 1e-6, step = 0.01),
                  num("pool_size1", "Individuals per pool", 50, min = 1, step = 1),
                  num("depth1", "Read depth per pool at a site", 60, min = 1, step = 5),
                  num("n_pools1", "Replicate pools per group", 1, min = 1, step = 1),
                  num("ploidy", "Ploidy", 2, min = 1, step = 1),
                  num("overdispersion", "Overdispersion (>= 1)", 1, min = 1, step = 0.1),
                  sel("target", "Inference target", c("Population" = "population", "The pools themselves" = "pool"), "population")),
    nvars = c("Read depth per pool" = "depth1", "Individuals per pool" = "pool_size1", "Replicate pools per group" = "n_pools1"),
    effect = "p2", default_alpha = 0.05,
    args = function(p, n = NULL, nvar = "depth1") {
      a <- list(p1 = p$p1, p2 = p$p2, pool_size1 = p$pool_size1, depth1 = p$depth1,
                n_pools1 = p$n_pools1, ploidy = p$ploidy, overdispersion = p$overdispersion,
                alpha = p$alpha, target = p$target)
      if (!is.null(n)) a[[nvar]] <- n
      a
    },
    interval = function(p) if (p$p2 > p$p1) c(p$p1 + 1e-6, 1 - 1e-9) else c(1e-9, p$p1 - 1e-6),
    curve = c(5, 400)),
  casecontrol = list(
    title = "Case-control", fun = "ps_casecontrol", sim = "casecontrol",
    note = "Allelic test, HWE, multiplicative model. alpha = 5e-8 is the human genome-wide convention; derive yours from the effective number of tests in your marker set.",
    params = list(num("maf", "Risk-allele frequency in controls", 0.2, min = 1e-6, max = 1 - 1e-6, step = 0.01),
                  num("odds_ratio", "Allelic odds ratio", 1.3, min = 1e-3, step = 0.05),
                  num("n_cases", "Cases", 2000, min = 2, step = 100),
                  num("ratio", "Controls per case", 1, min = 0.05, step = 0.5)),
    nvars = c("Cases" = "n_cases"), effect = "odds_ratio", default_alpha = 5e-8,
    args = function(p, n = NULL, nvar = "n_cases") {
      nc <- if (is.null(n)) p$n_cases else n
      list(maf = p$maf, odds_ratio = p$odds_ratio, n_cases = nc, n_controls = ceiling(nc * p$ratio), alpha = p$alpha)
    },
    interval = function(p) if (p$odds_ratio > 1) c(1 + 1e-6, 50) else c(1e-3, 1 - 1e-6),
    curve = c(100, 20000)),
  qtl = list(
    title = "QTL / quantitative trait", fun = "ps_qtl", sim = "qtl",
    note = "Additive SNP effect in trait-SD units, linear regression, HWE. Same alpha caveat as case-control.",
    params = list(num("maf", "Effect-allele frequency", 0.3, min = 1e-6, max = 1 - 1e-6, step = 0.01),
                  num("beta_sd", "Per-allele effect (trait SD)", 0.1, step = 0.01),
                  num("n", "Individuals", 5000, min = 3, step = 100)),
    nvars = c("Individuals" = "n"), effect = "beta_sd", default_alpha = 5e-8,
    args = function(p, n = NULL, nvar = "n") {
      list(maf = p$maf, beta_sd = p$beta_sd, n = if (is.null(n)) p$n else n, alpha = p$alpha)
    },
    interval = function(p) sign(p$beta_sd) * c(1e-5, 0.999 / sqrt(2 * p$maf * (1 - p$maf))),
    curve = c(100, 50000))
)

# ---------------------------------------------------------------------------
# design module
# ---------------------------------------------------------------------------
design_ui <- function(id, spec) {
  ns <- NS(id)
  layout_sidebar(
    sidebar = sidebar(width = 330,
      p(class = "small text-muted", "Starting values are examples. Replace them with estimates from your own data or the literature."),
      lapply(spec$params, function(p) param_input(ns, p)),
      numericInput(ns("alpha"), "alpha (two-sided, per test)", spec$default_alpha, min = 1e-15, max = 0.5),
      numericInput(ns("target_power"), "Target power", 0.8, min = 0.05, max = 0.999, step = 0.05),
      if (length(spec$nvars) > 1) selectInput(ns("nvar"), "Sample-size variable", spec$nvars),
      hr(),
      p(strong("Curve range")),
      fluidRow(column(6, numericInput(ns("cmin"), "from", spec$curve[1], min = 1)),
               column(6, numericInput(ns("cmax"), "to", spec$curve[2], min = 2))),
      checkboxInput(ns("logx"), "Log x-axis", spec$curve[2] / spec$curve[1] > 50)
    ),
    p(class = "text-muted", spec$note),
    uiOutput(ns("boxes")),
    card(full_screen = TRUE, card_header("Power, Type S and exaggeration across sample size"),
         plotOutput(ns("curve"), height = "300px")),
    navset_card_tab(
      nav_panel("Result table", tableOutput(ns("table"))),
      nav_panel("Curve data", tableOutput(ns("curve_table"))),
      nav_panel("Simulation check",
        p("Simulates the test a user would run, with estimated nuisance parameters, and compares it with the analytic values."),
        fluidRow(column(4, numericInput(ns("nsims"), "Simulated studies", if (spec$sim == "qtl") 1000 else 10000, min = 100, step = 1000)),
                 column(4, numericInput(ns("seed"), "Seed", 1, step = 1)),
                 column(4, br(), actionButton(ns("run_sim"), "Run simulation", class = "btn-primary"))),
        tableOutput(ns("sim"))),
      nav_panel("R code", verbatimTextOutput(ns("code")))
    )
  )
}

design_server <- function(id, spec) {
  moduleServer(id, function(input, output, session) {
    pars <- reactive({
      p <- lapply(c(vapply(spec$params, `[[`, "", "id"), "alpha", "target_power"), function(k) input[[k]])
      names(p) <- c(vapply(spec$params, `[[`, "", "id"), "alpha", "target_power")
      req(all(vapply(p, function(v) !is.null(v) && !(is.numeric(v) && is.na(v)), TRUE)))
      p
    })
    nvar <- reactive(if (length(spec$nvars) > 1) req(input$nvar) else unname(spec$nvars[1]))
    fun <- function() get(spec$fun)

    current <- reactive({
      a <- spec$args(pars())
      res <- tryCatch(do.call(fun(), a), error = function(e) e)
      validate(need(!inherits(res, "error"), if (inherits(res, "error")) conditionMessage(res)))
      list(args = a, res = res)
    })
    n_req <- reactive({
      p <- pars(); nv <- nvar()
      f <- function(n) tryCatch(do.call(fun(), spec$args(p, n, nv))$power, error = function(e) 0)
      solve_n(f, p$target_power, lo = if (nv %in% c("n_pools1", "pool_size1", "depth1")) 1 else 2, hi = 1e6)
    })
    mde <- reactive({
      p <- pars()
      g <- function(e) { q <- p; q[[spec$effect]] <- e; do.call(fun(), spec$args(q))$power - p$target_power }
      iv <- spec$interval(p)
      r <- tryCatch(uniroot(g, sort(iv), tol = 1e-9)$root, error = function(e) NA_real_)
      r
    })
    curve <- reactive({
      p <- pars(); nv <- nvar()
      req(input$cmin, input$cmax, input$cmax > input$cmin)
      ns <- unique(round(if (isTRUE(input$logx)) exp(seq(log(input$cmin), log(input$cmax), length.out = 80))
                         else seq(input$cmin, input$cmax, length.out = 80)))
      rbindlist(lapply(ns, function(n) {
        r <- tryCatch(do.call(fun(), spec$args(p, n, nv)), error = function(e) NULL)
        if (is.null(r)) NULL else r[, c(nv, "power", "type_s", "exaggeration"), with = FALSE]
      }))
    })

    output$boxes <- renderUI({
      r <- current()$res; p <- pars(); nr <- n_req(); m <- mde()
      nl <- names(spec$nvars)[spec$nvars == nvar()]
      layout_column_wrap(width = 1/4, fill = FALSE,
        value_box("Power", fmt(r$power), p(class = "small", "at the current design")),
        value_box("Type S rate", fmt(r$type_s), p(class = "small", "P(wrong sign | significant)")),
        value_box("Exaggeration ratio", fmt(r$exaggeration), p(class = "small", sprintf("min. %s for any significant result", fmt(r$min_sig_ratio)))),
        value_box(sprintf("%s for power %s", nl, fmt(p$target_power, 2)),
                  if (is.na(nr)) "not reachable" else format(nr, big.mark = ","),
                  if (is.na(nr)) p(class = "small", sprintf("power at %s = 1e6 is %s; other design variables limit it",
                                                         nvar(), fmt(do.call(fun(), spec$args(p, 1e6, nvar()))$power))),
                  p(class = "small", sprintf("MDE at current design: %s = %s", spec$effect, fmt(m))))
      )
    })
    output$curve <- renderPlot({
      d <- curve(); req(nrow(d) > 2)
      ps_plot(d, nvar(), target_power = pars()$target_power, log_x = isTRUE(input$logx)) +
        ggplot2::labs(x = names(spec$nvars)[spec$nvars == nvar()])
    }, res = 96)
    output$table <- renderTable({
      r <- copy(current()$res)
      drop <- intersect(c("design", "scale"), names(r))
      r <- r[, setdiff(names(r), drop), with = FALSE]
      data.table(quantity = names(r), value = vapply(r, function(v) if (is.numeric(v)) fmt(v, 5) else as.character(v), ""))
    })
    output$curve_table <- renderTable(curve(), digits = 5)
    output$code <- renderText({
      p <- pars(); a <- current()$args; nv <- nvar()
      a2 <- spec$args(p, NULL, nv); a2[[nv]] <- NULL
      a2 <- a2[!names(a2) %in% c("n2", "n_controls")]
      link <- if ("ratio" %in% names(p) && !isTRUE(p$paired)) {
        sprintf(", link = list(%s = %s)", if (spec$fun == "ps_casecontrol") "n_controls" else "n2", p$ratio)
      } else ""
      paste0(
        "library(powerscope)\n\n# current design\n", r_call(spec$fun, a), "\n\n",
        "# smallest ", nv, " for target power\n",
        "ps_sample_size(", spec$fun, ", \"", nv, "\", target_power = ", p$target_power, ",\n  ",
        sub("^[^(]*\\(", "", sub("\\)$", "", r_call(spec$fun, a2))), link, ")\n\n",
        "# check by simulation\nps_simulate(\"", spec$sim, "\", ",
        paste(deparse(a, width.cutoff = 500L), collapse = ""),
        ", n_sims = 1e4, seed = 1)\n")
    })
    sim_res <- eventReactive(input$run_sim, {
      a <- current()$args; r <- current()$res
      if (spec$sim == "poolseq" && (a$overdispersion != 1 || a$target != "population"))
        return(data.table(note = "Simulation supports target = population and overdispersion = 1 only."))
      if (spec$sim == "means" && isTRUE(a$paired))
        return(data.table(note = "Simulation check covers the two-sample t-test only."))
      if (spec$sim == "props" && a$method != "pooled")
        return(data.table(note = "Simulation check covers the pooled z-test only."))
      s <- withProgress(message = "Simulating", value = 0.5,
                        ps_simulate(spec$sim, a, n_sims = input$nsims, seed = input$seed))
      data.table(quantity = c("power", "Type S rate", "exaggeration ratio"),
                 analytic = c(fmt(r$power, 4), fmt(r$type_s, 3), fmt(r$exaggeration, 4)),
                 simulated = c(fmt(s$sim_power, 4), fmt(s$sim_type_s, 3), fmt(s$sim_exaggeration, 4)),
                 `Monte Carlo SE` = c(fmt(s$mc_se_power, 2), "", ""),
                 `simulated studies` = c(format(s$n_sims, big.mark = ","), "", ""))
    })
    output$sim <- renderTable(sim_res())
  })
}

# ---------------------------------------------------------------------------
# RNA-seq gene-table module (FDR across genes)
# ---------------------------------------------------------------------------
genes_ui <- function(id) {
  ns <- NS(id)
  layout_sidebar(
    sidebar = sidebar(width = 330,
      fileInput(ns("file"), "Gene table (CSV/TSV)", accept = c(".csv", ".tsv", ".txt", ".gz")),
      p(class = "small text-muted", "Columns: mean_count, dispersion (optional gene_id). From DESeq2: mcols(dds)$baseMean and mcols(dds)$dispersion."),
      actionButton(ns("demo"), "Load simulated example table"),
      hr(),
      numericInput(ns("log2fc"), "True |log2 FC| of DE genes", 1, step = 0.1),
      numericInput(ns("n1"), "Replicates per group", 4, min = 2, step = 1),
      numericInput(ns("depth"), "Depth multiplier vs. current", 1, min = 0.01, step = 0.25),
      numericInput(ns("q"), "Target FDR (BH)", 0.05, min = 1e-4, max = 0.5, step = 0.01),
      numericInput(ns("pi1"), "Fraction of genes truly DE", 0.1, min = 1e-4, max = 0.99, step = 0.01),
      p(class = "small text-muted", "The fraction of truly DE genes must be estimated for your system (e.g. from a pilot or comparable study); results depend strongly on it."),
      numericInput(ns("nmax"), "Curve: max replicates", 20, min = 3, step = 1)
    ),
    uiOutput(ns("source")),
    uiOutput(ns("boxes")),
    layout_column_wrap(width = 1/2,
      card(full_screen = TRUE, card_header("Expected discoveries by replicates per group"), plotOutput(ns("curve"), height = "280px")),
      card(full_screen = TRUE, card_header("Per-gene power at the current design"), plotOutput(ns("hist"), height = "280px"))),
    navset_card_tab(
      nav_panel("Summary", tableOutput(ns("summary"))),
      nav_panel("Curve data", tableOutput(ns("curve_table"))),
      nav_panel("Per-gene results", p("First 200 rows; download for all."), tableOutput(ns("per_gene")),
                downloadButton(ns("dl"), "Download per-gene TSV")),
      nav_panel("R code", verbatimTextOutput(ns("code"))))
  )
}

genes_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    tab <- reactiveVal(NULL); src <- reactiveVal(NULL)
    observeEvent(input$file, {
      d <- tryCatch(fread(input$file$datapath), error = function(e) NULL)
      if (is.null(d) || !all(c("mean_count", "dispersion") %in% names(d))) {
        showNotification("File must contain columns mean_count and dispersion.", type = "error"); return()
      }
      tab(d); src(sprintf("Uploaded file %s: %s rows.", input$file$name, format(nrow(d), big.mark = ",")))
    })
    observeEvent(input$demo, {
      set.seed(1)
      tab(data.table(gene_id = sprintf("g%05d", 1:5000), mean_count = rlnorm(5000, log(100), 1.5),
                     dispersion = rlnorm(5000, log(0.1), 0.6)))
      src("Simulated example: 5,000 genes, log-normal means (median 100) and dispersions (median 0.1). Not real data.")
    })
    output$source <- renderUI({
      if (is.null(tab())) p("Upload a gene table or load the simulated example.")
      else p(class = "text-muted", src())
    })
    fit <- reactive({
      req(tab(), input$log2fc != 0, input$n1 >= 2)
      tryCatch(ps_rnaseq_genes(tab(), log2fc = input$log2fc, n1 = input$n1, depth_factor = input$depth,
                               q = input$q, pi1 = input$pi1), error = function(e) e)
    })
    output$boxes <- renderUI({
      f <- fit()
      validate(need(!inherits(f, "error"), if (inherits(f, "error")) conditionMessage(f)))
      s <- f$summary
      layout_column_wrap(width = 1/4, fill = FALSE,
        value_box("Mean power (DE genes)", fmt(s$mean_power)),
        value_box("Expected true discoveries", fmt(s$expected_true_discoveries, 4)),
        value_box("Per-gene alpha for BH FDR", fmt(s$alpha_per_gene)),
        value_box("Exaggeration among discoveries", fmt(s$exaggeration_among_discoveries),
                  p(class = "small", "raw fold change / true, expected")))
    })
    curve <- reactive({
      req(tab(), input$nmax >= 3)
      rbindlist(lapply(2:input$nmax, function(n) {
        r <- tryCatch(ps_rnaseq_genes(tab(), input$log2fc, n, depth_factor = input$depth, q = input$q, pi1 = input$pi1)$summary,
                      error = function(e) NULL)
        if (is.null(r)) data.table(n1 = n, mean_power = 0, expected_true_discoveries = 0, alpha_per_gene = NA_real_,
                                   exaggeration_among_discoveries = NA_real_)
        else r[, .(n1, mean_power, expected_true_discoveries, alpha_per_gene, exaggeration_among_discoveries)]
      }))
    })
    output$curve <- renderPlot({
      d <- curve(); m1 <- nrow(tab()) * input$pi1
      ggplot(d, aes(n1, expected_true_discoveries)) +
        geom_hline(yintercept = m1, colour = "grey55", linetype = 2, linewidth = 0.4) +
        geom_line(colour = "#2a78d6", linewidth = 0.7) + geom_point(colour = "#2a78d6", size = 2) +
        annotate("text", x = min(d$n1), y = m1, label = sprintf("truly DE genes = %s", format(round(m1), big.mark = ",")),
                 hjust = 0, vjust = -0.5, size = 3.3, colour = "grey30") +
        labs(x = "Replicates per group", y = "Expected true discoveries") +
        theme_minimal(base_size = 11) + theme(panel.grid.minor = element_blank())
    }, res = 96)
    output$hist <- renderPlot({
      f <- fit(); req(!inherits(f, "error"))
      ggplot(f$per_gene, aes(power)) +
        geom_histogram(binwidth = 0.025, boundary = 0, fill = "#2a78d6", colour = "white", linewidth = 0.3) +
        labs(x = "Power at the per-gene alpha (if the gene were DE)", y = "Genes") +
        theme_minimal(base_size = 11) + theme(panel.grid.minor = element_blank())
    }, res = 96)
    output$summary <- renderTable({
      f <- fit(); req(!inherits(f, "error"))
      s <- f$summary
      data.table(quantity = names(s), value = vapply(s, function(v) fmt(as.numeric(v), 5), ""))
    })
    output$curve_table <- renderTable(curve(), digits = 5)
    output$per_gene <- renderTable({ f <- fit(); req(!inherits(f, "error")); head(f$per_gene, 200) }, digits = 4)
    output$dl <- downloadHandler(filename = function() "powerscope_per_gene.tsv",
                                 content = function(file) fwrite(fit()$per_gene, file, sep = "\t"))
    output$code <- renderText(paste0(
      "library(powerscope)\nlibrary(data.table)\n",
      "genes <- fread(\"your_gene_table.tsv\")   # columns mean_count, dispersion\n",
      sprintf("res <- ps_rnaseq_genes(genes, log2fc = %s, n1 = %s, depth_factor = %s, q = %s, pi1 = %s)\n",
              input$log2fc, input$n1, input$depth, input$q, input$pi1),
      "res$summary\n"))
  })
}

# ---------------------------------------------------------------------------
# overview: the retrodesign calculator
# ---------------------------------------------------------------------------
overview_ui <- function(id) {
  ns <- NS(id)
  layout_sidebar(
    sidebar = sidebar(width = 330,
      numericInput(ns("effect"), "True effect", 2, step = 0.5),
      numericInput(ns("se"), "Standard error of the estimate", 8.1, min = 1e-9, step = 0.5),
      numericInput(ns("alpha"), "alpha (two-sided)", 0.05, min = 1e-15, max = 0.5),
      numericInput(ns("df"), "df (blank = normal)", NA, min = 1),
      p(class = "small text-muted", "Defaults reproduce Gelman, Hill & Vehtari (2020), Regression and Other Stories, Figure 16.1.")),
    p("Any design reduces to a true effect and the standard error of its estimate. Given those, the calculator reports power, ",
      "the probability that a significant estimate has the wrong sign (Type S), and the expected overestimation of a significant ",
      "estimate (Type M, exaggeration ratio). The other tabs compute the standard error from design parameters."),
    uiOutput(ns("boxes")),
    card(full_screen = TRUE, card_header("Sampling distribution of the estimate"), plotOutput(ns("dens"), height = "320px")),
    card(card_header("Methods and assumptions"),
      tags$ul(
        tags$li("Engine: Gelman & Carlin (2014), Perspect Psychol Sci 9:641-651."),
        tags$li("t-tests and proportions are exact. RNA-seq, pool-seq and association use Wald tests with asymptotic normal estimates; each tab has a simulation check of the test actually run."),
        tags$li("FDR across genes: per-gene alpha from Jung (2005), Bioinformatics 21:3097-3104."),
        tags$li("RNA-seq variance: Hart et al. (2013), J Comput Biol 20:970-978. Pool-seq two-stage sampling: Gautier et al. (2013), Mol Ecol 22:3766-3779."),
        tags$li("All biological inputs (effects, SDs, dispersions, frequencies, overdispersion, fraction DE) must be estimated from your data or the literature."),
        tags$li(HTML("R package and source: <a href='https://github.com/Somatogenomics/powerscope' target='_blank'>github.com/Somatogenomics/powerscope</a>"))))
  )
}

overview_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    res <- reactive({
      req(input$effect, input$se, input$alpha); validate(need(input$effect != 0, "effect must be nonzero"))
      df <- if (is.na(input$df)) Inf else input$df
      ps_retrodesign(input$effect, input$se, input$alpha, df)
    })
    output$boxes <- renderUI({
      r <- res()
      layout_column_wrap(width = 1/4, fill = FALSE,
        value_box("Power", fmt(r$power)),
        value_box("Type S rate", fmt(r$type_s)),
        value_box("Exaggeration ratio", fmt(r$exaggeration)),
        value_box("Significance cutoff", paste0("|est| >= ", fmt(r$sig_cutoff)),
                  p(class = "small", sprintf("= %s x true effect", fmt(r$min_sig_ratio)))))
    })
    output$dens <- renderPlot({
      r <- res(); D <- r$effect; s <- r$se; cut <- r$sig_cutoff
      df <- r$df
      xs <- seq(min(-cut, D) - 3 * s, max(cut, D) + 3 * s, length.out = 1500)
      d <- data.table(x = xs, y = if (is.infinite(df)) dnorm(xs, D, s) else dt((xs - D) / s, df) / s)
      ggplot(d, aes(x, y)) +
        geom_area(data = d[x >= cut], fill = "#2a78d6", alpha = 0.45) +
        geom_area(data = d[x <= -cut], fill = "#2a78d6", alpha = 0.45) +
        geom_line(linewidth = 0.6) +
        geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.3) +
        geom_vline(xintercept = D, linewidth = 0.8) +
        annotate("text", x = D, y = max(d$y) * 1.04, label = sprintf("true effect = %s", fmt(D)), size = 3.4) +
        labs(x = "Estimated effect", y = "Density",
             caption = "Shaded: estimates declared significant. Thick line: true effect.") +
        theme_minimal(base_size = 11) + theme(panel.grid.minor = element_blank())
    }, res = 96)
  })
}

# ---------------------------------------------------------------------------
ui <- page_navbar(
  title = "powerscope",
  theme = bs_theme(version = 5, primary = "#2a78d6", base_font = font_collection("system-ui", "sans-serif")),
  nav_panel("Overview", overview_ui("ov")),
  nav_panel(specs$means$title, design_ui("means", specs$means)),
  nav_panel(specs$props$title, design_ui("props", specs$props)),
  nav_menu("RNA-seq",
    nav_panel("One gene", design_ui("rnaseq", specs$rnaseq)),
    nav_panel("All genes (FDR)", genes_ui("genes"))),
  nav_panel(specs$poolseq$title, design_ui("poolseq", specs$poolseq)),
  nav_menu("Association",
    nav_panel("Case-control", design_ui("casecontrol", specs$casecontrol)),
    nav_panel("QTL", design_ui("qtl", specs$qtl)))
)

server <- function(input, output, session) {
  overview_server("ov")
  for (nm in names(specs)) local({ n <- nm; design_server(n, specs[[n]]) })
  genes_server("genes")
}

shinyApp(ui, server)
