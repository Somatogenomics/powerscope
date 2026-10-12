#!/usr/bin/env Rscript
# Server-side tests of the Shiny app modules (no browser needed).
# Run from the package root after tools/build_app.R:  Rscript tools/test_app.R
suppressPackageStartupMessages(library(shiny))
setwd("app")
src <- readLines("app.R")
src <- src[seq_len(max(grep("^shinyApp\\(", src)) - 1)]   # everything except the launch
eval(parse(text = src), envir = globalenv())

ok <- function(cond, msg) { cat(if (cond) "PASS" else "FAIL", msg, "\n"); if (!cond) quit(status = 1) }
defaults <- function(spec) {
  v <- lapply(spec$params, function(p) switch(p$type, num = p$value, chk = p$value, sel = p$selected))
  names(v) <- vapply(spec$params, `[[`, "", "id")
  c(v, list(alpha = spec$default_alpha, target_power = 0.8, nvar = unname(spec$nvars[1]),
            cmin = spec$curve[1], cmax = spec$curve[2], logx = FALSE, nsims = 2000, seed = 1))
}

for (nm in names(specs)) {
  spec <- specs[[nm]]
  testServer(design_server, args = list(spec = spec), {
    do.call(session$setInputs, defaults(spec))
    r <- current()$res
    pkg <- do.call(get(spec$fun), spec$args(defaults(spec)))
    ok(isTRUE(all.equal(r$power, pkg$power)), sprintf("[%s] app power equals package power (%.4f)", nm, r$power))
    n <- n_req()
    if (!is.na(n)) {
      p <- defaults(spec); nv <- unname(spec$nvars[1])
      pw  <- do.call(get(spec$fun), spec$args(p, n, nv))$power
      pw1 <- if (n > 2) do.call(get(spec$fun), spec$args(p, n - 1, nv))$power else 0
      ok(pw >= 0.8 && pw1 < 0.8, sprintf("[%s] required %s = %d (power %.3f; at n-1 %.3f)", nm, nv, n, pw, pw1))
    } else ok(TRUE, sprintf("[%s] target not reachable within 1e6", nm))
    m <- mde()
    q <- defaults(spec); q[[spec$effect]] <- m
    ok(is.na(m) || abs(do.call(get(spec$fun), spec$args(q))$power - 0.8) < 1e-6,
       sprintf("[%s] MDE %s = %s gives power 0.8", nm, spec$effect, format(m, digits = 5)))
    ok(nrow(curve()) > 10, sprintf("[%s] curve has %d points", nm, nrow(curve())))
    ok(nchar(output$code) > 50, sprintf("[%s] R code rendered", nm))
    ok(grepl("Power", as.character(output$boxes$html)), sprintf("[%s] value boxes rendered", nm))
    session$setInputs(run_sim = 1)
    s <- sim_res()
    ok("simulated" %in% names(s), sprintf("[%s] simulation check ran: analytic %s vs simulated %s",
                                          nm, s$analytic[1], s$simulated[1]))
  })
}

testServer(overview_server, {
  session$setInputs(effect = 2, se = 8.1, alpha = 0.05, df = NA)
  ok(abs(res()$exaggeration - 9.5366) < 1e-3, "[overview] Figure 16.1 exaggeration 9.54")
  ok(grepl("0.0570", as.character(output$boxes$html)), "[overview] power box shows 0.0570")
})

testServer(genes_server, {
  session$setInputs(log2fc = 1, n1 = 4, depth = 1, q = 0.05, pi1 = 0.1, nmax = 8, demo = 1)
  s <- fit()$summary
  ok(abs(s$expected_false_discoveries / (s$expected_false_discoveries + s$expected_true_discoveries) - 0.05) < 1e-6,
     "[genes] FDR equals 0.05 at the solved per-gene alpha")
  ok(nrow(curve()) == 7, "[genes] curve over 2..8 replicates")
  ok(all(diff(curve()$expected_true_discoveries) > 0), "[genes] expected discoveries increase with replicates")
})
cat("all app tests passed\n")
