# Internal helpers shared by the design functions.

# Build the full crossing of all vector arguments (NULLs dropped).
.ps_grid <- function(args) {
  args <- args[!vapply(args, is.null, logical(1))]
  lens <- vapply(args, length, integer(1))
  if (any(lens == 0)) stop("zero-length argument: ", paste(names(args)[lens == 0], collapse = ", "))
  do.call(data.table::CJ, c(args, list(sorted = FALSE)))
}

# Attach engine output to a design grid.
# grid: data.table of inputs; effect/se/df: vectors aligned with grid rows.
.ps_finish <- function(grid, design, scale, effect, se, df, alpha) {
  eng <- ps_retrodesign(effect = effect, se = se, alpha = alpha, df = df)
  eng[, alpha := NULL]
  out <- cbind(data.table::data.table(design = design, scale = scale), grid, eng)
  data.table::setattr(out, "class", c("ps_result", class(out)))
  out[]
}

.check_pos <- function(x, name, allow_zero = FALSE) {
  bad <- if (allow_zero) any(!is.finite(x) | x < 0) else any(!is.finite(x) | x <= 0)
  if (bad) stop("`", name, "` must be ", if (allow_zero) ">= 0" else "> 0", call. = FALSE)
}

.check_prob <- function(x, name) {
  if (any(!is.finite(x) | x <= 0 | x >= 1)) stop("`", name, "` must be in (0, 1)", call. = FALSE)
}

#' @export
print.ps_result <- function(x, ...) {
  if (all(c("design", "scale") %in% names(x)) && nrow(x)) {
    cat("<powerscope result:", x$design[1], "| estimate scale:", x$scale[1], ">\n")
  }
  y <- data.table::copy(x)
  data.table::setattr(y, "class", c("data.table", "data.frame"))
  print(y, ...)
  invisible(x)
}
