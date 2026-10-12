.ps_palette <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6")

#' Plot power, Type S and exaggeration against a design parameter
#'
#' Three panels on separate y-scales (no shared axis): power, Type S rate and
#' exaggeration ratio. A second varying parameter, if any, is shown as
#' separate lines in a fixed colour order (at most 7; more are refused).
#'
#' @param x A `ps_result` from any design function.
#' @param xvar Column to put on the x-axis (e.g. `"n1"`).
#' @param group Optional column distinguishing lines.
#' @param target_power Optional reference line on the power panel.
#' @param log_x Logical; log10 x-axis.
#' @return A ggplot object.
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   ps_plot(ps_means(delta = c(0.3, 0.5), sd = 1, n1 = 4:80), "n1", group = "delta",
#'           target_power = 0.8)
#' }
#' @export
ps_plot <- function(x, xvar, group = NULL, target_power = NULL, log_x = FALSE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required for ps_plot()")
  d <- data.table::as.data.table(x)
  keep <- c(xvar, group, "power", "type_s", "exaggeration")
  if (!all(keep %in% names(d))) stop("missing columns: ", paste(setdiff(keep, names(d)), collapse = ", "))
  d <- d[, keep, with = FALSE]
  d[, grp := if (is.null(group)) "all" else factor(format(get(group), digits = 4))]
  if (data.table::uniqueN(d$grp) > length(.ps_palette)) stop("at most 7 lines; facet or subset instead")
  long <- data.table::melt(d, id.vars = c(xvar, "grp"),
                           measure.vars = c("power", "type_s", "exaggeration"))
  labs <- c(power = "Power", type_s = "Type S rate (wrong sign | significant)",
            exaggeration = "Exaggeration ratio E(|est| | sig) / |true|")
  long[, variable := factor(labs[as.character(variable)], levels = labs)]
  p <- ggplot2::ggplot(long, ggplot2::aes(x = .data[[xvar]], y = value, colour = grp)) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::facet_wrap(~variable, scales = "free_y", nrow = 1) +
    ggplot2::scale_colour_manual(values = .ps_palette, name = group) +
    ggplot2::labs(x = xvar, y = NULL) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   strip.text = ggplot2::element_text(hjust = 0, face = "bold"),
                   legend.position = if (is.null(group)) "none" else "top")
  ref <- data.table::data.table(variable = factor(labs["exaggeration"], levels = labs), y = 1)
  p <- p + ggplot2::geom_hline(data = ref, ggplot2::aes(yintercept = y), colour = "grey55",
                               linetype = 2, linewidth = 0.4)
  if (!is.null(target_power)) {
    tp <- data.table::data.table(variable = factor(labs["power"], levels = labs), y = target_power)
    p <- p + ggplot2::geom_hline(data = tp, ggplot2::aes(yintercept = y), colour = "grey55",
                                 linetype = 2, linewidth = 0.4)
  }
  if (log_x) p <- p + ggplot2::scale_x_log10()
  p
}
