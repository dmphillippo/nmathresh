#' Producing two-dimensional invariant regions
#'
#' This function produces two-dimensional threshold lines and invariant regions,
#' as shown by Phillippo \emph{et al.} (2016).
#'
#' @param thresh A \code{thresh} object, as produced by
#'   \code{\link{nma_thresh}}.
#' @param idx Integer specifying the index (with respect to \code{thresh$thresholds}) of the first data point to consider adjusting. Will be shown on the x axis.
#' @param idy Integer specifying the index (with respect to \code{thresh$thresholds}) of the second data point to consider adjusting. Will be shown on the y axis.
#' @param xlab Character string giving the label for the x axis.
#' @param ylab Character string giving the label for the y axis.
#' @param xlim Numeric vector of length 2, giving the x axis limits.
#' @param ylim Numeric vector of length 2, giving the y axis limits.
#' @param breaks Numeric vector giving position of tick marks on the x and y axes. Calculated automatically by default.
#' @param xbreaks Numeric vector giving position of tick marks on the x axis. Equal to \code{breaks} by default, if set this overrides any value given to \code{breaks}.
#' @param ybreaks Numeric vector giving position of tick marks on the y axis. Equal to \code{breaks} by default, if set this overrides any value given to \code{breaks}.
#' @param fill Fill colour for invariant region. Defaults to a nice shade of blue \code{rgb(.72, .80, .93)}.
#' @param lwd Line width for threshold lines. Default 1.
#' @param fontsize Font size for labels. Default 12.
#'
#' @return
#' @export
#'
#' @examples
#'
thresh_2d <- function(thresh, idx, idy,
                      xlab = "", ylab = "",
                      xlim = NULL, ylim = NULL,
                      breaks = NULL, xbreaks = breaks, ybreaks = breaks,
                      fill = rgb(.72, .80, .93),
                      lwd = 1,
                      fontsize = 12){
  NULL
}
