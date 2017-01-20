#' Producing two-dimensional invariant regions
#'
#' This function produces two-dimensional threshold lines and invariant regions,
#' as shown by Phillippo \emph{et al.} (2016).
#'
#' @param thresh A \code{thresh} object, as produced by
#'   \code{\link{nma_thresh}}.
#' @param idx Integer specifying the index (with respect to
#'   \code{thresh$thresholds}) of the first data point to consider adjusting.
#'   Will be shown on the x axis.
#' @param idy Integer specifying the index (with respect to
#'   \code{thresh$thresholds}) of the second data point to consider adjusting.
#'   Will be shown on the y axis.
#' @param xlab Character string giving the label for the x axis.
#' @param ylab Character string giving the label for the y axis.
#' @param xlim Numeric vector of length 2, giving the x axis limits.
#' @param ylim Numeric vector of length 2, giving the y axis limits.
#' @param breaks Numeric vector giving position of tick marks on the x and y
#'   axes. Calculated automatically by default.
#' @param xbreaks Numeric vector giving position of tick marks on the x axis.
#'   Equal to \code{breaks} by default, if set this overrides any value given to
#'   \code{breaks}.
#' @param ybreaks Numeric vector giving position of tick marks on the y axis.
#'   Equal to \code{breaks} by default, if set this overrides any value given to
#'   \code{breaks}.
#' @param fill Fill colour for invariant region. Defaults to a nice shade of
#'   blue \code{rgb(.72, .80, .93, .7)}.
#' @param lwd Line width for threshold lines. Default 1.
#' @param fontsize Font size for labels. Default 12.
#'
#' @import ggplot2
#'
#' @return A \code{ggplot} object containing the 2D threshold plot, which is
#'   returned invisibly and plotted (unless assigned).
#' @export
#'
#' @examples
#'
thresh_2d <- function(thresh, idx, idy,
                      xlab = "Adjustment to data point 1", ylab = "Adjustment to data point 2",
                      xlim = NULL, ylim = NULL,
                      breaks = ggplot2::waiver(), xbreaks = breaks, ybreaks = breaks,
                      fill = rgb(.72, .80, .93, .7),
                      lwd = 1,
                      fontsize = 12){

  # Number of treatments
  K <- nrow(thresh$Ukstar) + 1

  # Derive intercept and gradient of threshold lines using Ukstar matrix
  linedat <- data.frame(intercept = thresh$Ukstar[, idy],
                        gradient = - thresh$Ukstar[, idy] / thresh$Ukstar[, idx])


  # Intercepts array [line1, line2, x/y]
  l1s <- rep(1:(K-1), each = K-1)
  l2s <- rep(1:(K-1), times = K-1)

  intarray <- array(NA, dim = c(K-1 , K-1, 2))

  intarray[,,1] <- (linedat[l2s, "intercept"] - linedat[l1s, "intercept"]) /
                      (linedat[l1s, "gradient"] - linedat[l2s, "gradient"])

  diag(intarray[,,1]) <- NA   # intercept of line with itself

  intarray[,,2] <- linedat[, "gradient"] * intarray[,,1] + linedat[, "intercept"]

  # Intercept array to data frame [x,y]
  uptri <- upper.tri(intarray[,,1])
  intdat <- data.frame(x = intarray[,,1][uptri],
                       y = intarray[,,2][uptri],
                       l1 = matrix(rep(1:(K-1), times = K-1), nrow = K-1)[uptri],
                       l2 = matrix(rep(1:(K-1), each = K-1), nrow = K-1)[uptri])

  # Calculate invariant region
  # For each threshold line:
  #  1. Check which side of the line the origin lies
  #  2. Exclude all intercept points which lie the opposite side
  #  3. Repeat 1-2 for each line

  # Function to determine inclusion/exclusion of points
  IRvertex <- function(Mx,    # x value of test point
                       My,    # y value of test point
                       xint,  # x intercepts of all threshold lines
                       yint,  # y intercepts of all threshold lines
                       eps=1e-12){  # for testing equality

    # Quick parameter checks
    stopifnot(length(xint)==length(yint))
    stopifnot(length(Mx)==1)
    stopifnot(length(My)==1)

    # Is test point the same side of every line as the origin?
    # Ignore lines which a point lies on (within eps)
    all(ifelse(abs((xint - Mx)*(yint - My) - Mx*My) <= eps,
               TRUE,
               sign((xint - Mx)*(yint - My) - Mx*My) == sign(xint * yint)))
  }


  # For each intercept point, check inclusion/exclusion
  inIR <- mapply(IRvertex, Mx = intdat$x, My = intdat$y,
                 MoreArgs = list(xint = thresh$Ukstar[, idx],
                                 yint = thresh$Ukstar[, idy]))


  # Invariant region data
  IRdat <- intdat[inIR, ]

  # Sort into clockwise order
  IRdat$angle <- atan2(IRdat$y, IRdat$x)
  IRdat <- IRdat[order(IRdat$angle),]

  # Label data
  # labdat <- data.frame(
  #   x = c(-0.8, 0.03, -0.5),
  #   y = c(6.7, 9.5, -1.2),
  #   lab = paste0("paste(tilde(k),'* = ',", c(6, 5, 4),")")
  # )
  labdat <- list(x=NA_real_, y=NA_real_, lab=NA_character_)
  for (i in 1:(K-1)) {
    labdat$x[i] <- mean(IRdat[IRdat$l1 == i | IRdat$l2 == i, "x"])
    labdat$y[i] <- mean(IRdat[IRdat$l1 == i | IRdat$l2 == i, "y"])
    labdat$lab[i] <- paste0("paste(tilde(k),'* = ',", i + (i>=thresh$kstar), ")")
  }
  labdat <- as.data.frame(labdat)


  # Construct plot
  ggplot2::ggplot() +

    # Axes
    ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +

    # Threshold lines
    ggplot2::geom_abline(ggplot2::aes(intercept = intercept, slope = gradient),
                         data = linedat, colour = "grey60") +

    # Invariant region
    ggplot2::geom_polygon(ggplot2::aes(x = x, y = y), data = IRdat,
                          fill = fill, colour = "black") +
    #ggplot2::geom_point(ggplot2::aes(x = x, y = y), data = intdat, colour = ifelse(inIR, "red", "black")) +

    # Line labels
    ggplot2::geom_label(ggplot2::aes(x = x, y = y, label = lab),
                        data = labdat, parse = TRUE, na.rm = TRUE) +

    # Axis labels
    ggplot2::xlab(xlab) +
    ggplot2::ylab(ylab) +

    # Axis setup
    ggplot2::coord_cartesian(xlim = xlim, ylim = ylim) +
    ggplot2::scale_x_continuous(breaks = xbreaks) +
    ggplot2::scale_y_continuous(breaks = ybreaks) +

    # Theme
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
}
