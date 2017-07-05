#' Calculate thresholds and invariant intervals
#'
#' This function calculates decision-invariant bias-adjustment thresholds and
#' intervals for Bayesian network meta-analysis, as described by Phillippo
#' \emph{et al.} (under review). Thresholds are derived from the joint
#' posterior, and reflect the amount of change to a data point before the
#' treatment decision changes. Calculation is achieved using fast matrix
#' operations.
#'
#' @param mean.dk Posterior means of basic treatment parameters \eqn{d_k}.
#' @param lhood Likelihood (data) covariance matrix.
#' @param post Posterior covariance matrix (see details).
#' @param nmatype Character string, giving the type of NMA performed. One of
#'   "fixed" (fixed effects, the default) or "random" (random effects). May be
#'   abbreviated.
#' @param X [FE models only] Design matrix for basic treatment parameters.
#' @param mu.design [RE models only] Design matrix for any extra parameters.
#'   Defaults to NULL (no extra parameters).
#' @param delta.design [RE models only] Design matrix for delta, defaults to the
#'   \eqn{N \times N}{N x N} identity matrix.
#' @param opt.max Should the optimal decision be the maximal treatment effect
#'   (TRUE, default) or the minimum (FALSE).
#' @param trt.rank Rank of the treatment to derive thresholds for. Defaults to
#'   1, thresholds for the optimum treatment.
#' @param trt.bestn Derive thresholds for a decision over the best
#'   \code{trt.bestn} treatments. Defaults to 1, thresholds for the single
#'   optimum treatment only. Overrides \code{trt.rank}, with a warning.
#' @param trt.code Treatment codings of the reference trt and in the parameter
#'   vector \eqn{d_k}. Use if treatments re-labelled or re-ordered. Default is
#'   equivalent to 1:K.
#' @param trt.sub Only look at thresholds in this subset of treatments in
#'   trt.code, e.g. if some are excluded from the ranking. Default is equivalent
#'   to 1:K.
#'
#' @details This function provides bias-adjustment threshold analysis for both
#'   fixed and random effects NMA models, as described by Phillippo \emph{et
#'   al.} (under review). Parameters \code{mean.dk}, \code{lhood}, and
#'   \code{post} are always required, however there are differences in the
#'   specification of \code{post} and other required parameters and between the
#'   fixed and random effects cases:
#'
#'   \describe{ \item{\strong{FE models}}{The design matrix \code{X} for basic
#'   treatment parameters is required. The posterior covariance matrix specified
#'   in \code{post} should only refer to the basic treatment parameters.}
#'   \item{\strong{RE models}}{The design matrix \code{mu.design} for additional
#'   parameters (e.g. covariates) is required, as is the design matrix
#'   \code{delta.design} for random effects terms. The posterior covariance
#'   matrix specified in \code{post} should refer to the basic treatment
#'   parameters, RE terms, and additional parameters \emph{in that order}; i.e.
#'   \code{post} should be the posterior covariance matrix of the vector
#'   \eqn{(d^T, \delta^T, \mu^T)^T}.} }
#'
#' @section Model details: \strong{The FE NMA model}
#'
#'   The fixed effects NMA model is assumed to be of the form \describe{
#'   \item{Prior:}{\eqn{d \sim \mathrm{N} ( d_0, \Sigma_d )}{d ~ N(d_0,
#'   \Sigma_d)}}
#'   \item{Likelihood:}{\eqn{y|d \sim \mathrm{N} ( \delta, V )}{y|d ~ N(\delta,
#'   V)}} \item{FE model:}{\eqn{\delta = Xd + M\mu}} }
#'
#'   The additional parameters \eqn{\mu} may be given any sensible prior; they
#'   do not affect the threshold analysis in any way.
#'
#'   \strong{The RE NMA model}
#'
#'   The random effects NMA model is assumed to be of the form \describe{
#'   \item{Priors:}{\eqn{ d \sim \mathrm{N} ( d_0, \Sigma_d ), \quad \mu \sim
#'   \mathrm{N} ( \mu_0, \Sigma_\mu )}{d ~ N(d_0, \Sigma_d), \mu ~ N(\mu_0,
#'   \Sigma_\mu)}}
#'   \item{Likelihood:}{\eqn{y|d,\mu,\tau^2 \sim \mathrm{N} ( L\delta + M\mu, V
#'   )}{y|d,\mu,\tau^2 ~ N(L\delta + M\mu, V)}}
#'   \item{RE model:}{\eqn{\delta \sim \mathrm{N} ( Xd, \tau^2 )}{\delta ~ N(Xd,
#'   \tau^2)}} }
#'
#'   The between-study heterogeneity parameter \eqn{\tau^2} may be given any
#'   sensible prior. In the RE case, the threshold derivations make the
#'   approximation that \eqn{\tau^2} is fixed and known.
#'
#' @return An object of class \code{thresh}.
#' @seealso \code{\link{recon_vcov}}, \code{\link{thresh_forest}},
#'   \code{\link{thresh-class}}.
#' @export
#'
#' @examples
#' # Please see the vignette "Examples" for worked examples including use of
#' # this function, including more information on the brief code below.
#'
#' vignette("Examples", package = "nmathresh")
#'
#' ### Contrast level thresholds for Thrombolytic treatments NMA
#' K <- 6   # Number of treatments
#'
#' # Contrast design matrix is
#' X <- matrix(ncol = K-1, byrow = TRUE,
#'             c(1, 0, 0, 0, 0,
#'               0, 1, 0, 0, 0,
#'               0, 0, 1, 0, 0,
#'               0, 0, 0, 1, 0,
#'               0, -1, 1, 0, 0,
#'               0, -1, 0, 1, 0,
#'               0, -1, 0, 0, 1))
#'
#' # Reconstruct hypothetical likelihood covariance matrix using NNLS
#' lik.cov <- recon_vcov(Thrombo.post.cov, prior.prec = .0001, X = X)
#'
#' # Thresholds are then
#' thresh <- nma_thresh(mean.dk = Thrombo.post.summary$statistics[1:(K-1), "Mean"],
#'                      lhood = lik.cov,
#'                      post = Thrombo.post.cov,
#'                      nmatype = "fixed",
#'                      X = X,
#'                      opt.max = FALSE)
#'

nma_thresh <- function(mean.dk, lhood, post,
                       nmatype="fixed",
                       X=NULL,
                       mu.design=NULL, delta.design=diag(nrow=dim(lhood)),
                       opt.max=TRUE, trt.rank=1L, trt.bestn=1L,
                       trt.code=NULL, trt.sub=NULL) {


## Basic parameter checks --------------------------------------------------

  # Fixed or random effects
  tryCatch(isFE <- match.arg(tolower(nmatype), c("fixed","random")) == "fixed",
           error = function(err) {
             stop('nmatype should be one of "fixed" or "random"')
           })


  # Get number of data points N
  if (dim(lhood)[1] == dim(lhood)[2]) {
    N <- dim(lhood)[1]
    message("Likelihood for N = ", N, " data points.")
  } else stop("Likelihood covariance matrix lhood should be square.")

  # Get number of treatments
  K <- length(mean.dk) + 1
  message("Number of treatments is K = ", K, ".")

  # Check X matrix for FE model
  if (isFE) {
    if (is.null(X)) stop("Design matrix X must be provided for FE models.")
    else if (dim(X)[1] != N || dim(X)[2] != K-1) {
      stop("Design matrix X should be N x (K-1).")
    }
  }

  # Get number of extra parameters
  if (isFE || is.null(mu.design)) {
    m <- 0
  } else if (nrow(mu.design) != N) {
    stop("Number of rows in mu.design does not equal N.")
  } else {
    m <- ifelse(is.null(mu.design),0,ncol(mu.design))
    message("Number of extra parameters in m.design is m = ",m,".")
  }

  # Get number of delta parameters
  if (isFE) {
    n.delta <- 0
  } else if (nrow(delta.design) != N) {
    stop("Number of rows in delta.design does not equal N.")
  } else {
    n.delta <- ncol(delta.design)
    message("Number of delta parameters is n.delta = ",n.delta,".")
  }


  # Check posterior covariance matrix is n.delta+m+K-1 square
  if (nrow(post) != ncol(post)) {
    stop("Posterior covariance matrix should be square.")
  } else if (nrow(post) != n.delta + m + K-1) {
    if (isFE) stop("Posterior covariance matrix should be K-1 square.")
    else stop("Posterior covariance matrix should be n.delta+m+K-1 square.")
  }

  # Treatment rank
  if (length(trt.rank) > 1 | trt.rank != round(trt.rank)) {
    stop("trt.rank should be a single integer.")
  } else if (trt.rank < 1 | trt.rank > K) {
    stop("trt.rank should be between 1 and K.")
  }

  # Best n treatments, overrides trt.rank with a warning
  if (length(trt.bestn) > 1 | trt.bestn != round(trt.bestn)) {
    stop("trt.bestn should be a single integer.")
  } else if (trt.bestn < 1 | trt.bestn > K-1) {
    stop("trt.bestn should be between 1 and K-1.")
  } else if (trt.rank > 1 & trt.bestn > 1) {
    warning("trt.rank and trt.bestn specified - only trt.bestn will be used.")
    trt.rank <- 1L
  }

  # Note about recoded treatments
  if (is.null(trt.code)) {
    trt.code <- 1:K
  }
  else if(length(trt.code) != K) stop("trt.code should be of length K.")
  else {
    message("Using recoded treatments. Reference treatment is ", trt.code[1],
        ". Parameter vector is:\n",
        "\t", paste0("d[", trt.code[-1], "]", collapse=", ")
        )
  }

  # Treatment subset
  if (is.null(trt.sub)){
    trt.sub <- trt.code
  } else if(length(trt.sub)>K) stop("Length of trt.sub should be <= K.")
  else {
    message("Deriving thresholds on a subset of treatments:")
    message("\t", paste(trt.sub, collapse=", "))
  }

  trt.sub.internal <- which(trt.code %in% trt.sub)

  # Error if trt.rank > length(trt.sub)
  if(trt.rank > length(trt.sub)) {
    stop("trt.rank is larger than the length of trt.sub")
  }

  # Error if trt.bestn > length(trt.sub)
  if(trt.bestn >= length(trt.sub)) {
    stop("trt.bestn should be less than the length of trt.sub")
  }


## Pre-processing ----------------------------------------------------------

  # Create contrast "design" matrix
  D <- matrix(nrow=K*(K-1)/2, ncol=K-1)
  j <- 1
  for (i in 1:(K-1)) {
    if (i == 1) {
      D[j:(j+K-i-1),] <- diag(nrow=K-i)
    } else if (i == 2) {
       D[j:(j+K-i-1),] <- cbind(rep.int(-1,K-i), diag(nrow=K-i))
    } else {
      D[j:(j+K-i-1),] <- cbind(matrix(rep.int(0, (K-i)*(i-2)), ncol=i-2),
                               rep.int(-1, K-i),
                               diag(nrow=K-i))
    }
    j <- j+K-i
  }

  # Create vector of contrasts d_ab
  contr <- as.vector(D %*% mean.dk)

  # Augmented version of D to include a left column for d_1
  D.aug <- cbind(c(rep(-1, K-1), rep(0, (K-2)*(K-1)/2)), D)


## Derive influence matrix -------------------------------------------------

  ## FE models

  # If the likelihood covariance matrix was reconstructed to use in a
  # contrast-level analysis, e.g. using recon_vcov, there might be infinite
  # variances. We'll check that lhood is diagonal, and then handle these
  # correctly.

  if (isFE) {

    if (!Matrix::isDiagonal(lhood)) {
      inflmat <- post %*% crossprod(X, solve(lhood))
    } else {
      inflmat <- post %*% crossprod(X, diag(1/diag(lhood)))
    }

  } else {

  ## RE models
    Bstar <- post[1:(K-1), K:(K-1+n.delta)] # Posterior cov submatrix B*
    if (m > 0) {
      Cstar <- post[1:(K-1), (K+n.delta):(K-1+n.delta+m)] # Posterior cov submatrix C*
      inflmat <-
        (tcrossprod(Bstar, delta.design) + tcrossprod(Cstar, mu.design)) %*% solve(lhood)
    } else {
      inflmat <- Bstar %*% crossprod(delta.design, solve(lhood))
    }

  }


## Derive solution matrix U -------------------------------------------------

  threshmat <- sweep(1 / (D %*% inflmat),1, -contr, "*")

  # Add rownames
  contr.ab <- d_i2ab(1:(K*(K-1)/2), K)
  rownames(threshmat) <- paste0("d[", contr.ab$a, ",", contr.ab$b, "]")

  # Now we only need to look at contrasts involving the optimal treatment k*
  # Updated to handle trt.rank, to pick out other ranked treatments than the
  # optimal treatment k* in first place.
  # Updated to handle trt.sub, only look for k* in a subset of treatments.

  # Ignore treatments not in trt.sub
  mean.dk.subNA <- mean.dk
  mean.dk.subNA[!(1:(K-1) %in% (trt.sub.internal - 1))] <- NA

  # Get optimal treatment (or set of treatments)
  if (trt.bestn > 1) {
    kstar <- order(c(0, mean.dk.subNA), decreasing=opt.max)[1:trt.bestn]
  } else {
    kstar <- order(c(0, mean.dk.subNA), decreasing=opt.max)[trt.rank]
  }

  if (trt.rank == 1 & trt.bestn == 1) {
    message("Current optimal treatment is k* = ", trt.code[kstar], ".")
  } else if (trt.bestn == 1) {
    message("Current rank ", trt.rank, " treatment is k = ", trt.code[kstar], ".")
  } else {
    message("Current best ", trt.bestn, " treatments are k = ",
            paste(trt.code[kstar], collapse = ", "), ".")
  }

  # Pick out contrasts involving the optimal treatment(s)
  # For trt.bestn > 1, include contrasts between treatments in the set k* --
  # they will not affect the threshold, but they will inform the ranking.
  contr.kstar <- which(
    (contr.ab$a %in% kstar & contr.ab$b %in% trt.sub.internal) |
    (contr.ab$b %in% kstar & contr.ab$a %in% trt.sub.internal) )

  # So we look in the corresponding rows of the threshold matrix
  threshmat.kstar <- threshmat[contr.kstar, , drop = FALSE]


## Derive thresholds -------------------------------------------------------

    thresholds <- as.data.frame(
      do.call(rbind,
              apply(threshmat.kstar, 2,
                    get.int, kstar, trt.code, trt.sub
                    )
              )
      )


## Return thresh object ----------------------------------------------------
  return(structure(
    list(thresholds = thresholds,
         U = threshmat,
         Ukstar = threshmat.kstar,
         H = inflmat,
         kstar = trt.code[kstar],
         call = list(
           mean.dk = mean.dk,
           lhood = lhood,
           post = post,
           nmatype = nmatype,
           X = X,
           mu.design = mu.design,
           delta.design = delta.design,
           opt.max = opt.max,
           trt.rank = trt.rank,
           trt.bestn = trt.bestn,
           trt.code = trt.code,
           trt.sub = trt.sub
           )
         ),
    class="thresh")
    )

}



## Function get.int to return thresholds from U ----------------------------

# Return the positive and negative thresholds for each observation
# Define a function to do this which we can apply over the columns of U (observations)

#' Get thresholds from U matrix
#'
#' This function is intended for internal use only, and is called by
#' \code{nma_thresh} automatically.
#'
#' @param x Column of \eqn{U} matrix, containing all possible threshold
#'   solutions for a data point.
#' @param kstar Base-case optimal treatment.
#' @param trt.code Vector of (possibly recoded) treatments. See
#'   \code{nma_thresh} parameter of the same name.
#' @param trt.sub Vector of treatment indices, subset of trt.code, to consider.
#'   See \code{nma_thresh} parameter of the same name.
#'
#' @return Data frame of thresholds and new optimal treatments with columns
#'   \code{lo}, \code{lo.newkstar}, \code{hi}, and \code{hi.newkstar}.
#' @export
#'
get.int <- function(x, kstar, trt.code, trt.sub) {

  # Get internal recoded treatment list (possibly subset)
  trt.sub.internal <- which(trt.code %in% trt.sub)

  K <- length(trt.code)

  # Which rows have been passed in from a column of Ukstar?
  d_ab <- d_i2ab(1:(K*(K-1)/2), K)
  d_ab <- d_ab[(d_ab$a %in% kstar & d_ab$b %in% trt.sub.internal) |
                 (d_ab$b %in% kstar & d_ab$a %in% trt.sub.internal) , ]

  # Which contrasts are for thresholds (others are for switches within k*, if trt.bestn > 1)
  threshctr <- xor(d_ab$a %in% kstar, d_ab$b %in% kstar)

  xsub <- x[threshctr]

  # If both thresholds are infinite
  if (all(is.infinite(xsub))) {
    hi <- Inf
    lo <- -Inf
    hi.newkstar <- lo.newkstar <- NA_character_

    # If lower threshold is infinite
  } else if (all(xsub[!is.infinite(xsub)] > 0)) {
    hi <- min(xsub[!is.infinite(xsub)])
    hi.newkstar <- fnewkstar(x[x > 0 & !is.infinite(x)], d_ab[x > 0 & !is.infinite(x), ], hi, kstar, trt.code, trt.sub)

    lo <- -Inf
    lo.newkstar <- NA_character_

    # If upper threshold is infinite
  } else if (all(xsub[!is.infinite(xsub)] < 0)) {
    hi <- Inf
    hi.newkstar <- NA_character_
    lo <- max(xsub[!is.infinite(xsub)])
    lo.newkstar <- fnewkstar(x[x < 0 & !is.infinite(x)], d_ab[x < 0 & !is.infinite(x), ], lo, kstar, trt.code, trt.sub)

    # If neither threshold is infinite
  } else {
    hi <- min(xsub[xsub > 0 & !is.infinite(xsub)])
    hi.newkstar <-
      fnewkstar(x[x > 0 & !is.infinite(x)], d_ab[x > 0 & !is.infinite(x), ], hi, kstar, trt.code, trt.sub)
    lo <- max(xsub[xsub < 0 & !is.infinite(xsub)])
    lo.newkstar <-
      fnewkstar(x[x < 0 & !is.infinite(x)], d_ab[x < 0 & !is.infinite(x), ], lo, kstar, trt.code, trt.sub)
  }

  # To report new k* when trt.bestn >1, must translate into character strings
  tlo.newkstar <- paste0(lo.newkstar, collapse = ", ")
  thi.newkstar <- paste0(hi.newkstar, collapse = ", ")

  return(data.frame(lo=lo, lo.newkstar=tlo.newkstar, hi=hi, hi.newkstar=thi.newkstar))
}

# Internal function to calculate new k* from a set of all +ve or all -ve solutions
fnewkstar <- function(xsub, d_absub, thr, kstar, trt.code, trt.sub){

  newkstar <- kstar
  xsub <- abs(xsub)
  thr <- abs(thr)

  # Only concerned with elements up to the threshold
  ord <- order(xsub[xsub <= thr])
  xsub.lthr <- xsub[xsub <= thr][ord]
  d_absub.lthr <- d_absub[xsub <= thr, ][ord, ]

  nx <- length(xsub.lthr)
  for (i in 1:nx) {
    # Step along the solution vector, switching treatment ranks as we go
      temp <- newkstar
      temp[newkstar == d_absub.lthr[i, "a"]] <- d_absub.lthr[i, "b"]
      temp[newkstar == d_absub.lthr[i, "b"]] <- d_absub.lthr[i, "a"]
      newkstar <- temp
  }

  # Translate recoded treatments
  newkstar <- trt.code[newkstar]

  return(newkstar)
}
