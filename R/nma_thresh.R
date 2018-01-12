#' Calculate thresholds and invariant intervals
#'
#' This function calculates decision-invariant bias-adjustment thresholds and
#' intervals for Bayesian network meta-analysis, as described by Phillippo
#' \emph{et al.} (2017). Thresholds are derived from the joint
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
#' @param trt.code Treatment codings of the reference trt and in the parameter
#'   vector \eqn{d_k}. Use if treatments re-labelled or re-ordered. Default is
#'   equivalent to 1:K.
#' @param trt.sub Only look at thresholds in this subset of treatments in
#'   trt.code, e.g. if some are excluded from the ranking. Default is equivalent
#'   to 1:K.
#' @param mcid Minimal clinically important difference for the decision (when
#'   \code{mcid.type = 'decision'}) or for changing the decision (when
#'   \code{mcid.type = 'change'}). Defaults to 0, use the maximal efficacy
#'   decision rule.
#' @param mcid.type Default \code{'decision'}, the decision rule is based on
#'   MCID (see details). Otherwise \code{'change'}, use the maximum efficacy
#'   rule, but only consider changing the decision when the alternative
#'   treatment becomes more effective than the base case by \code{mcid} or more.
#'
#' @details This function provides bias-adjustment threshold analysis for both
#'   fixed and random effects NMA models, as described by Phillippo \emph{et
#'   al.} (2017). Parameters \code{mean.dk}, \code{lhood}, and
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
#' @section Decision rules:
#'
#'   The default decision rule is maximal efficacy; the optimal treatment is
#'   \eqn{ k^* = \mathrm{argmax}_k \mathbb{E}(d_{k})}{k* = argmax(E(d_k))}.
#'
#'   When \eqn{\epsilon} = \code{mcid} is greater than zero and
#'   \code{mcid.type = 'decision'}, the decision rule is no longer for a single
#'   best treatment, but is based on minimal clinically important difference. A
#'   treatment is in the optimal set if \eqn{\mathbb{E}(d_k) \ge
#'   \epsilon}{E(d_k) \ge \epsilon} and \eqn{\max_a \mathbb{E}(d_a) -
#'   \mathbb{E}(d_k) \le \epsilon}{max E(d_a) - E(d_k) \le \epsilon}.
#'
#'   When \code{mcid.type = 'change'}, the maximal efficacy rule is used, but
#'   thresholds are found for when a new treatment is better than the base-case
#'   optimal by at least \code{mcid}.
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
                       opt.max=TRUE, trt.rank=1, trt.code=NULL, trt.sub=NULL,
                       mcid=0, mcid.type='decision') {


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
    stop("trt.rank should be between 1 and K (number of trts).")
  }

  # Note about recoded treatments
  if (is.null(trt.code)) {
    trt.code <- 1:K
  }
  else if (length(trt.code) != K) stop("trt.code should be of length K.")
  else {
    message("Using recoded treatments. Reference treatment is ", trt.code[1],
        ". Parameter vector is:\n",
        "\t", paste0("d[", trt.code[-1], "]", collapse=", ")
        )
  }

  # Treatment subset
  if (is.null(trt.sub)){
    trt.sub <- trt.code
  } else if (length(trt.sub)>K) stop("Length of trt.sub should be <= K.")
  else {
    message("Deriving thresholds on a subset of treatments:")
    message("\t", paste(trt.sub, collapse=", "))
  }

  trt.sub.internal <- which(trt.code %in% trt.sub)

  # Error if trt.rank > length(trt.sub)
  if (trt.rank > length(trt.sub)) {
    stop("trt.rank is larger than the length of trt.sub")
  }

  # mcid should be a single non-negative numeric value
  if (!is.numeric(mcid) | length(mcid) != 1 | mcid < 0) {
    stop("mcid should be a single non-negative numeric value")
  }

  # Check mcid.type
  mcid.type <- match.arg(mcid.type, c('decision', 'change'))

  # Can't use mcid decision rule and treatment ranks
  if (mcid > 0 & mcid.type == 'decision' & trt.rank > 1) {
    stop("Can't use mcid decision rule and trt.rank at the same time.")
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

  if (mcid > 0 & mcid.type == 'change') {
    threshmat <- sweep(1 / (D %*% inflmat), 1, -contr - sign(contr)*mcid, "*")

    ## -- For mcid.type = "change" --
    # For mcid > 0, if a contrast is negative, we want a new decision when the
    # contrast is > +mcid. If a contrast is positive, we want a new decision
    # when the contrast is < -mcid. In other words, the contrast has to be
    # overturned by an extra mcid.new amount.
  } else {
    threshmat <- sweep(1 / (D %*% inflmat), 1, -contr + sign(contr)*mcid, "*")

    ## -- For mcid.type = "decision" --
    # And also for standard maximal efficacy rule, when mcid = 0 anyway.
    # For mcid > 0, if a contrast is negative, we want to know when the
    # contrast is = -mcid. If a contrast is positive, we want to know when
    # the contrast is = +mcid
  }

  # Now we only need to look at contrasts involving the optimal treatment k*
  # Updated to handle trt.rank, to pick out other ranked treatments than the
  # optimal treatment k* in first place.
  # Updated to handle trt.sub, only look for k* in a subset of treatments.

  mean.dk.subNA <- mean.dk
  mean.dk.subNA[!(1:(K-1) %in% (trt.sub.internal - 1))] <- NA

  if (opt.max){
    kstar <- order(c(0, mean.dk.subNA), decreasing=TRUE)[trt.rank]
  } else if (!opt.max){
    kstar <- order(c(0, mean.dk.subNA), decreasing=FALSE)[trt.rank]
  }

  if (trt.rank == 1) {
    message("Current optimal treatment is k* = ", trt.code[kstar], ".")
  } else {
    message("Current rank ", trt.rank, " treatment is k = ", trt.code[kstar], ".")
  }

  # And these contrasts have non-zero elements in the contrast design matrix D
  if (kstar > 1) {
    contr.kstar <- which(D[,kstar-1] != 0)
  } else {
    contr.kstar <- 1:(K-1)
  }

  # So we look in the corresponding rows of the threshold matrix
  threshmat.kstar <- threshmat[contr.kstar,]


## Derive thresholds -------------------------------------------------------

  # Only look in rows which correspond to contrasts with treatments in trt.sub
  contr.trt.sub <- trt.sub.internal[-which(trt.sub.internal == kstar)] -
    (trt.sub.internal[-which(trt.sub.internal == kstar)] >= kstar)*1

    thresholds <- as.data.frame(
      do.call(rbind,
              apply(threshmat.kstar[contr.trt.sub, , drop = FALSE], 2,
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
           trt.code = trt.code,
           trt.sub = trt.sub,
           mcid.new = mcid.new
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

  trt.sub.internal <- which(trt.code %in% trt.sub)

  # If both thresholds are infinite
  if (all(is.infinite(x))) {
    hi <- Inf
    lo <- -Inf
    hi.newkstar <- lo.newkstar <- NA

  # If lower threshold is infinite
  } else if (all(x[!is.infinite(x)] > 0)) {
    hi <- min(x[!is.infinite(x)])
    hi.newkstar <-
      trt.code[trt.sub.internal[which(x==hi) + (which(x==hi) >= which(trt.sub.internal==kstar))*1]]
    lo <- -Inf
    lo.newkstar <- NA

  # If upper threshold is infinite
  } else if (all(x[!is.infinite(x)] < 0)) {
    hi <- Inf
    hi.newkstar <- NA
    lo <- max(x[!is.infinite(x)])
    lo.newkstar <-
      trt.code[trt.sub.internal[which(x==lo) + (which(x==lo) >= which(trt.sub.internal==kstar))*1]]

  # If neither threshold is infinite
  } else {
    hi <- min(x[x>0 & !is.infinite(x)])
    hi.newkstar <-
      trt.code[trt.sub.internal[which(x==hi) + (which(x==hi) >= which(trt.sub.internal==kstar))*1]]
    lo <- max(x[x<0 & !is.infinite(x)])
    lo.newkstar <-
      trt.code[trt.sub.internal[which(x==lo) + (which(x==lo) >= which(trt.sub.internal==kstar))*1]]
  }
  return(data.frame(lo=lo,lo.newkstar=lo.newkstar,hi=hi,hi.newkstar=hi.newkstar))
}
