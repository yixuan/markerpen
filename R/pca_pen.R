# The old R implementation, now superseded by the C++ version below
# Kept here for record and for understanding the algorithm
#
# pca_pen = function(S, gr, lambda, gamma = 1.5, alpha = 0.01, maxit = 10, eps = 1e-4, verbose = 0)
# {
#     p = nrow(S)
#
#     # Ssub = S[gr, gr]
#     # e = RSpectra::eigs_sym(Ssub, 1)
#     # x = matrix(0, p, p)
#     # x[gr, gr] = tcrossprod(e$vectors)
#     # z1 = z2 = x
#
#     e = RSpectra::eigs_sym(S, 1)
#     x = z1 = z2 = tcrossprod(e$vectors)
#
#     # Time for Fantope projection
#     time_f = c()
#     # Time for penalty
#     time_p = c()
#     # Estimation error
#     err = c()
#
#     penmat = matrix(gamma * lambda, p, p)
#     penmat[gr, gr] = lambda
#     penmat[-gr, -gr] = gamma^2 * lambda
#
#     for(i in 1:maxit)
#     {
#         newx = (z1 + z2) / 2
#         resid = norm(newx - x, type = "F")
#         x = newx
#
#         t1 = Sys.time()
#         newz1 = z1 - x + prox_fantope(S - penmat, 2 * x - z1, alpha, 1, inc = 30)
#         t2 = Sys.time()
#         newz2 = z2 - x + pmax(2 * x - z2, 0)
#         t3 = Sys.time()
#
#         resid1 = norm(newz1 - z1, type = "F")
#         resid2 = norm(newz2 - z2, type = "F")
#         z1 = newz1
#         z2 = newz2
#         if(verbose > 0)
#             cat(sprintf("iter = %d, resid = %f, resid1 = %f, resid2 = %f\n", i, resid, resid1, resid2))
#
#         time_f = c(time_f, t2 - t1)
#         time_p = c(time_p, t3 - t2)
#         err = c(err, resid)
#         if(max(resid, resid1, resid2) < eps)
#             break
#     }
#     list(proj = x, z1 = z1, z2 = z2, time_f = time_f, time_p = time_p, time_t = time_f + time_p, error = err)
# }

##' Penalized Principal Component Analysis for Marker Gene Selection
##'
##' This function solves the optimization problem
##' \deqn{\min\quad-tr(SX) + \lambda p(X),}{min -tr(SX) + \lambda * p(X),}
##' \deqn{s.t.\quad O\preceq X \preceq I, X \ge 0, \text{and} tr(X)=1,}{s.t. O ≼ X ≼ I, X \ge 0, and tr(X) = 1,}
##' where \eqn{O\preceq X \preceq I}{O ≼ X ≼ I} means all eigenvalues of \eqn{X} are
##' between 0 and 1, \eqn{X \ge 0} means all elements of \eqn{X} are nonnegative,
##' and \eqn{p(X)} is a penalty function defined in the article
##' (see the \strong{References} section).
##'
##' @param S       The sample correlation matrix of gene expression.
##' @param gr      Indices of genes that are treated as markers in the prior information.
##' @param lambda  Tuning parameter to control the sparsity of eigenvectors.
##' @param w       Tuning parameter to control the weight on prior information.
##'                Larger \eqn{w} means genes not in the prior list are less likely
##'                to be selected as markers.
##' @param alpha   Step size of the optimization algorithm.
##' @param maxit   Maximum number of iterations.
##' @param eps     Tolerance parameter for convergence.
##' @param verbose Level of verbosity.
##'
##' @examples set.seed(123)
##' n = 200  # Sample size
##' p = 500  # Number of genes
##' s = 50   # Number of true signals
##'
##' # The first s genes are true markers, and others are noise
##' Sigma = matrix(0, p, p)
##' Sigma[1:s, 1:s] = 0.9
##' diag(Sigma) = 1
##'
##' # Simulate data from the covariance matrix
##' x = matrix(rnorm(n * p), n) %*% chol(Sigma)
##'
##' # Sample correlation matrix
##' S = cor(x)
##'
##' # Indices of prior marker genes
##' # Note that we have omitted 10 true markers, and included 10 false markers
##' gr = c(1:(s - 10), (s + 11):(s + 20))
##'
##' # Run the algorithm
##' res = pca_pen(S, gr, lambda = 0.1, verbose = 1)
##'
##' # See if we can recover the true correlation structure
##' image(res$projection, asp = 1)
##'
##' @references Qiu, Y., Wang, J., Lei, J., & Roeder, K. (2020).
##' Identification of cell-type-specific marker genes from co-expression patterns in tissue samples.
pca_pen = function(S, gr, lambda, w = 1.5, alpha = 0.01, maxit = 1000, eps = 1e-4, verbose = 0)
{
    p = nrow(S)

    # Ssub = S[gr, gr]
    # e = RSpectra::eigs_sym(Ssub, 1)
    # x = matrix(0, p, p)
    # x[gr, gr] = tcrossprod(e$vectors)
    # z1 = z2 = x

    e = RSpectra::eigs_sym(S, 1)
    x0 = tcrossprod(e$vectors)
    pca_pen_(
        S, gr, x0, lambda, w, alpha, maxit,
        fan_maxinc = 100, fan_maxiter = 10, eps = eps, verbose = verbose
    )
}
