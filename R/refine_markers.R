##' Marker Gene Selection via Penalized Principal Component Analysis
##'
##' This function refines a prior marker gene list by combining information from bulk
##' tissue data, based on the penalized principal component analysis. The current
##' implementation computes on one cell type at a time. To get marker genes for
##' multiple cell types, call this function iteratively.
##'
##' @param mat_exp The gene expression matrix in the original scale
##'                (not logarithm-transformed), with rows standing for observations and
##'                columns for genes. The matrix should include gene names as column names.
##' @param range   A character vector of gene names, representing the range of genes in which
##'                markers are sought.
##' @param markers A character vector of gene names giving the prior marker gene list.
##' @param lambda  A tuning parameter to control the number of selected marker genes. A larger
##'                value typically means a smaller number of genes.
##' @param w       Tuning parameter to control the weight on prior information.
##'                Larger \eqn{w} means genes not in the prior list are less likely
##'                to be selected as markers.
##' @param thresh  Below this threshold small factor loadings are treated as zeros.
##' @param alpha   Step size of the optimization algorithm.
##' @param maxit   Maximum number of iterations.
##' @param eps     Tolerance parameter for convergence.
##' @param verbose Level of verbosity.
refine_markers = function(
    mat_exp, range, markers, lambda, w = 1.5, thresh = 0.001,
    alpha = 0.01, maxit = 1000, eps = 1e-4, verbose = 0)
{
    # Restrict to the submatrix specified by range
    range = intersect(colnames(mat_exp), range)
    mat_range = mat_exp[, range]
    cor_range = cor(mat_range)
    markers = intersect(markers, range)
    pub_id = match(markers, range)

    spca = pca_pen(
        cor_range, pub_id, lambda = lambda, w = w,
        alpha = alpha, maxit = maxit, eps = eps, verbose = verbose
    )
    coefs = spca$evecs
    markers_coef = coefs[abs(coefs) > thresh]
    new_markers = range[abs(coefs) > thresh]
    ord = order(markers_coef, decreasing = TRUE)

    if(verbose > 0)
        cat(sprintf("original = %d, out = %d, in = %d, new = %d\n",
                    length(markers),
                    length(setdiff(markers, new_markers)),
                    length(setdiff(new_markers, markers)),
                    length(new_markers)))

    list(spca = spca, markers = new_markers[ord], markers_coef = markers_coef[ord])
}
