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
##'
##' @return A list containing the following components:
##' \describe{
##'   \item{spca}{The sparse PCA result as in \code{\link{pca_pen}()}.}
##'   \item{markers}{A character vector of selected markers genes.}
##'   \item{markers_coef}{The estimated factor loadings for the associated genes.}
##' }
##'
##' @examples # Data used in the vignette
##' load(system.file("examples", "gene_expr.RData", package = "markerpen"))
##' load(system.file("examples", "published_markers.RData", package = "markerpen"))
##' load(system.file("examples", "markers_range.RData", package = "markerpen"))
##'
##' # Get expression matrix - rows are observations, columns are genes
##' ind = match(rownames(dat), markerpen::gene_mapping$name)
##' ind = na.omit(ind)
##' ensembl = markerpen::gene_mapping$ensembl[ind]
##' mat_exp = t(dat[markerpen::gene_mapping$name[ind], ])
##' colnames(mat_exp) = ensembl
##'
##' # We compute the marker genes for two cell types here
##' # See the vignette for the full example
##'
##' # Markers for astrocytes
##' ast_re = refine_markers(mat_exp, markers_range$astrocytes, pub_markers$astrocytes,
##'                         lambda = 0.45, w = 1.5, maxit = 500, eps = 1e-3, verbose = 0)
##' # Remove selected markers from the expression matrix
##' mat_rest = mat_exp[, setdiff(colnames(mat_exp), ast_re$markers)]
##'
##' # Markers for oligodendrocytes
##' oli_re = refine_markers(mat_rest, markers_range$oligodendrocytes, pub_markers$oligodendrocytes,
##'                         lambda = 0.45, w = 1.5, maxit = 500, eps = 1e-3, verbose = 0)
##'
##' # Refined markers
##' markers_re = list(astrocytes       = ast_re$markers,
##'                   oligodendrocytes = oli_re$markers)
##'
##' # Post-process selected markers
##' # Pick the first 50 ordered markers
##' cor_markers = cor(mat_exp[, unlist(markers_re)])
##' markers_ord = sort_markers(cor_markers, markers_re)
##' markers_ord = lapply(markers_ord, head, n = 50)
##'
##' # Visualize the correlation matrix
##' image(cor(mat_exp[, unlist(markers_ord)]), asp = 1)
##'
##' @references Qiu, Y., Wang, J., Lei, J., & Roeder, K. (2020).
##' Identification of cell-type-specific marker genes from co-expression patterns in tissue samples.
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
