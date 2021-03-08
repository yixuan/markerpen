##' Post-processing Selected Marker Genes
##'
##' This function reorders the selected marker genes using information of the sample
##' correlation matrix.
##'
##' @param corr    The sample correlation matrix, whose row and column names are gene names.
##' @param markers A list of marker genes. Each component of the list is a vector of marker
##'                gene names corresponding to a cell type. All the gene names in this list
##'                must appear in the row/column names of \code{corr}.
##'
##' @return A list that has the same structure as the input \code{markers} argument, with
##'         the elements in each component reordered. See the example in
##'         \code{\link{refine_markers}()}.
sort_markers = function(corr, markers)
{
    gene_names = rownames(corr)
    markers_flat = unlist(markers)
    if(!all(markers_flat %in% gene_names))
        stop("all gene names in 'markers' must appear in the row/column names of 'corr'")

    markers_ord = markers
    for(i in 1:length(markers))
    {
        markeri = markers[[i]]
        ngene = length(markeri)
        score = numeric(ngene)
        for(j in 1:ngene)
        {
            corj = corr[markeri[j], ]
            in_gr = markeri[-j]
            out_gr = setdiff(markers_flat, markeri)
            score[j] = mean(abs(corj[in_gr])) / mean(abs(corj[out_gr]))
        }
        markers_ord[[i]] = markeri[order(score, decreasing = TRUE)]
    }
    markers_ord
}
