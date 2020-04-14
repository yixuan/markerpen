sort_markers = function(corr, markers)
{
    markers_flat = unlist(markers)
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
