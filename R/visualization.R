##' Visualization of Matrices
##'
##' \code{view_matrix()} visualizes a general matrix by mapping its elements to
##' colors. \code{view_evec()} visualizes a set of eigenvectors.
##'
##' @param mat,evec     The matrix to be visualized.
##' @param xlab         X-axis label.
##' @param ylab         Y-axis label.
##' @param legend_title Title of the colorbar legend.
##' @param asp          Aspect ratio of the plot.
##' @param bar_height   Height of the colorbar.
##' @param font_size    Base font size for the plot.
##'
##' @rdname visualization
##' @author Yixuan Qiu \url{https://statr.me}
##'
##' @examples
##' set.seed(123)
##' x = matrix(rnorm(200), 10, 20)
##' view_matrix(x)
##'
##' sigma1 = matrix(0.8, 20, 20) + 0.2 * diag(20)
##' sigma2 = matrix(0.6, 50, 50) + 0.4 * diag(50)
##' s1 = stats::rWishart(1, 100, sigma1)[, , 1]
##' s2 = stats::rWishart(1, 100, sigma2)[, , 1]
##' s = as.matrix(Matrix::bdiag(s1, s2))
##' view_matrix(s)
##'
##' v = eigen(s, symmetric = TRUE)$vectors[, 1:5]
##' view_evec(v)

# Visualization of a matrix by coloring its coefficients
view_matrix = function(mat, legend_title = "Coefficient", bar_height = 10, font_size = 20)
{
    mat = as.matrix(mat)
    lo = min(mat)
    hi = max(mat)
    # All data in the range [-r, r]
    r = max(abs(c(lo, hi)))

    # ggplot2 format
    gdat = data.frame(
        x = as.integer(col(mat)),
        y = as.integer(row(mat)),
        z = as.numeric(mat)
    )

    # Axis breaks
    breaks_x = pretty(gdat$x, n = 10)
    # Start from 1, not 0
    breaks_x[breaks_x == 0] = 1
    breaks_x = unique(breaks_x)

    breaks_y = pretty(gdat$y, n = 10)
    # Start from 1, not 0
    breaks_y[breaks_y == 0] = 1
    breaks_y = unique(breaks_y)

    # Map the color spectrum to [-r, r]
    ngrid = 1001
    col_pal = colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                                 "#FDDBC7", "#FFFFFF", "#D1E5F0",
                                 "#92C5DE", "#4393C3", "#2166AC", "#053061"))(ngrid)
    col_val = seq(-r, r, length.out = ngrid)
    lo_ind = findInterval(lo, col_val)
    hi_ind = findInterval(hi, col_val)
    colors = col_pal[lo_ind:hi_ind]

    ggplot(gdat, aes(x = x, y = y, fill = z)) +
        geom_tile() +
        scale_x_continuous("", breaks = breaks_x, expand = c(0, 0)) +
        scale_y_reverse("", breaks = breaks_y, expand = c(0, 0)) +
        scale_fill_gradientn(legend_title, colors = colors) +
        guides(fill = guide_colorbar(barheight = bar_height)) +
        coord_fixed() +
        theme_bw(base_size = font_size) +
        theme(axis.title = element_blank())
}

# Visualization of eigenvectors
view_evec = function(
    evecs,
    xlab = "Index of Variables", ylab = "Index of PCs", legend_title = "Factor\nLoading",
    asp = 0.2, bar_height = 6, font_size = 20
)
{
    v = as.matrix(evecs)
    lo = min(v)
    hi = max(v)
    # All values in the range [-r, r]
    r = max(abs(c(lo, hi)))

    # ggplot2 format
    gdat = data.frame(
        x = as.integer(row(v)),
        y = as.integer(col(v)),
        z = as.numeric(v)
    )

    # Map the color spectrum to [-r, r]
    ngrid = 1001
    col_pal = colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                                 "#FDDBC7", "#FFFFFF", "#D1E5F0",
                                 "#92C5DE", "#4393C3", "#2166AC", "#053061"))(ngrid)
    col_val = seq(-r, r, length.out = ngrid)
    lo_ind = findInterval(lo, col_val)
    hi_ind = findInterval(hi, col_val)
    colors = col_pal[lo_ind:hi_ind]

    ggplot(gdat, aes(x = x, y = y, fill = z)) +
        geom_tile() +
        scale_x_continuous(xlab, expand = c(0, 0)) +
        scale_y_reverse(ylab, expand = c(0, 0)) +
        scale_fill_gradientn(legend_title, colors = colors) +
        guides(fill = guide_colorbar(barheight = bar_height)) +
        theme_bw(base_size = font_size) +
        theme(aspect.ratio = asp)
}
