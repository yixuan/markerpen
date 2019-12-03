# Make row sum to 1
scale_row = function(x, rsum = NULL)
{
    if(is.null(rsum))
        rsum = rowSums(x)
    x / rsum
}

# Make column sum to 1
scale_col = function(x, csum = NULL)
{
    if(is.null(csum))
        csum = colSums(x)
    sweep(x, 2, csum, "/")
}

Deconv = R6Class("Deconv",
    public = list(
        initialize = function(dat_unnorm, markers, normalize = TRUE)
        {
            mat_unnorm = t(as.matrix(dat_unnorm %>% select(-gene_name)))
            # Normalize each observation if needed
            private$dat = if(normalize) scale_row(mat_unnorm) else mat_unnorm
            private$gene_name = dat_unnorm$gene_name

            # In case `markers` contains genes not in the data
            for(i in seq_along(markers))
            {
                markers[[i]] = intersect(markers[[i]], private$gene_name)
            }
            private$markers = markers

            ind = match(unlist(markers), private$gene_name)
            private$mat_exp = private$dat[, ind]
        },

        view_exp = function()
        {
            view_matrix(private$mat_exp)
        },

        view_F = function(from = 1, to = 20, ...)
        {
            nc = length(private$markers)
            view_matrix(t(private$mat_F[from:to, ]), legend_title = "Fraction", ...) +
                scale_y_continuous("", breaks = 1:nc,
                                   labels = substr(names(private$markers), 1, 3))
        },

        view_S = function()
        {
            print(private$mat_S)
        },

        view_G = function(...)
        {
            nc = length(private$markers)
            view_evec(private$mat_G, ...) +
                scale_y_reverse("", breaks = 1:nc,
                                labels = substr(names(private$markers), 1, 3))
        },

        view_SG = function(log = FALSE, ...)
        {
            sg = private$mat_G %*% Diagonal(x = private$mat_S)
            if(log)
                sg = log(1e6 * sg + 1)
            nc = length(private$markers)
            view_evec(sg, ...) +
                scale_y_reverse("", breaks = 1:nc,
                                labels = substr(names(private$markers), 1, 3))
        },

        get_F = function()
        {
            private$mat_F
        },

        get_S = function()
        {
            private$mat_S
        },

        get_G = function()
        {
            private$mat_G
        }
    ),

    private = list(
        dat       = NULL,
        gene_name = NULL,
        markers   = NULL,
        mat_exp   = NULL,

        mat_W = NULL,
        mat_F = NULL,
        mat_S = NULL,
        mat_G = NULL
    )
)

# Remove genes with expression level > 90% quantile or < 10% quantile
Deconv$set("public", "remove_extreme_genes", function(l = 0.1, u = 0.9) {
    # Mean expression level for each gene
    mean_exp = colMeans(private$mat_exp)
    qu = quantile(mean_exp, u)
    ql = quantile(mean_exp, l)
    ind = (mean_exp <= qu) & (mean_exp > ql)
    private$mat_exp = private$mat_exp[, ind]

    ind = relist(ind, private$markers)
    private$markers = mapply(`[`, private$markers, ind)

})

# Construct the (initial) signature matrix G
# Optionally select a subset of genes from each cell type
Deconv$set("public", "init_params", function(F0 = NULL, G0 = NULL, select = NULL) {
    X = private$mat_exp
    scl = sqrt(colSums(X^2))

    if(is.null(G0))
    {
        coefs = lapply(private$markers, function(marker) {
            ind = match(marker, private$gene_name)
            dat_marker = private$dat[, ind]
            coef = eigs_sym(cov(dat_marker), 1)$vectors
            rownames(coef) = marker
            abs(coef)
        })
        mat = bdiag(coefs)
        G0 = scale_col(as.matrix(mat))
    } else {
        G0 = scale_col(pmax(G0, 0))
    }

    if(is.null(F0))
    {
        # Xt = t(X)
        # rownames(Xt) = unlist(private$markers)
        # F0 = t(EstimateWeight(Xt, private$markers, method = "LM")$weight)
        # F0 = F0 / rowSums(F0)

        F0 = t(fcnnls(G0 / scl, t(X) / scl, pseudo = FALSE)$x)
        F0 = scale_row(pmax(F0, 0))
    } else {
        F0 = scale_row(pmax(F0, 0))
    }

    S0 = s_update(X, F0, G0)

    private$mat_W = as.matrix(1 - bdiag(lapply(private$markers, function(x) rep(1, length(x)))))
    private$mat_F = F0
    private$mat_S = S0
    private$mat_G = G0

    # Selection based on one-step correction
    if(!is.null(select))
    {
        prop = select

        # X ~= F * S * G
        G0 = g_update(X, F0, S0, 0, private$mat_W)
        rownames(G0) = unlist(private$markers)

        # Select `prop` proportion marker genes in each cell type for deconvolution
        nc = length(private$markers)
        markers_new = private$markers
        for(i in 1:nc)
        {
            Gsub = G0[markers_new[[i]], ]
            snr = Gsub[, i] / (mean(Gsub[, i]) + rowSums(Gsub[, -i]))
            # snr = Gsub[, i] / rowSums(Gsub)
            p0 = ceiling(length(markers_new[[i]]) * prop)
            markers_new[[i]] = names(sort(snr, decreasing = TRUE)[1:p0])
        }
        private$markers = markers_new

        ind = match(unlist(private$markers), private$gene_name)
        private$mat_exp = private$dat[, ind]

        # Recursively call this function to re-initialize (F, S, G)
        self$init_params(F0 = NULL, G0 = NULL, select = NULL)
    }
})

# Factorization
Deconv$set("public", "factorize", function(lambda = 0, niter = 100) {
    res = deconv_fsg(private$mat_exp, private$mat_W,
                     private$mat_F, private$mat_S, private$mat_G,
                     lambda, niter,
                     total_min = 0.95, total_max = 0.99, each_min = 0.01)
    private$mat_F = res$F
    private$mat_S = res$S
    private$mat_G = res$G

    invisible(res$loss)
})

Deconv$set("public", "factorize_bagging", function(lambda = 0, niter = 100, bagging_p = NULL, bagging_n = 10) {
    pm = min(sapply(private$markers, length))
    if(is.null(bagging_p))
        bagging_p = ceiling(pm / 2)
    bagging_p = min(bagging_p, pm)

    X = private$mat_exp
    W = private$mat_W

    res = vector("list", bagging_n)
    for(i in 1:bagging_n)
    {
        cat(sprintf("\n=============== Bagging %d ===============\n\n", i))

        subind = lapply(private$markers, function(x) {
            marker_sub = sample(x, bagging_p)
            ind = match(marker_sub, unlist(private$markers))
        })

        ind = unlist(subind)
        Xsub = X[, ind]
        Wsub = W[ind, ]
        F0 = private$mat_F
        G0 = private$mat_G[ind, ]
        G0 = sweep(G0, 2, colSums(G0), "/")
        S0 = s_update(Xsub, F0, G0)

        res[[i]] = deconv_fsg(Xsub, Wsub, F0, S0, G0, lambda, niter)
    }

    Fmat = Reduce("+", lapply(res, function(x) x$F)) / bagging_n

    private$mat_F = Fmat
    invisible(res)
})

