# F-update
# min  ||X - F * S * G'||^2
# s.t. F >= 0, rowsum(F) = 1
#
# Subproblem
# x is i-th row of X
# f is i-th row of F
# A = S * G'
# min  ||x - f * A||^2 = ||x' - A' * f'||^2
# s.t. f >= 0, c <= sum(f) <= 1
#
# Transform to a quadratic programming problem
# min  0.5 * u' * D * u - d'u
# s.t. u >= 0, c <= 1' * u <= 1
# where D = A * A', d = A * x'
f_sub = function(D, d, total_min = 0.95, total_max = 0.99, each_min = 0.01)
{
    nc = length(d)
    constr = rbind(rep(1, nc), rep(-1, nc), diag(nc))
    rhs = c(total_min, -total_max, rep(each_min, nc))
    sol = solve.QP(Dmat = D, dvec = d, Amat = t(constr), bvec = rhs, meq = 0)
    sol$solution
}

# S is a diagonal matrix, and the argument here is its diagonal elements
# X [n x p], S [nc x nc], G [p x nc]
f_update = function(X, S, G, total_min = 0.95, total_max = 0.99, each_min = 0.01)
{
    n = nrow(X)
    nc = length(S)

    scl = sqrt(colSums(X^2))
    X = sweep(X, 2, scl, "/")
    G = G / scl

    # B = AA' = SG'GS
    B = crossprod(G) * outer(S, S)

    # C = AX' = SG'X'
    C = t(X %*% G) * S

    Fmat = matrix(0, n, nc)
    for(i in 1:n)
    {
        Fmat[i, ] = f_sub(B, C[, i], total_min, total_max, each_min)
    }

    pmax(Fmat, 0)
}

# Test
if(FALSE)
{
    set.seed(123)
    n = 10
    p = 30
    nc = 5

    Fmat = matrix(runif(n * nc), n, nc)
    Frowsums = rowSums(Fmat)
    Fmat = Fmat / Frowsums

    S = runif(nc, 0, 10)

    G = matrix(runif(p * nc), p, nc)
    Gcolsums = colSums(G)
    G = sweep(G, 2, Gcolsums, "/")

    X = Fmat %*% diag(S) %*% t(G)
    res = f_update(X, S, G)
    max(abs(res - Fmat))
}



# S-update
# min  ||X - F * S * G'||^2
# s.t. S >= 0, S is diagonal
#
# Transform to a nonnegative LS problem
# min  ||vec(X) - A * s||^2
# s.t. s >= 0
# where A = (F[, 1] * G[, 1]', ..., F[, nc] * G[, nc]')
#
# X [n x p], F [n x nc], G [p x nc]
s_update = function(X, Fmat, G)
{
    n = nrow(X)
    p = ncol(X)
    nc = ncol(Fmat)

    scl = sqrt(colSums(X^2))
    X = sweep(X, 2, scl, "/")
    G = G / scl

    A = matrix(0, n * p, nc)
    for(i in 1:nc)
    {
        A[, i] = as.numeric(outer(Fmat[, i], G[, i]))
    }

    as.numeric(fcnnls(A, as.numeric(X), pseudo = FALSE)$x)
}

# Test
if(FALSE)
{
    set.seed(123)
    n = 10
    p = 30
    nc = 5

    Fmat = matrix(runif(n * nc), n, nc)
    Frowsums = rowSums(Fmat)
    Fmat = Fmat / Frowsums

    S = runif(nc, 0, 10)

    G = matrix(runif(p * nc), p, nc)
    Gcolsums = colSums(G)
    G = sweep(G, 2, Gcolsums, "/")

    X = Fmat %*% diag(S) %*% t(G)
    res = s_update(X, Fmat, G)
    max(abs(res - S))
}



# G-update
# min  ||X - F * S * G'||^2 + lambda * sum(W * G^2)
# s.t. G >= 0, colsum(F) = 1
#
# Transform to a quadratic programming problem
# min  0.5 * u' * D * u - d'u
# s.t. u >= 0, A * u = 1
# where D = H'H (x) I + lambda * diag(vec(W)), d = vec(X'H), A = I (x) 1'
#
# X [n x p], F [n x nc], S [nc x nc]
g_update = function(X, Fmat, S, lambda, W)
{
    n = nrow(X)
    p = ncol(X)
    nc = ncol(Fmat)

    scl = sqrt(colSums(X^2))
    X = sweep(X, 2, scl, "/")
    P = Diagonal(x = 1 / scl)

    # H = FS
    # H'H = SF'FS, X'H = X'FS
    hh = crossprod(Fmat) * outer(S, S)
    xh = P %*% t(S * crossprod(Fmat, X))

    d = as.numeric(xh)
    D = kronecker(hh, P^2) + lambda * Diagonal(x = as.numeric(W))
    A = rbind(kronecker(Diagonal(n = nc), matrix(1, 1, p)),
              Diagonal(n = p * nc))
    lb = c(rep(1, nc), rep(0, p * nc))
    ub = c(rep(1, nc), rep(Inf, p * nc))
    sol = solve_osqp(P = D, q = -d, A = A, l = lb, u = ub,
                     pars = osqpSettings(polish = TRUE, verbose = FALSE))
    if(sol$info$status_val != 1)
        stop(sprintf("rosqp: %s", sol$info$status))

    matrix(pmax(sol$x, 0), p, nc)
}

# Test
if(FALSE)
{
    set.seed(123)
    n = 10
    p = 30
    nc = 5

    Fmat = matrix(runif(n * nc), n, nc)
    Frowsums = rowSums(Fmat)
    Fmat = Fmat / Frowsums

    S = runif(nc, 0, 10)

    G = matrix(runif(p * nc), p, nc)
    Gcolsums = colSums(G)
    G = sweep(G, 2, Gcolsums, "/")

    W = matrix(rbinom(p * nc, 1, 0.5), p, nc)
    X = Fmat %*% diag(S) %*% t(G)
    res = g_update(X, Fmat, S, 0, W)
    max(abs(res - G))

    lambda = 0.01
    loss1 = lambda * sum(W * G^2)
    res = g_update(X, Fmat, S, lambda, W)
    loss2 = norm(X - Fmat %*% diag(S) %*% t(res), type = "F")^2 +
        lambda * sum(W * res^2)
    loss1
    loss2
}



loss_main = function(X, Fmat, S, G)
{
    scl = sqrt(colSums(X^2))
    X = sweep(X, 2, scl, "/")

    est = Fmat %*% diag(S) %*% t(G / scl)
    norm(X - est, type = "F")^2
}

deconv_fsg = function(X, W, F0, S0, G0, lambda, niter = 100, total_min = 0.95, total_max = 0.99, each_min = 0.01)
{
    obj1 = loss_main(X, F0, S0, G0)
    obj2 = sum(W * G0^2)
    obj = obj1 + lambda * obj2
    cat(sprintf("initial loss: obj1 = %f, obj2 = %f, obj = %f\n", obj1, obj2, obj))

    # Outer iteration
    Fmat = F0
    S = S0
    G = G0
    loss = c()
    for(i in 1:niter)
    {
        # F-update
        Fmat = f_update(X, S, G, total_min, total_max, each_min)

        # S-update
        S = s_update(X, Fmat, G)

        # G-update
        G = g_update(X, Fmat, S, lambda, W)

        # Objective function
        obj1 = loss_main(X, Fmat, S, G)
        obj2 = sum(W * G^2)
        obj = obj1 + lambda * obj2
        loss = c(loss, obj)
        cat(sprintf("i = %d, obj1 = %f, obj2 = %f, obj = %f\n", i, obj1, obj2, obj))
    }

    list(F = Fmat, S = S, G = G, loss = loss)
}



pgd_lg = function(X, W, F0, S0, G0, lambda, lr = 1e-3, niter = 100, eps = 1e-6)
{
    n = nrow(X)
    p = ncol(X)
    nc = length(S0)

    Y = log(X + eps)
    Fmat = F0
    S = S0
    G = G0

    loss = c()

    for(k in 1:niter)
    {
        GS = G %*% Diagonal(x = S)
        R = Fmat %*% t(GS)
        logR = log(R + eps)

        obj1 = norm(Y - logR, type = "F")^2
        obj2 = sum(W * G^2)
        obj = obj1 + lambda * obj2
        loss = c(loss, obj)
        cat(sprintf("k = %d, obj1 = %f, obj2 = %f, obj = %f\n", k, obj1, obj2, obj))

        # dR
        dR = 2 * (logR - Y) / (R + eps)

        # dF
        dF = as.matrix(dR %*% GS)

        # dS
        A = matrix(0, n * p, nc)
        for(i in 1:nc)
        {
            A[, i] = as.numeric(outer(Fmat[, i], G[, i]))
        }
        dS = as.numeric(crossprod(A, as.numeric(dR)))

        # dG
        dG = as.matrix(crossprod(dR, Fmat %*% Diagonal(x = S))) + 2 * lambda * G * W

        # Gradient descent
        cat(sprintf("===> ||dF|| = %f, ||dS|| = %f, ||dG|| = %f\n\n",
                    norm(dF, type = "F"),
                    sqrt(sum(dS^2)),
                    norm(dG, type = "F")))
        Fmat = Fmat - lr * dF
        S = S - lr * dS
        G = G - lr * dG

        # Projection
        for(i in 1:n)
        {
            Fmat[i, ] = proj_pos_simplex(Fmat[i, ], 1)
        }
        S = pmax(S, 0)
        for(i in 1:nc)
        {
            G[, i] = proj_pos_simplex(G[, i], 1)
        }
    }

    list(loss = loss, F = Fmat, S = S, G = G)
}



Deconv = R6Class("Deconv",
    public = list(
        initialize = function(dat_unnorm, markers, normalize = TRUE)
        {
            mat_unnorm = t(as.matrix(dat_unnorm %>% select(-gene_name)))
            if(normalize)
            {
                # Normalize each observation
                count_subj = rowSums(mat_unnorm)
                private$dat = mat_unnorm / count_subj
            } else {
                private$dat = mat_unnorm
            }

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

        view_F = function(n = 20, ...)
        {
            nc = length(private$markers)
            view_matrix(t(private$mat_F[1:n, ]), legend_title = "Fraction", ...) +
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

        view_SG = function(...)
        {
            nc = length(private$markers)
            view_evec(private$mat_G %*% Diagonal(x = private$mat_S), ...) +
                scale_y_reverse("", breaks = 1:nc,
                                labels = substr(names(private$markers), 1, 3))
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
Deconv$set("public", "init_params", function(select = NULL) {
    coefs = lapply(private$markers, function(marker) {
        ind = match(marker, private$gene_name)
        dat_marker = private$dat[, ind]
        coef = eigs_sym(cov(dat_marker), 1)$vectors
        rownames(coef) = marker
        abs(coef)
    })
    mat = bdiag(coefs)
    colnames(mat) = names(private$markers)
    rownames(mat) = unlist(lapply(coefs, rownames), use.names = FALSE)

    # Normalize G to make each column sum to 1
    G0 = sweep(as.matrix(mat), 2, colSums(mat), "/")

    # X ~= F * S * G
    X = private$mat_exp
    scl = sqrt(colSums(X^2))
    F0 = t(fcnnls(G0 / scl, t(X) / scl, pseudo = FALSE)$x)
    F0 = F0 / rowSums(F0)
    # Xt = t(X)
    # rownames(Xt) = unlist(private$markers)
    # F0 = t(EstimateWeight(Xt, private$markers, method = "LM")$weight)
    # F0 = F0 / rowSums(F0)
    S0 = s_update(X, F0, G0)

    private$mat_W = as.matrix(mat == 0)
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
        self$init_params(select = NULL)
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
