pca_pen_prox = function(S, gr, lambda, gamma = 1.5, alpha = 0.01, maxit = 10, eps = 1e-4, verbose = 0)
{
    p = nrow(S)

    # Ssub = S[gr, gr]
    # e = RSpectra::eigs_sym(Ssub, 1)
    # x = matrix(0, p, p)
    # x[gr, gr] = tcrossprod(e$vectors)
    # z1 = z2 = x

    e = RSpectra::eigs_sym(S, 1)
    x = z1 = z2 = tcrossprod(e$vectors)

    # Time for Fantope projection
    time_f = c()
    # Time for penalty
    time_p = c()
    # Estimation error
    err = c()

    penmat = matrix(gamma * lambda, p, p)
    penmat[gr, gr] = lambda
    penmat[-gr, -gr] = gamma^2 * lambda

    for(i in 1:maxit)
    {
        newx = (z1 + z2) / 2
        resid = norm(newx - x, type = "F")
        x = newx

        t1 = Sys.time()
        newz1 = z1 - x + prox_fantope(S - penmat, 2 * x - z1, alpha, 1, inc = 30)
        t2 = Sys.time()
        newz2 = z2 - x + pmax(2 * x - z2, 0)
        t3 = Sys.time()

        resid1 = norm(newz1 - z1, type = "F")
        resid2 = norm(newz2 - z2, type = "F")
        z1 = newz1
        z2 = newz2
        if(verbose > 0)
            cat(sprintf("iter = %d, resid = %f, resid1 = %f, resid2 = %f\n", i, resid, resid1, resid2))

        time_f = c(time_f, t2 - t1)
        time_p = c(time_p, t3 - t2)
        err = c(err, resid)
        if(max(resid, resid1, resid2) < eps)
            break
    }
    list(proj = x, z1 = z1, z2 = z2, time_f = time_f, time_p = time_p, time_t = time_f + time_p, error = err)
}

pca_pen = function(S, gr, lambda, gamma = 1.5, alpha = 0.01, maxit = 10, eps = 1e-4, verbose = 0)
{
    p = nrow(S)

    # Ssub = S[gr, gr]
    # e = RSpectra::eigs_sym(Ssub, 1)
    # x = matrix(0, p, p)
    # x[gr, gr] = tcrossprod(e$vectors)
    # z1 = z2 = x

    e = RSpectra::eigs_sym(S, 1)
    x0 = tcrossprod(e$vectors)
    pca_pen_(S, gr, x0, lambda, gamma, alpha, maxit, 100, 10, eps, verbose)
}
