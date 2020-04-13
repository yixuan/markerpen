# Project x = (x1, ..., xn) to {x: a <= sum(x) <= b, x >= 0}
proj_pos_simplex_ineq = function(x, a, b)
{
    n = length(x)
    constr = rbind(rep(1, n), rep(-1, n), diag(n))
    rhs = c(a, -b, rep(0, n))
    sol = quadprog::solve.QP(Dmat = diag(n), dvec = x, Amat = t(constr), bvec = rhs)
    sol$solution
}

proj_constr = function(x, gr, gr_weight)
{
    diagx = diag(x)
    diagx[gr] = proj_pos_simplex_ineq(diagx[gr], gr_weight, 1)
    w = sum(diagx[gr])
    diagx[-gr] = proj_pos_simplex(diagx[-gr], 1 - w)
    x = pmax(x, 0)
    diag(x) = diagx
    x
}

pca_pen_prox = function(S, gr, lambda, gr_weight = 0.8, maxit = 10, alpha = 0.01, eps = 1e-4, Pi = NULL)
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

    for(i in 1:maxit)
    {
        newx = (z1 + z2) / 2
        resid = norm(newx - x, type = "F")
        x = newx

        t1 = Sys.time()
        newz1 = z1 - x + prox_fantope(S - lambda, z2, alpha, 1, inc = 30)
        t2 = Sys.time()
        newz2 = z2 - x + proj_constr(z1, gr, gr_weight)
        t3 = Sys.time()

        resid1 = norm(newz1 - z1, type = "F")
        resid2 = norm(newz2 - z2, type = "F")
        z1 = newz1
        z2 = newz2
        cat(sprintf("iter = %d, resid = %f, resid1 = %f, resid2 = %f\n", i, resid, resid1, resid2))

        time_f = c(time_f, t2 - t1)
        time_p = c(time_p, t3 - t2)
        if(!is.null(Pi))
            err = c(err, norm(x - Pi, type = "F"))
        if(max(resid, resid1, resid2) < eps)
            break
    }
    list(proj = x, z1 = z1, z2 = z2, time_f = time_f, time_p = time_p, time_t = time_f + time_p, error = err)
}
