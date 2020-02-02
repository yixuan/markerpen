# Projection to Fantope
proj_fantope = function(x, d)
{
    p = nrow(x)
    e = eigen(x, symmetric = TRUE)
    gamma = e$values
    # Find a shift theta such that
    # sum(pmin(1, pmax(0, gamma - theta))) = d
    lb = gamma[p]
    ub = gamma[1]
    for(i in 1:100)
    {
        mid = (lb + ub) / 2
        gamma_new = pmin(1, pmax(0, gamma - mid))
        r = sum(gamma_new)
        if(r > d) { lb = mid } else { ub = mid }
        if(abs(r - d) < 1e-6)
            break
    }
    e$vectors %*% diag(gamma_new) %*% t(e$vectors)
}

# set.seed(123)
# x = matrix(rnorm(100^2, sd = 0.1), 100)
# S = crossprod(x)
# x = matrix(rnorm(100^2, sd = 0.1), 100)
# X = crossprod(x)
# alpha = 0.001
#
# r1 = proj_fantope(X + alpha * S, d = 5)
# r2 = prox_fantope(S, X, alpha, d = 5)
# max(abs(r1 - r2))
#
# head(eigen(r1, symmetric = TRUE)$values)
# head(eigen(r2, symmetric = TRUE)$values)



# Proximal-proximal gradient method
fps_prox_gr = function(S, d, gr, mu, lambda, gamma = 1.0, maxit = 10, alpha = 0.01, Pi = NULL)
{
    if(d != 1)
        stop("d must be 1")
    p = nrow(S)

    Ssub = S[gr, gr]
    e = RSpectra::eigs_sym(Ssub, d)
    x = matrix(0, p, p)
    x[gr, gr] = tcrossprod(e$vectors)
    z1 = z2 = x

    # e = RSpectra::eigs_sym(S, d)
    # x = z1 = z2 = tcrossprod(e$vectors)

    # Time for Fantope projection
    time_f = c()
    # Time for penalty
    time_p = c()
    # Estimation error
    err = c()

    penmat = matrix(0.5 * gamma * lambda, p, p)
    penmat[gr, gr] = -mu
    penmat[-gr, -gr] = -lambda

    for(i in 1:maxit)
    {
        newx = (z1 + z2) / 2
        resid = norm(newx - x, type = "F")
        x = newx

        t1 = Sys.time()
        # newz1 = z1 - x + proj_fantope(z2 + alpha * (S + penmat), d)
        newz1 = z1 - x + prox_fantope(S + penmat, z2, alpha, d, inc = 30)
        t2 = Sys.time()
        newz2 = z2 - x + pmax(z1, 0)
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
    }
    list(proj = x, z1 = z1, z2 = z2, time_f = time_f, time_p = time_p, time_t = time_f + time_p, error = err)
}



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

pca_pen_prox = function(S, gr, lambda, gr_weight = 0.8, maxit = 10, alpha = 0.01, Pi = NULL)
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
    }
    list(proj = x, z1 = z1, z2 = z2, time_f = time_f, time_p = time_p, time_t = time_f + time_p, error = err)
}

pca_pen = function(S, gr, lambda, gr_weight = 0.8, maxit = 10, rho = 1, dyn_rho = TRUE, Pi = NULL)
{
    p = nrow(S)

    e = RSpectra::eigs_sym(S, 1)
    x = tcrossprod(e$vectors)
    y = u = matrix(0, p, p)

    # Time for Fantope projection
    time_f = c()
    # Time for penalty
    time_p = c()
    # Estimation error
    err = c()

    for(i in 1:maxit)
    {
        t1 = Sys.time()
        # x = proj_fantope(y - u + (S - lambda) / rho, 1)
        x = prox_fantope(S - lambda, y - u, alpha = 1 / rho, d = 1, inc = 30)
        t2 = Sys.time()

        newy = proj_constr(x + u, gr, gr_weight)
        t3 = Sys.time()

        resid1 = norm(x - newy, type = "F")
        resid2 = rho * norm(newy - y, type = "F")

        cat(sprintf("iter = %d, resid1 = %f, resid2 = %f, rho = %f\n", i, resid1, resid2, rho))

        y = newy
        u = u + x - y
        if(isTRUE(dyn_rho))
        {
            if(resid1 > 10 * resid2)
                rho = rho * 2
            if(resid2 > 10 * resid1)
                rho = rho / 2
        }

        time_f = c(time_f, t2 - t1)
        time_p = c(time_p, t3 - t2)
        if(!is.null(Pi))
            err = c(err, norm(y - Pi, type = "F"))
    }
    list(proj = y, x = x, time_f = time_f, time_p = time_p, time_t = time_f + time_p, error = err)
}

