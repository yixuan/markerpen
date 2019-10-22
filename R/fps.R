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
    e = RSpectra::eigs_sym(S, d)
    x = z1 = z2 = tcrossprod(e$vectors)

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
