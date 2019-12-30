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

    as.numeric(fcnnls(A, as.numeric(X), pseudo = FALSE)$x) + 1e-5
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
