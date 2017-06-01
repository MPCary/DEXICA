# The below code is from the fastICA package.  It has been modified slightly
# to return the tolerance (lim[it]) reached at the last iteration and the number of
# iterations run (it)

dex.ica.R.def = function(X, n.comp, tol, fun, alpha, maxit, verbose, w.init) {
  if (verbose && fun == "logcosh")
    message("Deflation FastICA using logcosh approx. to neg-entropy function")
  if (verbose && fun == "exp")
    message("Deflation FastICA using exponential approx. to neg-entropy function")
  min.tol = 1
  p <- ncol(X)
  W <- matrix(0, n.comp, n.comp)
  for (i in 1:n.comp) {
    if (verbose)
      message("Component ", i)
    w <- matrix(w.init[i, ], n.comp, 1)
    if (i > 1) {
      t <- w
      t[1:length(t)] <- 0
      for (u in 1:(i - 1)) {
        k <- sum(w * W[u, ])
        t <- t + k * W[u, ]
      }
      w <- w - t
    }
    w <- w/sqrt(sum(w^2))
    lim <- rep(1000, maxit)
    it <- 1
    if (fun == "logcosh") {
      while (lim[it] > tol && it < maxit) {
        wx <- t(w) %*% X
        gwx <- tanh(alpha * wx)
        gwx <- matrix(gwx, n.comp, p, byrow = TRUE)
        xgwx <- X * gwx
        v1 <- apply(xgwx, 1, FUN = mean)
        g.wx <- alpha * (1 - (tanh(alpha * wx))^2)
        v2 <- mean(g.wx) * w
        w1 <- v1 - v2
        w1 <- matrix(w1, n.comp, 1)
        it <- it + 1
        if (i > 1) {
          t <- w1
          t[1:length(t)] <- 0
          for (u in 1:(i - 1)) {
            k <- sum(w1 * W[u, ])
            t <- t + k * W[u, ]
          }
          w1 <- w1 - t
        }
        w1 <- w1/sqrt(sum(w1^2))
        lim[it] <- Mod(Mod(sum((w1 * w))) - 1)
        if(lim[it] < min.tol) min.tol = lim[it]
        if (verbose)
          message("Iteration ", it - 1, " tol = ", format(lim[it]))
        w <- matrix(w1, n.comp, 1)
      }
    }
    if (fun == "exp") {
      while (lim[it] > tol && it < maxit) {
        wx <- t(w) %*% X
        gwx <- wx * exp(-(wx^2)/2)
        gwx <- matrix(gwx, n.comp, p, byrow = TRUE)
        xgwx <- X * gwx
        v1 <- apply(xgwx, 1, FUN = mean)
        g.wx <- (1 - wx^2) * exp(-(wx^2)/2)
        v2 <- mean(g.wx) * w
        w1 <- v1 - v2
        w1 <- matrix(w1, n.comp, 1)
        it <- it + 1
        if (i > 1) {
          t <- w1
          t[1:length(t)] <- 0
          for (u in 1:(i - 1)) {
            k <- sum(w1 * W[u, ])
            t <- t + k * W[u, ]
          }
          w1 <- w1 - t
        }
        w1 <- w1/sqrt(sum(w1^2))
        lim[it] <- Mod(Mod(sum((w1 * w))) - 1)
        if(lim[it] < min.tol) min.tol = lim[it]
        if (verbose)
          message("Iteration ", it - 1, " tol = ", format(lim[it]))
        w <- matrix(w1, n.comp, 1)
      }
    }
    W[i, ] <- w
  }
  return(list("W" = W, "iterations" = it, "actual.tol" = lim[it], "min.obs.tol" = min.tol))
}
