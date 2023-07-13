
Smooth <- function(y, X., B., B.X, M.inv, W.inv, Q., bn, xn, X.B, X, B, Q, d, n, power){

  ynorm2 <- sum(y^2)
  X.y <- X. %*% y
  B.y <- B. %*% y
  vec1 <- B.y - B.X%*%M.inv%*%X.y
  s <- (ynorm2 - t(X.y)%*%M.inv%*%X.y - t(vec1)%*%W.inv%*%vec1) %>% drop
  z <- diag(1/sqrt(d)) %*% Q. %*% vec1

  MSC <- numeric(bn)
  ALPHA <- matrix(0, bn, bn)
  BETA <- matrix(0, xn, bn)

  for(gamma in 1:bn)
  {
    z1 <- z[1:gamma]
    u <- sort(z1^2)
    ua <- c(0, cumsum(u))
    sa <- (s + ua) / (n - xn - gamma + (0:gamma))
    Ru <- c(u, Inf)
    Rl <- c(0, u)
    a <- which((Rl < sa) & (sa <= Ru))
    v <- (1 - (sa[a]/z1^2)) * (z1^2 > sa[a])
    Q1 <- Q[,1:gamma]
    d1 <- d[1:gamma]

    ALPHA[,gamma] <- Alpha.h <- Q1%*%diag(v/d1, gamma, gamma)%*%t(Q1)%*%vec1

    BETA[,gamma] <- Beta.h <- M.inv %*% (X.y - X.B%*%Alpha.h) %>% drop

    y.hat <- X %*% Beta.h + B %*% Alpha.h
    rss <- sum((y - y.hat)^2)/n
    df <- xn + gamma
    MSC[gamma] <- rss / ((1 - df/n)^power)
  } #end for gamma

  opt <- which.min(MSC)
  Beta.h <- BETA[,opt]
  Alpha.h <- ALPHA[,opt]
  y.hat <- X%*%Beta.h + B%*%Alpha.h

  return(
    list(
      y.hat=y.hat,
      Beta.hat=Beta.h,
      Alpha.hat=Alpha.h,
      eignum=opt
    )
  )
}
