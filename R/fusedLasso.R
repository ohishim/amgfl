
fl.cd <- function(y, CD, alpha, ADJA, W){

  y0 <- sqrt(sum(y^2))
  y <- y / y0

  CD1 <- sort(unique(CD))

  ADJA0 <- ADJA[order(ADJA[,1]),]
  CD2 <- ADJA0[,1]
  CD3 <- ADJA0[,2]

  n <- length(y)
  m <- length(CD1)

  N <- as.vector(tapply(y, CD, length))
  Y <- split(y, CD)
  Dcd <- split(CD3, CD2)

  Dn <- lapply(Dcd, function(x){which(CD1 %in% x)})

  M_l <- sapply(Dcd, function(x){length(x)})

  y_t <- t(y)
  y_bar <- mean(y)
  RR <- diag(N)
  RR_inv <- diag(1/N)
  Q2 <- sapply(Y,sum)

  MuLSE <- RR_inv %*% Q2

  syusoku <- 1/10000

  W_l <- unlist(lapply(W,function(x){sum(x)}))
  mu_inf <- mean(y)

  lam_max1 <- max(2*(Q2 - mu_inf*N)/W_l)
  lam_max2 <- max(2*(mu_inf*N - Q2)/W_l)

  lam_max <- max(lam_max1,lam_max2)

  Lambda <- rev(lam_max * ((3/4)^((1:100)-1)))

  EGCV <- numeric(100)
  MU <- matrix(0,m,100)

  Mu0 <- MuLSE

  for(lam_i in 1:100)
  {
    lambda <- Lambda[lam_i]

    Mu_bef <- Mu0
    Mu_aft <- Mu0
    cond1 <- 0

    while(cond1 == 0)
    {
      for(l in 1:m)
      {
        D_cd0 <- Dcd[[l]]
        m_l <- length(D_cd0)
        D_l <- Dn[[l]]

        q2_d <- Q2[l]

        t <- sort(Mu_aft[D_l])

        w0 <- W[[l]]
        w <- w0[order(Mu_aft[D_l])]

        if(m_l == 1)
        {
          w1 <- c(-w, w)
        } else
        {
          w1 <- c(
            -sum(w),
            lapply(1:(m_l-1),function(x){sum(w[1:x])-sum(w[(x+1):(m_l)])}) %>%
              unlist,
            sum(w)
          )
        }

        v <- (q2_d - (lambda * w1)/2) / N[l]

        R_vl <- c(-Inf,t)
        R_vu <- c(t,Inf)

        R_bl <- v[-1]
        R_bu <- v[-(m_l+1)]

        a_ast1 <- which((R_vl < v) & (v <= R_vu))
        a_ast2 <- which((R_bl <= t) & (t < R_bu))

        Mu_aft[l]<- ifelse(length(a_ast1)!=0, v[a_ast1],t[a_ast2])

      }

      cond3 <- 0
      while(cond3 == 0)
      {
        cond3 <- 1

        Mu_uni <- unique(Mu_aft)
        b <- length(Mu_uni)

        if(b == 1)
        {
          Mu_aft <- rep(y_bar,m)

          cond3 <- 1
        } else if(b < m)
        {
          E0 <- cbind(1:m, Mu_aft)

          for(l in 1:b)
          {
            E_l0 <- subset(E0,E0[,2]==Mu_uni[l])
            E_l <- E_l0[,1]
            b_l <- length(E_l)

            if((1 < b_l) && (b_l < m))
            {
              F_l <- lapply(1:b_l,function(x){
                Dn[[E_l[x]]] %>% setdiff(., E_l)
              }) %>% unlist
              c_l <- length(F_l)

              q1_f <- sum(N[E_l])
              q2_f <- sum(Q2[E_l])

              if(c_l == 0)
              {
                Mu_aft[E_l] <- q2_f / q1_f

              } else
              {
                t0 <- Mu_aft[F_l]
                t <-  sort(t0)

                w0 <- lapply(1:b_l,function(x){
                  W[[E_l[x]]][which(Dn[[E_l[x]]]%in%setdiff(Dn[[E_l[x]]],E_l))]
                }) %>% unlist
                w <- w0[order(t0)]

                if(c_l == 1)
                {
                  w1 <- c(-w, w)
                } else
                {
                  w1 <- c(
                    -sum(w),
                    lapply(1:(c_l-1),function(x){
                      sum(w[1:x])-sum(w[(x+1):(c_l)])
                    }) %>% unlist,
                    sum(w)
                  )
                }

                v <- (q2_f - (lambda * w1)/2)/q1_f

                R_vl <- c(-Inf,t)
                R_vu <- c(t,Inf)

                R_bl <- v[-1]
                R_bu <- v[-(c_l+1)]

                a_ast1 <- which((R_vl < v) & (v <= R_vu))
                a_ast2 <- which((R_bl <= t) & (t < R_bu))

                Mu_aft[E_l]<- ifelse(length(a_ast1)!=0, v[a_ast1],t[a_ast2])

                if(length(a_ast2) != 0)
                {
                  cond3 <- 0
                }

              }
            }
          }
        }
      }

      SA <- max((Mu_bef-Mu_aft)^2) / max(Mu_bef^2)

      if(SA < syusoku)
      {
        cond1 <- 1
        MU[,lam_i] <- Mu0 <- Mu_aft
      } else
      {
        cond1 <- 0
        Mu_bef <- Mu_aft
      }
    }

    sigh <- 1 - 2*sum(Q2*Mu_aft) + sum(N*(Mu_aft^2))
    df <- length(unique(Mu_aft))

    EGCV[lam_i] <- sigh / ((1-df/n)^(alpha))

  }

  opt <- which.min(EGCV)
  Mu_hat <- MU[,opt] * y0

  return(
    list(
      Mu.hat = Mu_hat,
      lam.hat = Lambda[opt]
    )
  )
}
