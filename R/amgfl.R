#' @title Trend estimations by semiparametric additive model
#'     with graph trend filtering via generalized fused Lasso (v0.7.0)
#' @description \code{amgfl} This function provides trend estimation results
#'     by semiparametric additive model with graph trend filtering via generalized fused Lasso
#'
#' @importFrom dplyr group_by mutate summarize
#' @importFrom ggplot2 aes element_blank facet_wrap geom_line geom_segment geom_sf
#'     ggplot scale_fill_gradientn theme xlab
#' @importFrom grDevices colorRampPalette
#' @importFrom magrittr %>% inset2 set_colnames set_names set_rownames subtract
#' @importFrom purrr exec map
#' @importFrom Rfast colmeans
#' @importFrom sf st_point st_sf st_sfc
#' @importFrom stats quantile sd
#'
#' @param dataset a list `list(y=y, X1=X1, X2=X2)` of which `y` is a vector of
#'     a response variable, `X1` is a matrix of explanatory variables, without
#'     an intercept, and with colnames, where the variables are fitted with a
#'     linear model, and `X2` is a matrix of explanatory variables with colnames,
#'     where the variables are fitted with an additive model;
#'     `X1` can be omitted and the intercept is automatically added into `X1`
#' @param space a list expressing spatial information;
#'     if a spatial division is given, `space` is given by the form of
#'     `list(area=area, adjacent=adjacent, Sp=Sp, div.pol=div.pol)`,
#'     where `area` is a vector of integers expressing sub-areas for each individual,
#'     `adjacent` is a matrix with two columns expressing adjacent relationships among sub-areas,
#'     `Sp` is an optional parameter which is a matrix with two columns expressing
#'     coordinates of observed points, and `div.pol` is an optional parameter which
#'     is an sfc_POLYGON object expressing a spatial division;
#'     if not, i.e., a spatial division is given by `amgfl`, `space` is given by the form of
#'     `list(Sp=Sp, m=m, spB=spB, subregion=subregion, seed=seed)`;
#'     where `Sp` is the same of above but must be required,
#'     `m` is an optional parameter expressing the number of sub-areas,
#'     `spB` is an optional parameter which is an sfc_POLYGON object expressing
#'     a boundary of the space subjected to analysis, `subregion` is an optional
#'     parameter which is a vector expressing subregions for each individual,
#'     and `seed` is an optional parameter which is the first parameter of `set.seed`;
#'     if `subregion` is given, `m` is a vector of the numbers of sub-areas for
#'     each subregion and `spB` is an sfc_POLYGON object expressing subregion-wise
#'     boundary and must be given
#' @param nknots the number of knots
#' @param scale1 if `TRUE`, `scale` is applied to `X1`
#' @param scale2 if `TRUE`, `scale` is applied to `X2`
#' @param power1 a positive value adjusting the strength of the penalty in EGCV for
#'     selecting tuning parameter in generalized fused Lasso
#' @param power2 a positive value adjusting the strength of the penalty in EGCV for
#'     selecting the number of eigenvalues in nonparametric smoothing; default is `log(length(y))`
#' @param maxit an upper bound for iteration
#' @param tol a tolerance for convergence
#' @param progress if `TRUE`, estimation progress is displayed
#' @param complete if `TRUE`, the message is displayed when the process terminates
#'
#' @return a list having the following elements:
#' \item{mu}{a vector of the estimates of spatial effects}
#'
#' \item{beta}{a vector of the estimates of coefficients for `X`}
#'
#' \item{alpha}{a vector of the estimates of coefficients for `B`}
#'
#' \item{fitted.values}{a vector of the fitted values}
#'
#' \item{OLS}{results of the ordinary least square (OLS) estimation given by the form of
#'     a list `list(mu, beta, alpha, fitted.values)`,
#'     where `mu`, `beta`, and `alpha` are the OLS estimates
#'     and `fitted.values` is a vector of fitted values based on the OLS estimates}
#'
#' \item{summary}{a summary of estimation results given by the form of
#'     a data.frame `data.frame(lambda, eigenvalues, division, group, runtime, iteration, convergence)`,
#'     where `lambda` is the optimal tuning parameter,
#'     `eigenvalues` is the optimal number of eigenvalues,
#'     `division` is the number of sub-areas,
#'     `group` is the number of estimated clusters,
#'     `runtime` is computational time (seconds),
#'     `iteration` is the number of iterations until convergence,
#'     and `convergence`: whether the algorithm converges}
#'
#' \item{dataset}{a list of dataset `list(y, X1, X2)`}
#'
#' \item{space}{a list expressing a spatial information}
#'
#' \item{par}{a list of the entered parameters}
#'
#' \item{scale}{a list of attributes of `X1` and `X2` given by applying `scale`}
#'
#' \item{X}{a centralized matrix for which coefficients are estimated without the penalty}
#'
#' \item{B}{a centralized matrix for which coefficients are estimated with the penalty}
#'
#' \item{center}{a list of vectors using for centering of `X` and `B`}
#'
#' \item{knots}{a matrix of knots}
#'
#' \item{trends}{a list having figures of trends for each variable of `X2`}
#'
#' \item{trends.all}{a figure consists of all figures of `trends`}
#'
#' \item{division}{a figure expressing the spatial division}
#'
#' \item{cluster}{a figure expressing `mu` and the estimated cluster}
#'
#' @export
#' @examples
#' #amgfl(dataset, space, nknots)

amgfl <- function(
    dataset, space, nknots=NULL,
    scale1=FALSE, scale2=FALSE, power1=2, power2=NULL,
    maxit=100, tol=1e-8, progress=FALSE, complete=TRUE
){

  T1 <- proc.time()[3]

  y <- dataset$y
  X1 <- dataset$X1
  X2 <- dataset$X2

  if(is.null(y)){stop("`dataset$y` must be given")}
  if(is.null(X2)){stop("`dataset$X2` must be given")}

  n <- length(y)
  is.X1 <- !is.null(X1)

  if(is.null(power2)){power2 <- log(n)}

  if(is.null(space$adjacent))
  {
    if(is.null(space$Sp)){stop("`space$Sp` must be given")}
    if(is.null(space$seed)){space$seed <- 123}

    div <- divsp(space$Sp, space$m, space$spB, space$subregion, space$seed)

    space$area <- div$area
    space$adjacent <- div$adjacent
    space$`div.pol` <- div$div.pol$geometry
  } else
  {
    if(is.null(space$area)){stop("`space$area` must be given")}
  }

  indiv.order <- order(space$area)
  input.order <- (1:n)[indiv.order] %>% order
  area <- space$area[indiv.order]
  adjacent <- space$adjacent
  div.pol <- space$`div.pol`

  space$area <- area

  if(!is.null(space$Sp)){space$Sp <- space$Sp[indiv.order,]}

  y <- y[indiv.order]

  if(is.X1){X1 <- X1[indiv.order,,drop=F]}
  X2 <- X2[indiv.order,,drop=F]

  ##############################################################################
  ###　  preparation 1
  ##############################################################################

  V1 <- colnames(X1)
  V2 <- colnames(X2)

  if(is.X1 & is.null(V1)){stop("column names of `dataset$X1` must be given")}
  if(is.null(V2)){stop("column names of `dataset$X2` must be given")}

  .X1 <- X1
  .X2 <- X2
  scaleX1 <- scaleX2 <- NULL

  if(scale1 & is.X1)
  {
    X1 <- scale(.X1)
    scaleX1 <- list(
      center = attr(X1, "scaled:center"),
      scale = attr(X1, "scaled:scale")
    )
  }

  if(scale2)
  {
    X2 <- scale(.X2)
    scaleX2 <- list(
      center = attr(X2, "scaled:center"),
      scale = attr(X2, "scaled:scale")
    )
  }

  CDlist <- table(area) %>% data.frame

  R <- do.call(
    rbind,
    lapply(1:NROW(CDlist), function(x){
      M <- matrix(0, CDlist$Freq[x], NROW(CDlist))
      M[,x] <-1
      return(M)
    })
  )

  CD1 <- sort(unique(area))

  ADJA0 <- adjacent[order(adjacent[,1]),]
  CD2 <- ADJA0[,1]
  CD3 <- ADJA0[,2]

  Dcd <- split(CD3, CD2)
  Dn <- lapply(Dcd,function(x){which(CD1 %in% x)})

  fx <- function(x){cbind(x, x^2, x^3)}

  fb <- function(x, knots, nknots){
    out <- lapply(1:nknots, function(j){
      ((x - knots[j])^3) * (x > knots[j])
    }) %>% do.call(cbind, .)
    return(out)
  }

  k1 <- ifelse(is.X1, ncol(X1), 0)
  k2 <- ncol(X2)
  m <- ncol(R)

  nknots.max <- apply(X2, 2, function(x){unique(x) %>% length}) %>% min %>%
    subtract(4) %>% c(floor((n - (k1+3*k2) - m)/k2)) %>% min

  if(is.null(nknots))
  {
    nknots <- c(10, nknots.max-2) %>% min
  } else
  {
    if(nknots > nknots.max)
    {
      stop(paste0("`nknots` must be smaller than ", nknots.max, " at least"))
    }
  }

  KNOTS <- map(1:k2, ~{
    quantile(X2[,.x], probs = seq(0, 1, (1/(nknots+1))))[-c(1, nknots+2)]
  }) %>% exec(rbind, !!!.) %>% set_rownames(V2)

  ##############################################################################
  ###   calculate initial values
  ##############################################################################

  .X <- map(1:k2, ~fx(X2[,.x])) %>% exec(cbind, !!!.) %>%
    set_colnames(map(V2, ~paste0(.x, 1:3)) %>% unlist)

  if(is.X1)
  {
    .X <- cbind(
      X1 %>% set_colnames(V1),
      .X
    )
  }

  X <- cbind(
    "icpt" = 1,
    scale(.X, scale=F)
  )

  .B <- map(1:k2, ~fb(X2[,.x], KNOTS[.x,], nknots)) %>%
    exec(cbind, !!!.) %>%
    set_colnames(
      map(V2, ~paste0(.x, 1:nknots)) %>% unlist
    )

  B <- scale(.B, scale=F)

  xn <- ncol(X)
  bn <- ncol(B)

  B. <- t(B)
  B.B <- B. %*% B
  X. <- t(X)
  M <- X. %*% X

  M.inv <- try(solve(M), silent=TRUE)
  if("try-error" %in% class(M.inv)){stop("cubic spline probably generated too large values or caused a rank deficient")}

  X.B <- X. %*% B
  B.X <- t(X.B)
  W <- B.B - B.X%*%M.inv%*%X.B
  W.inv <- solve(W)
  eigW <- eigen(W)
  d <- eigW$values

  if(any(d < 0)){stop("`nknots` is probably too large")}

  Q <- eigW$vectors
  Q. <- t(Q)

  res0 <- Smooth(y, X., B., B.X, M.inv, W.inv, Q., bn, xn, X.B, X, B, Q, d, n, power2)

  Mu0 <- res0$Beta.hat[1] %>% rep(ncol(R))
  Beta0 <- res0$Beta.hat[-1]
  Alpha0 <- res0$Alpha.hat

  ##############################################################################
  ###   preparation 2
  ##############################################################################

  X <- X[,-1]
  xn <- ncol(X)

  X. <- t(X)
  M <- X. %*% X
  M.inv <- solve(M)
  X.B <- X. %*% B
  B.X <- t(X.B)
  W <- B.B - B.X%*%M.inv%*%X.B
  W.inv <- solve(W)
  eigW <- eigen(W)
  d <- eigW$values
  Q <- eigW$vectors
  Q. <- t(Q)

  Z <- cbind(X, B, R); Z. <- t(Z)

  LSE <- try(solve(Z.%*%Z) %*% Z. %*% y, silent=TRUE)
  if("try-error" %in% class(LSE)){stop("`nknots` is probably too large")}

  Beta.LSE <- LSE[1:xn]
  Alpha.LSE <- LSE[(xn+1):(xn+bn)]
  Mu.LSE <- LSE[(xn+bn+1):ncol(Z)]

  Weight  <- lapply(1:NROW(CDlist), function(x){2/abs(Mu.LSE[x]- Mu.LSE[Dn[[x]]])})

  ##############################################################################
  ###   estimation
  ##############################################################################

  IND <- 0
  it <- 0

  Mu.bef <- Mu.aft <- Mu0
  Beta.bef <- Beta.aft <- Beta0
  Alpha.bef <- Alpha.aft <- Alpha0

  while(IND == 0 & it < maxit)
  {
    it <- it + 1

    res1 <- fl.cd(y - X%*%Beta.aft - B%*%Alpha.aft, area, power1, ADJA0, Weight)
    Mu.aft <- res1$Mu.hat
    lam.h <- res1$lam.hat

    Rmu <- lapply(1:nrow(CDlist), function(j){rep(Mu.aft[j], CDlist$Freq[j])}) %>% unlist
    res2 <- Smooth(y - Rmu, X., B., B.X, M.inv, W.inv, Q., bn, xn, X.B, X, B, Q, d, n, power2)
    Beta.aft <- res2$Beta.hat
    Alpha.aft <- res2$Alpha.hat
    eignum <- res2$eignum

    sa1 <- max((Mu.aft - Mu.bef)^2) / max(Mu.bef^2)
    sa2 <- max((Beta.aft - Beta.bef)^2) / max(Beta.bef^2)

    if(any(Alpha.bef != 0))
    {
      sa3 <- max((Alpha.aft - Alpha.bef)^2) / max(Alpha.bef^2)
    } else
    {
      sa3 <- 0
    }

    if(all(c(sa1, sa2, sa3) < tol))
    {
      IND <- 1
    } else
    {
      Mu.bef <- Mu.aft
      Beta.bef <- Beta.aft
      Alpha.bef <- Alpha.aft
    }

    if(progress){print(paste(it, sa1, sa2, sa3, sep="_"))}
  }

  T2 <- proc.time()[3]

  runtime <- (T2 - T1) %>% set_names(NULL)
  convergence <- IND == 1

  if(complete)
  {
    message(paste0("complete: ", it, " iterations; ", round(runtime, 2), " sec."))

    if(!convergence){
      message("The algorithm did not converge")
    }
  }

  #############################################################################
  ###   summarize
  #############################################################################

  Mu.hat <- Mu.aft
  Beta.hat <- Beta.aft %>% set_names(colnames(X))
  Alpha.hat <- Alpha.aft %>% set_names(colnames(B))

  out <- list(
    mu = Mu.hat,
    beta = Beta.hat,
    alpha = Alpha.hat,
    fitted.values = drop(
      X%*%Beta.aft + B%*%Alpha.aft + Rmu
    )[input.order],
    OLS = list(
      mu = Mu.LSE,
      beta = Beta.LSE,
      alpha = Alpha.LSE,
      fitted.values = drop(
        X%*%Beta.LSE + B%*%Alpha.LSE + unlist(
          lapply(1:nrow(CDlist), function(j){rep(Mu.LSE[j], CDlist$Freq[j])})
        )
      )[input.order]
    ),
    summary = data.frame(
      lambda = lam.h,
      eigenvalues = eignum,
      division = CD1 %>% length,
      group = Mu.aft %>% unique %>% length,
      runtime = runtime,
      iteration = it,
      convergence = convergence
    ),
    dataset = dataset,
    space = space %>% inset2("Sp", space$Sp[input.order,]) %>%
      inset2("area", space$area[input.order]),
    par = list(
      nknots=nknots, scale1=scale1, scale2=scale2, power1=power1, power2=power2,
      maxit=maxit, tol=tol, progress=progress, complete=complete
    ),
    scale = list(
      X1 = scaleX1,
      X2 = scaleX2
    ),
    X = X[input.order,],
    B = B[input.order,],
    center = list(
      X = colMeans(.X),
      B = attributes(B)$`scaled:center`
    ),
    knots = KNOTS
  )

  ##############################################################################
  ###   trends
  ##############################################################################

  C <- apply(.X2, 2, range)
  C.div <- map(1:k2, ~seq(
    C[1, .x], C[2, .x], length=100
  )) %>% exec(cbind, !!!.) %>% set_colnames(V2)

  if(scale2)
  {
    C.div <- scale(
      C.div,
      center = colmeans(.X2),
      scale = scale(.X2, scale=F) %>% apply(2, sd)
    )
  }

  H <- map(1:k2, ~cbind(
    fx(C.div[,.x]) %>%
      scale(center=.X[,paste0(V2[.x], 1:3)] %>% colmeans(), scale=F),
    fb(C.div[,.x], KNOTS[.x,], nknots) %>%
      scale(center=.B[,paste0(V2[.x], 1:nknots)] %>% colmeans(), scale=F)
  ))

  ggD <- map(1:k2, ~{
    data.frame(
      x = seq(C[1, .x], C[2, .x], length=100),
      y = H[[.x]] %*% c(
        Beta.hat[paste0(V2[.x], 1:3)],
        Alpha.hat[paste0(V2[.x], 1:nknots)]
      ) %>% drop
    ) %>% unique %>% mutate(v = V2[.x])
  }) %>% exec(rbind, !!!.) %>% mutate(v = factor(v, levels=V2))

  ggDseg0 <- split(ggD$y, ggD$v) %>% map(~c(min(.x), max(.x) - min(.x))) %>%
    exec(rbind, !!!.) %>% as.data.frame %>% set_colnames(c("m", "r"))

  ggDseg <- V2 %>% map(~{
    data.frame(
      x = .X2[, .x],
      y = ggDseg0[.x, 1],
      v = .x
    ) %>% mutate(yend = y - ggDseg0[.x, 2]/50)
  }) %>% exec(rbind, !!!.) %>% mutate(v = factor(v, levels=V2))

  out$trends <- map(1:k2, ~{
    xlab <- V2[.x]
    ggD %>% dplyr::filter(v == xlab) %>%
      ggplot() +
      geom_segment(
        data = dplyr::filter(ggDseg, v == xlab),
        aes(x=x, y=y, xend=x, yend=yend), linewidth=0.1
      ) +
      geom_line(aes(x=x, y=y)) +
      facet_wrap(~v) +
      theme(axis.title = element_blank())
  }) %>% set_names(V2)

  out$trends.all <- ggplot(ggD) +
    geom_segment(
      data = ggDseg,
      aes(x=x, y=y, xend=x, yend=yend), linewidth=0.1
    ) +
    geom_line(aes(x=x, y=y)) +
    theme(axis.title = element_blank()) +
    facet_wrap(.~v, scales="free")

  ##############################################################################
  ###   spatial effects
  ##############################################################################

  if(!is.null(div.pol))
  {
    divD <- st_sf(
      area = CD1,
      geometry = div.pol
    )

    jet.colors <- colorRampPalette(
      c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000")
    )

    if(!is.null(space$Sp))
    {
      sp <- st_sfc(
        lapply(1:n, function(i){st_point(space$Sp[i,])})
      )
    }

    if(exists("div"))
    {
      out$division <- div$fig.div
    } else
    {
      out$division <- ggplot(div.pol) + geom_sf()

      if(exists("sp"))
      {
        out$division <- out$division +
          geom_sf(data=sp, size=0.8, alpha=0.1)
      }
    }

    out$cluster <- mutate(divD, mu=Mu.hat - mean(Rmu)) %>%
      group_by(mu) %>% summarize() %>%
      ggplot() +
      geom_sf(aes(fill=mu)) +
      scale_fill_gradientn(colours = jet.colors(10000)) +
      theme(legend.title = element_blank())

    if(exists("sp"))
    {
      out$cluster <- out$cluster +
        geom_sf(data=sp, size=0.8, alpha=0.1)
    }
  }

  ##############################################################################
  ###   output
  ##############################################################################

  return(out)
}
