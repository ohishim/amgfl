#' @title Simultaneous confidence intervals via the bootstrap
#' @description \code{amgfl.CI} This function provides simultaneous confidence intervals
#'     via the bootstrap
#'
#' @importFrom dplyr group_by mutate summarize
#' @importFrom ggplot2 aes element_blank facet_wrap geom_line geom_ribbon
#'     geom_segment geom_sf ggplot scale_fill_gradientn theme
#' @importFrom magrittr %>% divide_by multiply_by raise_to_power set_colnames set_names
#' @importFrom purrr invoke map map_dbl
#' @importFrom sf st_sf
#'
#' @param res the output of `amgfl`
#' @param N the number of bootstrap sampling
#' @param Ndiv the number of divisions for intervals
#' @param .alpha calculate (1-`.alpha`)-simultaneous confidence intervals
#' @param seed the first parameter of `set.seed`
#' @param progress if TRUE, the progress of bootstrap is displayed
#'
#' @return a list having the following elements:
#' \item{trends}{a list having figures of trends with simultaneous confidence
#'     for each variable of `X2`}
#'
#' \item{trends.all}{a figure consists of all figures of `trends`}
#'
#' \item{cluster}{a figure consists of three panels;
#'     the left expresses the lower of simultaneous confidence intervals;
#'     the middle expresses `mu` and the estimated cluster;
#'     the right expresses the upper of simultaneous confidence intervals}
#'
#' @export
#' @examples
#' #amgfl.CI(res)

amgfl.CI <- function(
    res, N=100, Ndiv=100, .alpha=0.05, seed=NULL, progress=FALSE
){

  ##############################################################################
  ###   preparation
  ##############################################################################

  for(x in c("y", "X1", "X2"))
  {
    assign(x, res$dataset[[x]])
  }

  space <- res$space
  div.pol <- res$space$div.pol

  for(x in c("nknots", "scale1", "scale2", "power1", "power2", "maxit", "tol"))
  {
    assign(x, res$par[[x]])
  }

  if(!is.null(seed)){set.seed(seed)}

  Fit <- res$fitted.values
  ybar <- mean(y)
  mubar <- res$mu[space$area] %>% mean
  e <- y - Fit - (ybar - mubar)

  n <- length(y)

  V2 <- colnames(X2)
  k2 <- length(V2)
  m <- space$area %>% unique %>% length

  fx <- function(x){cbind(x, x^2, x^3)}

  fb <- function(x, knots, nknots){
    out <- lapply(1:nknots, function(j){
      ((x - knots[j])^3) * (x > knots[j])
    }) %>% do.call(cbind, .)
    return(out)
  }

  ##############################################################################
  ###   bootstrap
  ##############################################################################

  XI <- map(1:k2, ~matrix(0, N, 3+nknots)) %>% set_names(V2)
  MU <- matrix(0, N, m)

  for(bs.i in 1:N)
  {
    De <- sample(e, n, replace=T)
    .y <- Fit + (ybar - mubar) + De

    .res <- amgfl(
      list(y=.y, X1=X1, X2=X2), space, nknots=nknots, scale1=scale1, scale2=scale2,
      power1=power1, power2=power2, maxit=maxit, tol=tol, progress=F, complete=F
    )

    for(j in 1:k2)
    {
      XI[[j]][bs.i,] <- c(
        .res$beta[paste0(V2[j], 1:3)],
        .res$alpha[paste0(V2[j], 1:nknots)]
      ) %>% set_names(NULL)
    }

    MU[bs.i, ] <- .res$mu

    if(progress)
    {
      print(paste0("bootstrap: ", bs.i))
    }
  } #end for bs.i

  ##############################################################################
  ###   calculate confidence intervals
  ##############################################################################

  XIbar <- map(XI, ~Rfast::colmeans(.x)) %>% invoke(rbind, .)
  XI.C <- map(XI, ~scale(.x, scale=F))
  S <- map(XI.C, ~{t(.x)%*%.x / N})

  Mubar <- Rfast::colmeans(MU)
  Mu.C <- scale(MU, scale=F)
  s <- Mu.C %>% raise_to_power(2) %>% Rfast::colmeans()

  C <- apply(X2, 2, range)
  C.div <- map(1:k2, ~seq(
    C[1, .x], C[2, .x], length=Ndiv
  )) %>% invoke(cbind, .) %>% set_colnames(V2)

  if(scale2)
  {
    C.div <- scale(
      C.div,
      center = Rfast::colmeans(X2),
      scale = scale(X2, scale=F) %>% apply(2, sd)
    )
  }

  H <- map(1:k2, ~cbind(
    fx(C.div[,.x]) %>%
      scale(center=res$center$X[paste0(V2[.x], 1:3)], scale=F),
    fb(C.div[,.x], res$knots[.x,], nknots) %>%
      scale(center=res$center$B[paste0(V2[.x], 1:nknots)], scale=F)
  ))

  hSh <- map(1:k2, ~{
    j <- .x
    return(
      map_dbl(1:Ndiv, ~{
        h <- H[[j]][.x,]

        return(
          drop(S[[j]]%*%h) %>% multiply_by(h) %>% sum %>% sqrt
        )
      })
    )
  }) %>% invoke(cbind, .) %>% set_colnames(V2)

  .T <- map(1:k2, ~{
    j <- .x

    return(
      map(1:Ndiv, ~{
        h <- H[[j]][.x,]

        return(
          XI.C[[j]] %>% scale(center=F, scale=1/h) %>% Rfast::rowsums() %>% abs %>%
            divide_by(hSh[.x, j])
        )
      }) %>% invoke(rbind, .) %>% Rfast::colMaxs(value=T)
    )
  }) %>% invoke(cbind, .) %>% set_colnames(V2)

  .t <- abs(Mu.C) %>% divide_by(sqrt(s)) %>% Rfast::rowMaxs(value=T)

  z <- map_dbl(1:k2, ~quantile(.T[,.x], 1 - (.alpha/2)))
  z0 <- quantile(.t, 1 - (.alpha/2))

  ##############################################################################
  ###   trends
  ##############################################################################

  ggD <- map(1:k2, ~{
    data.frame(
      x = seq(C[1, .x], C[2, .x], length=Ndiv),
      y = H[[.x]] %*% c(
        res$beta[paste0(V2[.x], 1:3)],
        res$alpha[paste0(V2[.x], 1:nknots)]
      ) %>% drop,
      z = z[.x]*hSh[,.x]
    ) %>% unique %>% mutate(v = V2[.x])
  }) %>% invoke(rbind, .) %>% mutate(v = factor(v, levels=V2))

  ggDseg0 <- split(ggD, ggD$v) %>% map(~{
    y <- .x$y; z <- .x$z
    y1 <- y - z; y2 <- y + z
    return(c(min(y1), max(y2) - min(y1)))
  }) %>%
    invoke(rbind, .) %>% as.data.frame %>% set_colnames(c("m", "r"))

  ggDseg <- V2 %>% map(~{
    data.frame(
      x = res$dataset$X2[, .x],
      y = ggDseg0[.x, 1],
      v = .x
    ) %>% mutate(yend = y - ggDseg0[.x, 2]/50)
  }) %>% invoke(rbind, .) %>% mutate(v = factor(v, levels=V2))

  trends <- map(1:k2, ~{
    .xlab <- V2[.x]
    ggD %>% dplyr::filter(v == .xlab) %>%
      ggplot() +
      geom_segment(
        data = dplyr::filter(ggDseg, v == .xlab),
        aes(x=x, y=y, xend=x, yend=yend), linewidth=0.1
      ) +
      geom_ribbon(aes(x=x, ymin=y-z, ymax=y+z), color="grey80", fill="grey80") +
      geom_line(aes(x=x, y=y)) +
      facet_wrap(~v) +
      theme(axis.title = element_blank())
  }) %>% set_names(V2)

  trends.all <- ggplot(ggD) +
    geom_segment(
      data = ggDseg,
      aes(x=x, y=y, xend=x, yend=yend), linewidth=0.1
    ) +
    geom_ribbon(aes(x=x, ymin=y-z, ymax=y+z), color="grey80", fill="grey80") +
    geom_line(aes(x=x, y=y)) +
    theme(axis.title = element_blank()) +
    facet_wrap(.~v, scales="free")

  ##############################################################################
  ###   spatial effects
  ##############################################################################

  divD0 <- st_sf(
    area = space$area %>% unique,
    geometry = div.pol,
    mu = res$mu
  )

  divD <- rbind(
    divD0 %>% mutate(mu = mu - z0*sqrt(s)) %>%
      group_by(mu) %>% summarize() %>% mutate(key = "lower"),
    divD0 %>%
      group_by(mu) %>% summarize() %>% mutate(key = "estimates"),
    divD0 %>% mutate(mu = mu + z0*sqrt(s)) %>%
      group_by(mu) %>% summarize() %>% mutate(key = "upper")
  ) %>% mutate(
    key = factor(key, levels=c("lower", "estimates", "upper")),
    mu = mu - mean(res$mu[space$area])
  )

  jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

  cluster <- ggplot(divD) +
    geom_sf(aes(fill=mu)) +
    scale_fill_gradientn(colours = jet.colors(10000)) +
    facet_wrap(~key) +
    theme(legend.title = element_blank())

  ##############################################################################
  ###   output
  ##############################################################################

  out <- list(
    trends = trends,
    trends.all = trends.all,
    cluster = cluster
  )

  return(out)
}
