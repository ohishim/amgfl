#' @title Estimation for future observations (v0.0.2)
#' @description \code{amgfl.pred} This function provides predictive values for future observations
#'
#' @importFrom magrittr %>%
#' @importFrom purrr exec map map_dbl
#' @importFrom sf st_point st_intersects
#'
#' @param res the output of `amgfl`
#' @param X1 future observation: a matrix of explanatory variables, without an intercept, and with colnames,
#'     where the variables are fitted with a linear model;
#'     the intercept is automatically added
#' @param X2 future observation: a matrix of explanatory variables with colnames,
#'     where the variables are fitted with an additive model
#' @param area future observation: a vector of integers expressing sub-areas for each individual;
#'     if `Sp` is given, `area` can be omitted
#' @param Sp future observation: a matrix with two columns expressing coordinates of observed points;
#'     if `area` is given, `Sp` can be omitted
#'
#' @return a vector of predictive values
#'
#' @export
#' @examples
#' #amgfl.pred(res, X1, X2, area, Sp)

amgfl.pred <- function(res, X1=NULL, X2, area=NULL, Sp=NULL){

  is.X1 <- !is.null(X1)

  .X1 <- X1
  .X2 <- X2

  nknots <- res$par$nknots
  scale1 <- res$par$scale1
  scale2 <- res$par$scale2
  KNOTS <- res$knots

  fx <- function(x){cbind(x, x^2, x^3)}

  fb <- function(x, knots, nknots){
    out <- lapply(1:nknots, function(j){
      ((x - knots[j])^3) * (x > knots[j])
    }) %>% do.call(cbind, .)
    return(out)
  }

  if(is.X1 & scale1){X1 <- scale(.X1, center=res$scale$X1$center, scale=res$scale$X1$scale)}
  if(scale2){X2 <- scale(.X2, center=res$scale$X2$center, scale=res$scale$X2$scale)}

  k2 <- ncol(X2)

  X <- map(1:k2, ~fx(X2[,.x])) %>% exec(cbind, !!!.)

  if(is.X1)
  {
    X <- cbind(X1, X)
  }

  X <- scale(X, center=res$center$X, scale=F)

  B <- map(1:k2, ~fb(X2[,.x], KNOTS[.x,], nknots)) %>%
    exec(cbind, !!!.) %>% scale(center=res$center$B, scale=F)

  if(!is.null(area))
  {
    Mu <- res$mu[area]
  } else
  {
    div.pol <- res$space$div.pol
    n <- nrow(X2)
    Sp <- map(1:n, ~st_point(Sp[.x,])) %>% st_sfc
    Mu <- st_intersects(Sp, div.pol) %>%
      map_dbl(~.x[1]) %>% res$mu[.]
  }

  return(drop(
    X%*%res$beta + B%*%res$alpha + Mu
  ))
}
