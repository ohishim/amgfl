
# v0.4.1

#' @importFrom magrittr %>% divide_by
#' @importFrom purrr array_tree map exec map_dbl
#' @importFrom sf st_agr st_collection_extract st_convex_hull st_intersection
#'     st_intersects st_multipoint st_point st_sfc st_touches st_voronoi
#' @importFrom stats kmeans

divsp <- function(Sp, m=NULL, spB=NULL, subregion=NULL, seed=123){

  set.seed(seed)

  is.subregion <- !is.null(subregion)

  if(is.null(m))
  {
    if(is.subregion)
    {
      m <- table(subregion) %>% divide_by(100) %>% ceiling
    } else
    {
      m <- nrow(Sp) %>% divide_by(100) %>% ceiling
    }
  }

  ##############################################################################
  ###   spatial division
  ##############################################################################

  if(is.subregion)
  {
    .Sp <- as.data.frame(Sp) %>% split(subregion)
    q <- length(spB)

    div <- map(1:q, ~{
      if(m[.x] == 1)
      {
        return(spB[.x])
      } else
      {
        KM <- kmeans(.Sp[[.x]], m[.x])

        return(
          KM$centers %>% st_multipoint %>% st_voronoi(spB[.x]) %>%
            st_collection_extract(type="POLYGON") %>% st_intersection(spB[.x])
        )
      }
    }) %>% exec(c, !!!.) %>%
      st_sf(
        area = 1:sum(m),
        geometry = .
      )
  } else
  {
    ##############################################################################
    ###   k-means
    ##############################################################################

    KM <- kmeans(Sp, m)

    ##############################################################################
    ###   voronoi division
    ##############################################################################

    if(is.null(spB)){spB <- st_multipoint(Sp) %>% st_convex_hull %>% st_sfc}

    div <- KM$centers %>% st_multipoint %>% st_voronoi(spB) %>%
      st_collection_extract(type="POLYGON") %>% st_intersection(spB) %>% st_sfc %>%
      st_sf(
        area = 1:m,
        geometry = .
      )
  }

  sp <- array_tree(Sp, 1) %>% map(st_point) %>% st_sfc

  fig.div <- ggplot() +
    geom_sf(data=div) +
    geom_sf(data=sp, size=0.8, alpha=0.1)

  ##############################################################################
  ###   adjacent information
  ##############################################################################

  adj <- st_touches(div)
  adjD <- lapply(1:sum(m), function(j){
    data.frame(area=j, adj=adj[[j]])
  }) %>% exec(rbind, !!!.)

  ##############################################################################
  ###   output
  ##############################################################################

  return(list(
    area = st_intersects(sp, div) %>% map_dbl(~.x[1]),
    adjacent = adjD,
    div.pol = div,
    fig.div = fig.div
  ))
}
