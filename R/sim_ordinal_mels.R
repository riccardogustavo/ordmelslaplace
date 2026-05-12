#' @export
sim_ordinal_mels <- function( n_clusters = 50, ni = 30,
                             xmean = -0.3, xsd = 1, b = c( -0.2, 0.8 ),
                             y_star_sd = 0.5, cutoffs = c( -2.3, -0.7, 1.3 ),
                             b_scale = -0.4 ,
                             covm = matrix( c( 1, -0.1, -0.15, -0.1, 0.6, -0.2, -0.15, -0.2, 0.5 ), nrow = 3 )

) {
  n <- ni * n_clusters
  len_cutoff <- length( cutoffs )

  X <- numeric( 0 )
  Z <- numeric( 0 )
  y_star <- numeric( 0 )

  ranef <- mvtnorm::rmvnorm(n = n_clusters, mean = rep( 0, 3 ), sigma = covm )
  ranef_loc <- ranef[ , 1:2 ]
  ranef_scale <- ranef[ , 3, drop = F ]
  id <- numeric(0)
  for ( cluster in 1:n_clusters ) {
    id <- c( id, rep( cluster, times = ni ) )
    xi <- rnorm( n = ni, mean = xmean, sd = xsd )
    Xi_loc <- matrix( c( rep( 1, times = ni ), xi ), ncol = 2, nrow = ni, byrow = F )
    Zi <- Xi_loc
    y_star_sd <- y_star_sd * exp( -0.5 * ( xi * b_scale ) )
    res_star_i <- rnorm( n = ni, mean = 0, sd = y_star_sd )
    y_star <- c( y_star, Xi_loc %*% b + Zi %*% ranef_loc[ cluster,  ] + res_star_i )
    X <- rbind( X, Xi_loc )
  }

  y <- rep( NA, times = n )
  for ( i in 1:len_cutoff ) {
    y[ is.na( y ) ] <- ifelse( y_star[ is.na( y ) ] <= cutoffs[ i ], i, y[ is.na( y )] )
  }
  y <- ifelse( is.na( y ), len_cutoff + 1, y )

  print( table( y ) )
  barplot(height = table( y ), names.arg = 1:( len_cutoff + 1 ) )

  df <- data.frame(
    id = id,
    y = y,
    x = X[ , 2 ]
    )

  return( df )

}

