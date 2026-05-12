
make_parm_list <- function( parm, parm_table ) {
  
  # --------------------------------------
  # unpack parameters
  # --------------------------------------
  
  thresh_idx <- parm_table[ parm_table$type == "thresh", ]$idx
  thresh <- parm[ thresh_idx ]
  n_thresh <- length( thresh )
  
  b_loc_idx <- parm_table[ parm_table$type == "b_loc", ]$idx
  b_loc <- parm[ b_loc_idx ]
  
  b_scale_idx <- parm_table[ parm_table$type == "b_scale", ]$idx
  b_scale <- parm[ b_scale_idx ]
  
  Sigma_u_dim <- max( parm_table[ parm_table$type == "Sigma_u", ]$pos1 )
  L <- matrix( 0, ncol = Sigma_u_dim, nrow = Sigma_u_dim )
  
  for ( i in parm_table[ parm_table$type == "Sigma_u", ]$idx ) {
    L[ parm_table$pos1[ i ], parm_table$pos2[ i ] ] <- parm[ i ]
  }
  
  
  # --------------------------------------
  # parametrize
  # --------------------------------------
  
  # sum-exp - parametrization for thresholds
  thresh_raw <- thresh
  for ( i in 2:n_thresh ) {
    thresh[ i ] <- thresh[ i - 1 ] + exp( thresh[ i ] )
  }
  
  # Log-Cholesky parametrization for RE-Covariances
  log_L <- L
  diag( L ) <- exp( diag( L ) )
  Sigma_u <- L %*% t( L )
  inv_Sigma_u <- solve( Sigma_u )
  
  return( list( thresh = thresh,
                thresh_raw = thresh_raw,
                b_loc = b_loc,  
                b_scale = b_scale,
                log_L = log_L,
                L = L,
                Sigma_u = Sigma_u,
                inv_Sigma_u = inv_Sigma_u ) )
  
}
