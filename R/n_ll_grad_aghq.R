
n_ll_grad_aghq <- function( parm, agh_list, n_clusters,
                          X, Xis, Zis, W, Wis, yis, clust_ind,
                          parm_table,
                          latent_error_dist ) {


  # --------------------------------------
  # Preparations
  # --------------------------------------

  n_par <- length( parm )

  pl <- make_parm_list( parm = parm, parm_table = parm_table )

  b_loc <- pl$b_loc
  b_scale <- pl$b_scale
  thresh <- pl$thresh
  thresh_raw <- pl$thresh_raw
  L <- pl$L
  Sigma_u <- pl$Sigma_u
  inv_Sigma_u <- pl$inv_Sigma_u

  # get grid, points and log-weights from AGH-List
  n_grid <- agh_list$n_grid
  pts <- agh_list$pts
  log_wgh <- agh_list$log_wgh

  # initialize ll und gradient
  ll <- 0
  g <- numeric( n_par )

  Xb <- X %*% b_loc
  Wb <- W %*% b_scale


  # --------------------------------------
  # Calculations
  # --------------------------------------

  # make List with derivatives of Sigma_u
  grad_Su <- make_grad_Su_list( L = L, parm_table = parm_table )
  n_Su <- length( grad_Su )

  # iterate over clusters
  for ( i in 1:n_clusters ) {

    Xi <- Xis[[ i ]]
    Xbi <- Xb[ clust_ind[[ i ]], ]
    Zi <- Zis[[ i ]]
    Wi <- Wis[[ i ]]
    Wbi <- Wb[ clust_ind[[ i ]], ]
    yi <- yis[[ i ]]

    #mod <- mod_list$modes[[ i ]]
    #ninv_Hess_mod <- mod_list$ninv_mod_Hess[[ i ]]

    # calc Cholesky of negative inverse Hessian of modes
    #chol_ninv_Hess_mod <- t( chol( ninv_Hess_mod ) )


    # initialize vector / matrix for ll / gradients at grid-spots
    ll_is <- rep( NA, times = n_grid )
    g_i <- matrix(NA, nrow = n_grid, ncol = n_par )

    # iterate over grid of AGH-Points
    for ( j in 1:n_grid ) {

      # adjust points
      d <- agh_list$adj_pts[[ i ]][ j, ]

      # ll = f(y) + Phi(d)
      # calculate score and ll for f(y)
      g_fy <- ll_score_fy( d = d, yi = yi, Xi = Xi, Zi = Zi, Wi = Wi, Xbi = Xbi,
                           Wbi = Wbi,
                           thresh = thresh, thresh_raw = thresh_raw, b_loc = b_loc,
                           b_scale = b_scale,
                           latent_error_dist = latent_error_dist )

      # calculate score and ll for Phi(d)
      g_Phid <- ll_score_Phid( d = d, Sigma_u = Sigma_u,
                               inv_Sigma_u = inv_Sigma_u,
                               grad_Su = grad_Su, n_Su = n_Su )


      ll_is[ j ] <- g_fy$ll + g_Phid$ll

      # add normal distribution term for enumerator
      ll_is[ j ] <- ll_is[ j ] + ( 0.5 * sum( pts[ j, ]^2 ) )

      # add log-weights
      ll_is[ j ] <- ll_is[ j ] + log_wgh[ j ]

      # put them together as one gradient vector
      g_i[ j, ] <- c( g_fy$grad, g_Phid$grad )

    }

    # log-sum-exp Trick
    a <- max( ll_is )
    b <- exp( ll_is - a )
    c <- sum( b )

    ll_i <- a + log ( c )
    #ll_i <- ll_i + ( 0.5 * log( det( ninv_Hess_mod ) ) )

    ll <- ll + ll_i

    for ( j in 1:n_par ) {

      g[ j ] <- g[ j ] + ( sum( b * g_i[ ,j ] ) / c )

    }

  }

  return( list( objective = -ll, gradient = -g ) )
}
