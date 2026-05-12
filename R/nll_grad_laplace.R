#' Function to compute Laplace-Approximation and gradient
#'
#'
nll_grad_laplace <- function(
    parm,
    mod_list,
    yis, Zis, Wis, Xis,
    X, W,
    n_clusters,
    categ,
    clust_ind,
    parm_table,
    latent_error_dist,
    nlminb_ctrl
) {

  # find modes for all clusters
  # (this should later be integrated in Loop over all
  # clusters, so for example make_parm_List is not called twice)
  mod_list <- find_modi_numeric(
    parm = parm,
    yis = yis,
    Xis = Xis,
    Zis = Zis,
    Wis = Wis,
    n_clusters = n_clusters,
    categ = categ,
    parm_table = parm_table,
    latent_error_dist = latent_error_dist,
    nlminb_ctrl = nlminb_ctrl
  )

  pl <- make_parm_list( parm = parm, parm_table = parm_table )

  b_loc <- pl$b_loc
  b_scale <- pl$b_scale
  thresh <- pl$thresh
  thresh_raw <- pl$thresh_raw
  L <- pl$L
  Sigma_u <- pl$Sigma_u
  inv_Sigma_u <- pl$inv_Sigma_u

  Xb <- X %*% b_loc
  Wb <- W %*% b_scale

  ll <- 0
  g <- numeric( length( parm ) )

  grad_Su <- make_grad_Su_list( L = L, parm_table = parm_table )
  n_Su <- length( grad_Su )

  for ( i in 1:n_clusters ) {

    Xi <- Xis[[ i ]]
    Xbi <- Xb[ clust_ind[[ i ]], ]
    Zi <- Zis[[ i ]]
    Wi <- Wis[[ i ]]
    Wbi <- Wb[ clust_ind[[ i ]], ]
    yi <- yis[[ i ]]

    # Calculate Gradient of u(theta)
    # ACHTUNG Vorzeichen
    d_u_theta <- -1 * numDeriv::jacobian(
      func = just_n_score_modi,
      x = parm,
      u = mod_list$modes[[ i ]],
      yi = yi,
      Xi = Xi,
      Zi = Zi,
      Wi = Wi,
      parm_table = parm_table,
      categ = categ,
      latent_error_dist = latent_error_dist
    )
    grad_u <- mod_list$ninv_mod_Hess[[ i ]] %*% d_u_theta

    # calculate gradient of log(det(omega))
    # wrapper für dlogdet(omega)/du und für dlogdet(omega)/dtheta bauen
    # dann dlogdet(omega)/du * grad_u + dlogdet(omega)/dtheta berechnen


    g_fy <- ll_score_fy( d = mod_list$modes[[ i ]],
                         yi = yi, Xi = Xi, Zi = Zi, Wi = Wi, Xbi = Xbi,
                         Wbi = Wbi,
                         thresh = thresh, thresh_raw = thresh_raw, b_loc = b_loc,
                         b_scale = b_scale,
                         latent_error_dist = latent_error_dist
                         )

    g_Phid <- ll_score_Phid(
      d = mod_list$modes[[ i ]],
      Sigma_u = Sigma_u,
      inv_Sigma_u = inv_Sigma_u,
      grad_Su = grad_Su, n_Su = n_Su
    )

    ll_i <- g_fy$ll + g_Phid$ll

    ll_i <- ll_i + 0.5 * log( det( mod_list$ninv_mod_Hess[[ i ]] ) )
    ll_i <- ll_i + ( ncol( Sigma_u ) / 2 ) * log( 2 * pi )

    ll <- ll + ll_i

    ll

    g_i <- c( g_fy$grad, g_Phid$grad )

    g <- g + g_i



  }


  return( list(
    objective = -ll,
    gradient = -g
  ))

}



  # ?? mit find_modi_numeric

