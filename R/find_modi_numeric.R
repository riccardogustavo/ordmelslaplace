
find_modi_numeric <- function(  parm, yis, Xis, Zis, Wis, n_clusters,
                                categ, parm_table,
                                latent_error_dist,
                                nlminb_ctrl ) {
  
  # --------------------------------------
  # Preparations
  # --------------------------------------
  
  # unpack and parametrize parameters
  pl <- make_parm_list( parm = parm, parm_table = parm_table )
  
  thresh <- pl$thresh
  b_loc <- pl$b_loc
  b_scale <- pl$b_scale
  inv_Sigma_u <- pl$inv_Sigma_u
  n_modi <- ncol( inv_Sigma_u )
  
  ### Rich: maybe algo works better without using modes as starting values
  # initialize list to store modes in
  modes <- vector ( "list", n_clusters )
  
  # initialize list to store negative inverses of Hessians
  ninv_mod_Hess <- vector( "list", length = n_clusters )
  
  # initialize counter for number of clusters that did not converge
  clusters_not_converged <- 0
  
  # --------------------------------------
  # Estimation
  # --------------------------------------
  
  # iterate over clusters
  for ( i in 1:n_clusters ) {
    
    ### Richard: Maybe the algorithm works better, if no start-values for modi
    ### are provided
    # # set modi to last mode-estimate. In first iteration, set to zero
    # modi <- modes[[ i ]]
    # if ( is.null( modi ) ) {
    #   modi <- rep( 0, times = n_modi )
    # }
    modi <- numeric( n_modi )

    Xi <- Xis[[ i ]]
    Zi <- Zis[[ i ]]
    Wi <- Wis[[ i ]]
    yi <- yis[[ i ]]
    
    # make functions for nlminb to minimize
    obj <- make_Opt_Fun( func = n_ll_score_modi, 
                         argslist = list( 
                            yi = yi, Xi = Xi, Zi = Zi, Wi = Wi,
                            thresh = thresh, b_loc = b_loc,
                            b_scale = b_scale, inv_Sigma_u = inv_Sigma_u,
                            categ = categ, 
                            latent_error_dist = latent_error_dist )
                   )
    
    # estimate optfun using nlminb
    est <- suppressWarnings( 
        stats::nlminb( start = modi, objective = obj$fn, gradient = obj$gr, 
                       control = nlminb_ctrl
                       )
      )
    
    # monitor number of non-convergent clusters
    clusters_not_converged <- clusters_not_converged + est$convergence
    
    # generate negative inverse Hessian and add to list
    Hessian_i <- ll_score_Hess_modi(parm = est$par,
                                               yi = yi, Xi = Xi, Zi = Zi, Wi = Wi,
                                               thresh = thresh, b_loc = b_loc,
                                               b_scale = b_scale, inv_Sigma_u = inv_Sigma_u,
                                               categ = categ, 
                                               latent_error_dist = latent_error_dist )$Hessian
    ninv_mod_Hess[[ i ]] <- solve( -Hessian_i )
    
    # also add modi to list
    modes[[ i ]] <- est$par
    
  }
  
  # return new modes, negative inverse Hessians, number of non-convergend clusters
  return( list( modes = modes, ninv_mod_Hess = ninv_mod_Hess,
                clusters_not_converged = clusters_not_converged ) )
  
}
  


