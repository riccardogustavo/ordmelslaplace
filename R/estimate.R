
estimate <- function( start, yis, X, Xis, Zis, W, Wis, clust_ind, n_thresh, n_clusters,
                      categ,
                      dimension,
                      parm_table,
                      latent_error_dist,
                      SEs,
                      starting_values,
                      nlminb_ctrl ) {


  # --------------------------------------
  # Preparations
  # --------------------------------------

  VarCov <- NULL

  # set starting values
  if ( !is.null( starting_values ) ) {
    if ( length( starting_values ) != nrow( parm_table ) ) {
      stop("length of startvalues must match the number of parameters
           that need to be estimated.")
    }
    parm <- starting_values
  } else {
    parm <- c( rep( -1.5, n_thresh ), rep( 0, nrow( parm_table ) - n_thresh ) )
  }

  # initialize ll
  ll_old <- -Inf



  ### Rich: maybe algo works better without starting values for modes
  # make intial list for modes
  #modes <- vector ( "list", n_clusters )


  # --------------------------------------
  # Estimation with nlminb
  # --------------------------------------




  # did nlminb even converge on the last iteration?
  if ( est$convergence != 0 ) {
    warning( "nlminb did not converge: \n", est$message )
  }

  if ( mod_list$clusters_not_converged > 0 ) {
    warning( "nlminb did not converge estimating the modes of ",
             mod_list$clusters_not_converged, " clusters" )
  }

  if ( SEs ) {
    cat( "finished estimating parameters... now numerically estimating
         standard errors")
    Hess_parm <- numDeriv::jacobian( func = n_ll_grad_aghq, x = parm1,
                        agh_list = agh_list, mod_list = mod_list,
                        n_clusters = n_clusters,
                        Xis = Xis, Zis = Zis, Wis = Wis, yis = yis,
                        latent_error_dist = latent_error_dist
                        )
    # Hessian is already negative hessian
    VarCov <- solve( Hess_parm )

  }

  return( list( parm = parm1, ll = ll, VarCov = VarCov ) )
}
