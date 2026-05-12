
estimate <- function( start, yis, X, Xis, Zis, W, Wis, clust_ind, n_thresh, n_clusters,
                      categ,
                      nPoints, dimension,
                      parm_table,
                      maxiter,
                      conv_tol,
                      latent_error_dist,
                      SEs,
                      starting_values,
                      trace,
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

  # make aghq points
  agh_list <- make_AGH_list( dimension = dimension, nPoints = nPoints )


  ### Rich: maybe algo works better without starting values for modes
  # make intial list for modes
  #modes <- vector ( "list", n_clusters )


  # --------------------------------------
  # Actual Estimation Loop
  # --------------------------------------

  iter <- 0

  repeat {

    # Stop if maximum number of iteration is reached
    if ( iter > maxiter ) {
      stop( "Maximal number of iterations reached without convergence! ")
    }

    # Step 1:
    start_time_mod <- Sys.time()
    # estimate random-effects modi ( + Hessians ), given parm
    mod_list <- find_modi_numeric( parm = parm, yis = yis,
                      Xis = Xis, Zis = Zis, Wis = Wis,
                      n_clusters = n_clusters, categ = categ,
                      parm_table = parm_table,
                      latent_error_dist = latent_error_dist,
                      nlminb_ctrl = nlminb_ctrl )
    ### Rich: maybe algo works better without starting values for modes
    #modes <- mod_list$modes
    end_time_mod <- Sys.time()
    time_mod <- as.numeric( end_time_mod - start_time_mod, unit = "secs" )

    # Step 2:
    start_time_parm <- Sys.time()
    # adjust quadrature_points
    agh_list$adj_pts <- Map(
      function( modes, ninv_hess ) {
        Lt <- chol( ninv_hess )
        sweep( agh_list$pts %*% Lt, 2, modes, "+" )
      },
      mod_list$modes,
      mod_list$ninv_mod_Hess
    )
    # prepare functions to optimize in nlminb
    obj <- make_Opt_Fun( func = n_ll_grad_aghq,
                   argslist = list(
                     agh_list = agh_list,
                     n_clusters = n_clusters, clust_ind = clust_ind,
                     X = X, Xis = Xis, Zis = Zis, W = W, Wis = Wis, yis = yis,
                     parm_table = parm_table,
                     latent_error_dist = latent_error_dist
                   ))

    # estimate parm given the random effects ( + Hessians )
    est <- suppressWarnings(

        stats::nlminb( start = parm, objective = obj$fn,
                       gradient = obj$gr, control = nlminb_ctrl
                       )
      )

    end_time_parm <- Sys.time()
    time_parm <- as.numeric( end_time_parm - start_time_parm, unit = "secs" )


    parm1 <- est$par
    ll <- -est$objective

    # Check for convergence
    if ( max( abs( parm1 - parm ) ) < conv_tol ) {
      if ( max( abs( ll_old - ll ) ) < conv_tol ) {
        if ( trace ) {
          cat( "Converged!\n\n" )
        }
        break
      }
    }

    if ( trace ) {

      cat_trace_msg( iter = iter, nlminb_msg = est$message,
                     cls_not_cnvg = mod_list$clusters_not_converged,
                     time_mod = time_mod,
                     time_parm = time_parm,
                     parm_old = parm, parm_new = parm1,
                     parm_change = max( abs( parm1 - parm ) ),
                     ll_old = ll_old, ll_new = ll )
    }

    # update values
    iter <- iter + 1
    parm <- parm1
    ll_old <- ll



  }

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
