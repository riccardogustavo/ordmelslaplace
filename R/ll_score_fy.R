
ll_score_fy <- function( d, yi, Xi, Xbi, Zi, Wi, Wbi, thresh, thresh_raw,
                             latent_error_dist ) {


  # --------------------------------------
  # Preparations
  # --------------------------------------

  #extract Info and parameters
  ni <- length( yi )
  n_thresh <- length( thresh )

  n_re <- length( d )
  n_m_loc <- ncol( Zi )
  m_loc <- d[ 1:n_m_loc ]
  m_scale <- d[ n_re ]


  # initialize ll and score
  ll_i <- 0
  g_i <-  numeric( n_thresh + ncol(Xi) + ncol(Wi) )


  # --------------------------------------
  # Calculations
  # --------------------------------------

  for ( j in 1:ni ) {

    Zij <- Zi[ j, ]
    Xij <- Xi[ j, ]
    Wij <- Wi[ j, ]
    yij <- yi[ j ]

    ## beachte:
    # O ist abhängig von m_loc ist abhängig von parm
    # b ist abhängig von m_scale ist abhängig von parm

    O <- Xbi[ j, ] + ( Zij %*% m_loc )


    log_b <- -0.5 * ( Wbi[ j, ] + m_scale )
    b <- exp( max( log_b, -700 ) )

    # ++++++++++++++++++++++++++++++++++++
    # Reference-Category
    # ++++++++++++++++++++++++++++++++++++

    if ( yij > n_thresh ) {

      # Likelihood
      a <- thresh[ yij - 1 ] - O
      #b <- exp( -0.5 * ( ( Wij %*% b_scale ) + m_scale ) )
      eta_ijm1 <- a * b

      F_ijm1 <- switch( latent_error_dist,
                        "logit" = plogis( eta_ijm1 ),
                        "probit" = pnorm( eta_ijm1 )
      )
      F_ijm1 <- as.numeric( F_ijm1 )

      #pij <- 1 - F_ijm1
      #more stable then pij <- 1 - Fijm1 :
      pij <- switch( latent_error_dist,
                     "logit" = plogis( -eta_ijm1 ),
                     "probit" = pnorm( eta_ijm1, lower.tail = F )
      )

      # Gradient for b_loc
      Fp_ijm1 <- switch( latent_error_dist,
                         "logit" = F_ijm1 * ( 1 - F_ijm1 ),
                         "probit" = dnorm( eta_ijm1 )
      )
      Fp_ijm1 <- as.numeric( Fp_ijm1 )

      g_t <- c( 1, exp( thresh_raw[ 2:( yij - 1 ) ] ) )
      g_t <- g_t * as.numeric( b )

      g_b_loc <- -Xij * as.numeric( b )

      # Gradient for b_scale
      g_b_scale <- ( -0.5 * Wij ) * eta_ijm1

      g_ij <- c(g_t, g_b_loc, g_b_scale ) * Fp_ijm1
      # negative, because its the reference-categorie
      g_ij <- -g_ij


    # ++++++++++++++++++++++++++++++++++++
    # First / Lowest Category
    # ++++++++++++++++++++++++++++++++++++

    } else if ( yij == 1 ) {

      a <- thresh[ yij ] - O
      #b <- exp( -0.5 * ( ( Wij %*% b_scale ) + m_scale ) )
      eta_ij <- a * b

      F_ij <- switch( latent_error_dist,
                      "logit" = plogis( eta_ij ),
                      "probit" = pnorm( eta_ij )
      )
      F_ij <- as.numeric( F_ij )

      pij <- F_ij

      Fp_ij <- switch( latent_error_dist,
                       "logit" = F_ij * ( 1 - F_ij ),
                       "probit" = dnorm( eta_ij )
      )
      Fp_ij <- as.numeric( Fp_ij )

      g_t <- numeric( n_thresh )
      g_t[ 1 ] <- 1
      g_t <- g_t * as.numeric( b )

      g_b_loc <- -Xij * as.numeric( b )

      # Gradient for b_scale
      g_b_scale <- ( -0.5 * Wij ) * eta_ij

      g_ij <- c(g_t, g_b_loc, g_b_scale ) * Fp_ij



    # ++++++++++++++++++++++++++++++++++++
    # middlish category
    # ++++++++++++++++++++++++++++++++++++

    } else {

      a <- thresh[ yij ] - O
      #b <- exp( -0.5 * ( ( Wij %*% b_scale ) + m_scale ) )
      eta_ij <- a * b

      F_ij <- switch( latent_error_dist,
                      "logit" = plogis( eta_ij ),
                      "probit" = pnorm( eta_ij )
      )
      F_ij <- as.numeric( F_ij )

      Fp_ij <- switch( latent_error_dist,
                       "logit" = F_ij * ( 1 - F_ij ),
                       "probit" = dnorm( eta_ij )
      )
      Fp_ij <- as.numeric( Fp_ij )

      g_t <- numeric( n_thresh )
      g_t[ 1:yij ] <- c( 1, exp( thresh_raw[ 2:( yij ) ] ) )
      g_t <- g_t * as.numeric( b )


      g_b_loc <- -Xij * as.numeric( b )

      # Gradient for b_scale
      g_b_scale <- ( -0.5 * Wij ) * eta_ij

      g_ij <- c(g_t, g_b_loc, g_b_scale ) * Fp_ij


      a_m1 <- thresh[ yij - 1 ] - O
      #b_m1 <- exp( -0.5 * ( ( Wij %*% b_scale ) + m_scale ) )
      eta_ijm1 <- a_m1 * b

      F_ijm1 <- switch( latent_error_dist,
                        "logit" = plogis( eta_ijm1 ),
                        "probit" = pnorm( eta_ijm1 )
      )
      F_ijm1 <- as.numeric( F_ijm1 )

      Fp_ijm1 <- switch( latent_error_dist,
                         "logit" = F_ijm1 * ( 1 - F_ijm1 ),
                         "probit" = dnorm( eta_ijm1 )
      )
      Fp_ijm1 <- as.numeric( Fp_ijm1 )

      g_t_m1 <- numeric( n_thresh )
      if (yij > 2) {
        g_t_m1[ 1:( yij - 1 ) ] <- c( 1, exp( thresh_raw[ 2:( yij - 1 ) ] ) )
      } else {
        g_t_m1[1] <- 1
      }
      g_t_m1 <- g_t_m1 * as.numeric( b )

      #g_b_loc <- -Xij * as.numeric( b )

      # Gradient for b_scale
      g_b_scale_m1 <- ( -0.5 * Wij ) * eta_ijm1

      g_ij_m1 <- c(g_t_m1, g_b_loc, g_b_scale_m1 ) * Fp_ijm1

      pij <- F_ij - F_ijm1
      g_ij <- g_ij - g_ij_m1

    }

    pij <- max( pij, 1e-300 )
    ll_i <- ll_i + log( pij )
    g_ij <- g_ij / pij
    g_i <- g_i + g_ij


  }

  if ( any( is.na( g_i ) ) ) {
    stop(
      "NAN-Gradient evaluation occured while estimating parameters (step 2)"
      )
  }

  return( list( ll = ll_i, grad = g_i ) )

}
