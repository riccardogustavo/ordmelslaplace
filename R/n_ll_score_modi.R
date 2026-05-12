n_ll_score_modi <- function( parm,
                          yi, Xi, Zi, Wi, thresh, b_loc,
                          b_scale, inv_Sigma_u, categ,
                          latent_error_dist ) {


  # --------------------------------------
  # Preparations
  # --------------------------------------

  # extract Info and parameters
  ni <- length( yi )
  n_re <- length( parm )
  n_thresh <- length( thresh )
  n_m_loc <- ncol( Zi )
  m_loc <- parm[ 1:n_m_loc ]
  m_scale <- parm[ n_re ]


  # initialize ll and score
  ll_i <- 0
  score_i <-  numeric( n_re )

  # --------------------------------------
  # Calculations
  # --------------------------------------
  # Rich: this Code could surely be more concise
  # e. g. b is the same for everything, deal with this later
  # there are also some inconsistensies in namings of variables

  # iterate over observations (of specific cluster)
  for ( j in 1:ni ) {

    Zij <- Zi[ j, ]
    Xij <- Xi[ j, ]
    Wij <- Wi[ j, ]
    yij <- yi[ j ]

    O <- ( Xij %*% b_loc ) + ( Zij %*% m_loc )

    #index <- which( categ == yij )
    index <- yij

    log_b <- -0.5 * ( ( Wij %*% b_scale ) + m_scale )
    b <- exp( max( log_b, -700 ) )

    # ++++++++++++++++++++++++++++++++++++
    # Reference-Category
    # ++++++++++++++++++++++++++++++++++++

    if ( index > n_thresh ) {

      a <- thresh[ index - 1 ] - O
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

      Fp_ijm1 <- switch( latent_error_dist,
                         "logit" = F_ijm1 * ( 1 - F_ijm1 ),
                         "probit" = dnorm( eta_ijm1 )
      )
      Fp_ijm1 <- as.numeric( Fp_ijm1 )

      g_loc <- -Zij * as.numeric( b )
      g_scale <- -0.5 * eta_ijm1
      g <- c( g_loc, g_scale )
      score_ij <- -( Fp_ijm1 * g )


    # ++++++++++++++++++++++++++++++++++++
    # First / Lowest Category
    # ++++++++++++++++++++++++++++++++++++

    } else if ( index == 1 ) {

      a <- thresh[ index ] - O
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

      g_loc <- -Zij * as.numeric( b )
      g_scale <- -0.5 * eta_ij
      g <- c( g_loc, g_scale )
      score_ij <- Fp_ij * g


    # ++++++++++++++++++++++++++++++++++++
    # middlish category
    # ++++++++++++++++++++++++++++++++++++

    } else {

      a <- thresh[ index ] - O
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

      g_loc <- -Zij * as.numeric( b )
      g_scale <- -0.5 * eta_ij
      g <- c( g_loc, g_scale )
      s <- Fp_ij * g

      a_m1 <- thresh[ index - 1 ] - O
      #b <- exp( -0.5 * ( ( Wij %*% b_scale ) + m_scale ) )
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

      g_loc_m1 <- -Zij * as.numeric( b )
      g_scale_m1 <- -0.5 * eta_ijm1
      g_m1 <- c( g_loc_m1, g_scale_m1 )
      s_m1 <- Fp_ijm1 * g_m1

      pij <- F_ij - F_ijm1
      score_ij <- s - s_m1

    }

    # ++++++++++++++++++++++++++++++++++++
    # Last calculations independent from Category
    # ++++++++++++++++++++++++++++++++++++
    pij <- max( pij, 1e-300 )
    ll_i <- ll_i + log( pij )

    score_ij <- score_ij / pij
    score_i <- score_i + score_ij


  }

  score_i <- score_i - ( inv_Sigma_u %*% parm )

  ll_i <-  ll_i - ( 0.5 * ( t( parm ) %*% inv_Sigma_u %*% parm ) )

  # leave out constants for log-Likelihoods
  #ll_i <- ll_i + ( 0.5 * ( log( det( inv_Sigma_u ) ) ) )
  #ll_i <- ll_i - ( ( n_re / 2 ) * log( 2 * pi ) )

  if ( any( is.na( score_i ) ) ) {
    stop(
      "NAN-Gradient evaluation occured while estimating modes (step 1)"
    )
  }

  return( list( objective = -ll_i, gradient = -score_i ) )

}
