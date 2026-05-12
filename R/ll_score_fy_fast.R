
ll_score_fy <- function( d, yi, Xi, Xbi, Zi, Wi, Wbi, thresh, thresh_raw, b_loc,
                             b_scale, categ,
                             latent_error_dist ) {


  # --------------------------------------
  # Preparations
  # --------------------------------------

  # #extract Info and parameters

  n_b_loc <- length( b_loc )
  n_b_scale <- length( b_scale )
  n_thresh <- length( thresh )
  ni <- length( yi )
  n_re <- length( d )
  n_m_loc <- ncol( Zi )
  m_loc <- d[ 1:n_m_loc ]
  m_scale <- d[ n_re ]


  # # initialize ll and score
  # ll_i <- 0
  # g_i <-  numeric( n_thresh + length( b_loc ) + length( b_scale ) )


  # --------------------------------------
  # Calculations
  # --------------------------------------

  # Likelihood

  O <- Xbi + ( Zi %*% m_loc )
  O <- as.vector( O )

  log_b <- -0.5 * ( Wbi + m_scale )
  b <- exp( pmax( log_b, -700 ) )


  # thresh - O berechnen
  A <- outer( O, thresh, function( o, t ) t - o )
  eta <- sweep( A, 1, b, "*" )

  if ( latent_error_dist == "logit" ) {

    F_ <- plogis( eta )
    Fp_ <- F_ * ( 1 - F_ )

  } else if ( latent_error_dist == "probit" ) {

    F_ <- pnorm( eta )
    Fp_ <- dnorm( eta )

  }

  pij <- numeric( ni )

  ref <- yi > n_thresh
  first <- yi == 1
  mid <- !( ref | first )

  ### maybe more stable, but slower:
  # if ( any( ref ) ) {
  #   if ( latent_error_dist == "logit" ) {
  #     pij[ ref ] <- plogis( -eta[ ref, n_thresh ])
  #   } else if ( latent_error_dist == "probit" ) {
  #     pij[ ref ] <- pnorm( eta[ ref, n_thresh ], lower.tail = F )
  #   }
  # }

  #pij[ ref ] <- 1 - F_[ ref, n_thresh ]



  # Threshholds
  dt <- c( 1, exp( thresh_raw[ 2:n_thresh ] ) )


  G_thresh <- matrix( 0, ni, n_thresh )
  G_b_loc <- matrix( 0, ni, n_b_loc )
  G_b_scale <- matrix( 0, ni, n_b_scale )


  if ( any( first ) ) {

    rows <- which( first )

    # Likelihood
    pij[ first ] <- F_[ rows, 1 ]

    fp <- Fp_[ rows, 1 ]
    bb <- b[ rows ]

    # Thresholds
    G_thresh[ rows, 1 ] <- bb * fp

    # b_loc
    G_b_loc[ rows, ] <- ( -Xi[ rows, , drop = F ] * bb ) * fp

    # b_scale
    eta1 <- eta[ rows, 1 ]
    G_b_scale[ rows, ] <- (( -0.5 * Wi[ rows, , drop = F ] ) * eta1 ) * fp

  }

  if ( any( ref ) ) {

    rows <- which( ref )
    n_ref <- length( rows )

    # Likelihood
    #pij[ ref ] <- 1 - F_[ ref, n_thresh ]
    if (latent_error_dist == "logit") {
      pij[ ref ] <- plogis( -eta[ ref, n_thresh ] )
    } else if ( latent_error_dist == "probit" ) {
      pij[ref] <- pnorm( eta[ ref, n_thresh ], lower.tail = FALSE )
    }

    fp <- Fp_[ rows, n_thresh ]
    bb <- b[ rows ]

    # Threshholds
    Gt <- matrix( dt, nrow = n_ref, ncol = n_thresh, byrow = T )
    #Gt <- matrix( rep( dt, each = length( rows ) ), nrow = length( rows ) )
    Gt <- sweep( Gt, 1, bb * fp, "*" )
    G_thresh[ rows, ] <- -Gt

    # b_loc
    G_b_loc[ rows, ] <- -( ( -Xi[ rows, , drop = F ] ) * bb ) * fp

    # b_scale
    eta1 <- eta[ rows, n_thresh ]
    G_b_scale[ rows, ] <- -(( -0.5 * Wi[ rows, , drop = F ] ) * eta1 ) * fp


  }

  if ( any( mid ) ) {

    rows <- which( mid )
    cols <- yi[ mid ]
    n_mid <- length( rows )

    #Likelihood
    pij[ mid ] <- F_[ cbind( rows, cols ) ] - F_[ cbind( rows, cols - 1 ) ]

    #fp_1 <- Fp_[ rows, cols ]
    bb <- b[ rows ]

    eta_1 <- eta[ cbind( rows, cols ) ]
    eta_m1 <- eta[ cbind( rows, cols - 1 ) ]

    Fp_1 <- Fp_[ cbind( rows, cols ) ]
    Fp_m1 <- Fp_[ cbind( rows, cols - 1 ) ]

    # Threshholds
    # get boolean-indices for active parts of threshholds in categ and categ - 1
    sl <- seq_len( n_thresh )
    m_1 <- outer( cols, sl, ">=" )
    m_m1 <- outer( cols - 1, sl, ">=" )

    Gt <- matrix( dt, nrow = n_mid, ncol = n_thresh, byrow = T )
    Gt_1 <- m_1 * Gt
    Gt_m1 <- m_m1 * Gt

    Gt_1 <- sweep( Gt_1, 1, bb * Fp_1, "*" )
    Gt_m1 <- sweep( Gt_m1, 1, bb * Fp_m1, "*" )

    G_thresh[ rows, ] <- Gt_1 - Gt_m1


    # b_loc
    # -Xi * bb * Fp_1 - -Xi * bb * Fp_m1
    # ausklammern -Xi * bb * (Fp_1 - Fp_m1)
    G_b_loc[ rows, ] <- sweep( -Xi[ rows, , drop = F ], 1, bb * (Fp_1 - Fp_m1 ), "*" )

    # b_scale
    # Ausklammern:
    # Fp_1 * -0.5 * Wij * eta_1 - Fp_m1 * - 0.5 * Wij * eta_m1
    # = -0.5 * Wij * (Fp_1 * eta_1 - Fp_m1 * eta_m1)
    c <- -0.5 * ( ( Fp_1 * eta_1 ) - ( Fp_m1 * eta_m1 ) )
    G_b_scale[ rows, ] <- sweep( Wi[ rows, , drop = F ], 1, c, "*" )

  }



  # Gradient noch durhc pij teilen
  # --> nach der Summierung??

  pij <- pmax( pij, 1e-300 )
  ll_i <- sum( log( pij ) )

  g_t <- colSums( sweep( G_thresh, 1, pij, "/" ) )
  g_b_loc <- colSums( sweep( G_b_loc, 1, pij, "/" ) )
  g_b_scale <- colSums( sweep( G_b_scale, 1, pij, "/" ) )
  g_i <- c( g_t, g_b_loc, g_b_scale )


  if ( any( is.na( g_i ) ) ) {
    stop(
      "NAN-Gradient evaluation occured while estimating parameters (step 2)"
      )
  }

  return( list( ll = ll_i, grad = g_i ) )

}
