
ll_score_Hess_modi <- function( parm,
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


  # initialize Hessian, score and ll
  Hessian_i <- matrix( 0, ncol = n_re, nrow = n_re )
  score_i <- numeric( n_re )
  ll_i <- 0

  # --------------------------------------
  # Calculations
  # --------------------------------------

  for ( j in 1:ni ) {

    Zij <- Zi[ j, ]
    Xij <- Xi[ j, ]
    Wij <- Wi[ j, ]
    yij <- yi[ j ]

    O <- ( Xij %*% b_loc ) + ( Zij %*% m_loc )

    #index <- which( categ == yij )
    index <- yij

    # ++++++++++++++++++++++++++++++++++++
    # Reference-Category
    # ++++++++++++++++++++++++++++++++++++

    if ( index > n_thresh ) {

      ll_s_H_ij <- ll_s_H_mod_refcat( index = index, thresh = thresh, O = O,
                                     Zij = Zij, Wij = Wij, m_loc = m_loc ,
                                     m_scale = m_scale, b_scale = b_scale,
                                     n_re = n_re, n_m_loc = n_m_loc,
                                     latent_error_dist = latent_error_dist )


    # ++++++++++++++++++++++++++++++++++++
    # First / Lowest Category
    # ++++++++++++++++++++++++++++++++++++

    } else if ( index == 1 ) {

      ll_s_H_ij <- ll_s_H_mod_firstcat( index = index, thresh = thresh, O = O,
                                       Zij = Zij, Wij = Wij, m_loc = m_loc ,
                                       m_scale = m_scale, b_scale = b_scale,
                                       n_re = n_re, n_m_loc = n_m_loc,
                                       latent_error_dist = latent_error_dist )


    # ++++++++++++++++++++++++++++++++++++
    # middlish category
    # ++++++++++++++++++++++++++++++++++++

    } else {

      ll_s_H_ij <- ll_s_H_mod_middlecat( index = index, thresh = thresh, O = O,
                                         Zij = Zij, Wij = Wij, m_loc = m_loc ,
                                         m_scale = m_scale, b_scale = b_scale,
                                         n_re = n_re, n_m_loc = n_m_loc,
                                         latent_error_dist = latent_error_dist )

    }

    # Add up contributions of clusters

    Hessian_i <- Hessian_i + ll_s_H_ij$Hessian_ij

    score_i <- score_i + ll_s_H_ij$score_ij

    ll_i <- ll_i + ll_s_H_ij$ll_ij

  }

  score_i <- score_i - ( inv_Sigma_u %*% parm )

  Hessian_i <- Hessian_i - inv_Sigma_u

  ll_i <-  ll_i - ( 0.5 * ( t( parm ) %*% inv_Sigma_u %*% parm ) )
  ## leave out constants for log-Likelihoods
  #ll_i <- ll_i + ( 0.5 * ( log( det( inv_Sigma_u ) ) ) )
  #ll_i <- ll_i - ( ( n_re / 2 ) * log( 2 * pi ) )

  return( list( objective = ll_i, gradient = score_i, Hessian = Hessian_i ) )

}
