
ll_s_H_mod_middlecat <- function( index, thresh, O, Zij, Wij, m_loc, m_scale, b_scale, 
                              n_re, n_m_loc, latent_error_dist ) {
  
  ### First part
  # Likelihood
  a <- thresh[ index ] - O
  b <- exp( -0.5 * ( ( Wij %*% b_scale ) + m_scale ) )
  eta_ij <- a * b
  
  F_ij <- switch( latent_error_dist, 
                  "logit" = plogis( eta_ij ),
                  "probit" = pnorm( eta_ij )
  )
  F_ij <- as.numeric( F_ij )

  
  # score
  Fp_ij <- switch( latent_error_dist,
                   "logit" = F_ij * ( 1 - F_ij ),
                   "probit" = dnorm( eta_ij )
  )
  Fp_ij <- as.numeric( Fp_ij )
  
  g_loc <- -Zij * as.numeric( b )
  g_scale <- -0.5 * eta_ij
  g <- c( g_loc, g_scale )
  s <- Fp_ij * g
  
  # Hessian
  Fpp_ij <- switch( latent_error_dist, 
                    "logit" = Fp_ij * ( 1 - ( 2 * F_ij ) ),
                    "probit" = -eta_ij * dnorm( eta_ij )
  )
  Fpp_ij <- as.numeric( Fpp_ij )
  
  G <- matrix( 0, nrow = n_re, ncol = n_re )
  G[ n_m_loc + 1, 1:n_m_loc ] <- -0.5 * t( g_loc )
  G[ 1:n_m_loc, n_m_loc +1 ] <- -0.5 * g_loc
  G[ n_m_loc + 1, n_m_loc + 1 ] <- -0.5 * g_scale
  
  h2_1 <- Fpp_ij * tcrossprod( g )
  h2_1 <- h2_1 + ( Fp_ij * G )
  
  ### Second part
  # Likelihood
  a_m1 <- thresh[ index - 1 ] - O
  b_m1 <- exp( -0.5 * ( ( Wij %*% b_scale ) + m_scale ) )
  eta_ijm1 <- a_m1 * b_m1
  
  F_ijm1 <- switch( latent_error_dist, 
                    "logit" = plogis( eta_ijm1 ),
                    "probit" = pnorm( eta_ijm1 )
  )
  F_ijm1 <- as.numeric( F_ijm1 )

  
  # score
  Fp_ijm1 <- switch( latent_error_dist,
                     "logit" = F_ijm1 * ( 1 - F_ijm1 ),
                     "probit" = dnorm( eta_ijm1 )
  )
  Fp_ijm1 <- as.numeric( Fp_ijm1 )
  
  g_loc_m1 <- -Zij * as.numeric( b_m1 )
  g_scale_m1 <- -0.5 * eta_ijm1
  g_m1 <- c( g_loc_m1, g_scale_m1 )
  s_m1 <- Fp_ijm1 * g_m1

  
  # Hessian
  Fpp_ijm1 <- switch( latent_error_dist, 
                      "logit" = Fp_ijm1 * ( 1 - ( 2 * F_ijm1 ) ),
                      "probit" = -eta_ijm1 * dnorm( eta_ijm1 )
  )
  Fpp_ijm1 <- as.numeric( Fpp_ijm1 )
  
  G_m1 <- matrix( 0, nrow = n_re, ncol = n_re )
  G_m1[ n_m_loc + 1, 1:n_m_loc ] <- -0.5 * t( g_loc_m1 )
  G_m1[ 1:n_m_loc, n_m_loc +1 ] <- -0.5 * g_loc_m1
  G_m1[ n_m_loc + 1, n_m_loc + 1 ] <- -0.5 * g_scale_m1
  
  h2_2 <- Fpp_ijm1 * tcrossprod( g_m1 )
  h2_2 <- h2_2 + ( Fp_ijm1 * G_m1 )
  
  
  # put them together
  pij <- F_ij - F_ijm1
  
  score_ij <- s - s_m1

  h1 <- tcrossprod( score_ij ) / -( pij^2 )
  h2 <- ( h2_1 - h2_2 ) / pij
  Hessian_ij <- h1 + h2
  
  score_ij <- score_ij / pij  
  
  ll_ij <- log( pij )
  
  return( list( ll_ij = ll_ij, score_ij = score_ij, Hessian_ij = Hessian_ij ) )

}
