
ll_s_H_mod_refcat <- function( index, thresh, O, Zij, Wij, m_loc, m_scale, b_scale, 
                               n_re, n_m_loc, latent_error_dist ) {
  
  # Likelihood
  a <- thresh[ index - 1 ] - O
  b <- exp( -0.5 * ( ( Wij %*% b_scale ) + m_scale ) )
  eta_ijm1 <- a * b
  
  F_ijm1 <- switch( latent_error_dist, 
                   "logit" = plogis( eta_ijm1 ),
                   "probit" = pnorm( eta_ijm1 )
                   )
  F_ijm1 <- as.numeric( F_ijm1 )
  
  pij <- 1 - F_ijm1
  
  # score
  Fp_ijm1 <- switch( latent_error_dist,
                     "logit" = F_ijm1 * ( 1 - F_ijm1 ),
                     "probit" = dnorm( eta_ijm1 )
                     )
  Fp_ijm1 <- as.numeric( Fp_ijm1 )
  
  g_loc <- -Zij * as.numeric( b )
  g_scale <- -0.5 * eta_ijm1
  g <- c( g_loc, g_scale )
  score_ij <- -( Fp_ijm1 * g )
  
  # Hessian
  Fpp_ijm1 <- switch( latent_error_dist, 
                      "logit" = Fp_ijm1 * ( 1 - ( 2 * F_ijm1 ) ),
                      "probit" = -eta_ijm1 * dnorm( eta_ijm1 )
                      )
  Fpp_ijm1 <- as.numeric( Fpp_ijm1 )
  
  G <- matrix( 0, nrow = n_re, ncol = n_re )
  G[ n_m_loc + 1, 1:n_m_loc ] <- -0.5 * t( g_loc )
  G[ 1:n_m_loc, n_m_loc +1 ] <- -0.5 * g_loc
  G[ n_m_loc + 1, n_m_loc + 1 ] <- -0.5 * g_scale
  
  h1 <- -tcrossprod( score_ij ) / ( pij^2 )
  h2 <- ( Fpp_ijm1 / pij ) * tcrossprod( g )
  h3 <- ( Fp_ijm1 / pij ) * G
  
  Hessian_ij <- h1 - h2 - h3
  
  score_ij <- score_ij / pij
  
  ll_ij <- log( pij )
  
  return( list( ll_ij = ll_ij, score_ij = score_ij, Hessian_ij = Hessian_ij ) )
  
}
