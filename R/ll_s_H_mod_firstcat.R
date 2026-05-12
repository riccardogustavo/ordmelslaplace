
ll_s_H_mod_firstcat <- function(  index, thresh, O, Zij, Wij, m_loc, m_scale, b_scale, 
                                  n_re, n_m_loc, latent_error_dist ) {
  
  # Likelihood
  a <- thresh[ index ] - O
  b <- exp( -0.5 * ( ( Wij %*% b_scale ) + m_scale ) )
  eta_ij <- a * b
  
  F_ij <- switch( latent_error_dist, 
                    "logit" = plogis( eta_ij ),
                    "probit" = pnorm( eta_ij )
  )
  F_ij <- as.numeric( F_ij )
  
  pij <- F_ij
  
  # score
  Fp_ij <- switch( latent_error_dist,
                     "logit" = F_ij * ( 1 - F_ij ),
                     "probit" = dnorm( eta_ij )
  )
  Fp_ij <- as.numeric( Fp_ij )
  
  g_loc <- -Zij * as.numeric( b )
  g_scale <- -0.5 * eta_ij
  g <- c( g_loc, g_scale )
  score_ij <- Fp_ij * g
  
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
  
  h1 <- -tcrossprod( score_ij ) / ( pij^2 )
  h2 <- ( Fpp_ij / pij ) * tcrossprod( g )
  h3 <- ( Fp_ij / pij ) * G
  
  Hessian_ij <- h1 + h2 + h3
  
  score_ij <- score_ij / pij
  
  ll_ij <- log( pij )
  
  return( list( ll_ij = ll_ij, score_ij = score_ij, Hessian_ij = Hessian_ij ) )
  
}
