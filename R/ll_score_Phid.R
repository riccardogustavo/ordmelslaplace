
ll_score_Phid <- function( d, Sigma_u, inv_Sigma_u, grad_Su, n_Su ) {
  
  # Likelihood
  # notice that -0.5 * log(det(Sigma_u)) will be added later
  # notice that -dim/2 * log( 2*pi ) always disappears in AGHQ
  a <- ( t( d ) %*% inv_Sigma_u %*% d )
  b <- log( det( Sigma_u ) )
  
  ll <- -0.5 * ( a + b )
  
  # Gradient for all covariance-parameters
  
  grad <- numeric( n_Su )
  
  for ( i in 1:n_Su ) {
    
    A <- t( d ) %*% ( -inv_Sigma_u %*% grad_Su[[ i ]] %*% inv_Sigma_u ) %*% d
    B <- sum( inv_Sigma_u * t( grad_Su[[ i ]] ) )
    
    grad[ i ] <- -0.5 * ( A + B )
    
  }
  
  return( list( ll = ll, grad = grad ) )
  
}
