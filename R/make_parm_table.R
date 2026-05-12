make_parm_table <- function( X, W, Z , n_thresh ) {
  
  n_b_loc <- ncol( X )
  n_re_loc <- ncol( Z )
  n_b_scale <- ncol( W )
  
  # calc number of covariances, random scale effect included
  n_re <- n_re_loc + 1
  n_Sigma_u <- ( n_re * ( n_re + 1 ) ) / 2

  # create type-column
  type <- c( rep( "thresh", n_thresh ) )
  type <- c( type, rep( "b_loc", n_b_loc ) )
  type <- c( type, rep( "b_scale", n_b_scale ) )
  type <- c( type, rep( "Sigma_u", n_Sigma_u ) )
  
  # create pos1 and pos2 columns
  pos1 <- c( seq_len( n_thresh ), seq_len( n_b_loc ), seq_len( n_b_scale ) )
  pos2 <- c( rep( 1, n_thresh ), rep( 1, n_b_loc ), rep( 1, n_b_scale ) )
  
  for ( i in 1:n_re ) {
    
    pos1 <- c( pos1, i:n_re )    
    pos2 <- c( pos2, rep( i, times = n_re + 1 - i) )

  }
  
  idx <- 1:length( type )
  
  parm_table <- data.frame(
    type = type,
    pos1 = pos1,
    pos2 = pos2,
    idx = idx
  )
  
  return( parm_table )
  
}
