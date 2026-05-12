
make_grad_Su_list <- function( L, parm_table ) {
  
  # this function expects that L is a log-Cholesky Matrix,
  # meaning diag(L) is already the exp
  
  pt_Su <- parm_table[ parm_table$type == "Sigma_u", ]
  n_Sigm_u <- nrow( pt_Su )
  
  grad_Su <- vector( "list", n_Sigm_u )
  
  for ( i in 1:n_Sigm_u ) {
    
    p1 <- pt_Su$pos1[ i ]
    p2 <- pt_Su$pos2[ i ]
    g_L <- 0 * L
    
    if ( p1 == p2 ) {
      
      g_L[ p1, p2 ] <- L[ p1, p2 ]
      
    } else {
      
      g_L[ p1, p2 ] <- 1
      
    }
    
    grad_Su[[ i ]] <- ( g_L %*% t( L ) ) + ( L %*% t( g_L ) )
  }
  return( grad_Su )
}
