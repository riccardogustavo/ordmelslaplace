
make_AGH_list <- function( dimension, nPoints ) {
  
  ### 1. Ein Gitter an statischen Quadraturpunkten erstellen ###
  
  # Gridmatrix über die möglichen Kombinationen an Quadraturpunkten erstellen
  idx <- as.matrix( expand.grid( rep( list( 1:nPoints ), dimension ) ) )
  
  n_grid <- nPoints^dimension
  
  # Einmal Quadraturpunkte und Gewichte erstellen (dise sind für jede Dimension gleich)
  val <- fastGHQuad::gaussHermiteData( nPoints )
  
  # Gitter an Quadraturpunkten erstellen
  pts <- matrix( val$x[ idx ], nrow = nrow( idx ), ncol = dimension )
  
  # Ersten Schritt der Umskalierung vornehmen, also mit sqrt(2) multiplizieren
  pts <- sqrt( 2 ) * pts
  
  # Gitter an Gewichten erstellen
  wgh <- matrix( val$w[ idx ], nrow = nrow( idx ), ncol = dimension )
  
  # driekt log der Gewichte errechnen und aufsummieren
  log_wgh <- rowSums( log( wgh ) ) - ( dimension / 2 ) * log( pi )
  
  return( 
    list(
      n_grid = n_grid,
      pts = pts,
      log_wgh = log_wgh
    )
    )
  
}
