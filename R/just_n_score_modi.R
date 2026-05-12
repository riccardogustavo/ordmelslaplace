## wrapper function to compute gradient with numDeriv

just_n_score_modi <- function( parm, u,
                             yi, Xi, Zi, Wi, categ,
                             parm_table,
                             latent_error_dist ) {

  # thresh noch
  pl <- make_parm_list( parm = parm, parm_table = parm_table )
  thresh <- pl$thresh
  b_loc <- pl$b_loc
  b_scale <- pl$b_scale
  inv_Sigma_u <- pl$inv_Sigma_u


  return(
    n_ll_score_modi(
      parm = u,
      yi = yi,
      Xi = Xi,
      Zi = Zi,
      Wi = Wi,
      thresh = thresh,
      b_loc = b_loc,
      b_scale = b_scale,
      inv_Sigma_u = inv_Sigma_u,
      categ = categ,
      latent_error_dist = latent_error_dist
    )$gradient
  )
}
