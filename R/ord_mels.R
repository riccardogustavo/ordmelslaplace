#' @export
ord_mels <- function( fixed,
                      random,
                      fixed_scale,
                      latent_error_dist = c("logit", "probit"),
                      data,
                      SEs = F,
                      starting_values = NULL,
                      trace = T,
                      nlminb_ctrl = list( eval.max = 500, iter.max = 350 )
                      ) {

  # --------------------------------------
  # Argument- Checks
  # --------------------------------------

  if ( missing( data ) || !is.data.frame( data ) ) {
    stop( "'data' must be a data.frame" )
  }

  if ( missing( latent_error_dist ) ) {
    stop( "latent_error_dist must be provided" )
  }

  latent_error_dist <- match.arg( latent_error_dist )

  # make sure random and fixed are valid formulas
  fixed <- as.formula( fixed )
  random <- as.formula( random )
  fixed_scale <- as.formula( fixed_scale )

  # --------------------------------------
  # NA- handling
  # --------------------------------------

  # get rid of all NAs
  all_vars <- c(
    all.vars( fixed ),
    all.vars( random ),
    all.vars( fixed_scale )
    )
  all_vars <- unique( all_vars )

  complete_cases <- complete.cases( data[ , all_vars ] )
  n_omitted <- sum( !complete_cases )

  if ( n_omitted > 0 ) {
    cat( n_omitted, " observations had missing data and were removed")
  }

  data <- data[ complete_cases, ]

  if ( nrow( data ) == 0 ) {
    stop( "all observations have missing values")
  }

  # --------------------------------------
  # Build target and model matrices
  # --------------------------------------

  # extract target and make sure there is only one target
  target <- fixed[[ 2 ]]

  if ( length( target ) > 1 ) {
    stop( "Only one target variable should be specified in 'fixed'" )
  }

  y <- data[ , as.character( target ) ]

  y <- as.numeric( y )

  categ <- sort( unique( y ), decreasing = F )
  n_categ <- length( categ )

  # make y the indexes of categ (1st, 2nd, ..., ) categ
  y <- match( y, categ )

  if ( n_categ < 3 ) {
    stop( "target must have at least three unique values")
  }

  if ( n_categ > 10 ) {
    warning( "Target Variable has more then 10 unique values,
         are you sure you want to treat it as ordinal?")
  }

  n_thresh <- n_categ - 1

  # extract Model matrix for predictors, no intercept
  X <- model.matrix( fixed, data = data )

  # remove Intercept
  X <- X[ , colnames( X ) != "(Intercept)"]
  X <- as.matrix( X )
  if ( ncol( X ) == 0 ) {
  stop("Atleast one Predictor must be included (
       instead of an Intercept, threshhold parameters will be estimated
  )")
  }

  # extract cluster ids and make sure only one is provided
  #id <- random[[ 2 ]][[ 3 ]]
  if ( length( random ) != 2 ) {
    stop( "'random' must be a formula of form ~ <random_effects> | <cluster_id>")
  }

  id <- random[[ 2 ]][[ 3 ]]

  if ( length( id ) > 1 ) {
    stop( "Only one cluster-id-Variable should be provided in 'random'" )
  }

  id <- data[ , as.character( id ) ]

  # create Model Matrix for location random effects
  Z <- model.matrix( as.formula( call( "~", random[[ 2 ]][[ 2 ]] ) ), data = data )


  # create Model Matrix for scale fixed effects
  W <- model.matrix( fixed_scale, data = data )

  # remove Intercept
  W <- W[ , colnames( W ) != "(Intercept)" ]
  W <- as.matrix( W )
  if ( ncol( X ) == 0 ) {
    stop("Atleast one Predictor must be included in 'fixed_scale' (
       instead of an Intercept, threshhold parameters will be estimated
  )")
  }

  # --------------------------------------
  # Last preparations
  # --------------------------------------

  # make parm_table
  parm_table <- make_parm_table( X = X, Z = Z, W = W, n_thresh = n_thresh )

  # make list of model matrices for clusters
  clusters <- unique( id )
  n_clusters <- length( clusters )

  dimension <- ncol( Z ) + 1

  Xis <- Zis <- Wis <- yis <- vector( "list", n_clusters )

  clust_ind <- vector( "list", n_clusters )

  for ( i in 1:n_clusters ) {

    clust_ind[[ i ]] <- id == clusters[ i ]

    Xis[[ i ]] <- X[ clust_ind[[ i ]], , drop = F ]
    Zis[[ i ]] <- Z[ clust_ind[[ i ]], , drop = F ]
    Wis[[ i ]] <- W[ clust_ind[[ i ]], , drop = F ]
    yis[[ i ]] <- y[ clust_ind[[ i ]] ]

  }


  # --------------------------------------
  # Estimate
  # --------------------------------------

  result <- estimate( starting_values = starting_values, yis = yis, Xis = Xis,
                      X = X, Zis = Zis, W = W, Wis = Wis, clust_ind = clust_ind,
                     parm_table = parm_table, n_clusters = n_clusters,
                     categ = categ, n_thresh = n_thresh,
                     dimension = dimension,
                     latent_error_dist = latent_error_dist,
                     SEs = SEs,
                     nlminb_ctrl = nlminb_ctrl )

  ll <- result$ll

  pl <- make_parm_list( parm = result$parm, parm_table = parm_table )

  return( list (
    thresh = pl$thresh,
    b_loc = pl$b_loc,
    b_scale = pl$b_scale,
    Sigma_u = pl$Sigma_u,
    VarCov = result$VarCov
  ))

}
