
make_Opt_Fun <- function( func, argslist ) {
  
  # the cache
  last_parm <- NULL
  last_val <- NULL
  
  # function to check if parm has changed
  # if parm has change, reevaluate
  # this function can always access the cache, because
  # its scope
  eval_cached <- function( parm ) {
    
    if ( !identical( parm, last_parm ) ) {
      last_parm <<- parm
      last_val <<- do.call( func, c( list( parm = parm ), argslist ) )
    }
    
    return ( last_val )
  }
  
  # returns a list of functions
  return( list (
    fn = function( parm ) eval_cached( parm )$objective,
    gr = function( parm ) eval_cached( parm )$gradient
  ))
  
}
