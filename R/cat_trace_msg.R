
cat_trace_msg <- function( iter,
                            nlminb_msg,
                            cls_not_cnvg,
                            time_mod,
                            time_parm,
                            parm_old,
                            parm_new,
                            parm_change,
                            ll_old,
                            ll_new
) {
  cat(
    "\n\nIteration number", iter, "completed.\n\n",
    
    "nlminb convergence (took", time_parm, "seconds): ", nlminb_msg, "\n",
    "Mode estimation (took", time_mod, "seconds): ", cls_not_cnvg, 
    " Clusters did not converge\n\n\n",
    
    "Old parameters were: " , parm_old, "\n",
    "New parameters are : ", parm_new, "\n\n",
    
    "Maximal change in parm was: ", parm_change, "\n\n\n",
    
    "Old Log-Likelihood was: ", ll_old, "\n",
    "New Log-Likelihood is : ", ll_new, "\n\n",
    
    "Change in Log-Likelihood was: ", ll_new - ll_old, "\n",
    "-----------------------------------------------------------------",
    "\n"
  )
}