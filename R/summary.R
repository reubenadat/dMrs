# ----------
# Summary
# ----------

#' @title opt_sum
#' @description Summarizes the solutions post-optimization
#' @param OPT Output list from \code{run_analyses()}
#' @export
opt_sum = function(OPT){
	
	if( is.null(OPT) || length(OPT) == 0 ){
		cat("No solutions\n")
		return(NULL)
	}
	
	res = sapply(OPT,function(xx){
		c(xx$copula,ifelse(xx$upKAPPA == 0,"weibull","expweibull"))
	})
	res = smart_df(t(res))
	names(res) = c("COPULA","DIST")
	
	res$IDX = seq(nrow(res))
	res = res[,c("IDX","COPULA","DIST")]
	
	EST = sapply(OPT,function(xx){
		EST = xx$RES$cout$EST
		names(EST) = xx$RES$cout$PARS
		EST
	})
	EST = smart_df(t(EST))
	EST$theta = round(EST$theta,3)
	res = smart_df(res,EST)
	
	res$LL = sapply(OPT,function(xx){
		xx$RES$LL
	})
	res$NP = sapply(OPT,function(xx){
		xx$RES$nparams
	})
	res$BIC = sapply(OPT,function(xx){
		xx$RES$BIC
	})
	
	res$POST = exp(res$BIC - Rcpp_logSumExp(res$BIC))
	res$POST = round(res$POST,8)
	res
	
	return(res)
	
}

##