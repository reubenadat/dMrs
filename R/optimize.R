# ----------
# Optimization functions
# ----------

opt_replicate = function(REP,param_grid,theta,
	upKAPPA,ncores = 1,gTHRES = 8e-2,verb,PLOT){

	if(FALSE){
		REP 				= REP
		theta 			= theta
		upKAPPA 		= upKAPPA
		gTHRES 			= gTHRES
		verb 				= verb
		ncores 			= ncores
		PLOT				= TRUE
		
	}
	
	verb 		= get_verb(verb = verb)
	PLOT		= get_PLOT(PLOT = PLOT)
	upTHETA = ifelse(is.na(theta),1,0)
	DATA 		= REP$DATA
	PARAMS 	= REP$PARAMS
	upPARS	= get_upPARS(upKAPPA = upKAPPA,
							THETA = theta)
	upTHETA	= upPARS[4]
	tmp_copula = ifelse(
		PARAMS$copula %in% c("Independent","Clayton"),
		"Clayton",PARAMS$copula)
	
	wrap_LL 	= function(PARS){
		# PARS = iPARS
		dMrs_cLL(XX = DATA$time,DELTA = DATA$delta,
			D2 = DATA$dens_t2,S2 = DATA$surv_t2,
			PARS = PARS,copula = tmp_copula,
			verb = !TRUE)
	}
	wrap_GRAD = function(PARS){
		# PARS = iPARS
		out = dMrs_cGRAD(XX = DATA$time,DELTA = DATA$delta,
			D2 = DATA$dens_t2,S2 = DATA$surv_t2,PARS = PARS,
			copula = tmp_copula,upPARS = upPARS)
		
		# out = grad(wrap_LL,PARS)
		c(out)
	}
	wrap_HESS = function(PARS){
		# PARS = iPARS
		dMrs_cHESS(XX = DATA$time,DELTA = DATA$delta,
			D2 = DATA$dens_t2,S2 = DATA$surv_t2,PARS = PARS,
			copula = tmp_copula,upPARS = upPARS)
	}
	
	# Estimate parameters
	param_grid = sort(unique(c(0,param_grid)))
	if( verb ) cat(sprintf("%s: Grid search for initial parameters ...\n",date()))
	if( upPARS[3] == 0 ){
		unc_KAPPA = 0
	} else {
		unc_KAPPA = param_grid
	}
	if( upTHETA == 0 ){
		if( tmp_copula == "Clayton" ) 	log_THETA = log(theta)
		if( tmp_copula == "Gumbel" )	log_THETA = log(theta - 1)
	} else if( upTHETA == 1 ){
		# log_THETA = c(-Inf,param_grid)
		log_THETA = param_grid
		log_THETA = sort(unique(log_THETA))
	}
	
	
	gout = dMrs_GRID(XX = DATA$time,DELTA = DATA$delta,
		D2 = DATA$dens_t2,S2 = DATA$surv_t2,
		log_ALPHA = param_grid,log_LAMBDA = param_grid,
		unc_KAPPA = unc_KAPPA,log_THETA = log_THETA,
		copula = tmp_copula,verb = verb,ncores = ncores)
	# str(gout)
	colnames(gout$DAT) = c("log_alpha1","log_lambda1",
		"unc_kappa1","log_theta","LL")
	gout$DAT[1:5,]
	gprof = get_PROFILE(GRID = gout$DAT)
	# dim(gprof); head(gprof)
	
	gopt = get_LOCAL_OPTS(GPROF = gprof,verb = verb)
	# dim(gopt); head(gopt)
	
	if( verb && PLOT ){
		cat(sprintf("%s: Two-dimensional plot of the profile log likelihood ...\n",date()))
		sub_gopt = gopt[order(-gopt$LL),]
		sub_gopt = sub_gopt[seq(min(c(3,nrow(sub_gopt)))),,drop = FALSE]
		print(plot_LL(GPROF = gprof,
			GOPT = sub_gopt,
			COPULA = PARAMS$copula))
	}
	
	# gopt = gprof
	gopt$fin_LL = NA; gopt$nGRAD = NA; gopt$fin_logA = NA
	gopt$fin_logL = NA; gopt$fin_logK = NA; gopt$fin_logT = NA
	# gopt
	
	if( verb ) cat(sprintf("%s: BFGS optimization ...\n",date()))
	nn = nrow(gopt)
	for(ii in seq(nn)){
		# ii = 1
		# if( verb ) cat(".")
		if( verb ) smart_progress(ii = ii,nn = nn,
			iter = 5,iter2 = 2e2)
		iPARS = as.numeric(gopt[ii,1:4]); iPARS
		dMrs_BFGS(XX = DATA$time,DELTA = DATA$delta,
			D2 = DATA$dens_t2,S2 = DATA$surv_t2,
			PARS = iPARS,copula = tmp_copula,upPARS = upPARS,
			max_iter = 2e2,eps = 1e-6,verb = FALSE)
		if(FALSE){
		
		upPARS = c(1,1,0,0)
		dMrs_BFGS(XX = DATA$time,DELTA = DATA$delta,
			D2 = DATA$dens_t2,S2 = DATA$surv_t2,
			PARS = iPARS,copula = tmp_copula,upPARS = upPARS,
			max_iter = 2e2,eps = 1e-6,verb = !FALSE)
		
		}
		tmp_LL = wrap_LL(PARS = iPARS)
		tmp_GR = wrap_GRAD(PARS = iPARS)
		
		gopt$nGRAD[ii] = Rcpp_norm(tmp_GR)
		gopt$fin_logA[ii] = iPARS[1]
		gopt$fin_logL[ii] = iPARS[2]
		gopt$fin_logK[ii] = iPARS[3]
		gopt$fin_logT[ii] = iPARS[4]
		gopt$fin_LL[ii] = round(tmp_LL,3)
	}
	if( verb ){
		cat("\n")
		print(head(gopt))
	}
	gopt = gopt[gopt$fin_LL != -999,,drop = FALSE]
	gopt_pre = gopt
	
	# Remove non-local optimum solutions
	gopt = gopt[which(gopt$nGRAD < gTHRES),]
	if( nrow(gopt) == 0 ){
		if( verb ) cat(sprintf("%s: No converged solutions! ...\n",date()))
		return(list(GRID = gout$DAT,GPROF = gprof,GOPT = gopt,
			GOPT_PRE = gopt_pre,
			out = NULL,cout = NULL,LL = NULL,GRAD = NULL,
			HESS = NULL,COVAR = NULL))
	}
	
	# Order remaining, remove duplicates
	gopt = gopt[order(-gopt$fin_LL),]
	gopt = gopt[!duplicated(round(gopt[,
		sprintf("fin_%s",c("LL","logA",
		"logL","logK","logT"))],3)),]
	rownames(gopt) = NULL
	# gopt
	
	# Estimate covariance
	if( verb ) cat(sprintf("%s: Get covariance ...\n",date()))
	gopt$neg_var = NA
	for(jj in seq(nrow(gopt))){
		# jj = 1
		iPARS = as.numeric(gopt[jj,
			paste0("fin_log",c("A","L","K","T"))]); iPARS
		hess = wrap_HESS(PARS = iPARS); hess
		nz = which(diag(hess) != 0)
		# print(hess)
		if( rcond(hess[nz,nz,drop = FALSE]) == 0 ){
			next
			# stop("Error with hessian 1: rcond = 0\n")
		}
		covar = matrix(0,4,4)
		covar[nz,nz] = solve(-hess[nz,nz,drop = FALSE])
		# covar
		gopt$neg_var[jj] = any(diag(covar) < 0)
	}
	gopt = gopt[!is.na(gopt$neg_var),,drop = FALSE]
	gopt = gopt[which(gopt$neg_var == FALSE),]
	if( nrow(gopt) == 0 ){
		if( verb ) cat(sprintf("%s: Negative variance(s)! ...\n",date()))
		return(list(GRID = gout$DAT,GPROF = gprof,GOPT = NULL,
			GOPT_PRE = gopt_pre,
			out = NULL,cout = NULL,LL = NULL,GRAD = NULL,
			HESS = NULL,COVAR = NULL))
	}
	gopt = gopt[which.max(gopt$fin_LL),]
	
	gopt$log_alpha1 = gopt$fin_logA
	gopt$log_lambda1 = gopt$fin_logL
	gopt$unc_kappa1 = gopt$fin_logK
	gopt$log_theta = gopt$fin_logT
	gopt
	
	iPARS = as.numeric(gopt[1,
		paste0("fin_log",c("A","L","K","T"))]); iPARS
	if(FALSE){
		iPARS = true_PARS
	}
	hess = wrap_HESS(PARS = iPARS); hess
	nz = which(diag(hess) != 0)
	if( rcond(hess[nz,nz,drop = FALSE]) == 0 ){
		print(hess)
		cat("Error with hessian 1: rcond = 0\n")
		return(list(GRID = gout$DAT,GPROF = gprof,GOPT = NULL,
			GOPT_PRE = gopt_pre,
			out = NULL,cout = NULL,LL = NULL,GRAD = NULL,
			HESS = NULL,COVAR = NULL))
	}
	covar = matrix(0,4,4)
	covar[nz,nz] = solve(-hess[nz,nz,drop = FALSE])
	covar
	
	out = smart_df(PARS = c("log_alpha1","log_lambda1",
		"unc_kappa1","log_theta"),
		EST = c(iPARS),SE = c(sqrt(diag(covar))))
	z_alpha = qnorm(0.975)
	out$lowCI 	= out$EST - z_alpha * out$SE
	out$highCI 	= out$EST + z_alpha * out$SE
	out
	
	EST = exp(out$EST)
	if( PARAMS$copula == "Gumbel" ) EST[4] = EST[4] + 1
	
	cout = smart_df(PARS = c("alpha1","lambda1","kappa1","theta"),
		EST = EST,SE = NA)
	
	# Delta method: Calculate nabla g()
	nab_par	= rep(NA,4)
	nab_par[1:3] = cout$EST[1:3] # alpha, lambda, kappa
	nab_par[4] = ifelse(tmp_copula == "Clayton",
		cout$EST[4],exp(out$EST[4]))
	nab_par = diag(nab_par)
	
	cout$SE 			= sqrt(diag(nab_par %*% covar %*% nab_par))
	cout$lowCI 		= cout$EST - z_alpha * cout$SE
	cout$highCI 	= cout$EST + z_alpha * cout$SE
	cout$lowCI_2 	= cout$EST * exp(-z_alpha * out$SE)
	cout$highCI_2 = cout$EST * exp(z_alpha * out$SE)
	cout
	
	LL 		= wrap_LL(iPARS)
	GRAD 	= wrap_GRAD(iPARS)
	HESS 	= wrap_HESS(iPARS)
	
	# est. alpha1 + lambda1
	nparams = 2 								
	# est. kappa
	nparams = nparams + upKAPPA	
	# est. theta
	nparams = nparams + ifelse(
		( tmp_copula == "Clayton" & cout$EST[4] == 0 )
		| ( PARAMS$copula == "Gumbel" & cout$EST[4] == 1 ),
		0,1)
	
	NN = nrow(DATA)
	BIC = 2 * LL - nparams * log(NN)
	
	
	names(GRAD) 		= out$PARS
	dimnames(HESS) 	= list(out$PARS,out$PARS)
	dimnames(covar) = list(out$PARS,out$PARS)
	
	return(list(GRID = gout$DAT,GPROF = gprof,
		GOPT = gopt,GOPT_PRE = gopt_pre,
		out = out,cout = cout,LL = LL,GRAD = GRAD,
		HESS = HESS,COVAR = covar,nparams = nparams,
		BIC = BIC))
	
}

run_analysis = function(DATA,theta,upKAPPA,
	gTHRES,copula,param_grid,vec_time,verb,PLOT,ncores = 1){
	
	if(FALSE){
		DATA 				= one_rep$DATA
		theta 			= 0
		upKAPPA 		= c(0,1)[1]
		gTHRES 			= 8e-2
		copula 			= c("Clayton","Gumbel")[1]
		param_grid 	= param_grid
		vec_time 		= vec_time
		
		verb = TRUE; ncores = 1
		
	}
	
	verb = get_verb(verb = verb)
	PLOT = get_PLOT(PLOT = PLOT)
	req_names = c("time","delta","dens_t2","surv_t2")
	if( !all(req_names %in% names(DATA)) ){
		miss_names = req_names[!(req_names %in% names(DATA))]
		stop(sprintf("DATA colnames missing: %s",
			paste(miss_names,collapse = ", ")))
	}
	
	REP = list(DATA = DATA,PARAMS = smart_df(copula = copula))
	opt_out = opt_replicate(REP = REP,param_grid = param_grid,
		theta = theta,upKAPPA = upKAPPA,gTHRES = gTHRES,
		verb = verb,PLOT = PLOT,ncores = ncores)
	if( is.null(opt_out$out) ){
		return(list(upKAPPA = upKAPPA,
			copula = copula,RES = opt_out,
			PRED = NULL))
	}
	names(opt_out)
	opt_out[c("out","cout")]
	
	# Estimate covariance
	iPARS = as.numeric(opt_out$out$EST)
	iPARS
	
	est_alpha1 	= exp(opt_out$out$EST[which(opt_out$out$PARS == "log_alpha1")])
	est_lambda1 = exp(opt_out$out$EST[which(opt_out$out$PARS == "log_lambda1")])
	# est_kappa1 	= 1/(1 + exp(-opt_out$out$EST[which(opt_out$out$PARS == "logit_kappa1")]))
	est_kappa1 	= exp(opt_out$out$EST[which(opt_out$out$PARS == "unc_kappa1")])
	est_theta		= exp(opt_out$out$EST[which(opt_out$out$PARS == "log_theta")])
	if( copula == "Gumbel" ) est_theta = est_theta + 1
	
	est_SE = function(tt){
		if(FALSE){
			ii = 2; tt = pred$time[ii]
		}
		
		if( tt == 0 ) return(0)
		XDL = tt / est_lambda1
		XDLa = XDL^est_alpha1
		neXDLa = exp(-XDLa)
		nab_surv = c(
			-est_alpha1 * est_kappa1 * 
				(1 - neXDLa)^(est_kappa1-1) *
				neXDLa * XDLa * log(XDL),
			est_alpha1 * est_kappa1 *
				(1 - neXDLa)^(est_kappa1-1) *
				neXDLa * XDLa,
			-est_kappa1 * (1 - neXDLa)^est_kappa1 *
				log(1 - neXDLa),
			0
			)
		nab_surv[3] = nab_surv[3] * upKAPPA
		var_err = as.numeric(t(nab_surv) %*% opt_out$COVAR %*% nab_surv)
		# if( var_err < 0 ) var_err = 0
		sqrt(var_err)
	}
	est_AA = function(tt){
		# tt = 0.278
		if(tt == 0) return(0)
		XDL = tt / est_lambda1
		XDLa = XDL^est_alpha1
		neXDLa = exp(-XDLa)
		tmp_surv = 1 - (1 - neXDLa)^est_kappa1
		if( tmp_surv == 1 ) return(0)
		part_L_surv = 1 / ( tmp_surv * log(tmp_surv) )
		# Partial survival wrt constrained parameters
		nab_surv = c(
			-est_alpha1 * est_kappa1 * 
				(1 - neXDLa)^(est_kappa1-1) *
				neXDLa * XDLa * log(XDL),
			est_alpha1 * est_kappa1 *
				(1 - neXDLa)^(est_kappa1-1) *
				neXDLa * XDLa,
			-est_kappa1 * (1 - neXDLa)^est_kappa1 *
				log(1 - neXDLa),
			0
			)
		nab_surv[3] = nab_surv[3] * upKAPPA
		nab_vec = nab_surv * part_L_surv
		var_err = as.numeric(t(nab_vec) %*% opt_out$COVAR %*% nab_vec)
		# if( var_err < 0 ) var_err = 0
		sqrt(var_err)
	}
	
	pred = smart_df(time = vec_time)
	pred$log_F1 = est_kappa1 * log(1 - exp(-(pred$time / est_lambda1)^est_alpha1))
	pred$F1 = exp(pred$log_F1)
	pred$log_surv = round(log(1 - pred$F1),4)
	pred$surv = exp(pred$log_surv)
	pred$SE = sapply(pred$time,function(tt) est_SE(tt = tt))
	pred$low_surv = pred$surv - 1.96 * pred$SE
	pred$high_surv = pred$surv + 1.96 * pred$SE
	pred$AA = 1.96 * sapply(pred$time,function(tt) est_AA(tt = tt))
	pred$log_low_surv = exp(pred$AA) * pred$log_surv
	pred$log_high_surv = exp(-pred$AA) * pred$log_surv
	pred$low_surv2 = exp(pred$log_low_surv)
	pred$high_surv2 = exp(pred$log_high_surv)
	# round(pred[1:10,],5)
	
	return(list(upKAPPA = upKAPPA,copula = copula,
		RES = opt_out,PRED = pred))
	
}

#' @title run_analyses
#' @description This function performs a full analysis of 
#'	an inputted dataframe. The user may specify one of two 
#'	copulas, a \code{theta} value, a parametric grid to 
#'	search over, and a vector of times for predicting survival.
#' @param DATA A data.frame containing column names \code{time},
#'	\code{delta}, \code{dens_t2}, and \code{surv_t2}.
#' @param THETAs A vector of theta values to explore and optimize over.
#' @param upKAPPA An integer value taking values 0 or 1. If set 
#'	to 1, the exponentiated Weibull distribution is assumed. Otherwise,
#'	the Weibull distribution is assumed and optimized over. If
#'	undefined, the optimization will search over both distributions.
#' @param gTHRES A numeric threshold on the L2 norm of 
#'	the gradient evaluated at the MLE.
#' @param param_grid Vector of values spanning possible 
#'	log(alpha1), log(lambda1), log(kappa1), unconstrained 
#'	theta parameters
#' @inheritParams sim_replicate
#' @param COPULAS If undefined, will optimize over all copulas. 
#'	Otherwise set to 'Independent', 'Clayton' or 'Gumbel'
#' @param vec_time Vector of times to calculate predicted survival
#'	on the same scale as times provided in the \code{DATA} data.frame.
#' @param ncores An integer for the number of threads
#' @param PLOT A logical variable, set to \code{TRUE} by default to
#'	show the two-dimensional heatmap of the profile likelihood if 
#'	\code{verb = TRUE}.
#' @export
run_analyses = function(DATA,THETAs = NULL,upKAPPA,
	gTHRES = 8e-2,COPULAS,param_grid,vec_time,ncores = 1,
	verb,PLOT){
	
	if(FALSE){
		DATA 				= one_rep$DATA
		DATA				= rd
		
		THETAs 			= NULL
		upKAPPA 		= NULL
		gTHRES 			= 8e-2
		COPULAS			= c("Clayton","Gumbel")[2]
		param_grid 	= seq(-3,5,0.5)
		vec_time 		= round(seq(0,max(c(100,max(DATA$time))),
										length.out = 100),2)
		ncores			= 1
		verb 				= TRUE
		
	}
	
	COPULAS 		= get_copula(copula = COPULAS,PO = FALSE)
	verb 				= get_verb(verb = verb)
	PLOT				= get_PLOT(PLOT = PLOT)
	vec_time		= get_vecTIME(TIME = DATA$time,vec_time = vec_time)
	upKAPPA			= get_upKAPPA(upKAPPA = upKAPPA)
	if( missing(param_grid) ) param_grid = get_parGRID()
	
	# Multiple copulas
	if( length(COPULAS) > 1 || length(upKAPPA) > 1 ){
		out = list()
		cnt = 1
		for(copula in COPULAS){
		for(up_kappa in upKAPPA){
			if( verb ) cat(sprintf("%s: Try copula = %s, dist = %s, ...\n",
				date(),copula,ifelse(up_kappa == 0,"Weibull","Exp-Weibull")))
			
			tmp_out = run_analyses(DATA = DATA,THETAs = THETAs,
				upKAPPA = up_kappa,gTHRES = gTHRES,COPULAS = copula,
				param_grid = param_grid,vec_time = vec_time,
				ncores = ncores,verb = verb,PLOT = PLOT)
			
			nn_tmp_out = length(tmp_out)
			
			if( nn_tmp_out == 0 ) next
			
			for(cnt2 in seq(nn_tmp_out)){
				tmp_out2 = tmp_out[[cnt2]]
				if( is.null(tmp_out2$RES$cout) ) next
				out[[cnt]] = tmp_out2
				cnt = cnt + 1
			}
			
		}}
		return(out)
	}
	
	# Once one copula is supplied, determine THETAs
	THETAs = get_THETAs(THETAs = THETAs,COPULA = COPULAS)
	
	# Store output
	out = list()
	
	if( COPULAS == "Independent" ){
		if( verb ) cat(sprintf("\n%s: Assume independence ...\n",date()))
		tmp_out = run_analysis(DATA = DATA,
			theta = 0,upKAPPA = upKAPPA,gTHRES = gTHRES,
			copula = COPULAS,param_grid = param_grid,
			vec_time = vec_time,verb = verb,PLOT = PLOT,
			ncores = ncores)
		if( !is.null(tmp_out$RES$cout) )
			out[[length(out) + 1]] = tmp_out
		
		return(out)
	}
	
	# Assuming dependence
	if( verb ) cat(sprintf("\n%s: Assume dependence: ...\n",date()))
	
	## Estimate theta
	if( is.null(THETAs) ){
		if( verb ) cat(sprintf("%s: theta = MLE ...\n",date()))
		tmp_out = run_analysis(DATA = DATA,
			theta = NA,upKAPPA = upKAPPA,gTHRES = gTHRES,
			copula = COPULAS,param_grid = param_grid,
			vec_time = vec_time,verb = verb,PLOT = PLOT,
			ncores = ncores)
		if( !is.null(tmp_out$RES$cout) )
			out[[length(out) + 1]] = tmp_out
		
		return(out)
	}
	
	## Fixed theta
	THETAs = sort(unique(THETAs))
	cnt = 1
	for(theta in THETAs){
		# theta = THETAs[1]; theta
		if( verb ) cat(sprintf("\n%s: theta = %s ...\n",date(),theta))
		tmp_name = sprintf("theta = %s",round(theta,4))
		tmp_out = run_analysis(DATA = DATA,
			theta = theta,upKAPPA = upKAPPA,
			gTHRES = gTHRES,copula = COPULAS,
			param_grid = param_grid,
			vec_time = vec_time,verb = verb,
			PLOT = PLOT,ncores = ncores)
		out[[cnt]] = tmp_out
		cnt = cnt + 1
	}
	
	return(out)
	
}


###