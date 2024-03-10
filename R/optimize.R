# ----------
# Optimization functions
# ----------
get_PROFILE = function(GRID,COPULA,PLOT = FALSE){
	if(FALSE){
		GRID = gout
		PLOT = TRUE
		
	}
	
	GRID = smart_df(GRID)
	names(GRID) = c("log_alpha1","log_lambda1",
		"unc_kappa1","log_theta","LL")
	GRID[1:5,]
	
	PARS = colnames(GRID)[1:4]; PARS
	PROF = list()
	for(PAR in PARS){
		# PAR = PARS[4]; PAR
		
		xx = sort(unique(GRID[[PAR]])); xx
		yy = sapply(xx,function(zz){
			# zz = xx[1]; zz
			if( is.infinite(zz) ){
				idx = which(GRID[[PAR]] == zz); idx
			} else {
				idx = which(abs(GRID[[PAR]] - zz) < 1e-5)
			}
			max(GRID$LL[idx])
		},USE.NAMES = FALSE)
		tmp_df = smart_df(xx,yy)
		names(tmp_df) = c(PAR,"LL")
		tmp_df
		PROF[[PAR]] = tmp_df
		rm(tmp_df)
		
	}
	
	if( PLOT ){
		par(mfrow = c(2,2),mar = c(4.5,4,2,2),oma = c(0,0,2,0))
		for(PAR in PARS){
			# PAR = PARS[2]; PAR
			tmp_df = PROF[[PAR]]
			# qnt = quantile(tmp_df$LL,0.15); qnt
			# tmp_df = tmp_df[which(tmp_df$LL >= qnt),]
			# tmp_df
			
			if( nrow(tmp_df) == 1 ){
				plot(0,0,type = "n",xlab = PAR,
					ylab = "Prof.LL")
				next
			}
			
			plot(tmp_df,xlab = PAR,
				ylab = "Prof.LL",
				type = "b",pch = 16)
		}
		
		DIST = ifelse(all(GRID$unc_kappa1 == 0),"Weibull","Exp-Weibull")
		main = sprintf("Copula = %s; DIST = %s",COPULA,DIST)
		mtext(main,outer = TRUE,cex = 1.3)
		par(mfrow = c(1,1),mar = c(5,4,4,2)+0.1,oma = rep(0,4))
	}
	
	# Get local opts
	LOCAL = list()
	for(PAR in PARS){
		# PAR = PARS[4]; PAR
		
		tmp_df = PROF[[PAR]]
		tmp_df$maxLL = NA
		nn = nrow(tmp_df)
		
		if( nn == 1 ){
			tmp_df$maxLL = tmp_df$LL
			tmp_df = smart_rmcols(tmp_df,"maxLL")
			LOCAL[[PAR]] = tmp_df
			next
		}
		
		for(ii in seq(nn)){
			# ii = 1
			jj = ii
			old_LL = tmp_df$LL[jj]
			while(TRUE){
				
				# Shift
				if( jj == 1 ){
					kk = jj + 1
				} else if( jj == nn ){
					kk = jj - 1
				} else {
					kk = jj + c(-1,1)
				}
				
				vec_LL = tmp_df$LL[kk]
				tmp_LL = max(vec_LL)
				idx = which(vec_LL == tmp_LL); idx
				kk = kk[idx]
				kk
				
				if( tmp_LL > old_LL ){
					old_LL = tmp_LL
					jj = kk
					next
				}
				
				break
			}
			
			# idx = which(tmp_df$LL == old_LL)
			tmp_df$maxLL[ii] = old_LL
			
		}
		
		tmp_df = tmp_df[which(tmp_df$LL == tmp_df$maxLL),]
		tmp_df = smart_rmcols(tmp_df,"maxLL")
		LOCAL[[PAR]] = tmp_df
		
	}
	
	LOCAL
	
	res = c()
	for(PAR in PARS){
		# PAR = PARS[4]; PAR
		LOCAL[[PAR]]
		nlocal = nrow(LOCAL[[PAR]]); nlocal
	for(jj in seq(nlocal)){
		# jj = 1
		value = LOCAL[[PAR]][[PAR]][jj]; value
		if( is.infinite(value) ){
			idx = which(GRID[[PAR]] == value
				& GRID$LL == LOCAL[[PAR]]$LL[jj])
		} else {
			idx = which(abs(GRID[[PAR]] - value) < 1e-5
				& GRID$LL == LOCAL[[PAR]]$LL[jj])
		}
		length(idx)
		# if( length(idx) != 1 ) stop("Check this")
		res = rbind(res,GRID[idx,])
	}}
	res = res[!duplicated(res),,drop = FALSE]
	rownames(res) = NULL
	res
	
	return(res)
}
opt_replicate = function(DATA,COPULA,param_grid,theta,
	upKAPPA,ncores = 1,gTHRES = 1e-1,max_iter = 2e2,verb,PLOT){

	if(FALSE){
		DATA 				= one_rep$DATA
		param_grid	= list(A = seq(0,0.3,0.02),
										L = seq(1,2,0.05),
										K = 0,
										T = seq(1,3,0.1))
		theta 			= ifelse(COPULA == "Independent",0,NA)
		upKAPPA 		= ifelse(DIST == "weibull",0,1)
		gTHRES 			= 1e-1
		verb 				= TRUE
		ncores 			= 1
		PLOT				= TRUE
		
	}
	
	verb 		= get_verb(verb = verb)
	PLOT		= get_PLOT(PLOT = PLOT)
	upTHETA = ifelse(is.na(theta),1,0)
	upPARS	= get_upPARS(upKAPPA = upKAPPA,
							THETA = theta)
	upTHETA	= upPARS[4]
	
	wrap_LL 	= function(PARS){
		# PARS = iPARS
		dMrs_cLL(XX = DATA$time,DELTA = DATA$delta,
			log_D2 = DATA$log_dens_t2,log_F2 = DATA$log_cdf_t2,
			PARS = PARS,copula = COPULA,
			verb = !TRUE)
	}
	wrap_GRAD = function(PARS){
		# PARS = iPARS
		out = dMrs_cGRAD(XX = DATA$time,DELTA = DATA$delta,
			log_D2 = DATA$log_dens_t2,log_F2 = DATA$log_cdf_t2,
			PARS = PARS,copula = COPULA,upPARS = upPARS)
		
		# out = grad(wrap_LL,PARS)
		c(out)
	}
	wrap_HESS = function(PARS){
		# PARS = iPARS
		dMrs_cHESS(XX = DATA$time,DELTA = DATA$delta,
			log_D2 = DATA$log_dens_t2,log_F2 = DATA$log_cdf_t2,
			PARS = PARS,copula = COPULA,upPARS = upPARS)
	}
	
	# Estimate parameters
	if( is(param_grid,"numeric") ){
		tmp_list = list()
		for(nm in c("A","L","K","T")){
			tmp_list[[nm]] = param_grid
		}
		param_grid = tmp_list; rm(tmp_list)
	}
	for(nm in names(param_grid)){
		tmp_vec = param_grid[[nm]]
		if( min(tmp_vec) > 0 || max(tmp_vec) < 0 )
			next
		param_grid[[nm]] = sort(unique(c(0,param_grid[[nm]])))
	}
	if( verb ) cat(sprintf("%s: Grid search for initial parameters ...\n",date()))
	if( upPARS[3] == 0 ){
		unc_KAPPA = 0
	} else {
		unc_KAPPA = param_grid$K
	}
	if( upTHETA == 0 ){
		if( COPULA == "Independent" ) 	log_THETA = log(theta)
		if( COPULA == "Clayton" ) 			log_THETA = log(theta)
		if( COPULA == "Gumbel" )				log_THETA = log(theta - 1)
	} else if( upTHETA == 1 ){
		# log_THETA = c(-Inf,param_grid)
		log_THETA = param_grid$T
		log_THETA = sort(unique(log_THETA))
	}
	
	if( verb ) cat(sprintf("%s: Start GRID\n",date()))
	gout = dMrs_GRID(XX = DATA$time,DELTA = DATA$delta,
		log_D2 = DATA$log_dens_t2,log_F2 = DATA$log_cdf_t2,
		log_ALPHA = param_grid$A,log_LAMBDA = param_grid$L,
		unc_KAPPA = unc_KAPPA,log_THETA = log_THETA,
		copula = COPULA,verb = verb,ncores = ncores)
	if( verb ) cat(sprintf("%s: End GRID\n",date()))
	# str(gout)
	colnames(gout) = c("log_alpha1","log_lambda1",
		"unc_kappa1","log_theta","LL")
	gout = smart_df(gout)
	dim(gout); head(gout)
	gout = gout[which(gout[,"LL"] != -999),,drop = FALSE]
	dim(gout)
	gout[1:5,]
	chk_NA = any(is.na(gout[,"LL"])); chk_NA
	if( chk_NA ) stop("NAs in LL grid, debug this")
	
	if( nrow(gout) == 0 ){
		if( verb ) cat(sprintf("%s: No valid grid points! ...\n",date()))
		return(list(GRID = gout,GPROF = NULL,GOPT = NULL,
			GOPT_PRE = NULL,out = NULL,cout = NULL,
			LL = NULL,GRAD = NULL,HESS = NULL,COVAR = NULL))
	}
	
	gprof = get_PROFILE(GRID = gout,
		COPULA = COPULA,PLOT = PLOT)
	gprof
	
	gopt = gprof
	gopt$fin_logA = NA; gopt$fin_logL = NA; 
	gopt$fin_logK = NA; gopt$fin_logT = NA;
	gopt$fin_LL = NA; gopt$nGRAD = NA; 
	# gopt
	
	nn = nrow(gopt)
	if( verb ) cat(sprintf("%s: NR optimization w/ %s profile point(s) ...\n",date(),nn))
	DIST = ifelse(upKAPPA == 0,"weibull","expweibull")
	for(ii in seq(nn)){
		# ii = 3
		if( verb ){
			cat(sprintf("#---# Copula=%s, Dist=%s, Prof. point=%s out of %s\n",
				COPULA,DIST,ii,nn))
		}
		if( !is.na(gopt$fin_LL[ii]) ) next
		
		# Init pars
		iPARS = as.numeric(gopt[ii,1:4]); iPARS
		
		# Run Newton Raphson
		dMrs_NR(XX = DATA$time,DELTA = DATA$delta,
			log_D2 = DATA$log_dens_t2,log_F2 = DATA$log_cdf_t2,
			PARS = iPARS,copula = COPULA,upPARS = upPARS,
			max_iter = max_iter,eps = 5e-2,verb = verb)
		
		tmp_LL = wrap_LL(PARS = iPARS)
		tmp_GR = wrap_GRAD(PARS = iPARS)
		if(FALSE){
			tmp_LL
			tmp_GR
			tmp_HE = wrap_HESS(PARS = iPARS)
			tmp_HE
			
			covar = matrix(0,4,4)
			nz = which(diag(tmp_HE) != 0); nz
			smart_solve(tmp_HE[nz,nz])
			covar[nz,nz] = solve(-tmp_HE[nz,nz])
			covar
			
		}
		
		gopt$nGRAD[ii] = Rcpp_norm(tmp_GR)
		gopt$fin_logA[ii] = iPARS[1]
		gopt$fin_logL[ii] = iPARS[2]
		gopt$fin_logK[ii] = iPARS[3]
		gopt$fin_logT[ii] = iPARS[4]
		gopt$fin_LL[ii] = round(tmp_LL,3)
	}
	gopt = gopt[order(-gopt$fin_LL),]
	if( verb ){
		# cat("\n")
		print(head(gopt))
	}
	gopt = gopt[gopt$fin_LL != -999 & gopt$nGRAD > 0,,drop = FALSE]
	rownames(gopt) = NULL
	gopt_pre = gopt
	gopt_pre
	
	# Remove non-local optimum solutions
	gopt = gopt[which(gopt$nGRAD < gTHRES),]
	if( nrow(gopt) == 0 ){
		if( verb ) cat(sprintf("%s: No converged solutions! ...\n",date()))
		return(list(GRID = gout,GPROF = gprof,GOPT = gopt,
			GOPT_PRE = gopt_pre,
			out = NULL,cout = NULL,LL = NULL,GRAD = NULL,
			HESS = NULL,COVAR = NULL))
	}
	
	# Order remaining, remove duplicates
	gopt = gopt[order(-gopt$fin_LL),]
	gopt = gopt[!duplicated(round(gopt[,
		sprintf("fin_%s",c("LL","logA",
		"logL","logK","logT"))],2)),]
	rownames(gopt) = NULL
	gopt
	
	# Estimate covariance
	if( verb ) cat(sprintf("%s: Get covariance ...\n",date()))
	gopt$neg_var = NA
	for(jj in seq(nrow(gopt))){
		# jj = 1
		iPARS = as.numeric(gopt[jj,
			paste0("fin_log",c("A","L","K","T"))]); iPARS
		hess = wrap_HESS(PARS = iPARS); hess
		nz = which(diag(hess) != 0)
		if( length(nz) == 0 ) next
		# print(hess)
		if( rcond(hess[nz,nz,drop = FALSE]) == 0 ) next
		covar = matrix(0,4,4)
		tmp_mat = smart_solve(MAT = -hess[nz,nz,drop = FALSE])
		if( is.null(tmp_mat) ) next
		covar[nz,nz] = tmp_mat
		# covar
		gopt$neg_var[jj] = any(diag(covar) < 0)
	}
	gopt = gopt[!is.na(gopt$neg_var),,drop = FALSE]
	gopt = gopt[which(gopt$neg_var == FALSE),]
	if( nrow(gopt) == 0 ){
		if( verb ) cat(sprintf("%s: Negative variance(s)! ...\n",date()))
		return(list(GRID = gout,GPROF = gprof,GOPT = NULL,
			GOPT_PRE = gopt_pre,
			out = NULL,cout = NULL,LL = NULL,GRAD = NULL,
			HESS = NULL,COVAR = NULL))
	}
	
	if( verb ) print(gopt)
	
	# Grab best solution
	gopt = gopt[which.max(gopt$fin_LL),]
	
	iPARS = as.numeric(gopt[1,paste0("fin_log",c("A","L","K","T"))])
	# iPARS = as.numeric(gopt[2,paste0("fin_log",c("A","L","K","T"))])
	iPARS
	if(FALSE){
		iPARS = true_PARS
	}
	hess = wrap_HESS(PARS = iPARS); hess
	nz = which(diag(hess) != 0)
	if( rcond(hess[nz,nz,drop = FALSE]) == 0 ){
		print(hess)
		cat("Error with hessian 1: rcond = 0\n")
		return(list(GRID = gout,GPROF = gprof,GOPT = NULL,
			GOPT_PRE = gopt_pre,
			out = NULL,cout = NULL,LL = NULL,GRAD = NULL,
			HESS = NULL,COVAR = NULL))
	}
	covar = matrix(0,4,4)
	covar[nz,nz] = smart_solve(MAT = -hess[nz,nz,drop = FALSE])
	covar
	
	out = smart_df(PARS = c("log_alpha1","log_lambda1",
		"unc_kappa1","log_theta"),
		EST = c(iPARS),SE = c(sqrt(diag(covar))))
	z_alpha = qnorm(0.975)
	out$lowCI 	= out$EST - z_alpha * out$SE
	out$highCI 	= out$EST + z_alpha * out$SE
	out
	
	EST = exp(out$EST)
	if( COPULA == "Gumbel" ) EST[4] = EST[4] + 1
	
	cout = smart_df(PARS = c("alpha1","lambda1","kappa1","theta"),
		EST = EST,SE = NA)
	
	# Delta method: Calculate nabla g()
	nab_par	= rep(NA,4)
	nab_par[1:3] = cout$EST[1:3] # alpha, lambda, kappa
	nab_par[4] = ifelse(COPULA == "Clayton",
		cout$EST[4],exp(out$EST[4]))
	nab_par = diag(nab_par)
	
	cout$SE 		= sqrt(diag(nab_par %*% covar %*% nab_par))
	cout$lowCI 		= cout$EST - z_alpha * cout$SE
	cout$highCI 	= cout$EST + z_alpha * cout$SE
	cout$lowCI_2 	= cout$EST * exp(-z_alpha * out$SE)
	cout$highCI_2 = cout$EST * exp(z_alpha * out$SE)
	
	if( COPULA == "Gumbel" ){
		cout$lowCI_2[4] 	= exp(out$EST[4]) * exp(-z_alpha * out$SE[4]) + 1
		cout$highCI_2[4] 	= exp(out$EST[4]) * exp(z_alpha * out$SE[4]) + 1
	}
	
	if( verb ) print(cout)
	
	LL 		= wrap_LL(iPARS)
	GRAD 	= wrap_GRAD(iPARS)
	HESS 	= wrap_HESS(iPARS)
	
	LL = round(LL,3)
	
	nparams = 2 
		# alpha,lambda
	nparams = nparams + upKAPPA	
		# kappa(fixed or estimated)
	nparams = nparams + upTHETA
		# theta(fixed at the boundary or estimated)
	nparams
	
	NN = nrow(DATA)
	BIC = 2 * LL - nparams * log(NN)
	
	names(GRAD) 		= out$PARS
	dimnames(HESS) 	= list(out$PARS,out$PARS)
	dimnames(covar) = list(out$PARS,out$PARS)
	
	return(list(GRID = gout,GPROF = gprof,
		GOPT = gopt,GOPT_PRE = gopt_pre,
		out = out,cout = cout,LL = LL,GRAD = GRAD,
		HESS = HESS,COVAR = covar,nparams = nparams,
		BIC = BIC))
	
}
run_analysis = function(DATA,theta,upKAPPA,gTHRES,
	copula,param_grid,vec_time,max_iter = 2e2,verb,PLOT,
	ncores = 1){
	
	if(FALSE){
		DATA 			= one_rep$DATA
		theta 			= NA
		upKAPPA 		= ifelse(DIST == "weibull",0,1)
		gTHRES 			= 8e-2
		copula 			= COPULA
		param_grid 		= param_grid
		vec_time 		= vec_time
		
		verb = TRUE; ncores = 1
		
	}
	
	verb = get_verb(verb = verb)
	PLOT = get_PLOT(PLOT = PLOT)
	req_names = c("time","delta","log_dens_t2","log_cdf_t2")
	if( !all(req_names %in% names(DATA)) ){
		miss_names = req_names[!(req_names %in% names(DATA))]
		stop(sprintf("DATA colnames missing: %s",
			paste(miss_names,collapse = ", ")))
	}
	
	# REP = list(DATA = DATA,PARAMS = smart_df(copula = copula))
	opt_out = opt_replicate(DATA = DATA,COPULA = copula,
		param_grid = param_grid,theta = theta,upKAPPA = upKAPPA,
		gTHRES = gTHRES,max_iter = max_iter,verb = verb,
		PLOT = PLOT,ncores = ncores)
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
			ii = 40; tt = pred$time[ii]
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
#' @param DATA A data.frame containing column names 
#'	\code{time} (in years),
#'	\code{delta} (event indicator), 
#'	\code{log_dens_t2} (log-transformed population-based density), and 
#'	\code{log_cdf_t2} (log-transformed population-based cumulative density).
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
#' @param COPULAS If undefined, will optimize over all copulas. 
#'	Otherwise set to 'Independent', 'Clayton' or 'Gumbel'
#' @param vec_time Vector of times in years to calculate predicted survival.
#' @param max_iter Maximum Newton Raphson and Gradient Descent iterations
#'	to set.
#' @param ncores A positive integer for the number of threads to evaluate 
#'	log-likelihoods across the parameter grid.
#' @param PLOT A logical variable, set to \code{TRUE} by default to
#'	show the two-dimensional heatmap of the profile likelihood if 
#'	\code{verb = TRUE}.
#' @param verb Boolean value to display verbose information or not
#' @export
run_analyses = function(DATA,THETAs = NULL,upKAPPA,
	gTHRES = 1e-1,COPULAS,param_grid,vec_time,ncores = 1,
	max_iter = 2e2,verb,PLOT){
	
	stopifnot(is(ncores,"numeric"))
	stopifnot(ncores > 0 && round(ncores) == ncores)
	
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
			if( verb ) cat(sprintf("\n####\n%s: Try copula = %s, dist = %s, ...\n",
				date(),copula,ifelse(up_kappa == 0,"Weibull","Exp-Weibull")))
			
			tmp_out = run_analyses(DATA = DATA,THETAs = THETAs,
				upKAPPA = up_kappa,gTHRES = gTHRES,COPULAS = copula,
				param_grid = param_grid,vec_time = vec_time,
				ncores = ncores,max_iter = max_iter,verb = verb,
				PLOT = PLOT)
			
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
			vec_time = vec_time,max_iter = max_iter,
			verb = verb,PLOT = PLOT,ncores = ncores)
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
			vec_time = vec_time,
			max_iter = max_iter,verb = verb,
			PLOT = PLOT,ncores = ncores)
		out[[cnt]] = tmp_out
		cnt = cnt + 1
	}
	
	return(out)
	
}

###