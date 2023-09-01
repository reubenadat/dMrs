# ----------
# Optimization functions
# ----------

ref_LL = function(DATA,PARS,COPULA){
	# PARS = (ALPHA1,LAMBDA1,KAPPA1,THETA)
	if(FALSE){
		PARS = iPARS
		COPULA = PARAMS$copula; COPULA
		
	}
	
	LL = 0
	nn = nrow(DATA)
	error_num = -999
	
	ALPHA1 	= exp(PARS[1])
	LAMBDA1 = exp(PARS[2])
	KAPPA1 	= exp(PARS[3])
	THETA 	= exp(PARS[4])
	if( COPULA == "Gumbel" ) THETA = THETA + 1
	
	KAL = KAPPA1 * ALPHA1 / LAMBDA1
	
	for(ii in seq(nn)){
		# ii = 1
		TT_LL = DATA$time[ii] / LAMBDA1
		TLA = (TT_LL)^ALPHA1
		E_mTLA = exp(-TLA)
		CDF = (1 - E_mTLA)^KAPPA1
		PDF = KAL * (TT_LL)^(ALPHA1 - 1) * 
			(1 - E_mTLA)^(KAPPA1 - 1) * E_mTLA
		f1 = PDF
		f2 = DATA$dens_t2[ii]
		F1 = CDF
		if(F1 <= 0){
			print(sprintf("ii = %s: error1",ii))
			return(error_num)
		}
		F2 = 1 - DATA$surv_t2[ii]
		
		if( COPULA == "Independent" ){
			f_T1_T2 = f1 * F2 + f2 * F1
			F_T1_T2 = F1 * F2
		} else if( COPULA == "Clayton" ){
			F_T1_T2 = calc_copula(F1 = F1,F2 = F2,
				copula = COPULA,THETA = THETA)
			f_T1_T2 = ( (F1)^(-THETA) + (F2)^(-THETA) - 1 )^(-1/THETA - 1) *
				(f1 / F1^(THETA + 1) + f2 / F2^(THETA + 1))
			f_T1_T2 = calc_copula_dens(D1 = f1,D2 = f2,
				F1 = F1,F2 = F2,copula = COPULA,THETA = THETA,
				F_T1_T2 = F_T1_T2)
		} else if( COPULA == "Gumbel" ){
			nlog_u1 = -log(F1)
			nlog_u2 = -log(F2)
			F_T1_T2 = exp(-(nlog_u1^THETA + nlog_u2^THETA)^(1/THETA))
			f_T1_T2 = F_T1_T2 *
				( nlog_u1^THETA + nlog_u2^THETA )^(1 / THETA - 1) *
				( nlog_u1^(THETA-1) * f1 / F1 + nlog_u2^(THETA-1) * f2 / F2 )
		}
		
		if( DATA$delta[ii] == 1 ){
			tmp_num = f1 + f2 - f_T1_T2
			if( is.na(tmp_num) ){
				cat(sprintf("ii = %s, f1 = %s, f2 = %s, F_T1_T2 = %s\n",ii,f1,f2,F_T1_T2))
				break
			}
			
			if( tmp_num <= 0 ){
				print(sprintf("ii = %s: error2, f1=%s, f2=%s, f_T1_T2=%s, F1=%s, F2=%s, T1=%s, T2=%s",
					ii,f1,f2,f_T1_T2,F1,F2,( (F1)^(-THETA) + (F2)^(-THETA) - 1 )^(-1/THETA - 1),
					(f1 / F1^(THETA + 1) + f2 / F2^(THETA + 1))))
				return(error_num)
			}
			tmp_LL = log(tmp_num)
		} else {
			tmp_num = (1-F1) + (1-F2) - 1 + F_T1_T2
			if( tmp_num <= 0 ){
				print(sprintf("ii = %s: error3",ii))
				return(error_num)
			}
			tmp_LL = log(tmp_num)
		}
		
		LL = LL + tmp_LL
		
	}
	
	return(LL)
	
}
ref_LL_cpp = function(DATA,PARS,COPULA){
	if(FALSE){
		DATA 		= DATA
		PARS 		= iPARS
		COPULA 	= tmp_copula
		
	}
	
	# PARS = (ALPHA1,LAMBDA1,KAPPA1,THETA)
	
	LL = 0
	nn = nrow(DATA)
	error_num = -999
	
	ALPHA1 	= exp(PARS[1])
	LAMBDA1 = exp(PARS[2])
	KAPPA1 	= exp(PARS[3])
	THETA 	= exp(PARS[4])
	if( COPULA == "Gumbel" ) THETA = THETA + 1
	
	KAL = KAPPA1 * ALPHA1 / LAMBDA1
	PART1 = 0; PART2 = 0;
	
	for(ii in seq(nn)){
		# ii = 3623
		XDL = DATA$time[ii] / LAMBDA1
		enXDLa = exp(-(XDL)^ALPHA1)
		
		F1 = 1 - enXDLa
		if( KAPPA1 != 1.0 ) F1 = F1^KAPPA1
		S1 = 1.0 - F1
		D1 = KAL * (XDL)^(ALPHA1 - 1) * enXDLa
		if( KAPPA1 != 1.0 ) D1 = D1 * F1 / ( 1.0 - enXDLa )
		
		D2 = DATA$dens_t2[ii]
		S2 = DATA$surv_t2[ii]
		F2 = 1 - S2
		# if( F1 == 0.0 || F2 == 0.0 ) return(error_num)
		
		tmp_vec = calc_copula_CDF_PDF(D1 = D1,D2 = D2,
			F1 = F1,F2 = F2,copula = COPULA,THETA = THETA)
		tmp_vec = as.numeric(tmp_vec)
		if( any(is.na(tmp_vec)) ){
			tmp_vec
			print(sprintf("Check on ii = %s",ii))
			print(sprintf("D1=%s, D2=%s, F1=%s, F2=%s",D1,D2,F1,F2))
			return(error_num)
			
			F_T1_T2 = calc_copula(F1 = F1,F2 = F2,
				copula = COPULA,THETA = THETA); F_T1_T2
			
			calc_copula_dens(D1 = D1,D2 = D2,
				F1 = F1,F2 = F2,copula = COPULA,
				THETA = THETA,F_T1_T2 = F_T1_T2)
			
			
			# Check on calculation
			f_T1_T2 = F1^(-THETA) + F2^(-THETA) - 1.0
				f_T1_T2
				log(f_T1_T2)
			f_T1_T2 = f_T1_T2^(-1.0 / THETA - 1.0)
				f_T1_T2
				exp((-1/THETA - 1) * log(f_T1_T2))
			f_T1_T2 = f_T1_T2 * (D1 / F1^(THETA + 1.0) 
				+ D2 / F2^(THETA + 1.0)); f_T1_T2
			
			log_P = rep(NA,2)
			log_P[1] = log(D1) - (THETA + 1) * log(F1)
			log_P[2] = log(D2) - (THETA + 1) * log(F2)
			log_P
			
			log_TERM_1 = (-1/THETA-1) * log(F1^(-THETA) + F2^(-THETA) - 1); log_TERM_1
			log_TERM_2 = Rcpp_logSumExp(log_P); log_TERM_2
			exp(log_TERM_1 + log_TERM_2)
			
			
		}
		F1_F2 = tmp_vec[1]
		D1_D2 = tmp_vec[2]
		
		if( DATA$delta[ii] == 1 ){
			tmp_num = D1 + D2 - D1_D2
			
			if( tmp_num <= 0 ){
				print(sprintf("ii = %s: D1=%s, D2=%s, D1_D2=%s, F1=%s, F2=%s, T1=%s, T2=%s",
					ii,D1,D2,D1_D2,F1,F2,PART1,PART2))
				# print(sprintf("ii = %s",ii))
				return(error_num)
			}
			tmp_LL = log(tmp_num)
		} else {
			tmp_num = S1 + S2 - 1 + F1_F2
			if( tmp_num <= 0 ){
				print(sprintf("ii = %s",ii))
				return(error_num)
			}
			tmp_LL = log(tmp_num)
		}
		
		LL = LL + tmp_LL
		
	}
	
	LL
	
	return(LL)
}

#' @title calc_CDFs
#' @description Investigate the nature of the cumulative distribution 
#'	and copula functions for simulated or real data sets
#' @param DATA Dataframe containing observed time, second event time 
#'	survival probability
#' @param PARS Vector of parameter values
#' @param COPULA The copula considered by the user.
#' @export
calc_CDFs = function(DATA,PARS,COPULA){
	if(FALSE){
		DATA = one_rep$DATA
		# PARS = iPARS
		PARS = as.numeric(gout[1,1:4]); PARS
		COPULA
		
	}
	
	ALPHA1 	= exp(PARS[1])
	LAMBDA1 = exp(PARS[2])
	KAPPA1 	= exp(PARS[3])
	THETA 	= exp(PARS[4]) + 
							ifelse(COPULA %in% c("Independent","Clayton"),0,1)
	KAL = KAPPA1 * ALPHA1 / LAMBDA1
	
	TT_LL = DATA$time / LAMBDA1
	TLA = (TT_LL)^ALPHA1
	E_mTLA = exp(-TLA)
	CDF_1 = (1 - E_mTLA)^KAPPA1
	CDF_2 = 1 - DATA$surv_t2
	
	par(mfrow = c(3,2),mar = c(4.4,4.4,1,0.2))
	hist(CDF_1,breaks = 40,xlim = c(0,1))
	
	hist(CDF_2,breaks = 40,xlim = c(0,1))
	
	plot(DATA$time,CDF_1,col = "red",
		xlab = "Obs Time",ylab = "CDF")
	points(DATA$time,CDF_2,col = "blue")
	
	smoothScatter(CDF_1,CDF_2,
		xlim = c(0,1),ylim = c(0,1))
	
	F_T1_T2 = apply(smart_df(CDF_1,CDF_2),1,function(xx){
		calc_copula(F1 = xx[1],F2 = xx[2],
			copula = COPULA,THETA = THETA)
	})
	hist(F_T1_T2,breaks = 40,
		xlab = "Joint CDF Copula",
		xlim = c(0,1))
	
	if( all(c("T1","T2") %in% names(DATA)) ){
		smoothScatter(log(1 + DATA[,c("T1","T2")]),
			xlab = "log(1 + Time1)",
			ylab = "log(1 + Time2)",
			main = sprintf("Copula=%s, Theta=%s",COPULA,THETA))
	}
	
	par(mfrow = c(1,1),mar = c(5,4,4,2)+0.1)
	
	# Calc densities
	D1 = KAL * TLA^(ALPHA1 - 1) * E_mTLA
	if( KAPPA1 != 1.0 ) D1 = D1 * CDF_1 / (1 - E_mTLA)
	D2 = DATA$dens_t2
	
	out = smart_df(CDF_1 = CDF_1,
		CDF_2 = CDF_2,F_T1_T2 = F_T1_T2,
		D1 = D1,D2 = D2)
	
	out$D1_D2 = apply(out[,c("D1","D2","CDF_1","CDF_2","F_T1_T2")],
		1,function(xx){
		calc_copula_dens(D1 = xx[1],D2 = xx[2],
			F1 = xx[3],F2 = xx[4],copula = COPULA,
			THETA = THETA,F_T1_T2 = xx[5])
	})
	
	dim(out); head(out)
	
	return(out)
	
}

wrap_NR = function(DATA,PARS,COPULA,upPARS,verb = TRUE){
	
	dMrs_NR(XX = DATA$time,
		DELTA = DATA$delta,
		D2 = DATA$dens_t2,
		S2 = DATA$surv_t2,
		PARS = PARS,
		copula = COPULA,
		upPARS = upPARS,
		max_iter = 2e2,
		eps = 1e-7,
		verb = verb)
	
	out_PARS = PARS
	return(out_PARS)
	
}

opt_replicate = function(DATA,COPULA,param_grid,theta,
	upKAPPA,ncores = 1,gTHRES = 8e-2,verb,PLOT){

	if(FALSE){
		DATA 				= one_rep$DATA
		param_grid	= list(A = seq(0,1,0.05),
										L = seq(1,2,0.05),
										K = 0,
										T = seq(1,3,0.1))
		theta 			= ifelse(COPULA_2 == "Independent",0,NA)
		upKAPPA 		= ifelse(DIST_2 == "weibull",0,1)
		gTHRES 			= 8e-2
		verb 				= TRUE
		ncores 			= 1
		PLOT				= !TRUE
		
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
			D2 = DATA$dens_t2,S2 = DATA$surv_t2,
			PARS = PARS,copula = COPULA,
			verb = !TRUE)
	}
	wrap_GRAD = function(PARS){
		# PARS = iPARS
		out = dMrs_cGRAD(XX = DATA$time,DELTA = DATA$delta,
			D2 = DATA$dens_t2,S2 = DATA$surv_t2,PARS = PARS,
			copula = COPULA,upPARS = upPARS)
		
		# out = grad(wrap_LL,PARS)
		c(out)
	}
	wrap_HESS = function(PARS){
		# PARS = iPARS
		dMrs_cHESS(XX = DATA$time,DELTA = DATA$delta,
			D2 = DATA$dens_t2,S2 = DATA$surv_t2,PARS = PARS,
			copula = COPULA,upPARS = upPARS)
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
	
	gout = dMrs_GRID(XX = DATA$time,DELTA = DATA$delta,
		D2 = DATA$dens_t2,S2 = DATA$surv_t2,
		log_ALPHA = param_grid$A,log_LAMBDA = param_grid$L,
		unc_KAPPA = unc_KAPPA,log_THETA = log_THETA,
		copula = COPULA,verb = verb,ncores = ncores)
	# str(gout)
	colnames(gout) = c("log_alpha1","log_lambda1",
		"unc_kappa1","log_theta","LL")
	gout = smart_df(gout)
	dim(gout); head(gout)
	gout = gout[which(gout[,"LL"] != -999),,drop = FALSE]
	dim(gout)
	gout[1:5,]
	chk_NA = any(is.na(gout[,"LL"])); chk_NA
	if( chk_NA ){
		stop("NAs in LL grid, debug this")
		smart_table(!is.na(gout[,"LL"]))
		gout[is.na(gout[,"LL"]),]
		
		iPARS = c(-0.1,2.6,2.9,2.9)
		
		# Rcpp function
		dMrs_cLL(XX = DATA$time,DELTA = DATA$delta,
			D2 = DATA$dens_t2,S2 = DATA$surv_t2,
			PARS = iPARS,copula = COPULA,
			verb = TRUE)
		
		# R function
		ref_LL_cpp(DATA = DATA,PARS = iPARS,COPULA = COPULA)
		
		ref_LL(DATA = DATA,PARS = iPARS,COPULA = COPULA)
		
	}
	
	if( nrow(gout) == 0 ){
		if( verb ) cat(sprintf("%s: No valid grid points! ...\n",date()))
		return(list(GRID = gout,GPROF = NULL,GOPT = NULL,
			GOPT_PRE = NULL,out = NULL,cout = NULL,
			LL = NULL,GRAD = NULL,HESS = NULL,COVAR = NULL))
	}
	
	gprof = get_PROFILE(GRID = gout,PLOT = PLOT)
	# gprof = get_PROFILE(GRID = gout[which(gout$log_lambda1 > 2 & gout$log_theta > 1),],PLOT = verb)
	gprof
	
	gopt = gprof
	gopt$fin_logA = NA; gopt$fin_logL = NA; 
	gopt$fin_logK = NA; gopt$fin_logT = NA;
	gopt$fin_LL = NA; gopt$nGRAD = NA; 
	# gopt
	
	nn = nrow(gopt)
	if( verb ) cat(sprintf("%s: NR optimization w/ %s profile point(s) ...\n",date(),nn))
	for(ii in seq(nn)){
		# ii = 3
		# if( verb ) cat(".")
		if( verb ) smart_progress(ii = ii,nn = nn,
			iter = 1,iter2 = 5)
		if( !is.na(gopt$fin_LL[ii]) ) next
		
		iPARS = as.numeric(gopt[ii,1:4]); iPARS
		dMrs_NR(XX = DATA$time,DELTA = DATA$delta,
			D2 = DATA$dens_t2,S2 = DATA$surv_t2,
			PARS = iPARS,copula = COPULA,upPARS = upPARS,
			max_iter = 2e2,eps = 1e-6,verb = verb)
		
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
	
	cout$SE 			= sqrt(diag(nab_par %*% covar %*% nab_par))
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

run_analysis = function(DATA,theta,upKAPPA,
	gTHRES,copula,param_grid,vec_time,verb,PLOT,ncores = 1){
	
	if(FALSE){
		DATA 				= DATA
		theta 			= NA
		upKAPPA 		= upKAPPA
		gTHRES 			= 8e-2
		copula 			= c(COPULA,"Clayton","Gumbel")[1]
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
	
	# REP = list(DATA = DATA,PARAMS = smart_df(copula = copula))
	opt_out = opt_replicate(DATA = DATA,COPULA = copula,
		param_grid = param_grid,theta = theta,
		upKAPPA = upKAPPA,gTHRES = gTHRES,
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
		THETAs 			= NULL
		upKAPPA 		= ifelse(DIST == "weibull",0,1)
		COPULAS 		= COPULA
		param_grid = seq(-1,3,0.3)
		
		######
		vec_time 		= round(seq(0,max(c(100,max(DATA$time))),
										length.out = 100),2)
		gTHRES 			= 8e-2
		ncores			= 1
		verb 				= TRUE; PLOT = verb
		
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
			if( verb ) cat(sprintf("\n####\n%s: Try copula = %s, dist = %s, ...\n",
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

chk_profile_LL = function(GRID){
	
	if(FALSE){
		GRID = OPT[[3]]$RES$GRID
		
	}
	
	GRID = smart_df(GRID)
	# GRID = GRID[which(GRID$log_theta <= 1),]
	pars = colnames(GRID)[1:4]; pars
	
	# Get profile likelihood per param
	par(mfrow = c(2,2),mar = c(4.5,4,2,2))
	out = sapply(pars,function(xx){
		# xx = pars[1]; xx
		x1 = sort(unique(GRID[,xx]))
		y1 = sapply(x1,function(zz){
			# zz = x1[3]; zz
			max(GRID$LL[which(GRID[[xx]] == zz)],na.rm = TRUE)
		})
		
		dat = smart_df(x1 = x1,y1 = y1)
		plot(dat,xlab = xx,ylab = "Prof.LL",
			type = "b",pch = 16)
		abline(v = x1[which.max(y1)],lty = 2,lwd = 2,col = "red")
		out = c(VAL = x1[which.max(y1)],LL = max(dat$y1))
		out
	})
	out = smart_df(t(out))
	par(mfrow = c(1,1),mar = c(5,4,4,2)+0.1)
	
	return(out)
	
	# EXP: Get bivariate profile likelihood
	par_pair = pars[c(1,2)]; par_pair
	head(GRID)
	
	uPARs_1 = sort(unique(GRID[[par_pair[1]]]))
	uPARs_2 = sort(unique(GRID[[par_pair[2]]]))
	DAT = c()
	
	for(uPAR_1 in uPARs_1){
	for(uPAR_2 in uPARs_2){
		idx = which(GRID[[par_pair[1]]] == uPAR_1
			& GRID[[par_pair[2]]] == uPAR_2)
		length(idx)
		DAT = rbind(DAT,smart_df(uPAR_1 = uPAR_1,
			uPAR_2 = uPAR_2,LL = max(GRID$LL[idx])))
	}}
	
	names(DAT)[1:2] = c("xx","yy")
	dim(DAT); head(DAT)
	
	my_stepsize = diff(sort(unique(DAT[,1])))[1]; my_stepsize
	xrange = range(DAT[,1]); xrange
	yrange = range(DAT[,2]); yrange
	DAT = DAT[which(DAT[,1] >= xrange[1]
		& DAT[,1] <= xrange[2]
		& DAT[,2] >= yrange[1]
		& DAT[,2] <= yrange[2]),]
	# res$LL[is.na(res$LL)] = min(res$LL,na.rm = TRUE)
	nBREAKs = 10
	
	xbreaks = round(seq(xrange[1],xrange[2],length.out = nBREAKs),1)
	ybreaks = round(seq(yrange[1],yrange[2],length.out = nBREAKs),1)
	
	xx = yy = LL = NULL
	g1 = ggplot(data = DAT,mapping = aes(x = xx,
		y = yy,fill = LL)) +
		geom_tile(width = my_stepsize,height = my_stepsize) + 
		coord_cartesian(expand = 0,xlim = xrange,ylim = yrange) +
		scale_x_continuous(breaks = xbreaks) +
		scale_y_continuous(breaks = ybreaks) +
		labs(fill = "Profile LL") + xlab(par_pair[1]) +
		ylab(par_pair[2]) +
		scale_fill_viridis(discrete = FALSE)
	g1
	

}

###