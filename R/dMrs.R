
# ----------
# R functions: simulation, likelihood, gradient, optimization
# ----------
#' @title pexpweibull
#' @description CDF of exp-weibull
#' @param q quantile
#' @param lambda lambda parameter
#' @param alpha alpha parameter
#' @param kappa kappa parameter
#' @param log.p Boolean
#' @export
pexpweibull = function(q,lambda,alpha,kappa,log.p = FALSE){
	if( (!is.numeric(q)) || (!is.numeric(lambda)) 
		|| (!is.numeric(kappa)) || (!is.numeric(alpha)) )
		stop("non-numeric argument to mathematical function")
	
	if( (length(lambda)!= 1) || (length(kappa)!=1) || (length(alpha)!=1) )
		stop("Non-q parameters must be atomic")
	
	#if ((min(q) <= 0) || (lambda <= 0) || (kappa <= 0) || (alpha <= 0) ){ 
	# stop("Invalid arguments.  qp,lp,kp ap must be > 0 ")}
	log.cdf = kappa * pweibull(q,scale = lambda,shape = alpha,log.p = TRUE)
	
	# mpfr(ifelse(log.p,return(log.cdf),return(exp(log.cdf))),128L)
	ifelse(log.p,return(log.cdf),return(exp(log.cdf)))
	
}  

#' @title dexpweibull
#' @description PDF of exp-weibull
#' @inheritParams pexpweibull
#' @param x Quantile
#' @param log Boolean
#' @export
dexpweibull = function(x,lambda,alpha,kappa,log = FALSE){
	# t = 0.5; lambda = 1; alpha = 2; kappa = 2
	
	if( (!is.numeric(x)) || (!is.numeric(lambda)) 
		|| (!is.numeric(kappa)) || (!is.numeric(alpha)) )
		stop("non-numeric argument to mathematical function")
	
	if( (length(lambda)!= 1) || (length(kappa)!= 1)
		|| (length(alpha)!= 1) )
		stop("Non-x parameters must be atomic")
	
	log.pdf <- log(kappa) + (kappa - 1) * 
		pweibull(q = x,scale = lambda,shape = alpha,log.p = TRUE) +
		dweibull(x = x,scale = lambda,shape = alpha,log = TRUE)
	
	# mpfr(return(exp(log.pdf)), 128L)
	ifelse(log,return(log.pdf),return(exp(log.pdf)))
	
}

#' @title qexpweibull
#' @description Inv-quantile
#' @inheritParams pexpweibull
#' @param p Probability
#' @export
qexpweibull = function(p,lambda,alpha,kappa){
	# p = 0.5; lambda = 1; alpha = 2; kappa = 2
	# lambda * (-log(1-p^(1/kappa)))^(1/alpha)
	
	quant <- qweibull(p^(1/kappa),scale = lambda,shape = alpha)
	return(quant)
} 

sim_DATA = function(copula = "Clayton",dist1,n_obs,theta,
	alpha1,lambda1,kappa1,alpha2,lambda2){
	
	if(FALSE){
		copula = c("Clayton","Gumbel")[1]
		n_obs = 1e3; theta = 3/2
		alpha1 = 5; lambda1 = 1.2; kappa1 = 1
		alpha2 = 2; lambda2 = 2.1
		
	}
	
	if( copula == "Clayton" ){
		if( theta < 0 ) stop("theta must be >= 0")
		tmp_copula = claytonCopula(param = theta,
			dim = 2,use.indepC = c("message","TRUE","FALSE"))
	} else if( copula == "Gumbel" ){
		if( theta < 1 ) stop("theta must be >= 1")
		tmp_copula = gumbelCopula(param = theta,
			dim = 2,use.indepC =c("message","TRUE","FALSE"))
	} else {
		stop(sprintf("No code for copula = %s.",copula))
	}
	
	if( dist1 == "weibull" ){
		bvd = mvdc(copula = tmp_copula,
			margins = c("weibull","weibull"),
			paramMargins = list(
				list(shape = alpha1,scale = lambda1),
				list(shape = alpha2,scale = lambda2))
		)
	} else if( dist1 == "expweibull" ){
		
		bvd = mvdc(copula = tmp_copula,
			margins = c("expweibull","weibull"),
			paramMargins = list(
				list(alpha = alpha1,lambda = lambda1,kappa = kappa1),
				list(shape = alpha2,scale = lambda2))
			)
	} else {
		stop(sprintf("No code for dist1 = %s.",dist1))
	}
	
	all_data = rMvdc(n_obs,bvd)
	
	all_data
}
get_PROFILE = function(GRID){
	if(FALSE){
		GRID = gout$DAT
		
	}
	
	GRID = smart_df(GRID)
	names(GRID) = c("log_alpha1","log_lambda1",
		"unc_kappa1","log_theta","LL")
	GRID[1:5,]
	
	## New code, get profiled LL
	u_logA = sort(unique(GRID$log_alpha1))
	u_logL = sort(unique(GRID$log_lambda1))
	res = c()
	for(logA in u_logA){
	for(logL in u_logL){
		# logA = u_logA[1]; logL = u_logL[1]
		idx = which(GRID$log_alpha1 == logA
			& GRID$log_lambda1 == logL)
		length(idx)
		tmp_df = GRID[idx[which.max(GRID$LL[idx])],]
		tmp_df
		res = rbind(res,tmp_df)
	}}
	rownames(res) = NULL
	# dim(res); res[1:10,]
	
	return(res)
}
get_LOCAL_OPTS = function(GPROF,verb){
	if(FALSE){
		GPROF = gprof
		verb = verb
		
	}
	
	# Find all local optimal solutions
	verb = get_verb(verb = verb)
	if( verb ) cat(sprintf("%s: Find all local optimums\n",date()))
	res = c(); nn = nrow(GPROF)
	
	EPS = sort(unique(GPROF$log_alpha1))
	EPS = min(unique(diff(EPS)))
	EPS = c(EPS,sort(unique(GPROF$log_lambda1)))
	EPS = min(abs(diff(EPS))) + 0.1
	# EPS
	
	for(ii in seq(nn)){
		# ii = 1
		# GPROF[ii,]
		if( verb ) smart_progress(ii = ii,nn = nn,iter = 1e1,iter2 = 5e2)
		curr_PARS = as.numeric(GPROF[ii,c("log_alpha1","log_lambda1")])
		# curr_PARS
		curr_LL = GPROF$LL[ii]; # curr_LL
		while(TRUE){
			# Get neighborhood around curr_PARS
			curr_PARS; curr_LL
			idx = which(abs(GPROF$log_alpha1 - curr_PARS[1]) < EPS
				& abs(GPROF$log_lambda1 - curr_PARS[2]) < EPS)
			length(idx)
			GPROF[idx,]
			
			# Move to next optimal point
			tmp_df = GPROF[idx[which.max(GPROF$LL[idx])],]
			if( tmp_df$LL == curr_LL
				& tmp_df$log_alpha1 == curr_PARS[1]
				& tmp_df$log_lambda1 == curr_PARS[2] ){
				break
			}
			curr_PARS = as.numeric(tmp_df[,c("log_alpha1","log_lambda1")])
			curr_LL = tmp_df$LL
		}
		res = rbind(res,tmp_df)
		if( ii >= 5 ) res = unique(res)
	}
	res = res[order(res[,1],res[,2]),]
	rownames(res) = NULL
	dim(res); # res
	
	return(res)
}

#' @title plot_LL
#' @description Visualize grid of profile log likelihood
#' @param GPROF An outputted data.frame from \code{get_PROFILE()}.
#' @param GOPT An outputted data.frame from \code{get_LOCAL_OPTS()}.
#' @param nBREAKs A positive integer set to 5 by default to specify
#'	the number of breaks in axes ticks.
#' @param COPULA A string input, either "Clayton" or "Gumbel"
#' @export
plot_LL = function(GPROF,GOPT = NULL,nBREAKs = 5,COPULA){
	if(FALSE){
		GPROF 	= gprof
		GOPT 		= gopt
		nBREAKs = 5
		COPULA 	= PARAMS$copula
		
	}
	
	res = GPROF
	res$LL[res$LL == min(res$LL)] = NA
	low_quant = quantile(res$LL,0.3,na.rm = TRUE); low_quant
	table(res$LL >= low_quant)
	res$LL[res$LL < low_quant] = NA
	
	# Focus on tighter range of values
	xrange = range(res$log_alpha1[!is.na(res$LL)]); xrange
	yrange = range(res$log_lambda1[!is.na(res$LL)]); yrange
	res = res[which(res$log_alpha1 >= xrange[1]
		& res$log_alpha1 <= xrange[2]
		& res$log_lambda1 >= yrange[1]
		& res$log_lambda1 <= yrange[2]),]
	res$LL[is.na(res$LL)] = min(res$LL,na.rm = TRUE)
	
	xbreaks = round(seq(xrange[1],xrange[2],length.out = nBREAKs),1)
	ybreaks = round(seq(yrange[1],yrange[2],length.out = nBREAKs),1)
	
	LL = log_alpha1 = log_lambda1 = log_theta = NULL
	
	g1 = ggplot(data = res,aes(x = log_alpha1,y = log_lambda1)) +
		geom_tile(aes(fill = LL)) + 
		coord_cartesian(expand = 0,xlim = xrange,ylim = yrange) +
		scale_x_continuous(breaks = xbreaks) +
		scale_y_continuous(breaks = ybreaks) +
		xlab(expression(paste("log(",alpha[1],")"))) +
		ylab(expression(paste("log(",lambda[1],")"))) +
		labs(fill = "Profile LL") +
		scale_fill_viridis(discrete = FALSE)
	
	g1 = g1 + theme(text = element_text(size = 20),
			plot.title = element_text(hjust = 0.5),
			legend.key.height = unit(4,"line"),
			legend.key.width = unit(1,"cm"),
			legend.position = c("right","bottom")[1])
	
	# Add local opt solutions
	uniq_unc_kappa1 = unique(res$unc_kappa1)
	DIST = ifelse(length(uniq_unc_kappa1) != 1,
		"Multiple dists.",
		ifelse(uniq_unc_kappa1 == 0,"Weibull","Exp-Weibull"))
	TITLE = sprintf("%s, %s",COPULA,DIST)

	if( !is.null(GOPT) ){
		GOPT$log_theta = round(GOPT$log_theta,2)
		g1 = g1 + geom_point(data = GOPT,
			mapping = aes(x = log_alpha1,y = log_lambda1),
			color = "red",size = 5,shape = 18)
		g1 = g1 + geom_hline(yintercept = GOPT$log_lambda1,
			color = "red",linetype = 3,size = 0.75) +
			geom_vline(xintercept = GOPT$log_alpha1,
				color = "red",linetype = 3,size = 0.75)
		if( COPULA == "Clayton" ){
			g1 = g1 + geom_text(data = GOPT,
				aes(label = sprintf("log(\u03B8) = %s",log_theta)),
				vjust = -1,hjust = -0.1,size = 5,color = "red")
		} else if( COPULA == "Gumbel" ){
			g1 = g1 + geom_text(data = GOPT,
				aes(label = sprintf("log(\u03B8-1) = %s",log_theta)),
				vjust = -1,hjust = -0.1,size = 5,color = "red")
		}
		
	}
	
	g1 = g1 + ggtitle(TITLE)
	# g1
	
	return(g1)
}


# ----------
# Main functions
# ----------

#' @title refData_match
#' @description This function takes as input
#'	a working dataset of interest and a 
#'	reference dataset.
#' @param wDAT A working dataset data.frame containing required 
#'	columns \code{age}, \code{time}, \code{delta}, 
#'	\code{datediag_yr}, \code{dateEvent_yr}, and \code{sex} 
#'	corresponding to age, observed time, event status, 
#'	diagnosis year, observed year, and sex (coded 'female'/'male'), 
#'	respectively.
#' @param rDAT A reference dataset data.frame containing
#'	required columns \code{Year}, \code{Age}, \code{qx},
#'	and \code{sex} corresponding to reference year, 
#'	reference age, event hazard, and sex, respectively.
#' @export
refData_match = function(wDAT,rDAT){
	
	# wDAT = working dataset
	req_names1 = c("age","time","delta","datediag_yr","dateEvent_yr","sex")
	if( !all(req_names1 %in% names(wDAT)) ){
		miss_names = req_names1[!(req_names1 %in% names(wDAT))]
		stop(sprintf("Missing colnames in wDAT: %s",
			paste(miss_names,collapse = ", ")))
	}
	
	# rDAT = reference dataset
	req_names2 = c("Year","Age","qx","sex")
	# qx = hazard
	if( !all(req_names2 %in% names(rDAT)) ){
		miss_names = req_names2[!(req_names2 %in% names(rDAT))]
		stop(sprintf("Missing colnames in rDAT: %s",
			paste(miss_names,collapse = ", ")))
	}
	
	# Check missingness
	if( any(is.na(wDAT[,req_names1])) )
		stop("wDAT: Fix or remove rows of missing data")
	if( any(is.na(rDAT[,req_names2])) )
		stop("rDAT: Fix or remove rows of missing data")
	
	# Check variable classes
	if( !chk_int_vec(wDAT$age) ){
		stop("Set wDAT$age to integer")
	}
	if( !is(wDAT$time,"numeric") ){
		stop("Set wDAT$time to integer or numeric")
	}
	if( !chk_int_vec(wDAT$delta) ){
		stop("Set wDAT$delta to integers")
	}
	if( !chk_int_vec(wDAT$datediag_yr) ){
		stop("Set wDAT$datediag_yr to integers")
	}
	if( !chk_int_vec(wDAT$dateEvent_yr) ){
		stop("Set wDAT$dateEvent_yr to integers")
	}
	if( !is(wDAT$sex,"character") ){
		stop("Set wDAT$sex to character")
	}
	if( !all(wDAT$sex %in% c("female","male")) ){
		stop("wDAT$sex takes values 'female' and 'male'")
	}
	
	if( !chk_int_vec(rDAT$Year) ){
		stop("Set rDAT$Year to integers")
	}
	if( !chk_int_vec(rDAT$Age) ){
		stop("Set rDAT$Age to integers")
	}
	if( !is(rDAT$qx,"numeric") ){
		stop("Set rDAT$qx to numeric")
	}
	if( !is(rDAT$sex,"character") ){
		stop("Set rDAT$sex to character")
	}
	if( !all(rDAT$sex %in% c("female","male")) ){
		stop("rDAT$sex takes values 'female' and 'male'")
	}
	
	wDAT$sex2 = ifelse(wDAT$sex == "female",1,0)
	rDAT$sex2 = ifelse(rDAT$sex == "female",1,0)
	rDAT = rDAT[order(rDAT$Year,rDAT$Age,rDAT$sex),]
	
	# Do matching on age, sex, time of diagnosis, time of death
	cat("Get matching density and survival data ...\n")
	wDAT_vars = c("age","datediag_yr","dateEvent_yr","sex2")
	rDAT_vars = c("Year","Age","qx","sex2")
	OUT = dMrs_MATCH(wDAT = as.matrix(wDAT[,wDAT_vars]),
		rDAT = as.matrix(rDAT[,rDAT_vars]),ncores = 1,verb = TRUE)
	OUT = smart_df(OUT)
	names(OUT) = c("dens_t2","surv_t2")
	wDAT = cbind(wDAT,OUT)
	
	return(wDAT)
	
}

#' @title sim_replicate
#' @description This function simulates an analysis-ready dataset
#'	with arguments to specify the sample size and copula-related and 
#'	survival-related parameters.
#' @param copula A string input, either "Clayton" or "Gumbel"
#' @param dist1 A string input for distribution of work dataset 
#'	event times, either "weibull" or "expweibull" for Weibull with
#'	two parameters or exponentiated weibull.
#' @param NN Positive integer for simulated sample size
#' @param theta A positive value corresponding to the event time
#'	dependency depending on the copula. For Clayton, \code{theta} >= 0
#'	and for Gumbel, \code{theta} >= 1.
#' @param alpha1 First shape parameter for the first event.
#' @param lambda1 Scale parameter for the first event.
#' @param kappa1 Second shape parameter for the first
#'	event's exponentiated Weibull.
#' @param alpha2 Shape parameter for the second/reference event.
#' @param lambda2 Scale parameter for the second/reference event.
#' @param propC A numeric value between 0 and 1 for
#'	proportion of subjects with censored event times
#' @param verb Boolean, set to TRUE for verbose output.
#' @export
sim_replicate = function(copula,dist1,NN,theta,
	alpha1,lambda1,kappa1,alpha2,lambda2,propC,verb){
	
	copula = get_copula(copula = copula)
	verb = get_verb(verb = verb)
	if( verb ) cat(sprintf("%s: Simulate dataset ...\n",date()))
	
	# Simulate dataset
	kappa1 = ifelse(dist1 == "weibull",1,kappa1)
	dat = suppressMessages(sim_DATA(copula = copula,dist1 = dist1,
		n_obs = NN,theta = theta,alpha1 = alpha1,
		lambda1 = lambda1,kappa1 = kappa1,
		alpha2 = alpha2,lambda2 = lambda2))
	# names(dat); head(dat); dim(dat)
	dat = smart_df(dat); names(dat) = paste0("T",1:2); # head(dat)
	dat$T = apply(dat[,c("T1","T2")],1,function(xx) min(xx)); # head(dat)
	dat$D = apply(dat[,c("T1","T2")],1,function(xx) ifelse(xx[1] <= xx[2],1,0))
	# hist(dat$T,breaks=40,col="gray")
	# table(dat$D)
	
	# Tune prop subjects censored
	max_delta = max(dat$T) + 1; max_delta
	prop_event = 1 - propC; # prop_event
	width_delta = 1
	if( propC == 0 ){
		CC = rep(max_delta,NN)
	} else {
		while(TRUE){
			vec_delta = c(-1,1) * width_delta + max_delta; vec_delta
			# max_delta
			mat_CC = sapply(vec_delta,function(dd){
				runif(n = NN,0,dd)
			})
			vec_propO = apply(mat_CC,2,function(xx){
				mean(dat$T <= xx)
			})
			# print(vec_propO)
			
			if( vec_propO[1] <= prop_event && vec_propO[2] >= prop_event ){
				if( max(abs(vec_propO - prop_event)) < 1e-3 ){
					break
				} else {
					width_delta = width_delta * 0.8
				}
			} else if( vec_propO[1] > prop_event ){
				max_delta = max_delta * 0.8
			} else if( vec_propO[2] < prop_event ){
				max_delta = max_delta * 1.2
			}
			
		}
		CC = mat_CC[,1]
	}
	max_delta; width_delta
	dat$C = CC
	
	## Obtain observed time and event indicator
	dat$time = apply(dat[,c("T","C")],1,function(xx) min(xx))
	dat$delta = apply(dat[,c("T","C")],1,function(xx) ifelse(xx[1] <= xx[2],1,0))
	if( verb ) print(table(dat$delta) / NN)

	## Get dens_t2 and surv_t2, we assume alpha2/lambda2 are known and from Weibull
	dat$dens_t2 = dweibull(x = dat$time,shape = alpha2,scale = lambda2)
	dat$surv_t2 = 1 - pweibull(q = dat$time,shape = alpha2,scale = lambda2)
	
	params = smart_df(copula = copula,NN = NN,theta = theta,
		alpha1 = alpha1,lambda1 = lambda1,kappa1 = kappa1,
		alpha2 = alpha2,lambda2 = lambda2)
	
	list(DATA = dat,PARAMS = params)
	
}

opt_replicate = function(REP,param_grid,
	theta = NULL,upKAPPA,ncores = 1,gTHRES = 8e-2,verb){

	if(FALSE){
		REP 				= REP
		theta 			= theta
		upKAPPA 		= upKAPPA
		gTHRES 			= gTHRES
		verb 				= verb
		ncores 			= ncores
		
	}
	
	verb 		= get_verb(verb = verb)
	upTHETA = ifelse(is.null(theta),1,0)
	DATA 		= REP$DATA
	PARAMS 	= REP$PARAMS
	upPARS	= get_upPARS(upKAPPA = upKAPPA,
							THETA = theta)
	upTHETA	= upPARS[4]
	
	wrap_LL 	= function(PARS){
		# PARS = iPARS
		dMrs_cLL(XX = DATA$time,DELTA = DATA$delta,
			D2 = DATA$dens_t2,S2 = DATA$surv_t2,
			PARS = PARS,copula = PARAMS$copula,
			verb = !TRUE)
	}
	wrap_GRAD = function(PARS){
		# PARS = iPARS
		out = dMrs_cGRAD(XX = DATA$time,DELTA = DATA$delta,
			D2 = DATA$dens_t2,S2 = DATA$surv_t2,PARS = PARS,
			copula = PARAMS$copula,upPARS = upPARS)
		
		# out = grad(wrap_LL,PARS)
		c(out)
	}
	wrap_HESS = function(PARS){
		# PARS = iPARS
		dMrs_cHESS(XX = DATA$time,DELTA = DATA$delta,
			D2 = DATA$dens_t2,S2 = DATA$surv_t2,PARS = PARS,
			copula = PARAMS$copula,upPARS = upPARS)
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
		if( PARAMS$copula == "Clayton" ) 	log_THETA = log(theta)
		if( PARAMS$copula == "Gumbel" ) 	log_THETA = log(theta - 1)
	} else if( upTHETA == 1 ){
		log_THETA = c(-Inf,param_grid)
		log_THETA = sort(unique(log_THETA))
	}
	
	gout = dMrs_GRID(XX = DATA$time,DELTA = DATA$delta,
		D2 = DATA$dens_t2,S2 = DATA$surv_t2,
		log_ALPHA = param_grid,log_LAMBDA = param_grid,
		unc_KAPPA = unc_KAPPA,log_THETA = log_THETA,
		copula = PARAMS$copula,verb = verb,ncores = ncores)
	# str(gout)
	colnames(gout$DAT) = c("log_alpha1","log_lambda1",
		"unc_kappa1","log_theta","LL")
	gout$DAT[1:5,]
	gprof = get_PROFILE(GRID = gout$DAT)
	# dim(gprof); head(gprof)
	
	gopt = get_LOCAL_OPTS(GPROF = gprof,verb = verb)
	# dim(gopt); head(gopt)
	
	if( verb ){
		cat(sprintf("%s: Two-dimensional plot of the profile log likelihood ...\n",date()))
		print(plot_LL(GPROF = gprof,GOPT = gopt,
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
			PARS = iPARS,copula = PARAMS$copula,upPARS = upPARS,
			max_iter = 2e2,eps = 1e-6,verb = FALSE)
		if(FALSE){
		
		upPARS = c(1,1,0,0)
		dMrs_BFGS(XX = DATA$time,DELTA = DATA$delta,
			D2 = DATA$dens_t2,S2 = DATA$surv_t2,
			PARS = iPARS,copula = PARAMS$copula,upPARS = upPARS,
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
	nab_par[4] = ifelse(PARAMS$copula == "Clayton",
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
	nparams = 2 								# est. alpha1 + lambda1
	nparams = nparams + upKAPPA	# est. kappa
	nparams = nparams + upTHETA # est. theta
	NN = nrow(DATA)
	BIC = 2*LL - nparams * log(NN)
	
	
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
	gTHRES,copula,param_grid,vec_time,verb,ncores = 1){
	
	if(FALSE){
		DATA 				= one_rep$DATA
		theta 			= NULL
		upKAPPA 		= c(0,1)[2]
		gTHRES 			= 8e-2
		copula 			= c("Clayton","Gumbel")[1]
		param_grid 	= param_grid
		vec_time 		= vec_time
		
		verb = TRUE; ncores = 1
		
	}
	
	verb = get_verb(verb = verb)
	req_names = c("time","delta","dens_t2","surv_t2")
	if( !all(req_names %in% names(DATA)) ){
		miss_names = req_names[!(req_names %in% names(DATA))]
		stop(sprintf("DATA colnames missing: %s",
			paste(miss_names,collapse = ", ")))
	}
	
	REP = list(DATA = DATA,PARAMS = smart_df(copula = copula))
	opt_out = opt_replicate(REP = REP,param_grid = param_grid,
		theta = theta,upKAPPA = upKAPPA,gTHRES = gTHRES,
		verb = verb,ncores = ncores)
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
#'	the Weibull distribution is assumed and optimized over.
#' @param gTHRES A numeric threshold on the L2 norm of 
#'	the gradient evaluated at the MLE.
#' @param param_grid Vector of values spanning possible 
#'	log(alpha1), log(lambda1), log(kappa1), unconstrained 
#'	theta parameters
#' @inheritParams sim_replicate
#' @param COPULAS If \code{NULL}, will optimize over all copulas. 
#'	Otherwise set to 'Clayton' or 'Gumbel'
#' @param vec_time Vector of times to calculate predicted survival
#'	on the same scale as times provided in the \code{DATA} data.frame.
#' @param ncores An integer for the number of threads
#' @export
run_analyses = function(DATA,THETAs = NULL,upKAPPA = NULL,
	gTHRES = 8e-2,COPULAS = NULL,param_grid,vec_time,ncores = 1,verb){
	
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
	
	COPULAS 		= get_copula(copula = COPULAS)
	verb 				= get_verb(verb = verb)
	vec_time		= get_vecTIME(TIME = DATA$time,vec_time = vec_time)
	upKAPPA			= get_upKAPPA(upKAPPA = upKAPPA)
	if( missing(param_grid) ) param_grid = get_parGRID()
	
	# Multiple copula s
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
				ncores = ncores,verb = verb)
			
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
	
	if( COPULAS == "Clayton" ){
		if( !is.null(THETAs) && any(THETAs < 0) )
			THETAs = THETAs[THETAs >= 0]
	} else if( COPULAS == "Gumbel" ){
		if( !is.null(THETAs) && any(THETAs < 1) )
			THETAs = THETAs[THETAs >= 1]
	}
	
	out = list()
	
	if( is.null(THETAs) ){
		tmp_out = run_analysis(DATA = DATA,
			theta = THETAs,upKAPPA = upKAPPA,gTHRES = gTHRES,
			copula = COPULAS,param_grid = param_grid,
			vec_time = vec_time,verb = verb,ncores = ncores)
		out[[1]] = tmp_out
		return(out)
	}
	
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
			vec_time = vec_time,verb = verb,ncores = ncores)
		out[[cnt]] = tmp_out
		cnt = cnt + 1
	}
	
	return(out)
	
}

#' @title plot_SURVs
#' @param run_ANA The object outputed from \code{run_analyses}
#'	function.
#' @param MULTIPLE A boolean set to \code{TRUE} to display
#'	a survival plot per \code{theta} with accompanying maximum 
#'	log likelihood. Otherwise, the survival plots will be 
#'	overlayed with a legend.
#' @param ALPHA A numeric value between 0 and 1 to control the
#'	confidence band transparency.
#' @export
plot_SURVs = function(run_ANA,MULTIPLE,ALPHA = 0.5){
	
	if(FALSE){
		run_ANA 	= run_ana
		MULTIPLE 	= TRUE
		ALPHA 		= 0.5
		
	}
	
	if( is.null(run_ANA) || length(run_ANA) == 0 ){
		stop("No results!")
	}
	
	# Aggregate predictions
	res = c()
	for(vv in seq(length(run_ANA))){
		# vv = 1
		
		# Get MLEs, distribution, copula
		COPULA 		= run_ANA[[vv]]$copula
		cout 			= run_ANA[[vv]]$RES$cout; cout
		tmpALPHA	= cout$EST[which(cout$PARS == "alpha1")]; tmpALPHA
		tmpLAMBDA	= cout$EST[which(cout$PARS == "lambda1")]; tmpLAMBDA
		THETA 		= cout$EST[which(cout$PARS == "theta")]; THETA
		KAPPA		= cout$EST[which(cout$PARS == "kappa1")]; KAPPA
		DIST		= ifelse(KAPPA == 1,"Weibull","Exp-Weibull")
		LL			= run_ANA[[vv]]$RES$LL
		
		tmp_df = smart_df(COPULA = COPULA,ALPHA = tmpALPHA,LAMBDA = tmpLAMBDA,
			THETA = THETA,KAPPA = KAPPA,DIST = DIST,LL = LL,
			run_ANA[[vv]]$PRED[,c("time","surv","low_surv2","high_surv2")])
		
		names(tmp_df) = c("COPULA","ALPHA","LAMBDA",
			"THETA","KAPPA","DIST","LL",
			"time","surv","low_surv","high_surv")
		head(tmp_df)
		
		res = rbind(res,tmp_df)
	}
	res$LABEL = sprintf("%s + %s, \u03B1 = %s, \u03BB = %s, \n\u03B8 = %s, \u03BA = %s, LL = %s",
			res$COPULA,res$DIST,round(res$ALPHA,2),round(res$LAMBDA,2),
			round(res$THETA,2),round(res$KAPPA,2),round(res$LL,2))
	res = res[order(res$COPULA,res$DIST,res$THETA,res$KAPPA),]
	uLABEL = unique(res$LABEL)
	res$LABEL = factor(res$LABEL,levels = uLABEL)
	dim(res); res[1:5,]
	
	time = surv = low_surv = high_surv = LABEL = NULL
	
	if( MULTIPLE ){
		my_themes = theme(text = element_text(size = 28),
			axis.text.x = element_text(size = 20),
			legend.position = c("none","bottom")[1],
			panel.spacing = unit(1,"lines"),
			panel.background = element_blank(),
			panel.grid.major = element_line(colour = "grey50",
				size = 0.5,linetype = "dotted"),
			plot.title = element_text(hjust = 0.5),
			panel.border = element_rect(colour = "black",
				fill = NA,size = 1),
			strip.text.x = element_text(size = 15))
		
		gg = ggplot(data = res,mapping = aes(x = time,y = surv)) +
			geom_line(size = 1,aes(color = LABEL)) + 
			geom_ribbon(mapping = aes(ymin = low_surv,
				ymax = high_surv,fill = LABEL),alpha = ALPHA) +
			ylim(c(0,1)) + xlab("Time") + ylab("Survival Probability") +
			# ggtitle(my_title) + 
			facet_wrap(~ LABEL) + 
			my_themes
		# gg
		
	} else {
		my_themes = theme(text = element_text(size = 28),
			legend.position = c("none","bottom","right")[3],
			plot.title = element_text(hjust = 0.5),
			panel.background = element_blank(),
			panel.grid.major = element_line(colour = "grey50",
				size = 0.5,linetype = "dotted"),
			panel.border = element_rect(colour = "black",
				fill = NA,size = 1))
		
		gg = ggplot(data = res,
			mapping = aes(x = time,y = surv,group = LABEL,
				fill = LABEL)) +
			geom_line(size = 1.25,alpha = 1,
				aes(color = LABEL),show.legend = FALSE) +
			geom_ribbon(mapping = aes(ymin = low_surv,
				ymax = high_surv),alpha = ALPHA) +
			ylim(c(0,1)) + xlab("Time") + ylab("Survival Probability") +
			# ggtitle(my_title) + 
			labs(fill = "") + my_themes
	}
	
	return(gg)
}




###

