
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
get_LOCAL_OPTS = function(GPROF,show = TRUE){
	if(FALSE){
		GPROF = gprof
		rownames(GPROF) = NULL
		show = TRUE
		
	}
	
	# Find all local optimal solutions
	if( show ) cat(sprintf("%s: Find all local optimums\n",date()))
	res = c(); nn = nrow(GPROF)
	
	EPS = sort(unique(GPROF$log_alpha1))
	EPS = min(unique(diff(EPS)))
	EPS = c(EPS,sort(unique(GPROF$log_lambda1)))
	EPS = min(abs(diff(EPS))) + 0.1
	# EPS
	
	for(ii in seq(nn)){
		# ii = 1
		# GPROF[ii,]
		if( show ) smart_progress(ii = ii,nn = nn,iter = 1e1,iter2 = 5e2)
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
#' @param TITLE String for plot title
#' @export
plot_LL = function(GPROF,GOPT = NULL,nBREAKs = 5,TITLE = NULL){
	if(FALSE){
		GPROF = gprof
		GOPT = get_LOCAL_OPTS(GPROF = GPROF)
		
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
	LL = log_alpha1 = log_lambda1 = NULL
	
	g1 = ggplot(data = res,aes(x = log_alpha1,y = log_lambda1)) +
		geom_tile(aes(fill = LL)) + 
		coord_cartesian(expand = 0,xlim = xrange,ylim = yrange) +
		scale_x_continuous(breaks = xbreaks) +
		scale_y_continuous(breaks = ybreaks) +
		xlab(expression(paste("log(",alpha[1],")"))) +
		ylab(expression(paste("log(",lambda[1],")"))) +
		labs(fill = "Profile LL") +
		scale_fill_viridis(discrete = FALSE)
	
	if( !is.null(TITLE) ) g1 = g1 + ggtitle(TITLE)
	
	g1 = g1 + theme(text = element_text(size = 20),
			plot.title = element_text(hjust = 0.5),
			legend.key.height = unit(4,"line"),
			legend.key.width = unit(1,"cm"),
			legend.position = c("right","bottom")[1])
	
	# Add local opt solutions
	if( !is.null(GOPT) ){
		g1 = g1 + geom_point(data = GOPT,
			mapping = aes(x = log_alpha1,y = log_lambda1),
			color = "red",size = 5,shape = 18)
		g1 = g1 + geom_hline(yintercept = GOPT$log_lambda1,
			color = "red",linetype = 3,size = 0.75) +
			geom_vline(xintercept = GOPT$log_alpha1,
				color = "red",linetype = 3,size = 0.75)
	}
	
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
	if( class(wDAT$age) != "integer" ){
		stop("Set wDAT$age to integer")
	}
	if( !(class(wDAT$time) %in% c("integer","numeric")) ){
		stop("Set wDAT$time to integer or numeric")
	}
	if( !(class(wDAT$delta) %in% c("integer")) ){
		stop("Set wDAT$delta to integer")
	}
	if( !(class(wDAT$datediag_yr) %in% c("integer")) ){
		stop("Set wDAT$datediag_yr to integer")
	}
	if( !(class(wDAT$dateEvent_yr) %in% c("integer")) ){
		stop("Set wDAT$dateEvent_yr to integer")
	}
	if( !(class(wDAT$sex) %in% c("character")) ){
		stop("Set wDAT$sex to character")
	}
	if( !all(wDAT$sex %in% c("female","male")) ){
		stop("wDAT$sex takes values 'female' and 'male'")
	}
	
	if( !(class(rDAT$Year) %in% c("integer")) ){
		stop("Set rDAT$Year to integer")
	}
	if( !(class(rDAT$Age) %in% c("integer")) ){
		stop("Set rDAT$Age to integer")
	}
	if( !(class(rDAT$qx) %in% c("numeric")) ){
		stop("Set rDAT$qx to numeric")
	}
	if( !(class(rDAT$sex) %in% c("character")) ){
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
		rDAT = as.matrix(rDAT[,rDAT_vars]),ncores = 1,show = TRUE)
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
#' @param show Boolean, set to TRUE for verbose output.
#' @export
sim_replicate = function(copula,dist1,NN,theta,
	alpha1,lambda1,kappa1,alpha2,lambda2,propC,show){
	
	if( show ) cat(sprintf("%s: Simulate dataset ...\n",date()))
	
	# Simulate dataset
	kappa1 = ifelse(dist1 == "weibull",1,kappa1)
	dat = suppressMessages(sim_DATA(copula = copula,dist1 = dist1,
		n_obs = NN,theta = theta,alpha1 = alpha1,
		lambda1 = lambda1,kappa1 = kappa1,
		alpha2 = alpha2,lambda2 = lambda2))
	# names(dat); class(dat); head(dat); dim(dat)
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
	if( show ) print(table(dat$delta) / NN)

	## Get dens_t2 and surv_t2, we assume alpha2/lambda2 are known and from Weibull
	dat$dens_t2 = dweibull(x = dat$time,shape = alpha2,scale = lambda2)
	dat$surv_t2 = 1 - pweibull(q = dat$time,shape = alpha2,scale = lambda2)
	
	params = smart_df(copula = copula,NN = NN,theta = theta,
		alpha1 = alpha1,lambda1 = lambda1,kappa1 = kappa1,
		alpha2 = alpha2,lambda2 = lambda2)
	
	list(DATA = dat,PARAMS = params)
	
}

#' @title opt_replicate
#' @description This function optimizes the dataset generated by 
#'	\code{sim_replicate}. First this function runs a grid search
#'	for optimal initial parameter estimates and then proceeds to
#'	optimize using BFGS.
#' @param REP An object generated by \code{sim_replicate}.
#' @param param_grid Vector of values spanning possible 
#'	log(alpha1), log(lambda1), and log(kappa1) parameters
#' @inheritParams sim_replicate
#' @param upKAPPA An integer value taking values 0 or 1. If set 
#'	to 1, the exponentiated Weibull distribution is assumed. Otherwise,
#'	the Weibull distribution is assumed and optimized over.
#' @param gTHRES A numeric threshold on the L2 norm of 
#'	the gradient evaluated at the MLE.
#' @param ncores An integer for the number of threads
#' @export
opt_replicate = function(REP,param_grid,
	theta = NULL,upKAPPA,ncores = 1,gTHRES = 8e-2,show){
	
	if(FALSE){
		REP = one_rep; param_grid = seq(-2,3,0.5)
		theta = NULL; upKAPPA = upKAPPA;
		ncores = 1; show = FALSE
		
		# Test inputs
		REP = REP; param_grid = param_grid
		theta = theta; upKAPPA = upKAPPA; 
		show = show; ncores = ncores
		
	}
	
	upTHETA 	= ifelse(is.null(theta),1,0)
	DATA 			= REP$DATA
	PARAMS 		= REP$PARAMS
	upPARS		= rep(1,4)
	upPARS[3] = upKAPPA
	upPARS[4] = upTHETA
	
	wrap_LL 	= function(PARS){
		# PARS = iPARS
		dMrs_cLL(XX = DATA$time,DELTA = DATA$delta,
			D2 = DATA$dens_t2,S2 = DATA$surv_t2,
			PARS = PARS,copula = PARAMS$copula,
			show = !TRUE)
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
	if( show ) cat(sprintf("%s: Grid search for initial parameters ...\n",date()))
	if( upPARS[3] == 0 ){
		# logit_KAPPA = Inf
		unc_KAPPA = 0
	} else {
		# logit_KAPPA = param_grid
		unc_KAPPA = param_grid
	}
	# if( is.null(theta) && upTHETA == 0 ) stop("Specify theta or set upTHETA = 1!")
	if( upTHETA == 0 ){
		log_THETA = log(theta)
	} else if( upTHETA == 1 ){
		log_THETA = param_grid
	}
	
	gout = dMrs_GRID(XX = DATA$time,DELTA = DATA$delta,
		D2 = DATA$dens_t2,S2 = DATA$surv_t2,
		log_ALPHA = param_grid,log_LAMBDA = param_grid,
		unc_KAPPA = unc_KAPPA,log_THETA = log_THETA,
		copula = PARAMS$copula,show = show,ncores = ncores)
	# str(gout)
	colnames(gout$DAT) = c("log_alpha1","log_lambda1",
		"unc_kappa1","log_theta","LL")
	gout$DAT[1:5,]
	gprof = get_PROFILE(GRID = gout$DAT)
	gopt = get_LOCAL_OPTS(GPROF = gprof,show = show)
	dim(gopt); gopt
	
	if( show ){
		cat(sprintf("%s: Two-dimensional plot of the log likelihood ...\n",date()))
		print(plot_LL(GPROF = gprof,GOPT = gopt))
	}
	
	gopt$fin_LL = NA; gopt$nGRAD = NA; gopt$fin_logA = NA
	gopt$fin_logL = NA; gopt$fin_logK = NA; gopt$fin_logT = NA
	# gopt
	
	if( show ) cat(sprintf("%s: BFGS optimization ...\n",date()))
	for(ii in seq(nrow(gopt))){
		# ii = 1
		if( show ) cat(".")
		iPARS = as.numeric(gopt[ii,1:4]); iPARS
		dMrs_BFGS(XX = DATA$time,DELTA = DATA$delta,
			D2 = DATA$dens_t2,S2 = DATA$surv_t2,
			PARS = iPARS,copula = PARAMS$copula,upPARS = upPARS,
			max_iter = 2e2,eps = 1e-6,show = FALSE)
		tmp_LL = wrap_LL(PARS = iPARS)
		tmp_GR = wrap_GRAD(PARS = iPARS)
		
		gopt$nGRAD[ii] = Rcpp_norm(tmp_GR)
		gopt$fin_logA[ii] = iPARS[1]
		gopt$fin_logL[ii] = iPARS[2]
		gopt$fin_logK[ii] = iPARS[3]
		gopt$fin_logT[ii] = iPARS[4]
		gopt$fin_LL[ii] = round(tmp_LL,3)
	}
	if( show ){
		cat("\n")
		print(gopt)
	}
	gopt_pre = gopt
	
	# Remove non-local optimum solutions
	gopt = gopt[which(gopt$nGRAD < gTHRES),]
	if( nrow(gopt) == 0 ){
		if( show ) cat(sprintf("%s: No converged solutions! ...\n",date()))
		return(list(GRID = gout$DAT,GPROF = gprof,GOPT = gopt,
			GOPT_PRE = gopt_pre,
			out = NULL,cout = NULL,LL = NULL,GRAD = NULL,
			HESS = NULL,COVAR = NULL))
	}
	
	# Order remaining, remove duplicates
	gopt = gopt[order(-gopt$fin_LL),]
	rownames(gopt) = NULL
	gopt = gopt[!duplicated(round(gopt[,
		sprintf("fin_%s",c("LL","logA",
		"logL","logK","logT"))],3)),]
	# gopt
	
	# Estimate covariance
	if( show ) cat(sprintf("%s: Get covariance ...\n",date()))
	gopt$neg_var = NA
	for(jj in seq(nrow(gopt))){
		# jj = 1
		iPARS = as.numeric(gopt[jj,
			paste0("fin_log",c("A","L","K","T"))]); iPARS
		hess = wrap_HESS(PARS = iPARS); hess
		nz = which(diag(hess) != 0)
		# print(hess)
		if( rcond(hess[nz,nz]) == 0 ){
			next
			# stop("Error with hessian 1: rcond = 0\n")
		}
		covar = matrix(0,4,4)
		covar[nz,nz] = solve(-hess[nz,nz])
		# covar
		gopt$neg_var[jj] = any(diag(covar) < 0)
	}
	gopt = gopt[!is.na(gopt$neg_var),,drop = FALSE]
	gopt = gopt[which(gopt$neg_var == FALSE),]
	if( nrow(gopt) == 0 ){
		if( show ) cat(sprintf("%s: Negative variance(s)! ...\n",date()))
		return(list(GRID = gout$DAT,GPROF = gprof,GOPT = NULL,
			GOPT_PRE = gopt_pre,
			out = NULL,cout = NULL,LL = NULL,GRAD = NULL,
			HESS = NULL,COVAR = NULL))
	}
	gopt = gopt[which.max(gopt$fin_LL),]
	gopt
	
	iPARS = as.numeric(gopt[1,
		paste0("fin_log",c("A","L","K","T"))]); iPARS
	hess = wrap_HESS(PARS = iPARS); hess
	nz = which(diag(hess) != 0)
	if( rcond(hess[nz,nz]) == 0 ){
		print(hess)
		cat("Error with hessian 1: rcond = 0\n")
		return(list(GRID = gout$DAT,GPROF = gprof,GOPT = NULL,
			GOPT_PRE = gopt_pre,
			out = NULL,cout = NULL,LL = NULL,GRAD = NULL,
			HESS = NULL,COVAR = NULL))
	}
	covar = matrix(0,4,4)
	covar[nz,nz] = solve(-hess[nz,nz])
	covar
	
	out = smart_df(PARS = c("log_alpha1","log_lambda1",
		"unc_kappa1","log_theta"),
		EST = c(iPARS),SE = c(sqrt(diag(covar))))
	out$lowCI 	= out$EST - 1.96 * out$SE
	out$highCI 	= out$EST + 1.96 * out$SE
	out
	
	EST = exp(out$EST)
	if( PARAMS$copula == "Gumbel" ) EST[4] = EST[4] + 1
	
	cout = smart_df(PARS = c("alpha1","lambda1","kappa1","theta"),
		EST = EST,SE = NA)
	
	# Delta method: Calculate nabla g()
	# stop("Debug the SE delta method calculation!")
	nab_par	= rep(NA,4)
	nab_par[1:3] = cout$EST[1:3] # alpha, lambda, kappa
	nab_par[4] = ifelse(PARAMS$copula == "Clayton",
		cout$EST[4],exp(out$EST[4]))
	nab_par = diag(nab_par)
	
	cout$SE 			= sqrt(diag(nab_par %*% covar %*% nab_par))
	cout$lowCI 		= cout$EST - 1.96 * cout$SE
	cout$highCI 	= cout$EST + 1.96 * cout$SE
	cout$lowCI_2 	= cout$EST * exp(-1.96 * out$SE)
	cout$highCI_2 = cout$EST * exp(1.96 * out$SE)
	cout
	
	LL 		= wrap_LL(iPARS)
	GRAD 	= wrap_GRAD(iPARS)
	HESS 	= wrap_HESS(iPARS)
	
	names(GRAD) 		= out$PARS
	dimnames(HESS) 	= list(out$PARS,out$PARS)
	dimnames(covar) = list(out$PARS,out$PARS)
	
	return(list(GRID = gout$DAT,GPROF = gprof,GOPT = gopt,
		GOPT_PRE = gopt_pre,
		out = out,cout = cout,LL = LL,GRAD = GRAD,
		HESS = HESS,COVAR = covar))
}

#' @title full_sim
#' @description This function performs a full simulation after 
#'	supplying the corresponding copula, sample size, underlying
#'	parameters, number of simulation replicates, and grid of 
#'	values to search over.
#' @inheritParams sim_replicate
#' @inheritParams opt_replicate
#' @param RR Integer number of replicates
#' @export
full_sim = function(copula,dist1,NN,theta,
	alpha1,lambda1,kappa1,alpha2,lambda2,propC,RR,
	param_grid,upKAPPA,gTHRES = 8e-2,show,ncores = 1){
	
	if(FALSE){
		copula = copula; dist1 = dist1; NN = NN; theta = theta
		alpha1 = alpha1; lambda1 = lambda1; kappa1 = kappa1; alpha2 = alpha2
		lambda2 = lambda2; propC = propC; RR = 1e2; param_grid = seq(-1,3,0.5)
		upKAPPA = 0; show = TRUE
		
	}
	
	kappa1 = ifelse(dist1 == "weibull",1,kappa1)
	
	# set.seed(1); 
	res = c()
	for(rr in seq(RR)){
		# rr = 1
		if( show ) smart_progress(ii = rr,nn = RR)
		
		while(TRUE){
			one_rep = sim_replicate(copula = copula,dist1 = dist1,
				NN = NN,theta = theta,alpha1 = alpha1,lambda1 = lambda1,
				kappa1 = kappa1,alpha2 = alpha2,lambda2 = lambda2,
				propC = propC,show = FALSE)
			table(one_rep$DATA$D)
			one_opt = opt_replicate(REP = one_rep,
				param_grid = param_grid,theta = theta,
				upKAPPA = upKAPPA,gTHRES = gTHRES,
				show = FALSE,ncores = ncores)
			if( !is.null(one_opt$out) ) break
		}
		
		res = rbind(res,smart_df(RR = rr,one_opt$out))
		rm(one_rep)
	}

	# head(res)
	res$TRUTH = NA
		res$TRUTH[res$PARS == "log_alpha1"] = log(alpha1)
		res$TRUTH[res$PARS == "log_lambda1"] = log(lambda1)
		res$TRUTH[res$PARS == "log_kappa1"] = log(kappa1)
		res$TRUTH[res$PARS == "log_theta"] = log(theta)
	res$CP = ifelse(res$TRUTH >= res$lowCI & res$TRUTH <= res$highCI,1,0)
	res[1:10,]
	mean_PARS = sapply(unique(res$PARS),function(xx) 
		mean(res$EST[which(res$PARS == xx)]),USE.NAMES = FALSE)
	se_PARS = sapply(unique(res$PARS),function(xx) 
		sd(res$EST[which(res$PARS == xx)]),USE.NAMES = FALSE)
	see_PARS = sapply(unique(res$PARS),function(xx) 
		mean(res$SE[which(res$PARS == xx)]),USE.NAMES = FALSE)
	CP_PARS = sapply(unique(res$PARS),function(xx) 
		mean(res$CP[which(res$PARS == xx)]),USE.NAMES = FALSE)
	
	fin = smart_df(copula = copula,PARS = unique(res$PARS),
		NN = NN,TRUTH = log(c(alpha1,lambda1,kappa1,theta)),EST = mean_PARS)
	fin$BIAS = mean_PARS - fin$TRUTH
	fin$SE = se_PARS; fin$SEE = see_PARS
	fin$CP = CP_PARS
	# fin
	
	if( show ){
		par(mfrow = c(2,2),mar = c(4,4,1,0) + 0.2,
			cex.lab = 1.3,oma = c(0,0,2,0))
		for(vv in unique(res$PARS)){
			# vv = unique(res$PARS)[3]
			if( vv == "log_alpha1" ){
				vv2 = bquote(log(alpha[1]))
			} else if( vv == "log_lambda1" ){
				vv2 = bquote(log(lambda[1]))
			} else if( vv == "log_kappa1" ){
				vv2 = bquote(log(kappa[1]))
			} else {
				vv2 = bquote(log(theta))
			}
			tmp_vec = res$EST[res$PARS == vv]
			if( any(is.infinite(tmp_vec)) ) next
			if( sd(tmp_vec) == 0 ) next
			hist(tmp_vec,breaks = 30,main = "",
				xlab = vv2,freq = FALSE,col = "gray")
			true_value = ifelse(vv == "log_alpha1",log(alpha1),
				ifelse(vv == "log_lambda1",log(lambda1),
				ifelse(vv == "log_kappa1",log(kappa1),log(theta))))
			abline(v = true_value,lwd = 2,lty = 2,col = "green")
			abline(v = mean(res$EST[res$PARS==vv]),lwd = 3,lty = 3,col = "red")
			abline(v = median(res$EST[res$PARS==vv]),lwd = 2,lty = 3,col = "blue")
			legend("topleft",legend = c("true value","est mean","est median"),
				pch = rep(16,3),col = c("green","red","blue"),cex = 0.8)
		}
		mtext(bquote("Sim:"~.(copula)~"; N"==.(NN)~";"~theta==.(theta)~
			";"~alpha[1]==.(alpha1)~";"~lambda[1]==.(lambda1)~";"~kappa[1]==.(kappa1)),
			outer = TRUE,cex = 1.2)
		par(mfrow = c(1,1),mar = c(5,4,4,2) + 0.1,
			cex.lab = 1,oma = rep(0,4))
	}
	
	fin
}

run_analysis = function(DATA,theta,upKAPPA,
	gTHRES,copula,param_grid,vec_time,show,ncores = 1){
	
	if(FALSE){
		DATA = DATA; theta = theta; upKAPPA = upKAPPA
		copula = copula; param_grid = param_grid
		vec_time = vec_time; show = show
		
		
		# Test an issue
		DATA = rd
		theta = 4/3
		upKAPPA = 0
		copula = "Gumbel"
		param_grid = seq(-5,5,0.25)
		vec_time = vec_time
		show = TRUE
		
		ncores = 1
		
	}
	
	req_names = c("time","delta","dens_t2","surv_t2")
	if( !all(req_names %in% names(DATA)) ){
		miss_names = req_names[!(req_names %in% names(DATA))]
		stop(sprintf("DATA colnames missing: %s",
			paste(miss_names,collapse = ", ")))
	}
	
	REP = list(DATA = DATA,PARAMS = smart_df(copula = copula))
	opt_out = opt_replicate(REP = REP,param_grid = param_grid,
		theta = theta,upKAPPA = upKAPPA,gTHRES = gTHRES,
		show = show,ncores = ncores)
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
	
	list(upKAPPA = upKAPPA,copula = copula,RES = opt_out,PRED = pred)
}

#' @title run_analyses
#' @description This function performs a full analysis of 
#'	an inputted dataframe. The user may specify one of two 
#'	copulas, a \code{theta} value, a parametric grid to 
#'	search over, and a vector of times for predicting survival.
#' @param DATA A data.frame containing column names \code{time},
#'	\code{delta}, \code{dens_t2}, and \code{surv_t2}.
#' @param THETAs A vector of theta values to explore and optimize over.
#' @inheritParams sim_replicate
#' @inheritParams opt_replicate
#' @param vec_time Vector of times to calculate predicted survival
#'	on the same scale as times provided in the \code{DATA} data.frame.
#' @export
run_analyses = function(DATA,THETAs = NULL,upKAPPA,
	gTHRES = 8e-2,copula,param_grid = seq(-3,3,0.25),
	vec_time,ncores = 1,show){
	
	if(FALSE){
		DATA = rd
		THETAs = 0; upKAPPA = 1
		copula = my_copula; param_grid = seq(-2,5,0.4)
		vec_time = round(seq(0,max(c(100,max(DATA$time))),
			length.out = 50),2)
		show = TRUE
		
	}
	
	if( copula == "Clayton" ){
		if( !is.null(THETAs) && any(THETAs < 0) ) stop("THETAs should be >= 0")
	} else if( copula == "Gumbel" ){
		if( !is.null(THETAs) && any(THETAs < 1) ) stop("THETAs should be >= 1")
	} else {
		stop("not a coded copula!")
	}
	
	out = list()
	
	if( is.null(THETAs) ){
		out[["theta = MLE"]] = run_analysis(DATA = DATA,
			theta = THETAs,upKAPPA = upKAPPA,
			copula = copula,param_grid = param_grid,
			vec_time = vec_time,show = show,ncores = ncores)
		return(out)
	}
	
	THETAs = sort(THETAs)
	for(theta in THETAs){
		# theta = THETAs[1]; theta
		if( show ) cat(sprintf("\n%s: theta = %s ...\n",date(),theta))
		tmp_name = sprintf("theta = %s",round(theta,4))
		tmp_out = run_analysis(DATA = DATA,
			theta = theta,upKAPPA = upKAPPA,
			gTHRES = gTHRES,copula = copula,
			param_grid = param_grid,
			vec_time = vec_time,show = show,ncores = ncores)
		out[[tmp_name]] = tmp_out
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
#' @param LAB_THETA Boolean set to \code{TRUE} to replace
#'	theta values with character strings.
#' @param hide_LL Boolean set to \code{TRUE} to suppress
#'	plot title containing evaluated log likelihoods.
#' @export
plot_SURVs = function(run_ANA,MULTIPLE,ALPHA = 0.5,
	LAB_THETA = FALSE,hide_LL = TRUE){
	
	if(FALSE){
		run_ANA = RES[[nm]]
		MULTIPLE = TRUE
		ALPHA = 0.5
		LAB_THETA = TRUE
		
	}
	
	# Remove results without predictions
	run_ANA_2 = list()
	for(vv in names(run_ANA)){
		if( is.null(run_ANA[[vv]]$PRED) ) next
		run_ANA_2[[vv]] = run_ANA[[vv]]
	}
	run_ANA = run_ANA_2; rm(run_ANA_2)
	
	if( is.null(run_ANA) || length(run_ANA) == 0 ){
		stop("No results!")
	}
	
	my_copula = run_ANA[[1]]$copula
	my_upKAPPA = ifelse(run_ANA[[1]]$upKAPPA == 1,
		"Exponentiated Weibull","Weibull")
	my_title = sprintf("%s copula on %s distribution",
		my_copula,my_upKAPPA); my_title
	
	# Aggregate predictions
	res = c()
	for(vv in names(run_ANA)){
		# vv = names(run_ANA)[1]
		tmp_df = smart_df(theta = as.numeric(gsub("theta = ","",vv)),
			run_ANA[[vv]]$PRED[,c("time","surv","low_surv2","high_surv2")])
			# run_ANA[[vv]]$PRED[,c("time","surv","low_surv","high_surv")])			)
		
		names(tmp_df) = c("theta","time","surv","low_surv","high_surv")
		res = rbind(res,tmp_df)
	}
	# dim(res); res[1:5,]
	
	lev_theta = sort(unique(res$theta))
	# lev_theta = sprintf("theta = %s",
		# sort(as.numeric(gsub("theta = ","",lev_theta))))
	lev_theta
	
	# Get max log likelihood
	max_LL = c()
	for(vv in names(run_ANA)){
		max_LL = c(max_LL,run_ANA[[vv]]$RES$LL)
	}
	max_LL
	
	if( MULTIPLE ){
		if( LAB_THETA && length(lev_theta) == 4 ){
			text_theta = c("Independence","Mild","Strong","Very Strong")
			new_label = sprintf("%s; LL = %s",
				text_theta,round(max_LL,2))
			if( hide_LL ){
				new_label = sprintf("%s",text_theta)
			}
		} else {
			new_label = sprintf("\u03B8 = %s; LL = %s",
				as.numeric(gsub("theta = ","",lev_theta)),
				round(max_LL,2))
		}
		# new_label
	} else {
		# new_label = gsub("theta","\u03B8",lev_theta)
		new_label = sprintf("\u03B8 = %s",lev_theta)
	}
	new_label
	
	res$theta2 = factor(res$theta,levels = lev_theta,
		labels = new_label)
	# dim(res); res[1:5,]
	
	time = surv = theta2 = low_surv = high_surv = NULL
	
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
			geom_line(size = 1,aes(color = theta2)) + 
			geom_ribbon(mapping = aes(ymin = low_surv,
				ymax = high_surv,fill = theta2),alpha = ALPHA) +
			ylim(c(0,1)) + xlab("Time") + ylab("Survival Probability") +
			ggtitle(my_title) + facet_wrap(~ theta2) + 
			my_themes
	} else {
		my_themes = theme(text = element_text(size = 28),
			legend.position = c("none","bottom")[2],
			plot.title = element_text(hjust = 0.5),
			panel.background = element_blank(),
			panel.grid.major = element_line(colour = "grey50",
				size = 0.5,linetype = "dotted"),
			panel.border = element_rect(colour = "black",
				fill = NA,size = 1))
		
		gg = ggplot(data = res,
			mapping = aes(x = time,y = surv,group = theta2,
				fill = theta2)) +
			geom_line(size = 1.25,alpha = 1,
				aes(color = theta2),show.legend = FALSE) +
			geom_ribbon(mapping = aes(ymin = low_surv,
				ymax = high_surv),alpha = ALPHA) +
			ylim(c(0,1)) + xlab("Time") + ylab("Survival Probability") +
			ggtitle(my_title) + labs(fill = "") + my_themes
	}
	
	return(gg)
}

#' @title full_ana_opt
#' @description A comprehensive optimization over copulas
#'	and predicting survival probabilities.
#' @inheritParams run_analyses
#' @param max_year An integer number of years of survival 
#'	probability to calculate.
#' @export
full_ana_opt = function(DATA,max_year = 100,
	param_grid = seq(-3,5,0.5),show = TRUE){
	
	my_copulas	= c("Clayton","Gumbel")
	upKAPPAs		= c(0,1) # 0 = weibull, 1 = exp-weibull
	vec_time 		= seq(0,max_year)
	
	res = list()
	for(my_copula in my_copulas){
	for(upKAPPA in upKAPPAs){
		# my_copula = my_copulas[2]; upKAPPA = upKAPPAs[2]
		
		tmp_name = sprintf("copula = %s; upKAPPA = %s",
			my_copula,upKAPPA)
		
		if( !is.null(res[[tmp_name]]) ) next
		
		cat(sprintf("%s: %s ...\n",date(),tmp_name))
		if( my_copula == "Clayton" )
			THETAs = c(0,2/3,2,6)
		if( my_copula == "Gumbel" )
			THETAs = c(1,4/3,2,4)
		
		THETAs = round(THETAs,2)
		
		run_ana = run_analyses(DATA = DATA,
			THETAs = THETAs,upKAPPA = upKAPPA,
			gTHRES = 0.5,
			copula = my_copula,param_grid = param_grid,
			vec_time = vec_time,
			show = show)
		
		res[[tmp_name]] = run_ana
		rm(run_ana)
		
	}}
	
	return(res)
	
}

#' @title full_ana_est
#' @description A comprehensive summary of parameter estimates
#' 	in data.frame and outputted table format.
#' @param RES The output object from \code{full_ana_opt()}
#' @param COHORT A string to specify which dataset or cohort
#' 	is being analyzed.
#' @export
full_ana_est = function(RES,COHORT){
	
	## Provide Latex style estimates
	est = c()
	for(vv1 in names(RES)){
		# vv1 = names(RES)[1]; vv1
		THETAs = names(RES[[vv1]]); THETAs
		tmp_str = strsplit(vv1,"; ")[[1]]; tmp_str
		my_copula = gsub("copula = ","",tmp_str[1]); my_copula
		upKAPPA = as.numeric(gsub("upKAPPA = ","",tmp_str[2])); upKAPPA
	for(vv2 in THETAs){
		# vv2 = THETAs[1]; vv2
		# str(RES[[vv1]][[vv2]])
		tmp_est = RES[[vv1]][[vv2]]$RES$cout[,
			c("PARS","EST","lowCI_2","highCI_2")]
		if( is.null(tmp_est) ) next
		
		tmp_est[,-1] = round(tmp_est[,-1],2)
		if( upKAPPA == 0 ){
			tmp_est = tmp_est[!grepl("kappa1",tmp_est$PARS),]
		}
		idx_theta = which(tmp_est$PARS == "theta")
		tmp_est = tmp_est[-idx_theta,]
		# tmp_est
		tmp_est = smart_df(Copula = my_copula,
			Dist = ifelse(upKAPPA == 0,"Weibull","Exp. Weibull"),
			THETA = gsub("theta = ","",vv2),tmp_est)
		est = rbind(est,tmp_est); rm(tmp_est)
	}
	}
	est$CI = sprintf("(%s, %s)",est$lowCI_2,est$highCI_2)
	est = est[,c("Copula","Dist","THETA","PARS","EST","CI")]
	names(est) = c("Copula","Dist","Theta",
		"Parameter","Estimate","95% CI")
	est$Parameter[est$Parameter == "alpha1"] = "$\\alpha_1$"
	est$Parameter[est$Parameter == "lambda1"] = "$\\lambda_1$"
	est$Parameter[est$Parameter == "kappa1"] = "$\\kappa_1$"
	est$Theta = smart_digits(as.numeric(est$Theta),2)

	print_latex_table(DATA = est,my_align = "llccrr",
		repeat_VARS = c("Copula","Dist","Theta"),
		add_table = TRUE,
		caption = sprintf("Real data analysis results for %s.",COHORT))

	return(est)

}

#' @title full_ana_surv
#' @description Summarizes survival probabilities
#' 	in a Latex output and data.frame object
#' @inheritParams full_ana_est
#' @param myYEARS Defaults to \code{NULL} to list
#'	survivals per year. Otherwise supply an integer 
#'	vector of years of survival probabilies to output
#' @export
full_ana_surv = function(RES,myYEARS = NULL){
	if(FALSE){
		RES = res
		myYEARS = NULL
		
	}
	
	stab = c()
	all_nms = names(RES)
	vec_time = RES[[1]][[1]]$PRED$time
	
	for(tmp_name in all_nms){
		
		tmp_list 	= RES[[tmp_name]]
		upKAPPA 	= ifelse(grepl("upKAPPA = 0",tmp_name),0,1)
		my_copula = ifelse(grepl("Clayton",tmp_name),"Clayton","Gumbel")
		
		all_thetas = names(tmp_list)
		all_thetas
		cnt = 1
		sub_df = c()
		
	for(theta in all_thetas){
		# theta = all_thetas[1]; theta
		
		Dist = ifelse(upKAPPA == 0,"Weibull","Exp. Weibull")
		theta2 = strsplit(theta," = ")[[1]][2]
		if( my_copula == "Clayton" ){
			theta2 = ifelse(theta2 == 0,0,
				ifelse(theta2 == 0.67,0.25,
				ifelse(theta2 == 2,0.5,0.75)))
			theta2
		} else {
			theta2 = ifelse(theta2 == 1,0,
				ifelse(theta2 == 1.33,0.25,
				ifelse(theta2 == 2,0.5,0.75)))
		}
		
		Years = smart_digits(vec_time,2)
		chk_obj = is.null(tmp_list[[theta]]$PRED); chk_obj
		
		tmp_df = smart_df(Copula = my_copula,Dist = Dist,
			Years = Years)
		
		if( chk_obj ){
			tmp_df$OUT = ""
		} else {
			tmp_tab = tmp_list[[theta]]$PRED[,c("time","surv","SE")]
			tmp_tab$OUT = sprintf("%s (%s)",
				smart_digits(tmp_tab$surv,4),
				smart_digits(tmp_tab$SE,4))
			tmp_df$OUT = tmp_tab$OUT
		}
		
		tmp_df = name_change(tmp_df,"OUT",
			sprintf("theta:%s%%",round(theta2*100)))
		
		if( cnt == 1 ){
			sub_df = tmp_df
		} else {
			sub_df = cbind(sub_df,tmp_df[,4,drop = FALSE])
		}
		
		cnt = cnt + 1
	}
		
		stab = rbind(stab,sub_df)
	}
	
	# Subset specific years
	if( is.null(myYEARS) ){
		stab2 = stab
	} else {
		
		if( any(myYEARS != round(myYEARS)) ) 
			stop("myYEARS should be integer valued")
		
		stab2 = stab[which(as.numeric(stab$Years) %in% myYEARS),]
	}
	
	nms = names(stab2)
	names(stab2) = gsub("theta","$\\\\theta$",nms)
	stab2

	## Printing fancy latex
	print_latex_table(
		DATA = stab2,
		repeat_VARS = c("Copula","Dist"),
		add_table = TRUE,
		fontsize = "small",
		my_align = "llrcccc",
		caption = "Survival Results")
	
	return(stab2)
	
}

#' @title full_ana_plots
#' @description Visual summary of grid search, 
#'	survival curves.
#' @param curr_dir The directory where data analyses 
#'	are outputted.
#' @inheritParams full_ana_est
#' @param NAME Character string specifying a subdirectory
#'	of \code{curr_dir} to place plots for a specific 
#'	working dataset.
#' @export
full_ana_plots = function(curr_dir,RES,NAME){
	if(FALSE){
		RES = res1
		NAME = NAME
		
	}
	
	NAME_dir = file.path(curr_dir,NAME)
	smart_mkdir(curr_dir)
	smart_mkdir(NAME_dir)
	
	# Provide manuscript survival plots
	for(nm in names(RES)){
		# nm = names(RES)[1]; nm
		
		tmp_copula = gsub("copula = (.*); (.*)","\\1",nm); tmp_copula
		tmp_dist = gsub("(.*); upKAPPA = (.*)","\\2",nm)
		tmp_dist = ifelse(tmp_dist == 0,"Weibull","Exp-Weibull"); tmp_dist
		
		gg = plot_SURVs(run_ANA = RES[[nm]],
			MULTIPLE = TRUE,ALPHA = 0.5,
			LAB_THETA = TRUE)
		
		png_fn = file.path(NAME_dir,sprintf("%s_%s.png",
			tmp_copula,tmp_dist))
		
		ggsave(filename = png_fn,plot = gg,
			device = "png",width = 12,height = 8,
			units = "in")
		rm(gg)
		
	}

	# Provide manuscript heatmaps
	pdf_fn = file.path(NAME_dir,"heatmaps.pdf")
	pdf(file = pdf_fn,height = 9,width = 12)
	
	for(nm1 in names(RES)){
		# nm1 = names(RES)[1]; nm1
	for(nm2 in names(RES[[nm1]])){
		# nm2 = names(RES[[nm1]])[1]; nm2
		
		tmp_copula = gsub("copula = (.*); (.*)","\\1",nm1); tmp_copula
		tmp_dist = gsub("(.*); upKAPPA = (.*)","\\2",nm1)
		tmp_dist = ifelse(tmp_dist == 0,"Weibull","Exp-Weibull"); tmp_dist
		tmp_theta = nm2
		
		TITLE = sprintf("%s; %s; %s",tmp_copula,tmp_dist,tmp_theta)
		TITLE
		
		gg = plot_LL(GPROF = RES[[nm1]][[nm2]]$RES$GPROF,
			GOPT = RES[[nm1]][[nm2]]$RES$GOPT_PRE)
		print(gg)
		rm(gg)
		
	}}
	dev.off()
	
	return(NULL)
	
}




###

