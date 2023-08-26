
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

sim_DATA = function(copula,dist1,n_obs,theta,
	alpha1,lambda1,kappa1,alpha2,lambda2){
	
	if(FALSE){
		copula = c("Clayton","Gumbel")[1]
		n_obs = 1e3; theta = 3/2
		alpha1 = 5; lambda1 = 1.2; kappa1 = 1
		alpha2 = 2; lambda2 = 2.1
		
	}
	
	if( copula == "Independent" ){
		tmp_copula = claytonCopula(param = 0,
			dim = 2,use.indepC = c("message","TRUE","FALSE"))
	} else if( copula == "Clayton" ){
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
	
	return(all_data)
	
}
get_PROFILE = function(GRID,PLOT = FALSE){
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
		# PAR = PARS[1]; PAR
		
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
		par(mfrow = c(2,2),mar = c(4.5,4,2,2))
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
		par(mfrow = c(1,1),mar = c(5,4,4,2)+0.1)
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
	EPS = min(abs(diff(EPS))) * 1.1
	EPS
	
	for(ii in seq(nn)){
		# ii = 84
		# GPROF[ii,]
		if( verb ) smart_progress(ii = ii,nn = nn,iter = 5,iter2 = 2e2)
		curr_PARS = as.numeric(GPROF[ii,c("log_alpha1","log_lambda1")])
		# curr_PARS
		curr_LL = GPROF$LL[ii]; # curr_LL
		while(TRUE){
			# Get neighborhood around curr_PARS
			curr_PARS; curr_LL
			
			if(FALSE){
				# Check what set of parameters are 'close' to curr_PARS
				curr_kappa = GPROF$unc_kappa1[which(
					GPROF$log_alpha1 == curr_PARS[1]
					& GPROF$log_lambda1 == curr_PARS[2]
					)]; curr_kappa
				curr_theta = GPROF$log_theta[which(
					GPROF$log_alpha1 == curr_PARS[1]
					& GPROF$log_lambda1 == curr_PARS[2]
					)]; curr_theta
				
				GPROF[which(abs(GPROF$log_alpha1 - curr_PARS[1]) < EPS
					& abs(GPROF$log_lambda1 - curr_PARS[2]) < EPS
					& abs(GPROF$unc_kappa1 - curr_kappa) < EPS
					& abs(GPROF$log_theta - curr_theta) < EPS
					),]
				
			}
			
			# stop("Possible new idea for neighborhood local grid")
			
			idx = which(abs(GPROF$log_alpha1 - curr_PARS[1]) < EPS
				& abs(GPROF$log_lambda1 - curr_PARS[2]) < EPS)
			length(idx)
			GPROF[idx,]
			
			# Move to next optimal point
			tmp_df = GPROF[idx[which.max(GPROF$LL[idx])],]
			tmp_df
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
test_LOCAL_OPTS = function(GRID,PLOT = FALSE){
	if(FALSE){
		GRID = smart_df(gout$DAT)
		PLOT = TRUE
		
	}
	
	GRID = smart_df(GRID)
	if( any(is.na(GRID$LL)) ) 
		stop("NA/NaN for LL in GRID")
	
	PARS = names(GRID)[1:4]; PARS
	
	if( PLOT ){
		par(mfrow = c(2,2),mar = c(4,4,0.5,0.5))
	}
	
	out = c()
	for(PAR in PARS){
		# PAR = PARS[1]; PAR
		
		xx = sort(unique(GRID[[PAR]])); xx
		yy = sapply(xx,function(zz){
			max(GRID$LL[which(GRID[[PAR]] == zz)])
		})
		EPS = min(diff(xx)); EPS
		
		nvals = length(xx)
		
		# If no grid for this par, skip
		if( nvals == 1 ) next
		
		dat = smart_df(xx,yy,row = seq(nvals),max_idx = NA); # dat
		out_par = c()
		
		for(ii in seq(nvals)){
			# ii = 1
			
			if( !is.na(dat$max_idx[ii]) ) next
			
			curr_xx = dat$xx[ii]; curr_xx
			curr_LL = dat$yy[ii]; curr_LL
			idx_path = ii
			
			while(TRUE){
				idx = which( abs(dat$xx - curr_xx) <= EPS 
					& dat$yy >= curr_LL ); idx
				
				tmp_df = dat[idx,,drop = FALSE]
				tmp_df = tmp_df[order(-tmp_df$yy),,drop = FALSE]; tmp_df
				new_xx = tmp_df$xx[1]
				new_LL = tmp_df$yy[1]
				
				if( new_xx == curr_xx ){
					# stop("not ready yet")
					dat$max_idx[idx_path] = tmp_df$row[1]
					
					# Get data from GRID
					out_par = rbind(out_par,GRID[which(GRID[[PAR]] == new_xx
						& GRID$LL == new_LL),,drop = FALSE])
					
					break
				}
				
				curr_xx = new_xx
				curr_LL = new_LL
				idx_path = c(idx_path,tmp_df$row[1])
				
			}
			
		}
		
		out_par = unique(out_par)
		rownames(out_par) = NULL
		out_par
		
		out = rbind(out,out_par)
		out = unique(out)
		rownames(out) = NULL
		
		if( PLOT ){
			plot(dat[,c("xx","yy")],xlab = PAR,
				ylab = "Profile LL",type = "b",pch = 16)
			abline(v = out_par[[PAR]],lty = 2,lwd = 2,col = "red")
		}
		
	}
	
	if( PLOT ){
		par(mfrow = c(1,1),mar = c(5,4,4,2) + 0.1)
	}
	
	out
	return(out)
	
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
#' @param copula A string input, either "Independent", "Clayton" or "Gumbel"
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
	
	copula = get_copula(copula = copula,PO = TRUE)
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
	
	out = list(DATA = dat,PARAMS = params)
	
	return(out)
	
}


###

