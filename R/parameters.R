# ----------
# Parameters
# ----------

get_PLOT = function(PLOT){
	
	if( missing(PLOT) ){
		PLOT = TRUE
	}
	
	if( !is(PLOT,"logical") )
		stop("PLOT should be logical")
	
	return(PLOT)
	
}
get_copula = function(copula,PO){
	
	if( missing(PO) )
		stop("Set PO for get_copula")
	
	if( missing(copula) ){
		if( PO ){
			copula = make_menu(
				PROMPT = "Select an option",
				OPTS = c("Independent","Clayton","Gumbel"))
		} else {
			return(c("Independent","Clayton","Gumbel"))
		}
	}
	
	if( !is(copula,"character") ){
		return(get_copula(PO = PO))
	}
	
	copula = unique(copula)
	if( length(copula) != 1 ){
		return(get_copula(PO = PO))
	}
	
	if( !copula %in% c("Independent","Clayton","Gumbel") ){
		return(get_copula(PO = PO))
	}
	
	return(copula)
	
}
get_verb = function(verb){
	
	if( missing(verb) ){
		verb = make_menu(PROMPT = "Set verbose parameter",
			OPTS = c(TRUE,FALSE))
	}
	
	if( !class(verb) %in% c("logical","character") ){
		return(get_verb())
	}
	
	verb = unique(verb)
	if( length(verb) != 1 ){
		return(get_verb())
	}
	
	verb = as.logical(verb)
	if( is.na(verb) )
		return(get_verb())
	
	return(verb)
	
}
get_upPARS = function(upKAPPA,THETA){
	
	upPARS 		= rep(1,4)
	upPARS[3] = upKAPPA
	upPARS[4] = ifelse(is.na(THETA),1,0)
	
	return(upPARS)
	
}
get_upKAPPA = function(upKAPPA){
	
	if( missing(upKAPPA) )
		return(c(0,1))
	
	if( !all(upKAPPA %in% c(0,1)) )
		return(get_upKAPPA())
	
	return(upKAPPA)
	
}
get_parGRID = function(LOW,HIGH,STEP){
	
	if( missing(LOW) ){
		LOW = make_menu(PROMPT = "Enter a minimum value for the parameter grid (e.g. -2)",OPTS = NULL)
	}
	LOW = tryCatch(as.numeric(LOW),
		warning = function(ww){NULL},
		error = function(ee){NULL})
	if( is.null(LOW) ) 
		return(get_parGRID(HIGH = HIGH,STEP = STEP))
	
	if( missing(HIGH) ){
		HIGH = make_menu(PROMPT = "Enter a maximum value for the parameter grid (e.g. 3)",OPTS = NULL)
	}
	HIGH = tryCatch(as.numeric(HIGH),
		warning = function(ww){NULL},
		error = function(ee){NULL})
	if( is.null(HIGH) ) 
		return(get_parGRID(LOW = LOW,STEP = STEP))
	
	if( missing(STEP) ){
		STEP = make_menu(PROMPT = "Enter a stepsize for the parameter grid (e.g. 0.1)",OPTS = NULL)
	}
	STEP = tryCatch(as.numeric(STEP),
		warning = function(ww){NULL},
		error = function(ee){NULL})
	if( is.null(STEP) ) 
		return(get_parGRID(LOW = LOW,HIGH = HIGH))
	
	if( length(LOW) != 1 || length(HIGH) != 1 || length(STEP) != 1 ){
		return(get_parGRID())
	}
	
	param_grid = seq(from = LOW,to = HIGH,by = STEP)
	param_grid = sort(unique(param_grid))
	return(param_grid)
	
}
get_vecTIME = function(TIME,vec_time){
	
	if( missing(vec_time) ){
		vec_time = seq(0,round(max(TIME),1),
				length.out = 100)
	}
	
	vec_time = sort(unique(vec_time))
	
	return(vec_time)
	
}
get_THETAs = function(THETAs,COPULA){
	
	if( missing(THETAs) ){
		if( COPULA == "Independent" )
			THETAs = NULL
	}
	
	if( COPULA == "Clayton" ){
		if( !is.null(THETAs) && any(THETAs <= 0) )
			THETAs = THETAs[THETAs > 0]
	} else if( COPULA == "Gumbel" ){
		if( !is.null(THETAs) && any(THETAs <= 1) )
			THETAs = THETAs[THETAs > 1]
	}
	
	return(THETAs)
	
}
calc_uPARS = function(ALPHA1,LAMBDA1,KAPPA1,THETA,DIST,COPULA){
	
	PARS = rep(NA,4)
	names(PARS) = c("log_alpha1","log_lambda1","log_kappa1","unc_theta")
	
	PARS[1] = log(ALPHA1)
	PARS[2] = log(LAMBDA1)
	PARS[3] = ifelse(DIST == "weibull",0,log(KAPPA1))
	PARS[4] = ifelse(COPULA == "Independent",-Inf,
		ifelse(COPULA == "Clayton",log(THETA),
		log(THETA-1)))
	PARS
	
}

### EOF

