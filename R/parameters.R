# ----------
# Parameters
# ----------

get_copula = function(copula){
	
	if( missing(copula) ){
		copula = make_menu(
			PROMPT = "Select an option",
			OPTS = c("Clayton","Gumbel"))
	}
	
	if( is.null(copula) ){
		return(c("Clayton","Gumbel"))
	}
	
	if( !class(copula) %in% "character" ){
		return(get_copula())
	}
	
	copula = unique(copula)
	if( length(copula) != 1 ){
		return(get_copula())
	}
	
	if( !copula %in% c("Clayton","Gumbel") ){
		return(get_copula())
	}
	
	return(copula)
	
}
get_verb = function(verb){
	
	if( missing(verb) ){
		verb = make_menu(OPTS = c(TRUE,FALSE))
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
	upPARS[4] = ifelse(is.null(THETA),1,0)
	
	return(upPARS)
	
}
get_upKAPPA = function(upKAPPA){
	
	if( missing(upKAPPA) ){
		DIST = make_menu(PROMPT = "Select a distribution",
			OPTS = c("Weibull","Exp-Weibull"))
		upKAPPA = ifelse(DIST == "Weibull",0,1)
		return(upKAPPA)
	}
	
	if( is.null(upKAPPA) ){
		upKAPPA = c(0,1)
	}
	
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
		return(get_parGRID(LOW = LOW,HIGH = HIGH,STEP = STEP))
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

### EOF

