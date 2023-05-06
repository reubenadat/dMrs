# ----------
# Parameters
# ----------

get_copula = function(copula){
	
	if( missing(copula) ){
		copula = make_menu(
			PROMPT = "Select a copula",
			OPTS = c("Clayton","Gumbel"))
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

### EOF

