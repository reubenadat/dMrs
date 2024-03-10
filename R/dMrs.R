
# ----------
# Match functions
# ----------

#' @title refData_match
#' @description This function takes as input
#'	a working dataset of interest and a 
#'	reference dataset.
#' @param wDAT A working dataset data.frame containing 
#'	required columns \code{age}, \code{time}, \code{datediag_yr}, 
#'	and \code{sex} corresponding to age, observed time (in years), 
#'	diagnosis year, and sex (coded 'female'/'male'), respectively.
#' @param rDAT A reference dataset data.frame containing
#'	required columns \code{year}, \code{age}, \code{qx},
#'	and \code{sex} corresponding to reference year, 
#'	reference age, event's interval hazard, and sex, respectively.
#' @param ncores Integer number of parallel threads to decrease
#'	matching runtime.
#' @return A dataframe containing calculated log-transformed density and 
#'	log-transformed cumulative distribution
#' @export
refData_match = function(wDAT,rDAT,ncores = 1){
	
	# wDAT = working dataset
	req_names1 = c("age","time","datediag_yr","sex")
	nms = tolower(names(wDAT))
	names(wDAT) = nms
	if( !all(req_names1 %in% nms) ){
		miss_names = req_names1[!(req_names1 %in% nms)]
		stop(sprintf("Missing colnames in wDAT: %s",
			paste(miss_names,collapse = ", ")))
	}
	
	# rDAT = reference dataset
	req_names2 = c("year","age","qx","sex")
	nms = tolower(names(rDAT))
	names(rDAT) = nms
	if( !all(req_names2 %in% nms) ){
		miss_names = req_names2[!(req_names2 %in% nms)]
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
	if( !chk_int_vec(wDAT$datediag_yr) ){
		stop("Set wDAT$datediag_yr to integers")
	}
	if( !is(wDAT$sex,"character") ){
		stop("Set wDAT$sex to character")
	}
	if( !all(wDAT$sex %in% c("female","male")) ){
		stop("wDAT$sex takes values 'female' and 'male'")
	}
	
	if( !chk_int_vec(rDAT$year) ){
		stop("Set rDAT$year to integers")
	}
	if( !chk_int_vec(rDAT$age) ){
		stop("Set rDAT$age to integers")
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
	
	# Check if all values in wDAT are in rDAT
	if( !all(unique(wDAT$sex) %in% unique(rDAT$sex)) )
		stop("Not all wDAT sexes contained in rDAT")
	if( !all(unique(wDAT$age) %in% unique(rDAT$age)) )
		stop("Not all wDAT ages contained in rDAT")
	if( !all(unique(wDAT$datediag_yr) %in% unique(rDAT$year)) )
		stop("Not all wDAT years contained in rDAT")
	
	ncores = as.integer(ncores)
	if( !is(ncores,"integer") ) stop("ncores must be integer")
	if( ncores < 1 ) stop("ncores >= 1")
	
	wDAT$sex2 = ifelse(wDAT$sex == "female",1,0)
	rDAT$sex2 = ifelse(rDAT$sex == "female",1,0)
	rDAT = rDAT[order(rDAT$year,rDAT$age,rDAT$sex),]
	
	# Do matching on age, sex, time of diagnosis, time of death
	# cat("Get matching density and survival data ...\n")
	wDAT_vars = c("age","datediag_yr","sex2","time")
	rDAT_vars = c("year","age","qx","sex2")
	OUT = dMrs_MATCH(wDAT = as.matrix(wDAT[,wDAT_vars]),
		rDAT = as.matrix(rDAT[,rDAT_vars]),
		ncores = ncores,verb = TRUE)
	OUT = smart_df(OUT)
	names(OUT) = c("log_dens_t2","log_cdf_t2")
	
	return(OUT)
	
}

###

