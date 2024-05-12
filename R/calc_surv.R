#' @title calc_mDENS_mSURV
#' @description Once the optimization has finished and estimates
#'	obtained, the user may calculate net survival and net
#'	densities, strata-specific survivals and densities, and 
#'	observed survivals and observed densities for a set of times
#' @param vec_time A numeric vector of times in years to calculate 
#'	survival and density values for. 
#' @param run_ANA The list object outputted from \code{run_analyses()}
#' @inheritParams refData_match
#' @param SEXs A character vector of sexes to calculate quantities for.
#'	Values should be \code{'male'} or \code{'female'}. Note that input
#'	values will be subsetted to those included in \code{rDAT}
#' @param AGEs An integer vector of ages in years. Note that input 
#'	values will be subsetted to those included in \code{rDAT}.
#' @param YEARs An integer vector of calendar years diagnosed. 
#'	Note that input values will be subsetted to those included in 
#'	\code{rDAT}.
#' @export
calc_mDENS_mSURV = function(vec_time,run_ANA,rDAT,SEXs,AGEs,YEARs){
	
	nms = names(rDAT)
	names(rDAT) = tolower(nms)
	
	time = NULL
	tDAT = smart_df(time = vec_time)
	sDAT = smart_df(sex = SEXs)
	aDAT = smart_df(age = AGEs)
	yDAT = smart_df(datediag_yr = YEARs)
	
	# Perform cross join
	wDAT = sqldf("
	select 
		tt.*,ss.*,aa.*,yy.*
	from 
		tDAT as tt, 
		sDAT as ss, 
		aDAT as aa, 
		yDAT as yy
	")
	wDAT = wDAT[which(wDAT$sex %in% unique(rDAT$sex)
		& wDAT$age %in% unique(rDAT$age)
		& wDAT$datediag_yr %in% unique(rDAT$year)
		),]
	# dim(wDAT); head(wDAT)
	#wDAT$age = as.integer(wDAT$age)
	#wDAT$datediag_yr = as.integer(wDAT$datediag_yr)
	# str(wDAT)
	# str(rDAT)
	
	# Rcpp version of matching
	MET = refData_match(wDAT = wDAT,
		rDAT = rDAT,ncores = 1)
	
	wDAT$log_dens_t2 	= MET$log_dens_t2
	wDAT$log_cdf_t2 	= MET$log_cdf_t2
	
	fin = c()
	nsolu = length(run_ANA); nsolu
	
	for(solu in seq(nsolu)){
	
		# real data
		# solu = 1
		#run_ANA[[solu]]$RES$cout
		ePARS = run_ANA[[solu]]$RES$cout$EST
		LAM = ePARS[2]
		ALP = ePARS[1]
		KAP = ePARS[3]
		THETA = ePARS[4]
		tmp_COPULA = run_ANA[[solu]]$copula
		tmp_DIST = ifelse(KAP == 1,"Weibull","ExpWeibull")
		
		NN = nrow(wDAT)
		mat = smart_df(ID = seq(NN))
		mat$time = wDAT$time
		mat$CATEGORY = ""
		if( length(unique(wDAT$sex)) > 1 ) 
			mat$CATEGORY = sprintf("%sSex=%s;",mat$CATEGORY,wDAT$sex)
		if( length(unique(wDAT$age)) > 1 ) 
			mat$CATEGORY = sprintf("%sAge=%s;",mat$CATEGORY,wDAT$age)
		if( length(unique(wDAT$datediag_yr)) > 1 ) 
			mat$CATEGORY = sprintf("%sYear=%s;",mat$CATEGORY,wDAT$datediag_yr)
		mat$CONFIG = sprintf("%s,%s",tmp_COPULA,tmp_DIST)
		
		mat$log_CDF_1 = NA
		mat$log_CDF_2 = wDAT$log_cdf_t2
		mat$log_DEN_1 = NA
		mat$log_DEN_2 = wDAT$log_dens_t2
		mat$COP = NA
		mat$OFF = NA
		mat$mDENS = NA # mixture density
		mat$mSURV = NA # mixture survival
		mat$POST_1 = NA
		
		for(ID in seq(NN)){
			# ID = 1
			
			tmp_vals = calc_expweibull_logCDF_logPDF(
				XX = wDAT$time[ID],LAM = LAM,ALP = ALP,
				KAP = KAP)
			tmp_vals = c(tmp_vals); tmp_vals
			
			mat$log_CDF_1[ID] = tmp_vals[1]
			mat$log_DEN_1[ID] = tmp_vals[2]
			
			log_DENs = c(mat$log_DEN_1[ID],mat$log_DEN_2[ID])
			log_CDFs = c(mat$log_CDF_1[ID],mat$log_CDF_2[ID])
			
			out = calc_copula_CDF_OFF(log_DENs = log_DENs,
				log_CDFs = log_CDFs,copula = tmp_COPULA,
				THETA = THETA)
			out = c(out); out
			
			mat$COP[ID] = out[1]
			mat$OFF[ID] = out[2]
			
			DENs = exp(log_DENs); DENs
			mat$mDENS[ID] = DENs[1] + DENs[2] - mat$OFF[ID]
			
			CDFs = exp(log_CDFs); CDFs
			mat$mSURV[ID] = 1 - CDFs[1] - CDFs[2] + mat$COP[ID]
			
		}
		
		mat$DEN_1 = exp(mat$log_DEN_1)
		mat$DEN_2 = exp(mat$log_DEN_2)
		
		# warning("Label SURV_1 as disease-specific survival")
		mat$SURV_1 = 1 - exp(mat$log_CDF_1)
		mat$SURV_2 = 1 - exp(mat$log_CDF_2)
		dim(mat); # head(mat)
		
		fin = rbind(fin,mat)
		rm(mat)
		
	}
	# dim(fin); head(fin)
	
	# Convert fin from horizontal to vertical
	
	aa = sqldf("
	select distinct
		time,
		CATEGORY,
		'N/A' as CONFIG,
		'Individual' as Type,
		SURV_2 as Survival,
		DEN_2 as Density
	from
		fin
	where
		CONFIG = (select CONFIG from fin limit 1)
	
	union
	
	select
		time,
		'Net' as CATEGORY,
		CONFIG,
		'Individual' as Type,
		SURV_1 as Survival,
		DEN_1 as Density
	from
		fin
	
	union
	
	select
		time,
		CATEGORY,
		CONFIG,
		'Observed' as Type,
		mSURV as Survival,
		mDENS as Density
	from
		fin
	")
	aa$Type = factor(aa$Type,
		levels = c("Individual","Observed"))
	# dim(aa); head(aa)
	
	return(aa)
	
}

visual_surv_dens = function(SURV_DENS){
	#' @title visual_surv_dens
	#' @description This function
	#' @export
	
	# stop("Edit the labels CATEGORY and CONFIG, refresh the vignette and paper images")
	# CATEGORY to Strata
	# CONFIG to Model
	# Individual to Strata
	
	time = NULL
	Survival = NULL
	CATEGORY = NULL
	CONFIG = NULL
	Density = NULL
	
	aa = SURV_DENS
	
	# ggplot settings
	my_themes = theme(text = element_text(size = 28),
		axis.text.x = element_text(size = 20),
		legend.position = c("none","bottom","right")[3],
		legend.text = element_text(size = 10),
		legend.key.width = unit(2,"cm"),
		panel.spacing = unit(1,"lines"),
		panel.background = element_blank(),
		panel.grid.major = element_line(colour = "grey50",size = 0.5,linetype = "dotted"),
		plot.title = element_text(hjust = 0.5),
		panel.border = element_rect(colour = "black",fill = NA,size = 1),
		# strip.text.x = element_text(size = 12),
		strip.text = element_text(size = 25))
	
	# Survivals
	gg1 = ggplot(data = aa,mapping = aes(x = time,
			y = Survival,fill = CATEGORY)) +
		geom_line(linewidth = 1,aes(color = CATEGORY,linetype = CONFIG)) +
		facet_wrap(~Type,ncol = 1) +
		guides(color = guide_legend(override.aes = list(linewidth = 5)),
			linetype = guide_legend(override.aes = list(linewidth = 2))) +
		xlab("Time since Diagnosis (yrs)") +
		my_themes
	
	# Densities
	gg2 = ggplot(data = aa,mapping = aes(x = time,
			y = Density,fill = CATEGORY)) +
		geom_line(linewidth = 1,aes(color = CATEGORY,linetype = CONFIG)) +
		facet_wrap(~Type,ncol = 1) +
		guides(color = guide_legend(override.aes = list(linewidth = 5)),
			linetype = guide_legend(override.aes = list(linewidth = 2))) +
		xlab("Time since Diagnosis (yrs)") + 
		my_themes
	
	out = list(SURV = gg1,DENS = gg2)
	return(out)
	
}

##