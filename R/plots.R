# ----------
# Visualization functions
# ----------

#' @title plot_SURVs
#' @description Plot net survival probabilities per model
#' @param run_ANA The object outputted from \code{run_analyses}
#'	function.
#' @param MULTIPLE A boolean set to \code{TRUE} to display
#'	a survival plot per \code{theta} with accompanying maximum 
#'	log likelihood. Otherwise, the survival plots will be 
#'	overlayed with a legend.
#' @param ncol Integer number of columns of plots to display
#' @param ALPHA A numeric value between 0 and 1 to control the
#'	confidence band transparency.
#' @export
plot_SURVs = function(run_ANA,MULTIPLE,ncol = 1,ALPHA = 0.5){
	
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
		BIC			= run_ANA[[vv]]$RES$BIC
		
		tmp_df = smart_df(COPULA = COPULA,ALPHA = tmpALPHA,LAMBDA = tmpLAMBDA,
			THETA = THETA,KAPPA = KAPPA,DIST = DIST,LL = LL,BIC = BIC,
			run_ANA[[vv]]$PRED[,c("time","surv","low_surv2","high_surv2")])
		
		names(tmp_df) = c("COPULA","ALPHA","LAMBDA",
			"THETA","KAPPA","DIST","LL","BIC",
			"time","surv","low_surv","high_surv")
		head(tmp_df)
		
		res = rbind(res,tmp_df)
	}
	res$LABEL = sprintf("%s + %s, \u03B1 = %s, \u03BB = %s, \n\u03B8 = %s, \u03BA = %s, LL = %s, BIC = %s",
			res$COPULA,res$DIST,round(res$ALPHA,2),round(res$LAMBDA,2),
			round(res$THETA,2),round(res$KAPPA,2),round(res$LL,2),round(res$BIC,2))
	res = res[order(res$COPULA,res$DIST,res$THETA,res$KAPPA),]
	uLABEL = unique(res$LABEL)
	res$LABEL = factor(res$LABEL,levels = uLABEL)
	
	LABEL2 = NULL
	res$LABEL2 = sprintf("%s + %s",res$COPULA,res$DIST)
	dim(res); res[1:5,]
	
	time = surv = low_surv = high_surv = LABEL = NULL
	
	if( MULTIPLE ){
		my_themes = theme(text = element_text(size = 28),
			axis.text.x = element_text(size = 20),
			legend.position = c("none","bottom")[1],
			panel.spacing = unit(1,"lines"),
			panel.background = element_blank(),
			panel.grid.major = element_line(colour = "grey50",
				linewidth = 0.5,linetype = "dotted"),
			plot.title = element_text(hjust = 0.5),
			panel.border = element_rect(colour = "black",
				fill = NA,linewidth = 1),
			strip.text.x = element_text(size = 12))
		
		gg = ggplot(data = res,mapping = aes(x = time,y = surv)) +
			geom_line(linewidth = 1,aes(color = LABEL)) + 
			geom_ribbon(mapping = aes(ymin = low_surv,
				ymax = high_surv,fill = LABEL),alpha = ALPHA) +
			ylim(c(0,1)) + xlab("Time since Diagnosis (yrs)") + 
			ylab("Net Survival") +
			# ggtitle(my_title) + 
			facet_wrap(~ LABEL,ncol = ncol) + 
			my_themes
		# gg
		
	} else {
		my_themes = theme(text = element_text(size = 28),
			legend.position = c("none","bottom","right")[2],
			plot.title = element_text(hjust = 0.5),
			panel.background = element_blank(),
			panel.grid.major = element_line(colour = "grey50",
				linewidth = 0.5,linetype = "dotted"),
			panel.border = element_rect(colour = "black",
				fill = NA,linewidth = 1),
			legend.key.width = unit(0.85, "inch"),
			legend.key.size = unit(0.5, "inch"),
			legend.text = element_text(size = 20))
		
		gg = ggplot(data = res,
			mapping = aes(x = time,y = surv,group = LABEL2,
				fill = LABEL2)) +
			geom_line(size = 1.25,alpha = 1,
				aes(color = LABEL2),show.legend = FALSE) +
			geom_ribbon(mapping = aes(ymin = low_surv,
				ymax = high_surv),alpha = ALPHA) +
			ylim(c(0,1)) + xlab("Time since Diagnosis (yrs)") + 
			ylab("Net Survival") +
			# ggtitle(my_title) + 
			labs(fill = "") + my_themes
	}
	
	return(gg)
}

###
