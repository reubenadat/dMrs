# ----------
# Visualization functions
# ----------

#' @title plot_LL
#' @description Visualize grid of profile log likelihood
#' @param GPROF An outputted data.frame from \code{get_PROFILE()}.
#' @param GOPT An outputted data.frame from \code{get_LOCAL_OPTS()}.
#' @param nBREAKs A positive integer set to 5 by default to specify
#'	the number of breaks in axes ticks.
#' @param COPULA A string input, either "Clayton" or "Gumbel"
#' @param HJUST A numeric value for horizontal alignment of labeled
#'	ggplot points
#' @export
plot_LL = function(GPROF,GOPT = NULL,nBREAKs = 5,COPULA,
	HJUST = 1){
	
	if(FALSE){
		GPROF 	= gprof
		GOPT 		= gopt
		COPULA 	= PARAMS$copula
		
		nBREAKs = 5; HJUST = -0.1
		
	}
	
	res = GPROF
	res$LL[res$LL == min(res$LL)] = NA
	low_quant = quantile(res$LL,0.2,na.rm = TRUE); low_quant
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
	
	my_stepsize = diff(sort(unique(res$log_alpha1)))[1]
	
	g1 = ggplot(data = res,mapping = aes(x = log_alpha1,
		y = log_lambda1,fill = LL)) +
		geom_tile(width = my_stepsize,height = my_stepsize) + 
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
		if( COPULA %in% c("Independent","Clayton") ){
			g1 = g1 + geom_text(data = GOPT,
				aes(label = sprintf("log(\u03B8) = %s",log_theta)),
				vjust = -1,hjust = HJUST,size = 5,color = "red")
		} else if( COPULA == "Gumbel" ){
			g1 = g1 + geom_text(data = GOPT,
				aes(label = sprintf("log(\u03B8-1) = %s",log_theta)),
				vjust = -1,hjust = HJUST,size = 5,color = "red")
		}
		
	}
	
	g1 = g1 + ggtitle(TITLE)
	# g1
	
	return(g1)
}


#' @title plot_SURVs
#' @param run_ANA The object outputed from \code{run_analyses}
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
			strip.text.x = element_text(size = 12))
		
		gg = ggplot(data = res,mapping = aes(x = time,y = surv)) +
			geom_line(size = 1,aes(color = LABEL)) + 
			geom_ribbon(mapping = aes(ymin = low_surv,
				ymax = high_surv,fill = LABEL),alpha = ALPHA) +
			ylim(c(0,1)) + xlab("Time") + ylab("Survival Probability") +
			# ggtitle(my_title) + 
			facet_wrap(~ LABEL,ncol = ncol) + 
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
		# PARS = as.numeric(uPARS[1,1:4]); PARS
		PARS = TRUTH$uPARS
		PARS
		COPULA = COPULA; COPULA
		
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
	
	F_T1_T2 = apply(smart_df(CDF_1,CDF_2),1,function(xx){
		calc_copula(F1 = xx[1],F2 = xx[2],
			copula = COPULA,THETA = THETA)
	})
	
	# Calc densities
	D1 = dexpweibull(x = DATA$time,lambda = LAMBDA1,
		alpha = ALPHA1,kappa = KAPPA1,log = FALSE)
	D2 = DATA$dens_t2
	
	out = smart_df(time = DATA$time,
		CDF_1 = CDF_1,CDF_2 = CDF_2,
		D1 = D1,D2 = D2,F_T1_T2 = F_T1_T2)
	
	out$D1_D2 = apply(out[,c("D1","D2","CDF_1","CDF_2","F_T1_T2")],
		1,function(xx){
		calc_copula_offset(D1 = xx[1],D2 = xx[2],
			F1 = xx[3],F2 = xx[4],copula = COPULA,
			THETA = THETA,F_T1_T2 = xx[5])
	})
	
	out$H1 = out$D1 / (1 - out$CDF_1)
	out$H2 = out$D2 / (1 - out$CDF_2)
	
	dim(out); head(out)
	
	# Plot
	par(mfrow = c(4,3),mar = c(4.4,4.4,2,0.5))
	hist(out$CDF_1,main = "",xlab = "F1",
		breaks = 40,xlim = c(0,1))
	
	hist(out$CDF_2,main = "",xlab = "F2",
		breaks = 40,xlim = c(0,1))
	
	plot(DATA$time,out$CDF_1,col = "red",
		xlab = "Obs Time",ylab = "CDF",pch = 1,
		ylim = c(0,1))
	points(DATA$time,out$CDF_2,col = "blue",pch = 4)
	legend("bottomright",legend = sprintf("Event %s",c(1,2)),
		col = c("red","blue"),pch = rep(16,2),pt.cex = rep(2,2))
	
	# smoothScatter(out$CDF_1,out$CDF_2,
		# xlim = c(0,1),ylim = c(0,1),
		# xlab = "F1",ylab = "F2")
	
	max_yy = quantile(c(out$H1,out$H2),0.95) * 1.5
	plot(DATA$time,out$H1,col = "red",
		xlab = "Obs Time",ylab = "Hazard",
		pch = 1,ylim = c(0,max_yy))
	points(DATA$time,out$H2,col = "blue",
		pch = 4)
	points(DATA$time,0.5*(out$H1 + out$H2),
		col = "magenta",pch = 5)
	
	y_range = range(log(c(out$H1,out$H2,out$H1 + out$H2)))
	plot(DATA$time,log(out$H1),col = "red",
		xlab = "Obs Time",ylab = "Log(Hazard)",
		pch = 1,ylim = y_range)
	points(DATA$time,log(out$H2),col = "blue",
		pch = 4)
	points(DATA$time,log(0.5*(out$H1 + out$H2)),
		col = "magenta",pch = 5)
	
	hist(out$F_T1_T2,breaks = 40,
		xlab = "Joint CDF Copula",
		xlim = c(0,1),main = "")
	
	hist(out$D1_D2,breaks = 40,
		xlab = "Offset Copula",
		main = "")
	
	hist(log(out$D1),breaks = 40,
		xlab = "log(D1)",main = "")
	
	hist(log(out$D2),breaks = 40,
		xlab = "log(D2)",main = "")
	
	if( all(c("T1","T2") %in% names(DATA)) ){
		smoothScatter(log(1 + DATA[,c("T1","T2")]),
			xlab = "log(1 + Time1)",
			ylab = "log(1 + Time2)",
			main = sprintf("Copula=%s,Theta=%s",COPULA,round(THETA,3)))
	}
	
	par(mfrow = c(1,1),mar = c(5,4,4,2) + 0.1)
	
	return(out)
	
}


###
