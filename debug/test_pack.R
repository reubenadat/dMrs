# R script with R functions for dMrs

# ----------
# Shortcuts, Install/Load Libraries
# ----------
rm(list = ls())
setdirs = function(){
	
	where_file = gsub("\\\\","/",Sys.getenv("HOME"))
	clust_name = Sys.getenv("SLURM_CLUSTER_NAME")

	if( where_file == "C:/Users/Admin/Documents" ){
		desk_dir = "C:/Users/Admin/Desktop"
		git_dir = file.path(desk_dir,"github")
		curr_dir = file.path(desk_dir,"dMrs_sim")
	} else if( where_file == "/home/pllittle" && clust_name == "plough" ){
		git_dir = "~/github"
		curr_dir = "~/dMrs"
	} else if( where_file == "/Users/adat" ){
		git_dir = "/Users/adat/Github"
		curr_dir = "/Users/adat"
	}
	pack_dir = file.path(git_dir,"dMrs")
	
	# Package check
	req_packs = c("devtools","numDeriv","smarter",
		"ggplot2","viridis","relsurv","copula","dMrs")
	all_packs = as.character(installed.packages()[,1])
	for(pack in req_packs){
		
		if( pack %in% all_packs ){
			library(pack,character.only = TRUE)
			next
		}
		
		if( pack == "smarter" ){
			install_github(repo = "pllittle/smarter")
		} else if( pack == "dMrs" ){
			# install(pack_dir,build_vignettes = FALSE)
			Rcpp_fn = file.path(pack_dir,"src/dMrs.cpp")
			Rcpp::sourceCpp(file = Rcpp_fn,showOutput = TRUE)
			source(file.path(pack_dir,"R/smarter.R"))
			source(file.path(pack_dir,"R/parameters.R"))
			source(file.path(pack_dir,"R/plots.R"))
			source(file.path(pack_dir,"R/optimize.R"))
			source(file.path(pack_dir,"R/summary.R"))
			source(file.path(pack_dir,"R/dMrs.R"))
		} else {
			install.packages(pack)
		}
		
	}
	
	# Make directories
	smart_mkdir(curr_dir)
	
	# Output
	out = list(git_dir = git_dir,
		curr_dir = curr_dir,
		pack_dir = pack_dir)
	
	return(out)
	
}
my_dirs = setdirs()
my_dirs

# ----------
# Create package datasets
# ----------
if( TRUE ){

library(stringr)

usa_fn = file.path(my_dirs$pack_dir,
	"inst/extdata",
	"USA_2020_Males_Females")
usa = readLines(usa_fn,warn = FALSE)
usa = usa[usa != ""]
usa = stringr::str_squish(usa)
# usa = trimws(usa)
usa[1:5]

idx_male = grep("Males",usa); idx_male
idx_female = grep("Females",usa); idx_female

usa_men 	= usa[seq(idx_male + 1,idx_female - 1)]
usa_women = usa[seq(idx_female + 1,length(usa))]

usa_men_dat = t(sapply(usa_men,function(xx){
	strsplit(xx," ")[[1]]
},USE.NAMES = FALSE))
colnames(usa_men_dat) = usa_men_dat[1,]
usa_men_dat = usa_men_dat[-1,]
usa_men_dat = usa_men_dat[,c("Year","Age","qx")]
usa_men_dat = smart_df(usa_men_dat)
usa_men_dat$sex = "male"
dim(usa_men_dat); usa_men_dat[1:10,]

usa_women_dat = t(sapply(usa_women,function(xx){
	strsplit(xx," ")[[1]]
},USE.NAMES = FALSE))
colnames(usa_women_dat) = usa_women_dat[1,]
usa_women_dat = usa_women_dat[-1,]
usa_women_dat = usa_women_dat[,c("Year","Age","qx")]
usa_women_dat = smart_df(usa_women_dat)
usa_women_dat$sex = "female"
dim(usa_women_dat); usa_women_dat[1:10,]

usa = rbind(usa_men_dat,usa_women_dat)
dim(usa); usa[1:10,]

save(usa,file = file.path(my_dirs$pack_dir,"data/usa.rda"))

}

# ----------
# Test Rcpp functions
# ----------
if(FALSE){

# F1 = c(0,0.5)[1]; F1
# F1 = seq(0,1,length.out = 5e2)
ALPHA 	= 1.5
LAMBDA 	= 20
KAPPA 	= 2
TT			= seq(0,100,length.out = 5e2)
F1			= pexpweibull(q = TT,lambda = LAMBDA,alpha = ALPHA,kappa = KAPPA)
F2 = c(0,0.3)[1]; F2
COPULA = c("Independent","Clayton","Gumbel")[2]
THETA = 3

# Test derivatives of copula
LAMBDA2 = LAMBDA * 1.5
ALPHA2 	= ALPHA * 0.8
KAPPA2	= 1

ff = function(tt){
	F1 = pexpweibull(q = tt,lambda = LAMBDA,
		alpha = ALPHA,kappa = KAPPA)
	F2 = pexpweibull(q = tt,lambda = LAMBDA2,
		alpha = ALPHA2,kappa = KAPPA2)
	calc_copula(F1 = F1,F2 = F2,
		copula = COPULA,THETA = THETA)
}
tt = 5e3
ff(tt)
g1 = grad(ff,tt); g1
gg = function(tt){
	F1 = pexpweibull(q = tt,lambda = LAMBDA,
		alpha = ALPHA,kappa = KAPPA)
	D1 = dexpweibull(x = tt,lambda = LAMBDA,
		alpha = ALPHA,kappa = KAPPA)
	D2 = dexpweibull(x = tt,lambda = LAMBDA2,
		alpha = ALPHA2,kappa = KAPPA2)
	F2 = pexpweibull(q = tt,lambda = LAMBDA2,
		alpha = ALPHA2,kappa = KAPPA2)
	calc_copula_dens(D1 = D1,D2 = D2,
		F1 = F1,F2 = F2,
		copula = COPULA,THETA = THETA)
}
g2 = gg(tt); g2

abs(g1 - g2)


}

# ----------
# Simulation Ideas
# ----------
if(FALSE){

# Num replicates
RR = 5e2

# Sample size grid
NN = c(500,1000,2000)

# Set alpha1 = 5, lambda1 = 1.2, alpha2 = 2.1, lambda2 = 2.0
	# Pick based on real data

# Set kappa1 for Exp-Weibull to 2

# Degree of censoring: aim for 20% censoring

# Set grid size to stepsize = 0.1 or 0.25 for now

# Theta: Can be independent or moderate dependence
	# For Clayton: Theta = 0 or Theta = 2
	# For Gumbel: Theta = 1 or Theta = 2

# If we don't know the truth ...
	# How often is the true configuration 
	#		(copula and distribution) selected 
	#		based on max(LL)? 
	# Does this approach work in practice?
	# Maybe max(BIC) should max(LL) not work

# If we already know the configuration, then ...
	# Bias
	# SE
	# SEE
	# Coverage Probability for alpha = 0.05



}

# ----------
# Running code example
# ----------

if( TRUE ){ 	# Specify parameters/arguments

COPULAS = c("Independent","Clayton","Gumbel")
DISTs		= c("weibull","expweibull")
COPULA 	= COPULAS[3]
dist1 	= DISTs[1]
NN 			= 4e3
theta 	= ifelse(COPULA == "Independent",0,5)
if( COPULA == "Clayton" && theta < 0 ) stop("theta issue")
if( COPULA == "Gumbel" && theta < 1 ) stop("theta issue")
alpha1 	= 1.2
lambda1 = 4
kappa1 	= ifelse(dist1 == "weibull",1,1/2)
alpha2 	= 1.6
lambda2 = 5
propC 	= 0.2
uPARS = log(c(alpha1,lambda1,kappa1,theta))
if( COPULA == "Gumbel" ) uPARS[4] = log(theta - 1)
names(uPARS) = sprintf("log_%s",c("A","L","K","T"))
# true_PARS

cPARS = c(alpha1,lambda1,kappa1,theta)
names(cPARS) = c("alpha1","lambda1","kappa1","theta")
TRUTH = list(COPULA = COPULA,DIST = dist1,
	uPARS = uPARS,cPARS = cPARS)
print(TRUTH)

}
if( TRUE ){ 	# Simulate dataset

set.seed(2)
one_rep = sim_replicate(copula = TRUTH$COPULA,
	dist1 = TRUTH$DIST,NN = NN,theta = theta,
	alpha1 = alpha1,lambda1 = lambda1,
	kappa1 = kappa1,alpha2 = alpha2,
	lambda2 = lambda2,propC = propC,
	verb = TRUE)
#one_rep$PARAMS
#print(table(one_rep$DATA$D) / NN)
#print(table(one_rep$DATA$delta) / NN)

print(table(D = one_rep$DATA$D,
	Delta = one_rep$DATA$delta))

upKAPPA = ifelse(TRUTH$DIST == "weibull",0,1)

aa = calc_CDFs(DATA = one_rep$DATA,
	PARS = TRUTH$uPARS,COPULA = TRUTH$COPULA)

}
if( FALSE ){ 	# Optimize

print(TRUTH)

DATA = one_rep$DATA
log_time75 = as.numeric(log(quantile(DATA$time[DATA$delta == 1],0.75)))
log_time75

# Agnostic
param_grid = list(
	A = seq(-1,2,0.1),
	L = log_time75 + seq(-1,1,0.1),
	K = seq(-0.5,0.5,0.1),
	T = seq(-1,3,0.25)
	)
param_grid

# If we have the truth ...
ss = 0.1
param_grid = list(
	A = TRUTH$uPARS[1] + seq(-0.5,0.5,ss),
	L = TRUTH$uPARS[2] + seq(-0.5,0.5,ss),
	K = TRUTH$uPARS[3] + seq(-0.5,0.5,ss),
	T = TRUTH$uPARS[4] + seq(-0.5,0.5,ss))
if( is.infinite(TRUTH$uPARS[4]) ){
	param_grid$T = seq(-1,3,ss)
}
param_grid

sapply(param_grid,length)
nGRID = prod(sapply(param_grid,length)); nGRID

run_ana = run_analyses(
	DATA = one_rep$DATA,
	COPULAS = COPULA,
	# upKAPPA = upKAPPA,
	param_grid = param_grid,
	verb = TRUE,PLOT = TRUE)

class(run_ana)
length(run_ana)

OO = opt_sum(OPT = run_ana); OO
print(TRUTH)

solu = which(OO$COPULA == TRUTH$COPULA
	& OO$DIST == TRUTH$DIST)
solu
# run_ana[[solu]]$RES$GOPT_PRE
# run_ana[[solu]]$RES$out
run_ana[[solu]]$RES$cout
print(TRUTH)

GRID = smart_df(run_ana[[solu]]$RES$GRID)
dim(GRID); head(GRID)

out = get_PROFILE(GRID = GRID,PLOT = TRUE)
# out = get_PROFILE(GRID = GRID[which(GRID$log_theta < 5),],PLOT = TRUE)
out

iPARS = out$VAL
iPARS
upPARS = c(1,1,0,1)

wrap_NR(DATA = one_rep$DATA,
	PARS = iPARS,
	COPULA = COPULA,
	upPARS = upPARS)

# Check survival
plot_SURVs(run_ANA = run_ana,
	MULTIPLE = TRUE,ALPHA = 0.4)



}
if( TRUE ){ 	# Test optimization

my_dirs$rep_dir = file.path(my_dirs$curr_dir,"REPS")

tCOPULA = c("Independent","Clayton","Gumbel")[2]
tDIST		= c("weibull","expweibull")[2]
rr 			= 215
NN			= 2e4

repCDN_dir = file.path(my_dirs$rep_dir,
	sprintf("C.%s_D.%s",tCOPULA,tDIST))

rds_fn = file.path(repCDN_dir,sprintf("R.%s.rds",rr))
one_rep = readRDS(rds_fn)
one_rep$DATA = one_rep$DATA[seq(NN),]
smart_table(D = one_rep$DATA$D,Delta = one_rep$DATA$delta)
one_rep$PARAMS

uPARS = get_uPARS(PARAMS = one_rep$PARAMS)
uPARS

aa = calc_CDFs(DATA = one_rep$DATA,
	PARS = uPARS,COPULA = tCOPULA)

stepsize = 0.25
param_grid = list(
	A = uPARS[1] + seq(-0.5,0.5,stepsize),
	L = uPARS[2] + seq(-0.5,0.5,stepsize),
	K = uPARS[3] + seq(-0.5,0.5,stepsize),
	T = uPARS[4] + seq(-0.5,0.5,stepsize))

if( is.infinite(uPARS[4]) ){
	param_grid$T = seq(-1,1,stepsize)
}
param_grid

# Estimate assuming truth known
DIST 		= c(tDIST,"weibull")[1]
COPULA	= c(tCOPULA,"Clayton")[1]

run_ana = run_analyses(
	DATA = one_rep$DATA,
	COPULAS = COPULA,
	upKAPPA = ifelse(DIST == "weibull",0,1),
	param_grid = param_grid,
	# param_grid = seq(-1,3,0.05),
	verb = TRUE)


class(run_ana)
length(run_ana)

OO = opt_sum(OPT = run_ana); OO

solu = 2
GRID = smart_df(run_ana[[solu]]$RES$GRID); head(GRID)
out = get_PROFILE(GRID = GRID,PLOT = TRUE)

# Plot survival curves
plot_SURVs(run_ANA = run_ana,
	MULTIPLE = TRUE,ALPHA = 0.4)

# Check why no valid grid points
ref_LL(DATA = one_rep$DATA,
	PARS = uPARS,
	COPULA = COPULA)

ref_LL_cpp(DATA = one_rep$DATA,
	PARS = uPARS,
	COPULA = COPULA)

wrap_LL(DATA = one_rep$DATA,
	PARS = uPARS,verb = TRUE)

# Calculate subject's values
one_rep$PARAMS
KAPPA = one_rep$PARAMS$kappa1
ALPHA = one_rep$PARAMS$alpha1
LAMBDA = one_rep$PARAMS$lambda1
THETA = one_rep$PARAMS$theta
KAL = KAPPA * ALPHA / LAMBDA; KAL
ii = 2145

XDL = one_rep$DATA$time[ii] / LAMBDA
enXDLa = exp(-XDL^ALPHA)

F1 = 1.0 - enXDLa
if( KAPPA != 1.0 ) F1 = F1^KAPPA;
S1 = 1.0 - F1
D1 = KAL * XDL^(ALPHA - 1.0) * enXDLa
if( KAPPA != 1.0 ) D1 = D1 * F1 / ( 1.0 - enXDLa );
D2 = one_rep$DATA$dens_t2[ii]
F2 = 1 - one_rep$DATA$surv_t2[ii]
F_T1_T2 = calc_copula(F1 = F1,F2 = F2,
	copula = COPULA,THETA = THETA)
f_T1_T2 = calc_copula_dens(D1 = D1,D2 = D2,
	F1 = F1,F2 = F2,copula = COPULA,
	THETA = THETA,F_T1_T2 = F_T1_T2)
D1; D2; F1; F2; F_T1_T2; f_T1_T2

(F1^(-THETA) + F2^(-THETA) - 1)^(-1/THETA) * (D1/F1^(THETA+1) + D2/F2^(THETA+1))


}
if( FALSE ){ 	# Debug optimization

my_dirs$rep_dir = file.path(my_dirs$curr_dir,"REPS")
my_dirs$opt_dir = file.path(my_dirs$curr_dir,"OPTS")

COPULA 	= c("Independent","Clayton","Gumbel")[1]
DIST		= c("weibull","expweibull")[1]
rr 			= 4
NN			= 5e3

# Import rep
repCDN_dir = file.path(my_dirs$rep_dir,
	sprintf("C.%s_D.%s",COPULA,DIST))
rds_fn = file.path(repCDN_dir,sprintf("R.%s.rds",rr))
one_rep = readRDS(rds_fn)
one_rep$DATA = one_rep$DATA[seq(NN),]

# Import opt
tmp_name = sprintf("C.%s_D.%s_N.%s",COPULA,DIST,NN)
opt_fn = file.path(my_dirs$opt_dir,tmp_name,sprintf("R.%s.rds",rr))
run_ana = readRDS(opt_fn)

OO = opt_sum(OPT = run_ana); OO
	solu = 4
	names(run_ana[[solu]]$RES)
	GRID = smart_df(run_ana[[solu]]$RES$GRID); head(GRID)
	out = get_PROFILE(GRID = GRID,PLOT = TRUE)
	run_ana[[solu]]$RES$GOPT_PRE
	run_ana[[solu]]$RES$GOPT
	run_ana[[solu]]$RES$GRAD
	run_ana[[solu]]$RES$out
	run_ana[[solu]]$RES$cout

# Check gradient
COPULA_2 = run_ana[[solu]]$copula; COPULA_2
wrap_LL = function(PARS){
	# PARS = iPARS
	dMrs_cLL(XX = one_rep$DATA$time,
		DELTA = one_rep$DATA$delta,
		D2 = one_rep$DATA$dens_t2,
		S2 = one_rep$DATA$surv_t2,
		PARS = PARS,
		copula = COPULA_2,
		verb = FALSE)
}
wrap_GRAD = function(PARS){
	# PARS = iPARS
	out = dMrs_cGRAD(XX = one_rep$DATA$time,
		DELTA = one_rep$DATA$delta,
		D2 = one_rep$DATA$dens_t2,
		S2 = one_rep$DATA$surv_t2,
		PARS = PARS,
		copula = COPULA_2,
		upPARS = upPARS)
	
	# out = grad(wrap_LL,PARS)
	c(out)
}

iPARS = run_ana[[solu]]$RES$out$EST
iPARS
upPARS = rep(1,4)
upPARS[3] = ifelse(iPARS[3] == 0,0,1)
upPARS[4] = ifelse(iPARS[4] == -Inf,0,1)
upPARS
old_LL = wrap_LL(PARS = iPARS); old_LL
eps = 1e-10
test_GRAD = rep(0,4)
for(idx in seq(4)){
	# idx = 3
	if( upPARS[idx] == 0 ) next
	shift = rep(0,4); shift[idx] = eps
	new_LL = wrap_LL(PARS = iPARS + shift); new_LL
	diff_LL = new_LL - old_LL; diff_LL
	test_GRAD[idx] = diff_LL / eps; test_GRAD[idx]
}
test_GRAD
nGRAD = sqrt(sum(test_GRAD^2)); nGRAD

sqrt(sum(wrap_GRAD(PARS = iPARS)^2))


}

# Run full simulation
# fsim = full_sim(copula = copula,dist1 = dist1,NN = NN,theta = theta,
	# alpha1 = alpha1,lambda1 = lambda1,kappa1 = kappa1,alpha2 = alpha2,
	# lambda2 = lambda2,propC = propC,RR = 1e2,param_grid = seq(-1,3,0.5),
	# upKAPPA = 1,verb = TRUE)
# fsim

# Generate Figure 1 about comparing confidence intervals
gen_fig_CIs = function(theta = 3,lambda1 = 6,VERB = TRUE){
	if(FALSE){
		VERB = TRUE
		
	}
	
	# Simulate dataset
	COPULA 	= "Clayton"
	dist1 	= "weibull"
	upKAPPA = ifelse(dist1 == "weibull",0,1)
	NN			= 1e2
	# theta		= 3
	theta2	= ifelse(theta < 2,"Mild Dependency","High Dependency")
	if( COPULA == "Clayton" && theta < 0 ) stop("theta issue")
	if( COPULA == "Gumbel" && theta < 1 ) stop("theta issue")
	alpha1 	= 3
	# lambda1 = 6
	kappa1 	= ifelse(dist1 == "weibull",1,2)
	alpha2 	= 3
	lambda2 = 6
	propC 	= 0.4
	
	set.seed(1)
	one_rep = sim_replicate(copula = COPULA,dist1 = dist1,
		NN = NN,theta = theta,alpha1 = alpha1,lambda1 = lambda1,
		kappa1 = kappa1,alpha2 = alpha2,lambda2 = lambda2,
		propC = propC,verb = VERB)
	
	# Optimize
	param_grid 	= seq(-2,4,0.25)
	vec_time 		= seq(0,round(max(one_rep$DATA$time),0) + 2,
									length.out = 100)
	
	run_ana = run_analyses(DATA = one_rep$DATA,
		THETAs = theta,upKAPPA = upKAPPA,
		copula = COPULA,param_grid = param_grid,
		vec_time = vec_time,verb = VERB)
	names(run_ana)

	# Get/Prep Survival
	sdat = run_ana[[1]]$PRED
	# dim(sdat); sdat[1:5,]
	
	fdat1 = smart_df(sdat[,c("time","surv","low_surv","high_surv")],CI = 1)
	fdat2 = smart_df(sdat[,c("time","surv","low_surv2","high_surv2")],CI = 2)
	fdat2 = name_change(fdat2,"low_surv2","low_surv")
	fdat2 = name_change(fdat2,"high_surv2","high_surv")
	fdat 	= rbind(fdat1,fdat2)
	dim(fdat); fdat[1:3,]
	
	fdat$CI_2 = factor(fdat$CI,levels = c(1,2),
		labels = c("delta","log-log"))
	
	# Plot Survival
	my_themes = theme(text = element_text(size = 28),
		legend.position = c("none","bottom")[2],
		plot.title = element_text(hjust = 0.5),
		panel.background = element_blank(),
		panel.grid.major = element_line(colour = "grey50",
			size = 0.5,linetype = "dotted"),
		panel.border = element_rect(colour = "black",
			fill = NA,size = 1))
	
	if(FALSE){ #ggplot code
	gg = ggplot(data = fdat,
		mapping = aes(x = time,y = surv,
			group = CI_2,
			fill = CI_2,
			color = CI_2)) +
		geom_line(size = 1.3,alpha = 1,
			mapping = aes(color = CI_2
				#,linetype = CI_2
				),
			show.legend = FALSE) +
		geom_ribbon(alpha = 0.3,size = 0.8,
			linetype = "dashed",
			mapping = aes(ymin = low_surv,
				ymax = high_surv),
			show.legend = FALSE) +
		# facet_grid(~ CI_2) +
		# ylim(c(0,1)) + 
		# xlim(c(0,4)) +
		xlab("Time") + 
		ylab("Survival Probability") +
		ggtitle("Confidence Interval Types") + 
		labs(fill = "Type of CI") + 
		my_themes
	gg
	}
	
	plot(fdat1[,c("time","surv")],
		type = "l",lwd = 2,bty = "n",
		xlim = c(0,max(fdat$time)),
		ylim = c(min(fdat$low_surv),max(fdat$high_surv)),
		xaxs = "i",xlab = "Time (years)",
		ylab = "Survival Probability",
		main = sprintf("Types of Confidence Intervals\nwith %s",
			theta2),
		cex.lab = 1.4,las = 1,cex.main = 1.75)
	abline(h = c(0,1),lty = 2,lwd = 0.5)
	lines(fdat1[,c("time","low_surv")],
		col = "red",lwd = 2,lty = 2)
	lines(fdat1[,c("time","high_surv")],
		col = "red",lwd = 2,lty = 2)
	lines(fdat2[,c("time","low_surv")],
		col = "blue",lwd = 2,lty = 2)
	lines(fdat2[,c("time","high_surv")],
		col = "blue",lwd = 2,lty = 2)
	legend(0,0.4,
		legend = c("Delta","Log-Log"),
		col = c("red","blue"),
		# lwd = rep(4,2),
		pch = rep(15,2),
		cex = 1.75,pt.cex = 3,
		bty = "n")
	
	# Plot heatmap of grid
	
	
}

# smarter:::show_png()
png(file.path(my_dirs$curr_dir,"dMrs_0.png"),
	units = "px",height = 2000,width = 2000,res = 250,
	type = "cairo",pointsize = 20)
gen_fig_CIs(theta = 3,lambda1 = 6,VERB = FALSE)
dev.off()


# ----------
# Code for analyzing real data
# ----------

# Import data
wdat_fn = file.path(my_dirs$pack_dir,
	"inst/extdata","bc.txt")
wDAT = data.table::fread(wdat_fn,data.table = FALSE)
rdat_fn = file.path(my_dirs$pack_dir,
	"inst/extdata","frepop.txt")
rDAT = data.table::fread(rdat_fn,data.table = FALSE)

names(wDAT)[names(wDAT) == "Status"] = "delta"
names(wDAT)[names(wDAT) == "Age"] = "age"
names(wDAT)[names(wDAT) == "Time"] = "time"
wDAT[1:3,]

rDAT$sex = "female"
rDAT$Age = as.integer(gsub("\\+","",rDAT$Age))
rDAT = rDAT[order(rDAT$Year,rDAT$Age,rDAT$sex),]
rDAT = rDAT[,c("Year","Age","qx","sex")]
rDAT[1:3,]

# Run matching with reference dataset
rd = refData_match(wDAT = wDAT,rDAT = rDAT)
rd[1:5,]
table(rd$delta)

## Analyze
res1 = full_ana_opt(DATA = rd,max_year = 100,
	param_grid = seq(-3,5,0.5),verb = TRUE)
names(res1)

# Get Estimates
cohort = "French Breast Cancer Registry"
NAME 	 = "FPOP"
out_est = full_ana_est(RES = res1,COHORT = cohort)

# Get Survival Probabilities
out_surv = full_ana_surv(RES = res1,
	myYEARS = c(1,3,5,10))
out_surv

# Provide plots
full_ana_plots(curr_dir = my_dirs$curr_dir,
	RES = res1,
	NAME = NAME)


# ----------
# Run second real data analysis
# ----------
conv_date_to_year = function(DATE){
	DATE = as.character(DATE)
	out = as.integer(strsplit(DATE,"-")[[1]][1])
	return(out)
}
data(slopop)
data(rdata)

# Prep reference data: need cols Year, Age, qx, sex
str(slopop) # this is an array!
bb 		= as.array(slopop); dim(bb)
SEXs 	= dimnames(bb)[[3]]
YEARs = dimnames(bb)[[2]]
rDAT = c()
for(SEX in SEXs){
for(YEAR in YEARs){
	# SEX = SEXs[1]; YEAR = YEARs[1]
	
	tmp_vec = bb[,YEAR,SEX]
	
	tmp_df = smart_df(Year = YEAR,
		Age = names(tmp_vec),
		qx = as.numeric(tmp_vec),
		sex = SEX)
	rDAT = rbind(rDAT,tmp_df)
	rm(tmp_df)
	
}}
rDAT$Year = as.integer(rDAT$Year)
rDAT$Age = as.integer(rDAT$Age)
dim(rDAT); rDAT[1:5,]

# Prep working data: need cols age, time, delta, datediag_yr, dateEvent_yr, sex
wDAT = rdata
wDAT = name_change(wDAT,"cens","delta")
wDAT$sex[wDAT$sex == "1"] = "male"
wDAT$sex[wDAT$sex == "2"] = "female"
# note: wDAT$year is Julian date aka days since 1960-01-01
wDAT$datediag_yr = as.Date(wDAT$year)
wDAT$dateEvent_yr = wDAT$datediag_yr + wDAT$time
wDAT$datediag_yr = sapply(wDAT$datediag_yr,conv_date_to_year)
wDAT$dateEvent_yr = sapply(wDAT$dateEvent_yr,conv_date_to_year)
wDAT$time = wDAT$time / 365.25
wDAT[1:3,]

# Match data for input to dMrs
rd = refData_match(wDAT = wDAT,rDAT = rDAT)
rd[1:5,]
table(rd$delta)

cohort 	= "Slovenia Colon Cancer Registry"
NAME 		= "Slovenia"

# Analyze
run_ana = run_analyses(
	DATA = rd,
	THETAs = c(0,2/3,2,6),COPULAS = "Clayton",
	# THETAs = c(1,4/3,2,4),COPULAS = "Gumbel",
	# upKAPPA = c(0,1)[2],
	# COPULAS = c("Clayton","Gumbel")[2],upKAPPA = 1,
	param_grid = seq(-3,4,0.25),
	verb = TRUE)
length(run_ana)

# Check estimates

idx = 1
run_ana[[idx]]$copula
run_ana[[idx]]$RES[c("out","cout","LL")]
dim(run_ana[[idx]]$RES$GRID)

plot_LL(GPROF = run_ana[[idx]]$RES$GPROF,
	GOPT = run_ana[[idx]]$RES$GOPT,
	COPULA = run_ana[[idx]]$copula)

# Plot survival
plot_SURVs(run_ANA = run_ana,MULTIPLE = TRUE,ALPHA = 0.4)
plot_SURVs(run_ANA = run_ana,MULTIPLE = !TRUE,ALPHA = 0.4)

### OLDER analysis code below




res2 = full_ana_opt(DATA = rd,max_year = 100,
	param_grid = seq(-3,5,0.25),verb = TRUE)
names(res2)

# Get Estimates
out_est = full_ana_est(RES = res2,COHORT = cohort)

# Get Survival Probabilities
out_surv = full_ana_surv(RES = res2,
	myYEARS = c(1,3,5,10))
out_surv

# Provide plots
full_ana_plots(curr_dir = my_dirs$curr_dir,
	RES = res2,
	NAME = NAME)



###

