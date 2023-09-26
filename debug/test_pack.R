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
	pack_dir
	
	# Package check
	req_packs = c("devtools","pkgbuild","installr",
		"Rtools","numDeriv","smarter",
		"ggplot2","viridis","relsurv","copula","dMrs")
	for(pack in req_packs){
		# pack = req_packs[4]; pack
		
		chk_pack = tryCatch(find.package(pack),
			error = function(ee){NULL})
		chk_pack
		
		if( !is.null(chk_pack) ){
			library(pack,character.only = TRUE)
			next
		}
		
		if( pack == "smarter" ){
			install_github(repo = "pllittle/smarter")
		} else if( pack == "Rtools" ){
			aa1 = find_rtools(); aa1
			aa2 = grepl("Rtools",Sys.getenv("PATH")); aa2
			rtools4_dir = gsub("\\\\","/",Sys.getenv("RTOOLS40_HOME"))
			
			if( !aa1 || rtools4_dir == "" ){
				install.Rtools(check = TRUE,
					check_r_update = FALSE,
					GUI = FALSE)
				cat("Rerun setdirs\n")
				return(NULL)
			}
			
			rtools4_dir = file.path(rtools4_dir,"usr/bin")
			rtools4_dir
			
			if( !aa2 ){
				curr_PATH = Sys.getenv("PATH")
				new_PATH = sprintf("%s;%s",rtools4_dir,curr_PATH)
				Sys.setenv("PATH" = new_PATH)
			}
			
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

if(FALSE){

ALP = 2
LAM = 5
tt = seq(0,LAM * 2,length.out = 1e3)
PDF = dweibull(x = tt,shape = ALP,scale = LAM)
plot(tt,PDF,type = "l")

}
if(FALSE){ 		# Purposely shape bathtub

NN = 2e3
LAMs = 4
tt = seq(0.1,max(LAMs)*3,length.out = NN)

ALPs = c(0.5,1,1.2)
	ALPs = ALPs[ALPs > 1]
	# ALPs = ALPs[ALPs < 1]
	# ALPs = ALPs[ALPs == 1]
KAPs = c(0.1,0.5,1,3,5)
	# KAPs = KAPs[KAPs > 1]
	# KAPs = KAPs[KAPs < 1]
	# KAPs = KAPs[KAPs == 1]

len = length(ALPs) * length(LAMs) * length(KAPs)
HH = matrix(NA,nrow = NN,ncol = len)
colnames(HH) = seq(len)

ee = 1
for(ALP in ALPs){
for(LAM in LAMs){
for(KAP in KAPs){
	
	log_dd = dexpweibull(x = tt,lambda = LAM,
		alpha = ALP,kappa = KAP,log = TRUE)
	pp = pexpweibull(q = tt,lambda = LAM,
		alpha = ALP,kappa = KAP,log.p = FALSE)
	hh = exp(log_dd - log(1 - pp)) # hazard = dens / surv
	HH[,ee] = hh
	colnames(HH)[ee] = sprintf("ALP=%s;LAM=%s;KAP=%s",ALP,LAM,KAP)
	ee = ee + 1
	
}}}

LWD = 5
vec_colors = smart_colors(len)
for(ee in seq(len)){
	# ee = 1
	if( ee == 1 ){
		# max_yy = quantile(c(HH),0.95)
		plot(tt,HH[,ee],col = vec_colors[ee],
			ylim = c(0,1),type = "l",lwd = LWD)
	} else {
		lines(tt,HH[,ee],col = vec_colors[ee],
			lwd = LWD)
	}
}

legend("right",legend = colnames(HH),
	pch = rep(16,len),col = vec_colors,
	pt.cex = rep(2,len),cex = 0.75,bty = "n")
abline(v = LAMs,lty = 2)

# Takeaways to note:
#	To estimate kappa 'well', individuals from the 
#		first underlying event with earlier times
#		are more informative than those that occur 
#		later

}
if( TRUE ){ 	# Specify parameters/arguments

COPULAS = c("Independent","Clayton","Gumbel")
DISTs		= c("weibull","expweibull")
tCOPULA = COPULAS[3]
tDIST 	= DISTs[2]
NN 			= 5e3
theta 	= ifelse(tCOPULA == "Independent",0,5)
if( tCOPULA == "Clayton" && theta < 0 ) stop("theta issue")
if( tCOPULA == "Gumbel" && theta < 1 ) stop("theta issue")
alpha1 	= 1.2
lambda1 = 12
kappa1 	= ifelse(tDIST == "weibull",1,0.8)
alpha2 	= 3.0
# ALPHA: if < 1, events happen early
#		 if > 1, events happen later
lambda2 = 10
propC 	= 0.2
uPARS = log(c(alpha1,lambda1,kappa1,theta))
if( tCOPULA == "Gumbel" ) uPARS[4] = log(theta - 1)
names(uPARS) = sprintf("log_%s",c("A","L","K","T"))
# true_PARS

cPARS = c(alpha1,lambda1,kappa1,theta)
names(cPARS) = c("alpha1","lambda1","kappa1","theta")
TRUTH = list(COPULA = tCOPULA,DIST = tDIST,
	uPARS = uPARS,cPARS = cPARS)
# print(TRUTH)

}
if( TRUE ){ 	# Simulate dataset

set.seed(1)
one_rep = sim_replicate(copula = TRUTH$COPULA,
	dist1 = TRUTH$DIST,NN = NN,theta = theta,
	alpha1 = alpha1,lambda1 = lambda1,
	kappa1 = kappa1,alpha2 = alpha2,
	lambda2 = lambda2,propC = propC,
	verb = TRUE)

upKAPPA = ifelse(TRUTH$DIST == "weibull",0,1)

aa = calc_CDFs(DATA = one_rep$DATA,
	PARS = TRUTH$uPARS,COPULA = TRUTH$COPULA)
# head(aa)

print(table(D = one_rep$DATA$D,Delta = one_rep$DATA$delta))

# check num obs events
GROUP = bin_cont_var(VAR = one_rep$DATA$time,
	NUM_GROUPS = 4,binNUM = TRUE); # smart_table(GROUP)
print(smart_table(D = one_rep$DATA$D[one_rep$DATA$delta==1],
	G = GROUP[one_rep$DATA$delta==1]))


}
if( FALSE ){ 	# Optimize

one_rep$PARAMS

# If we make some constraints
len = 20
A_range = c(0.5,2)
L_range = exp(quantile(log(one_rep$DATA$time),c(0.65,0.95)))
K_range = c(0.5,2)
T_range = c(1,10)

A_ugrid = log(seq(A_range[1],A_range[2],length.out = len))
L_ugrid = log(seq(L_range[1],L_range[2],length.out = len))
K_ugrid = log(seq(K_range[1],K_range[2],length.out = len))
T_ugrid = log(seq(T_range[1],T_range[2],length.out = len))

param_grid = list(A = A_ugrid,
	L = L_ugrid,K = K_ugrid,T = T_ugrid)
# param_grid

sapply(param_grid,length)
nGRID = prod(sapply(param_grid,length)); nGRID

run_ana = run_analyses(
	DATA = one_rep$DATA,
	COPULAS = tCOPULA,
	# upKAPPA = upKAPPA,
	param_grid = param_grid,
	# param_grid = seq(-4,4,0.5),
	verb = TRUE,PLOT = TRUE)
length(run_ana)

OO = opt_sum(OPT = run_ana); OO
print(TRUTH)

# Check survival
plot_SURVs(run_ANA = run_ana,
	MULTIPLE = TRUE,ncol = 2,
	ALPHA = 0.4)

solu = 2
solu = which(OO$COPULA == TRUTH$COPULA
	& OO$DIST == TRUTH$DIST)
solu

uPARS = run_ana[[solu]]$RES$out$EST; # uPARS
run_ana[[solu]]$RES$GOPT_PRE
run_ana[[solu]]$RES$cout

# EST pars
aa = calc_CDFs(DATA = one_rep$DATA,
	PARS = uPARS,COPULA = run_ana[[solu]]$copula); # head(aa)

# True pars
aa = calc_CDFs(DATA = one_rep$DATA,
	PARS = TRUTH$uPARS,COPULA = TRUTH$COPULA)

out = get_PROFILE(
	GRID = run_ana[[solu]]$RES$GRID,
	COPULA = run_ana[[solu]]$copula,
	PLOT = TRUE); # out




}
if( TRUE ){ 	# Test optimization

my_dirs$rep_dir = file.path(my_dirs$curr_dir,"REPS")

tCOPULA = c("Independent","Clayton","Gumbel")[2]
tDIST		= c("weibull","expweibull")[1]
rr 			= 1
NN			= 1e4

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
head(aa)

stepsize = 0.2
bound = 1
param_grid = list(
	A = uPARS[1] + seq(-bound,bound,stepsize),
	L = uPARS[2] + seq(-bound,bound,stepsize),
	K = uPARS[3] + seq(-bound,bound,stepsize),
	T = uPARS[4] + seq(-bound,bound,stepsize))

if( is.infinite(uPARS[4]) ){
	# param_grid$T = seq(-bound,bound,stepsize)
	# param_grid$T = seq(-1,3,stepsize)
	param_grid$T = seq(4,9,0.5)
}
param_grid
prod(sapply(param_grid,length))

# Estimate assuming truth known
COPULA	= c(tCOPULA,"Independent","Clayton","Gumbel")[1]
DIST 		= c(tDIST,"weibull","expweibull")[1]

run_ana = run_analyses(
	DATA = one_rep$DATA,
	COPULAS = COPULA,
	upKAPPA = ifelse(DIST == "weibull",0,1),
	param_grid = param_grid,
	# param_grid = seq(-1,3,0.1),
	verb = TRUE)
if( length(run_ana) > 0 ){
	OO = opt_sum(OPT = run_ana)
	print(OO)
}

solu = which(OO$COPULA == tCOPULA & OO$DIST == tDIST); solu
solu = 1
GRID = smart_df(run_ana[[solu]]$RES$GRID); # head(GRID)
out = get_PROFILE(GRID = GRID,
	COPULA = run_ana[[solu]]$copula,
	PLOT = TRUE)
run_ana[[solu]]$RES$cout

# Check distribution of copula, any precision problems
aa = calc_CDFs(DATA = one_rep$DATA,
	PARS = run_ana[[solu]]$RES$out$EST,
	COPULA = COPULA)
head(aa)
smart_table(aa$CDF_1 %in% c(0,1))
smart_table(aa$CDF_2 %in% c(0,1))
smart_table(aa$F_T1_T2 == 0)
smart_table(aa$D1_D2 == 0)
smart_hist(aa$D1_D2)

# Plot survival curves
plot_SURVs(run_ANA = run_ana,
	MULTIPLE = TRUE,ALPHA = 0.4)

# Check gradient
uPARS = c(0.174083, 1.36877, 0, 5.72699)
upPARS = c(1,1,0,1)
# uPARS = wrap_NR(DATA = one_rep$DATA,PARS = uPARS,
	# COPULA = COPULA,upPARS = upPARS,mult = 5,
	# verb = TRUE); uPARS
aa = calc_CDFs(DATA = one_rep$DATA,
	PARS = uPARS,COPULA = COPULA)
head(aa)
smart_table(aa$CDF_2 %in% c(0,1))

stop("Double check how CDF and PDF copula are calculated, see if existing function exists")
if(FALSE){ # Check Gumbel

THETA = exp(uPARS[4]) + ifelse(tCOPULA == "Gumbel",1,0); THETA
ALPHA = exp(uPARS[1])
LAMBDA = exp(uPARS[2])

ALPHA_2 = 1.6; LAMBDA_2 = 5



head(aa[order(-aa$D2),])

tt = seq(0,100)
PDF = dweibull(x = tt,scale = LAMBDA,shape = ALPHA)
plot(tt,PDF,type = "l")

ii = 3
tt = one_rep$DATA$time[ii]
F1 = aa$CDF_1[ii]; #F1
F2 = aa$CDF_2[ii]; #F2
F_T1_T2 = aa$F_T1_T2[ii]; F_T1_T2
	# Use copula package functions to check calculation of copula CDF
	library(copula)
	# tmp_copula = gumbelCopula(param = THETA,dim = 2)
	tmp_copula = claytonCopula(param = THETA,dim = 2)
	bvd = mvdc(copula = tmp_copula,
		margins = c("weibull","weibull"),
		paramMargins = list(
			list(shape = ALPHA,scale = LAMBDA),
			list(shape = ALPHA_2,scale = LAMBDA_2))
		)
	pMvdc(x = c(tt,tt),mvdc = bvd) # we're good!!!

D1 = aa$D1[ii]; D1
D2 = aa$D2[ii]; D2
nlog_F1 = -log(F1)
nlog_F2 = -log(F2)
# D1_D2 = F_T1_T2 * 
	# ( (nlog_F1)^THETA + (nlog_F2)^THETA )^(1/THETA-1) *
	# ( (nlog_F1)^(THETA-1) * D1/F1 + (nlog_F2)^(THETA-1) * D2/F2 ); D1_D2
# log_D1_D2 = log(F_T1_T2) +
	# (1/THETA-1) * log( (nlog_F1)^THETA + (nlog_F2)^THETA ) +
	# log( (nlog_F1)^(THETA-1) * D1/F1 + (nlog_F2)^(THETA-1) * D2/F2 )
	# exp(log_D1_D2)
D1_D2 = aa$D1_D2[ii]; D1_D2
	dMvdc(x = c(tt,tt),mvdc = bvd)
	stop("The PDF of Gumbel")
	
	pweibull(q = tt,shape = ALPHA,scale = LAMBDA)
	

ALPHA2 = 1.6
LAMBDA2 = 5
D2 = dweibull(x = one_rep$DATA$time[ii],
	shape = ALPHA2,scale = LAMBDA2); D2
F2 = pweibull(q = one_rep$DATA$time[ii],
	shape = ALPHA2,scale = LAMBDA2); F2
one_rep$DATA[ii,]

}


smart_table(aa$CDF_1 %in% c(0,1))
smart_table(aa$F_T1_T2 %in% c(0,1))
smart_table(aa$D1 %in% c(0,1))
smart_table(aa$D1_D2 %in% c(0,1))

old_LL = wrapper_LL(DATA = one_rep$DATA,
	PARS = uPARS,COPULA = COPULA,verb = TRUE); old_LL
GRAD = wrapper_GRAD(DATA = one_rep$DATA,
	PARS = uPARS,COPULA = COPULA,
	upPARS = upPARS); GRAD
tGRAD = rep(0,4)
# ij = which(abs(GRAD) == max(abs(GRAD))); ij
ij = 4
tGRAD[ij] = GRAD[ij]
tGRAD = GRAD

for(cc in seq(0,40)){
	new_PARS = uPARS + 1/4^cc * tGRAD
	new_LL = wrapper_LL(DATA = one_rep$DATA,
		PARS = new_PARS,COPULA = COPULA,
		verb = TRUE); new_LL
	if( new_LL == -999 ) next
	
	if( new_LL > old_LL ){
		print("update")
		uPARS = new_PARS
		cat(sprintf("diff_LL = %s\n",new_LL - old_LL))
		break
	}
	
}
cc


shift = c(1e-5,5e-6,1e-6,1e-9,1e-11)[3]
tGRAD = rep(NA,4)

for(ij in seq(4)){
	# ij = 1
	if( upPARS[ij] == 0 ){
		tGRAD[ij] = 0
		next
	}
	new_PARS = uPARS
	new_PARS[ij] = uPARS[ij] + shift
	old_LL = wrapper_LL(DATA = one_rep$DATA,
		PARS = uPARS,COPULA = COPULA,verb = TRUE)
	new_LL = wrapper_LL(DATA = one_rep$DATA,
		PARS = new_PARS,COPULA = COPULA,verb = TRUE)
	tGRAD[ij] = (new_LL - old_LL) / shift
}
tGRAD
Rcpp_norm(tGRAD)

# Precision problem
XX = 1e-3
LAM = 5
ALP = 3

PL_log_CDF = function(XX,LAM,ALP){
	
	XDL = XX / LAM; XDL
	# exp(XDL^ALP)
	F1 = 1 - exp(-XDL^ALP); F1
	
	if( F1 <= 1e-10 ){
		log_F1 = ALP * log(XDL) + log(1 - 0.5 * XDL^ALP)
	} else {
		log_F1 = log(F1)
	}
	
	log_F1
	return(log_F1)
	
}

pweibull(q = XX,shape = ALP,scale = LAM,log.p = TRUE)
PL_log_CDF(XX = XX,LAM = LAM,ALP = ALP)


}
if( FALSE ){ 	# Debug optimization

my_dirs$rep_dir = file.path(my_dirs$curr_dir,"REPS")
my_dirs$opt_dir = file.path(my_dirs$curr_dir,"OPTS")

COPULA 	= c("Independent","Clayton","Gumbel")[1]
DIST		= c("weibull","expweibull")[1]
rr 			= 93
NN			= 2e4

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
	out = get_PROFILE(GRID = GRID,
		COPULA = run_ana[[solu]]$copula,
		PLOT = TRUE)
	run_ana[[solu]]$RES$GOPT_PRE
	run_ana[[solu]]$RES$GOPT
	run_ana[[solu]]$RES$GRAD
	run_ana[[solu]]$RES$out
	run_ana[[solu]]$RES$cout

# Check gradient
COPULA_2 = run_ana[[solu]]$copula; COPULA_2
iPARS = run_ana[[solu]]$RES$out$EST
iPARS

wrapper_GRAD(DATA = one_rep$DATA,
	PARS = iPARS,COPULA = COPULA_2,
	upPARS = rep(1,4))

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
if(FALSE){		# Test precision copula

aa = calc_CDFs(DATA = one_rep$DATA,
	PARS = run_ana[[solu]]$RES$out$EST,
	COPULA = COPULA)
head(aa)

idx 		= 1; aa[idx,]
F1 			= aa$CDF_1[idx]
F2 			= aa$CDF_2[idx]
ALPHA		= run_ana[[solu]]$RES$cout$EST[1]
LAMBDA	= run_ana[[solu]]$RES$cout$EST[2]
KAPPA		= run_ana[[solu]]$RES$cout$EST[3]
THETA 	= c(run_ana[[solu]]$RES$cout$EST[4],10)[1]
COPULA 	= "Clayton"

# Current Rcpp
F_T1_T2 = calc_copula(F1 = F1,F2 = F2,
	copula = COPULA,THETA = THETA); F_T1_T2

D1 = dexpweibull(x = one_rep$DATA$time[idx],
	lambda = LAMBDA,alpha = ALPHA,
	kappa = KAPPA); D1
D2 = one_rep$DATA$dens_t2[idx]; D2

calc_copula_dens(D1 = D1,D2 = D2,F1 = F1,F2 = F2,
	copula = COPULA,THETA = THETA,F_T1_T2 = F_T1_T2)

test_PDF = function(F1,F2,D1,D2,THETA,COPULA){
	
	if( COPULA == "Clayton" ){
		PDF_1 = (F1^(-THETA) + F2^(-THETA) - 1)^(-1/THETA-1)
		PDF_2 = (D1/F1^(THETA+1) + D2/F2^(THETA+1))
		PDF = PDF_1 * PDF_2
		
		if( PDF_1 == 0 || is.infinite(PDF_1) 
			|| PDF_2 == 0 || is.infinite(PDF_2) ){
			
			# stop("PDF_1 precision issue")
			log_CDFs = -THETA * log(c(F1,F2))
			log_mm = max(log_CDFs)
			log_PDF_1 = (-1/THETA - 1) *
				( log_mm + 
				log( sum(exp(log_CDFs - log_mm)) - 1 / exp(log_mm) ) )
			
			log_vec = rep(NA,2)
			log_vec[1] = log(D1) - (THETA+1)*log(F1)
			log_vec[2] = log(D2) - (THETA+1)*log(F2)
			log_PDF_2 = logSumExp(log_vec)
			
			log_PDF = log_PDF_1 + log_PDF_2
			log_PDF
			
			PDF = exp(log_PDF)
			PDF
		}
		
	} else if( COPULA == "Gumbel" ){
		stop("no code yet")
	}
	
	return(PDF)
}

test_PDF(F1 = F1,F2 = F2,D1 = D1,D2 = D2,
	THETA = THETA,COPULA = COPULA)

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

