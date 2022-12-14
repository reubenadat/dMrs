# R script with R functions for dMrs

# ----------
# Shortcuts, Install/Load Libraries
# ----------
rm(list=ls())
setdirs = function(){
	
	where_file = gsub("\\\\","/",Sys.getenv("HOME"))
	clust_name = Sys.getenv("SLURM_CLUSTER_NAME")

	if( where_file == "C:/Users/Admin/Documents" ){
		desk_dir = "C:/Users/Admin/Desktop"
		git_dir = file.path(desk_dir,"github")
		curr_dir = file.path(desk_dir,"dMrs")
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
		"ggplot2","viridis","relsurv","dMrs")
	all_packs = as.character(installed.packages()[,1])
	for(pack in req_packs){
		
		if( pack %in% all_packs ){
			library(pack,character.only = TRUE)
			next
		}
		
		if( pack == "smarter" ){
			install_github(repo = "pllittle/smarter")
		} else if( pack == "dMrs" ){
			# stop("Look at reubenadat/dMrs for installation")
			install(pack_dir,build_vignettes = FALSE)
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


# ----------
# Running code example
# ----------

# Specify parameters/arguments
COPULAS = c("Clayton","Gumbel")
DISTs		= c("weibull","expweibull")
COPULA 	= COPULAS[1]
dist1 	= DISTs[2]
NN 			= 2e3
theta 	= 4
if( COPULA == "Clayton" && theta < 0 ) stop("theta issue")
if( COPULA == "Gumbel" && theta < 1 ) stop("theta issue")
alpha1 	= 3
lambda1 = 4
kappa1 	= ifelse(dist1 == "weibull",1,2)
alpha2 	= 3
lambda2 = 6
propC 	= 0.1
true_PARS = log(c(alpha1,lambda1,kappa1,theta)); true_PARS

# Simulate dataset
set.seed(1)
one_rep = sim_replicate(copula = COPULA,dist1 = dist1,
	NN = NN,theta = theta,alpha1 = alpha1,lambda1 = lambda1,
	kappa1 = kappa1,alpha2 = alpha2,lambda2 = lambda2,
	propC = propC,show = TRUE)
one_rep$PARAMS
table(one_rep$DATA$D)
table(one_rep$DATA$delta)

# Run analyses
my_copula = COPULAS[1]; my_copula
if( my_copula == "Clayton" ) 	THETAs = c(0,2/3,2,4,6)
if( my_copula == "Gumbel" ) 	THETAs = c(1,4/3,2,4)
THETAs 	= round(THETAs,2); THETAs
upKAPPA = c(0,1)[2] # 0 = Weibull, 1 = Exp-Weibull

run_ana = run_analyses(DATA = one_rep$DATA,
	THETAs = THETAs,upKAPPA = upKAPPA,
	copula = my_copula,param_grid = seq(-2,4,0.25),
	vec_time = seq(0,round(max(one_rep$DATA$time),0)),
	show = TRUE)
names(run_ana)
run_ana[[1]]$RES[c("out","cout")]

plot_SURVs(run_ANA = run_ana,
	MULTIPLE = TRUE,ALPHA = 0.4)
plot_SURVs(run_ANA = run_ana,
	MULTIPLE = FALSE,ALPHA = 0.4)

## Test delta method and variance
tt = 5 # try time
tmp_surv = function(PARS){
	
	CDF = pexpweibull(q = tt,
		lambda = exp(PARS[2]),
		alpha = exp(PARS[1]),
		kappa = exp(PARS[3]))
	SURV = 1 - CDF
	return(SURV)
	
}
names(run_ana)
bb = run_ana[["theta = 4"]]
bb$RES$cout
bb$PRED[which(bb$PRED$time == tt),]
PARS = bb$RES$out$EST
PARS
tmp_surv(PARS = PARS)

nabla = numDeriv::grad(tmp_surv,PARS)
VV = t(nabla) %*% bb$RES$COVAR %*% nabla; VV
sqrt(VV)

# Run full simulation
# fsim = full_sim(copula = copula,dist1 = dist1,NN = NN,theta = theta,
	# alpha1 = alpha1,lambda1 = lambda1,kappa1 = kappa1,alpha2 = alpha2,
	# lambda2 = lambda2,propC = propC,RR = 1e2,param_grid = seq(-1,3,0.5),
	# upKAPPA = 1,show = TRUE)
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
		propC = propC,show = VERB)
	
	# Optimize
	param_grid 	= seq(-2,4,0.25)
	vec_time 		= seq(0,round(max(one_rep$DATA$time),0) + 2,
									length.out = 100)
	
	run_ana = run_analyses(DATA = one_rep$DATA,
		THETAs = theta,upKAPPA = upKAPPA,
		copula = COPULA,param_grid = param_grid,
		vec_time = vec_time,show = VERB)
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
	"inst/extdata","bc_git.txt")
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
	param_grid = seq(-3,5,0.5),show = TRUE)
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
res2 = full_ana_opt(DATA = rd,max_year = 100,
	param_grid = seq(-3,5,0.25),show = TRUE)
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

