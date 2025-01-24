# dMrs

## <a href="#TOC" name="TOC">Table of Contents</a>

* [Installation](#INSTALL)
* [Package Features](#FEAT)
* [Citation](#CITE)

## <a href="#INSTALL" name="INSTALL">Installation</a>

To install the package, run the following code for required dependencies and package.

```R
all_packs = as.character(installed.packages()[,1])
req_packs = c("devtools","Rcpp","RcppArmadillo",
	"Rmpfr","copula","ggplot2","viridis","smarter","dMrs")
run = 0

for(pack in req_packs){
	
	if( pack %in% all_packs ){
		library(pack,character.only = TRUE)
		next
	}
	
	if( pack == "smarter" ){
		devtools::install_github("pllittle/smarter")
		
	} else if( pack == "dMrs" ){
		devtools::install_github("reubenadat/dMrs")
		
	} else {
		install.packages(pack)
		
	}
	run = 1
	
}

if( run == 1 ) stop("Re-run this code")

```

Check out the vignette!

```R
vignette("dMrs")
```

## <a href="#FEAT" name="FEAT">Package Features</a>

The package's main functions include

1. Simulation: Varying underlying shape, scale, dependency parameters.
2. Optimization: Grid search and BFGS
3. Predicting survival
4. Matching working dataset subjects with reference dataset on age, sex, time period.

## <a href="#CITE" name="CITE">Citation</a>
1. Adatorwovor, R., Latouche, A., & Fine, J. P. (2022). A parametric approach to relaxing the independence assumption in relative survival analysis. The international journal of biostatistics, 18(2), 577-592.
2. Perme, M. P., Stare, J., & Estève, J. (2012). On estimation in relative survival. Biometrics, 68(1), 113-120.
3. Adatorwovor, R. & Little, P.L. (2024). OnlyA Dependent Competing Risks Model for Net Survival
Estimation. JASA, Under review.

