# Steps to create/check/install package from directory

rm(list = ls())
git_dir = sprintf("%s/../..",getwd()); git_dir
pack = "dMrs"
pack_dir = file.path(git_dir,pack)
if( !file.exists(pack_dir) ) q("no")

chk_pack = tryCatch(find.package(pack),
	error = function(ee){NULL}); chk_pack

if( !is.null(chk_pack) ){
	remove.packages(pack)
	q("no")
}

req_packs = c("Rcpp","devtools","usethis",
	"rmarkdown","knitr")
for(pp in req_packs){
	
	chk_pack = tryCatch(find.package(pp),
		error = function(ee){NULL},
		warning = function(ww){NULL})
	
	if( is.null(chk_pack) ) stop(sprintf("Install %s",pp))
	
	library(pp,character.only = TRUE)
}

compileAttributes(pkgdir = pack_dir)
document(pkg = pack_dir)
use_gpl3_license()
check_pandoc = pandoc_available(); check_pandoc
make_vign = check_pandoc

# Check: takes some time
chk = tryCatch(check(pkg = pack_dir,
	manual = TRUE,cran = TRUE,
	error_on = "note",
	vignettes = make_vign),
	error = function(ee){NULL},
	warning = function(ww){NULL})
chk

# Install
if( !is.null(chk) ){
	install(pack_dir,
		build_vignettes = make_vign,
		upgrade = FALSE)
}

###
