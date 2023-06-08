# Steps to create/check/install package from directory

rm(list = ls())
bb = strsplit(getwd(),"/")[[1]]
pack_dir = paste(bb[-length(bb)],collapse = "/")
pack = strsplit(pack_dir,"/")[[1]]
pack = pack[length(pack)]
pack

chk_pack = tryCatch(find.package(pack),
	error = function(ee){NULL}); chk_pack

if( !is.null(chk_pack) ){
	remove.packages(pack)
	q("no")
}

Rcpp::compileAttributes(pkgdir = pack_dir)
devtools::document(pkg = pack_dir)
usethis::use_gpl3_license()
# Sys.setenv("RSTUDIO_PANDOC" = "C:/Program Files/RStudio/bin/pandoc")
# check_pandoc = rmarkdown::pandoc_available(); check_pandoc
make_vign = FALSE

# Check: takes some time
devtools::check(pkg = pack_dir,manual = TRUE,
	cran = TRUE,error_on = c("warning","note")[2],
	vignettes = make_vign)

# Install
devtools::install(pack_dir,
	build_vignettes = make_vign,
	upgrade = FALSE)

###
