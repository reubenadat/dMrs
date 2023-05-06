# ----------
# Smarter functions
# ----------
smart_progress = function(ii,nn,string = ".",
	iter = 5,iter2 = 2e2,...){
	
	if(ii %% iter == 0)
		cat(string,...)
	
	if(ii %% iter2 == 0 || ii == nn)
		cat(sprintf("%s out of %s\n",ii,nn),...)
	
}
smart_table = function(...){
	table(...,useNA = 'ifany')
}
smart_df = function(...){
	data.frame(...,stringsAsFactors = FALSE)
}
smart_digits = function(x,digits = 2){
	sprintf(paste0("%.",digits,"f"),round(x,digits))
}
format_latex = function(INPUT){
	# INPUT = "optE_AIC%"
	
	if( length(grep("^\\$",INPUT)) == 1 && length(grep("\\$$",INPUT)) == 1 ){
		return(INPUT)
	}
	
	INPUT2 = gsub("%","\\\\%",INPUT)
	INPUT2 = gsub("_","\\\\_",INPUT2)
	INPUT2
}
clean_repeats = function(VEC){
	if(FALSE){
		VEC = c(rep("a",2),rep("b",2),"a","c")
		VEC
	}
	
	curr_string = NA
	for(ii in seq(length(VEC))){
		# ii = 1
		if( ii == 1 ){
			curr_string = VEC[ii]
		} else {
			if( VEC[ii] == curr_string ){
				VEC[ii] = ""
			} else {
				curr_string = VEC[ii]
			}
		}
	}
	
	VEC
}
print_latex_table = function(DATA,repeat_VARS = NULL,
	my_align = NULL,add_table = FALSE,fontsize = NULL,
	caption = NULL,label = NULL,midrule1 = NULL,
	latex_comment = NULL,...){
	
	orig_names = colnames(DATA)
	
	if( nrow(DATA) > 1 ){
		DATA = smart_df(apply(DATA,2,as.character))
	} else {
		DATA = smart_df(t(apply(DATA,2,as.character)))
	}
	
	if( !is.null(repeat_VARS) && length(repeat_VARS) > 0 ){
		# loop thru vector(column) to find repeats and replace with ""
		tmp_index = which(orig_names %in% repeat_VARS)
		DATA[,tmp_index] = apply(DATA[,tmp_index,drop=FALSE],2,clean_repeats)
	}
	
	prep_DATA = DATA
	
	cat("\n",...)
	
	if( !is.null(latex_comment) ){
		cat(sprintf("%% %s\n",latex_comment),...)
	}
	
	if( add_table ){
		cat(paste0("\\begin{table}[!htbp] \n\\centering\n"),...)
		if( !is.null(fontsize) )
			cat(paste0("\\",fontsize,"\n"),...)
		else
			cat(paste0("\\normalsize\n"),...)
		if( !is.null(caption) ){
			caption = gsub("\n","",caption)
			cat(paste0("\\caption{",caption,"}\n"),...)
		}
		if( !is.null(label) ) cat(paste0("\\label{tab:",label,"}\n"),...)
	}
	
	if( is.null(my_align) ){
		cat(paste0("\\begin{tabular}{l",
			paste(rep("c",ncol(prep_DATA)-1),collapse=""),"}\n"),...)
	} else {
		cat(paste0("\\begin{tabular}{",my_align,"}\n"),...)
	}
	
	cat("\\toprule\n",...)
	cat(paste0(paste(sapply(orig_names,format_latex),collapse=" & ")," \\\\\n"),...)
	
	if( is.null(midrule1) ){
		cat("\\midrule\n",...)
	} else {
		cat(paste0(midrule1,"\n"),...)
	}
	apply(prep_DATA,1,function(x)
		cat(paste0(paste(sapply(x,format_latex),
			collapse=" & ")," \\\\\n"),...))
	cat("\\bottomrule\n\\end{tabular}\n",...)

	if( add_table ){
		cat(paste0("\\end{table}\n"),...)
	}

	cat("\n",...)
	
}
smart_mkdir = function(input_dir){
	
	if( !file.exists(input_dir) || !dir.exists(input_dir) )
		dir.create(path = input_dir,recursive = TRUE)
	
}
name_change = function(DATA,ORIG_NAME,NEW_NAME){
	
	old_idx = which(colnames(DATA) == ORIG_NAME)
	new_idx = which(colnames(DATA) == NEW_NAME)
	if( length(new_idx) > 0 ){
		return(DATA)
	} else if( length(old_idx) > 0 ){
		colnames(DATA)[old_idx] = NEW_NAME
		return(DATA)
	} else {
		stop(sprintf("ORIG_NAME = %s missing",ORIG_NAME))
	}
	
}
make_menu = function(PROMPT,OPTS){
	
	if( missing(PROMPT) )
		PROMPT = "Select an option"
	if( missing(OPTS) )
		stop("Add a vector of options")
	
	INDENT = "   "
	cmd = sprintf("%s:",PROMPT)
	vec_seq = seq(length(OPTS))
	for( ii in vec_seq ){
		cmd = sprintf("%s\n%s%s) %s",cmd,INDENT,ii,OPTS[ii])
	}
	cmd = sprintf("%s\n%s> ",cmd,INDENT)
	
	while(TRUE){
		resp = readline(prompt = cmd)
		if( resp %in% vec_seq )
			break
		cat("Not an option, try again\n")
	}
	
	resp = as.integer(resp)
	return(OPTS[resp])
	
}

###

