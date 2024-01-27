#include <RcppArmadillo.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <Rmath.h>
// #include <cmath>

// [[Rcpp::depends("RcppArmadillo")]]

template<typename T>
void printR_obj(const T& obj){
	Rcpp::Rcout << obj << std::endl;	
}

const double LOG2PI = log(2*arma::datum::pi);

// --------------------
// Intermediate Functions
// --------------------

// [[Rcpp::export]]
double Rcpp_norm(const arma::vec& a){
	return arma::norm(a);
}

// [[Rcpp::export]]
double Rcpp_max_abs_diff(const arma::vec& aa,const arma::vec& bb){
	return arma::max(arma::abs(aa - bb));
}

// [[Rcpp::export]]
double Rcpp_logSumExp(const arma::vec& log_x){
	if( log_x.n_elem == 1 ){
		return log_x.at(0);
	} else {
		double max_val = max(log_x);
		arma::vec log_x_2 = log_x - max_val;
		return log(arma::sum(arma::exp(log_x_2))) + max_val;
	}
}

// [[Rcpp::export]]
void prt_vec(const arma::vec& aa){
	arma::uword ii, nn = aa.n_elem;
	
	Rcpp::Rcout << "(";
	for(ii = 0; ii < nn; ii++){
		Rcpp::Rcout << aa.at(ii);
		if( ii < nn - 1 ) Rcpp::Rcout << ", ";
	}
	Rcpp::Rcout << ")\n";
	
}

// [[Rcpp::export]]
double log_CDF_weibull(const double& XX,
	const double& LAM,const double& ALP){
	
	double log_CDF = R::pweibull(XX,ALP,LAM,1,1);
	
	if( log_CDF == 0.0 && XX != arma::datum::inf ){
		// Precision fix: -Surv_Weibull
		log_CDF = -1.0 * R::pweibull(XX,ALP,LAM,0,0);
	}
	
	return log_CDF;
}

// [[Rcpp::export]]
arma::vec calc_expweibull_logCDF_logPDF(const double& XX,
	const double& LAM,const double& ALP,const double& KAP){
	
	double log_cdf_weibull = log_CDF_weibull(XX,LAM,ALP);
	double log_cdf_expweibull = log_cdf_weibull;
	double log_pdf_expweibull = R::dweibull(XX,ALP,LAM,1);
	
	if( KAP != 1.0 ){
		log_cdf_expweibull *= KAP;
		log_pdf_expweibull += std::log(KAP) + (KAP - 1.0) * log_cdf_weibull;
	}
	
	arma::vec out = arma::zeros<arma::vec>(2);
	out.at(0) = log_cdf_expweibull;
	out.at(1) = log_pdf_expweibull;
	
	return out;
	
}

// [[Rcpp::export]]
double log_SURV_expweibull(const double& XX,
	const double& LAM,const double& ALP,const double& KAP){
	
	// Calc log_CDF_expweibull
	double log_CDF = KAP * log_CDF_weibull(XX,LAM,ALP);
	
	// Calc SURV and check for precision problem
	double log_SURV = std::log(1.0 - std::exp(log_CDF));
	if( log_SURV == -arma::datum::inf ){
		// aka log(KAP) + log_SURV_weibull
		log_SURV = std::log(KAP) + R::pweibull(XX,ALP,LAM,0,1);
	}
	
	return log_SURV;
}

// --------------------
// New Functions

// [[Rcpp::export]]
double calc_copula(const arma::vec& log_CDFs,
	const std::string& copula,const double& THETA){
	
	double cop_CDF = 0.0;
	
	// A general property of copulas, if any CDF is 0, copula is 0
	if( log_CDFs.has_inf() )
		return 0.0;
	
	if( copula == "Independent" ){
		
		cop_CDF = std::exp(arma::sum(log_CDFs));
		return cop_CDF;
		
	} else if( copula == "Clayton" ){
		
		// Precision method
		arma::vec log_vec = log_CDFs;
		log_vec *= -1.0 * THETA;
		double log_mm = arma::max(log_vec), log_COP;
		log_COP = -1.0 / THETA *
			( log_mm + 
			std::log( arma::sum(arma::exp(log_vec - log_mm)) - 
			1.0 / std::exp(log_mm)) );
		cop_CDF = std::exp(log_COP);
		return cop_CDF;
		
	} else if( copula == "Gumbel" ){
		
		arma::vec log_vec = THETA * arma::log(-log_CDFs);
		cop_CDF = std::exp(1.0 / THETA * Rcpp_logSumExp(log_vec));
		cop_CDF = std::exp(-cop_CDF);
		return cop_CDF;
		
	} else {
		Rcpp::stop("Not a valid copula!");
	}
	
}

// [[Rcpp::export]]
double calc_copula_offset(const arma::vec& log_DENs,
	const arma::vec& log_CDFs,const std::string& copula,
	const double& THETA,const double& cop_CDF){
	
	if( cop_CDF == 0.0 ) return 0.0; /*Assumed to be true no matter the copula ...*/
	
	double offset = 0.0;
	arma::vec log_vec = arma::zeros<arma::vec>(2);
	
	if( copula == "Independent" ){
		
		// offset = D1 * F2 + D2 * F1;
		
		/*
		Question: Should this be zero if any F1, F2,... are zero?
		*/
		
		log_vec.at(0) = log_DENs.at(0) + log_CDFs.at(1);
		log_vec.at(1) = log_DENs.at(1) + log_CDFs.at(0);
		offset = arma::sum(arma::exp(log_vec));
		return offset;
		
	} else if( copula == "Clayton" ){
		
		// if( cop_CDF == 0.0 ) return 0.0;
		
		log_vec = -1.0 * THETA * log_CDFs;
		double log_mm = arma::max(log_vec);
		
		offset = (-1.0 / THETA - 1.0) *
			( log_mm + 
			std::log( arma::sum(arma::exp(log_vec - log_mm)) - 
			1.0 / std::exp(log_mm)) );
		
		log_vec = log_DENs - (THETA + 1.0) * log_CDFs;
		offset += Rcpp_logSumExp(log_vec);
		
		offset = std::exp(offset);
		return offset;
		
	} else if( copula == "Gumbel" ){
		
		// if( cop_CDF == 0.0 ) return 0.0;
		
		double /*nlog_F1 = -std::log(F1),
			nlog_F2 = -std::log(F2),*/ tmp_lognum;
		
		/*
		double log_nlog_F1 = std::log(nlog_F1);
		double log_nlog_F2 = std::log(nlog_F2);
		*/
		
		offset = std::log(cop_CDF);
		
		/*
		log_vec.at(0) = THETA * log_nlog_F1;
		log_vec.at(1) = THETA * log_nlog_F2;
		*/
		arma::vec log_nlog_CDFs = arma::log(-1.0 * log_CDFs);
		
		log_vec = THETA * log_nlog_CDFs;
		tmp_lognum = Rcpp_logSumExp(log_vec);
		
		offset += (1.0 / THETA - 1.0) * tmp_lognum;
		
		/*
		log_vec.at(0) = (THETA - 1.0) * log_nlog_F1 + std::log(D1) + nlog_F1;
		log_vec.at(1) = (THETA - 1.0) * log_nlog_F2 + std::log(D2) + nlog_F2;
		*/
		log_vec = (THETA - 1.0) * log_nlog_CDFs + log_DENs - log_CDFs;
		offset += Rcpp_logSumExp(log_vec);
		
		offset = std::exp(offset);
		return offset;
		
	} else {
		Rcpp::stop("Not a coded copula!");
	}
	
}

// [[Rcpp::export]]
arma::vec calc_copula_CDF_OFF(const arma::vec& log_DENs,
	const arma::vec& log_CDFs,const std::string& copula,
	const double& THETA){
	
	double offset = 0.0, 
		cop_CDF = calc_copula(log_CDFs,copula,THETA);
	arma::vec out = arma::zeros<arma::vec>(2);
	
	out.at(0) = cop_CDF;
	
	offset = calc_copula_offset(log_DENs,log_CDFs,
		copula,THETA,cop_CDF);
	
	out.at(1) = offset;
	
	return out;
}

double dMrs_LL(const arma::vec& XX,const arma::uvec& DELTA,
	const arma::vec& log_D2,const arma::vec& log_F2,const double& THETA,
	const double& ALPHA,const double& LAMBDA,const double& KAPPA,
	const std::string& copula,const bool& verb = false){
	
	arma::uword ii, NN = XX.n_elem;
	double LL = 0.0, LL_ii, offset, 
		cop_CDF, AA, BB, error_num = -999.0;
	arma::vec out = arma::zeros<arma::vec>(2),
		log_DENs = out, log_CDFs = out,
		DENs = out,CDFs = out,
		vEW = out;
	
	for(ii = 0; ii < NN; ii++){
		
		vEW = calc_expweibull_logCDF_logPDF(
			XX.at(ii),LAMBDA,ALPHA,KAPPA);
		
		log_CDFs.at(0) = vEW.at(0);
		log_CDFs.at(1) = log_F2.at(ii);
		CDFs = arma::exp(log_CDFs);
		
		log_DENs.at(0) = vEW.at(1);
		log_DENs.at(1) = log_D2.at(ii);
		DENs = arma::exp(log_DENs);
		
		out = calc_copula_CDF_OFF(log_DENs,
			log_CDFs,copula,THETA);
		
		if( out.has_nan() ){
			// if( verb ){
				// Rcpp::Rcout << "ii = " << ii + 1 << "; ";
				// Rcpp::Rcout << "Delta = " << DELTA.at(ii) << "; ";
				// Rcpp::Rcout << "D1 = " << DENs.at(0) << "; ";
				// Rcpp::Rcout << "D2 = " << DENs.at(1) << "; ";
				// Rcpp::Rcout << "F1 = " << CDFs.at(0) << "; ";
				// Rcpp::Rcout << "F2 = " << CDFs.at(1) << "; ";
				// Rcpp::Rcout << "cop_CDF = " << out.at(0) << "; ";
				// Rcpp::Rcout << "D1_D2 = " << out.at(1) << "\n";
			// }
			return error_num;
		}
		
		cop_CDF = out.at(0);
		offset = out.at(1);
		
		if( DELTA.at(ii) == 1 ){
			AA = DENs.at(0) + DENs.at(1) - offset;
			if( AA <= 0.0 ) return error_num;
			LL_ii = std::log( AA );
			
		} else {
			BB = 1.0 - CDFs.at(0) - CDFs.at(1) + cop_CDF;
			if( BB <= 0.0 ) return error_num;
			LL_ii = std::log( BB );
			
		}
		
		LL += LL_ii;
	}
	
	return LL;
	
}

// [[Rcpp::export]]
double dMrs_cLL(const arma::vec& XX,const arma::uvec& DELTA,
	const arma::vec& log_D2,const arma::vec& log_F2,const arma::vec& PARS,
	const std::string& copula,const bool& verb = false){
	
	// For Clayton: 0 <= THETA < infinity
	// For Gumbel: 1 <= THETA < infinity
	arma::vec ePARS = arma::exp(PARS);
	double KAPPA = ePARS.at(2), THETA = ePARS.at(3);
	// double KAPPA = 1.0 / (1.0 + std::exp(-PARS.at(2)));
	
	if( copula == "Gumbel" )
		THETA += 1.0;
	
	return dMrs_LL(XX,DELTA,log_D2,log_F2,THETA,ePARS.at(0),
		ePARS.at(1),KAPPA,copula,verb);
}

arma::vec dMrs_GRAD(const arma::vec& XX,const arma::uvec& DELTA,
	const arma::vec& log_D2,const arma::vec& log_F2,const double& THETA,
	const double& ALPHA,const double& LAMBDA,const double& KAPPA,
	const std::string& copula,const arma::vec& upPARS,
	const double& shift = 5e-6){
	
	arma::uword ii;
	double old_LL, new_LL, error_num = -999.0, new_theta;
	arma::vec GRAD = arma::zeros<arma::vec>(4),
		old_PARS = GRAD, new_PARS = GRAD;
	old_PARS.at(0) = std::log(ALPHA);
	old_PARS.at(1) = std::log(LAMBDA);
	old_PARS.at(2) = std::log(KAPPA);
	
	if( copula == "Gumbel" ){
		old_PARS.at(3) = std::log(THETA - 1.0);
	} else { // Clayton and Independent
		old_PARS.at(3) = std::log(THETA);
	}
	
	old_LL = dMrs_LL(XX,DELTA,log_D2,log_F2,THETA,ALPHA,
		LAMBDA,KAPPA,copula,false);
	if( old_LL == error_num ){
		GRAD.fill(error_num);
		return GRAD;
	}
	
	for(ii = 0; ii < 4; ii++){
		new_PARS = old_PARS;
		if( upPARS.at(ii) == 0.0 ) continue;
		new_PARS.at(ii) += shift;
		
		new_theta = std::exp(new_PARS.at(3));
		if( copula == "Gumbel" ){
			new_theta += 1.0;
		}
		
		new_LL = dMrs_LL(XX,DELTA,log_D2,log_F2,new_theta,
			std::exp(new_PARS.at(0)),
			std::exp(new_PARS.at(1)),
			std::exp(new_PARS.at(2)),
			copula,false);
		if( new_LL == error_num ){
			GRAD.fill(error_num);
			return GRAD;
		}
		
		if( new_LL == old_LL ){
			return dMrs_GRAD(XX,DELTA,log_D2,log_F2,
				THETA,ALPHA,LAMBDA,KAPPA,copula,
				upPARS,2*shift);
		}
		
		GRAD.at(ii) = (new_LL - old_LL) / shift;
	}
	
	return GRAD;
}

// [[Rcpp::export]]
arma::vec dMrs_cGRAD(const arma::vec& XX,const arma::uvec& DELTA,
	const arma::vec& log_D2,const arma::vec& log_F2,const arma::vec& PARS,
	const std::string& copula,const arma::vec& upPARS){
	
	arma::vec ePARS = arma::exp(PARS);
	double THETA = ePARS.at(3),
		KAPPA = ePARS.at(2), shift = 5e-6;
	
	if( copula == "Gumbel" ){
		THETA += 1.0;
	}
	
	return dMrs_GRAD(XX,DELTA,log_D2,log_F2,THETA,ePARS.at(0),
		ePARS.at(1),KAPPA,copula,upPARS,shift);
}

arma::mat dMrs_HESS(const arma::vec& XX,const arma::uvec& DELTA,
	const arma::vec& log_D2,const arma::vec& log_F2,const double& THETA,
	const double& ALPHA,const double& LAMBDA,const double& KAPPA,
	const std::string& copula,const arma::vec& upPARS){
	
	arma::uword ii, np = 4;
	double shift = 5e-6, error_num = -999.0, new_theta;
	arma::vec old_GRAD = dMrs_GRAD(XX,DELTA,log_D2,log_F2,
		THETA,ALPHA,LAMBDA,KAPPA,copula,upPARS,shift),
		PARS = arma::zeros<arma::vec>(np),
		PARS_2 = PARS, tmp_vec = PARS;
	arma::mat HESS = arma::zeros<arma::mat>(np,np),
		I_np = arma::eye<arma::mat>(np,np);
	
	if( old_GRAD.at(0) == error_num ){
		HESS.fill(error_num);
		return HESS;
	}
	
	PARS.at(0) = std::log(ALPHA);
	PARS.at(1) = std::log(LAMBDA);
	PARS.at(2) = std::log(KAPPA);
	
	if( copula == "Gumbel" ){
		PARS.at(3) = std::log(THETA - 1.0);
	} else {
		PARS.at(3) = std::log(THETA);
	}
	
	for(ii = 0; ii < np; ii++){
		if( upPARS.at(ii) == 0.0 ) continue;
		PARS_2 = PARS + shift * I_np.col(ii) % upPARS;
		
		new_theta = std::exp(PARS_2.at(3));
		if( copula == "Gumbel" ){
			new_theta += 1.0;
		}
		
		tmp_vec = dMrs_GRAD(XX,DELTA,log_D2,log_F2,
			new_theta,std::exp(PARS_2.at(0)),
			std::exp(PARS_2.at(1)),std::exp(PARS_2.at(2)),
			copula,upPARS,shift);
		if( arma::any(tmp_vec == error_num) ){
			HESS.fill(error_num);
			return HESS;
		}
		tmp_vec = tmp_vec - old_GRAD;
		tmp_vec /= shift;
		HESS(arma::span(ii,np - 1),ii) = tmp_vec.subvec(ii,np - 1);
		if( ii < np - 1 ){
			HESS(ii,arma::span(ii + 1,np - 1)) = tmp_vec.subvec(ii + 1,np - 1).t();
		}
	}
	
	return HESS;
}

// [[Rcpp::export]]
arma::mat dMrs_cHESS(const arma::vec& XX,const arma::uvec& DELTA,
	const arma::vec& log_D2,const arma::vec& log_F2,const arma::vec& PARS,
	const std::string& copula,const arma::vec& upPARS){
	
	arma::vec ePARS = arma::exp(PARS);
	double THETA = ePARS.at(3);
	// double KAPPA = 1.0 / (1.0 + std::exp(-PARS.at(2)));
	double KAPPA = ePARS.at(2);
	
	if( copula == "Gumbel" ){
		THETA += 1.0;
	}
	
	return dMrs_HESS(XX,DELTA,log_D2,log_F2,THETA,ePARS.at(0),
		ePARS.at(1),KAPPA,copula,upPARS);
}

// [[Rcpp::export]]
void dMrs_NR(const arma::vec& XX,const arma::uvec& DELTA,
	const arma::vec& log_D2,const arma::vec& log_F2,arma::vec& PARS,
	const std::string& copula,const arma::vec& upPARS,
	const arma::uword& max_iter = 2e2,const double& eps = 5e-2,
	const arma::uword& mult = 5,const bool& verb = true){
	
	arma::uword iter = 0, jj, kk, uu,
		np = PARS.n_elem, reach = 0;
	
	// Initialize parameters
	if( verb ){
		Rcpp::Rcout << "iPARS = ";
		prt_vec(PARS);
	}
	arma::mat I_np = arma::eye<arma::mat>(np,np),
		HESS = I_np, iHESS = I_np;
	arma::vec xk = PARS, curr_xk = arma::zeros<arma::vec>(np),
		new_xk = curr_xk, GRAD = curr_xk, pk_NR = curr_xk,
		pk_GD = curr_xk;
	double error_num = -999.0, diff_LL = 0.0, 
		rcond_num, diff_PARS,
		orig_LL = dMrs_cLL(XX,DELTA,log_D2,log_F2,xk,copula,false),
		nGRAD, old_LL, new_LL;
	arma::uvec chk = arma::zeros<arma::uvec>(np),
		idx_fin = arma::find_finite(PARS);
	
	old_LL = orig_LL;
	while( iter < max_iter ){
		GRAD = dMrs_cGRAD(XX,DELTA,log_D2,log_F2,xk,copula,upPARS);
		if( GRAD.at(0) == error_num ){
			if( verb ) Rcpp::Rcout << "Invalid pars\n";
			return;
		}
		
		HESS = dMrs_cHESS(XX,DELTA,log_D2,log_F2,xk,copula,upPARS);
		if( HESS.at(0,0) == error_num ){
			if( verb ) Rcpp::Rcout << "Invalid pars\n";
			return;
		}
		
		chk.zeros();
		chk(arma::find(upPARS == 1.0 && HESS.diag() == 0.0)).ones();
		if( arma::any(chk == 1) ){
			if( verb ) Rcpp::Rcout << "Variance issue\n";
			return;
		}
		arma::uvec nz = arma::find(HESS.diag() != 0.0);
		rcond_num = arma::rcond(HESS.submat(nz,nz));
		if( rcond_num == 0.0 ){
			if( verb ) Rcpp::Rcout << "Variance issue, rcond\n";
			return;
		}
		iHESS.zeros();
		iHESS.submat(nz,nz) = arma::inv(-1.0 * HESS.submat(nz,nz));
		
		pk_NR = iHESS * GRAD;
		pk_NR /= std::max(1.0, Rcpp_norm(pk_NR));
		
		pk_GD = GRAD;
		pk_GD /= std::max(1.0, Rcpp_norm(pk_GD));
		
		uu = 0;
		for(kk = 0; kk < 1; kk++){
		for(jj = 0; jj <= 30; jj++){
			
			if( kk == 0 ) new_xk = xk + pk_NR / std::pow(4.0,jj);
			if( kk == 1 ) new_xk = xk + pk_GD / std::pow(4.0,jj);
			
			new_LL = dMrs_cLL(XX,DELTA,log_D2,log_F2,new_xk,copula,false);
			if( new_LL == error_num ) continue;
			if( new_LL <= old_LL ) continue;
			
			GRAD = dMrs_cGRAD(XX,DELTA,log_D2,log_F2,new_xk,copula,upPARS);
			if( arma::any(GRAD == error_num) ) continue;
			
			HESS = dMrs_cHESS(XX,DELTA,log_D2,log_F2,new_xk,copula,upPARS);
			if( HESS.at(0,0) == error_num ) continue;
			
			chk.zeros();
			chk(arma::find(upPARS == 1.0 && HESS.diag() == 0.0)).ones();
			if( arma::any(chk == 1) ) continue;
			
			arma::uvec nz = arma::find(HESS.diag() != 0.0);
			rcond_num = arma::rcond(HESS.submat(nz,nz));
			if( rcond_num == 0.0 ) continue;
			iHESS.zeros();
			iHESS.submat(nz,nz) = arma::inv(-1.0 * HESS.submat(nz,nz));
			if( arma::any(iHESS.diag() < 0.0) ) continue;
			
			uu = 1;
			break;
			
		}
			
			if( uu == 1 ) break;
			
		}
		
		if( uu == 0 ){
			if( verb ) Rcpp::Rcout << "No more update\n";
			break;
		}
		
		diff_LL = new_LL - old_LL;
		diff_PARS = Rcpp_norm(xk(idx_fin) - new_xk(idx_fin));
		nGRAD = Rcpp_norm(GRAD);
		if( diff_LL == 0.0 || diff_PARS == 0.0 ) break;
		
		if( diff_LL < eps * 1e-1 && diff_PARS < eps * 1e-1 ){
			if( nGRAD < eps ){
				reach++;
			} else {
				reach = 0; // resets if conditions not met
			}
			
			if( reach >= 15 ){
				if( verb ) Rcpp::Rcout << "Optimization criteria met\n";
				break;
			}
			
		} else {
			reach = 0;
		}
		
		if( verb ){
			if( (iter + 1) % mult == 0 ){
				Rcpp::Rcout << "iter=" << iter + 1 << "; LL=" << old_LL
					<< "; diff.LL=" << diff_LL << "; diff.PARS=" << diff_PARS
					<< "; nGRAD=" << nGRAD << "; meth=" << kk 
					<< "; reach=" << reach << "; PARS = ";
					prt_vec(new_xk);
			}
		}
		
		xk = new_xk;
		old_LL = new_LL;
		iter++;
		
	}
	
	old_LL = dMrs_cLL(XX,DELTA,log_D2,log_F2,xk,copula,false);
	if( old_LL > orig_LL ){ // Criteria for convergence
		PARS = xk;
	}

	if( verb ){
		Rcpp::Rcout << "####\nNum Iter = " << iter+1 << "\n";
		Rcpp::Rcout << "Params = "; prt_vec(xk);
		Rcpp::Rcout << "LL = " << old_LL << "\n";
		GRAD = dMrs_cGRAD(XX,DELTA,log_D2,log_F2,xk,copula,upPARS);
		HESS = dMrs_cHESS(XX,DELTA,log_D2,log_F2,PARS,copula,upPARS);
		arma::uvec nz = arma::find(HESS.diag() != 0.0);
		iHESS.submat(nz,nz) = arma::inv(-1.0 * HESS.submat(nz,nz));
		Rcpp::Rcout << "GRAD = "; prt_vec(GRAD);
		Rcpp::Rcout << "Convergence Indicators: \n"
			<< "   NormGrad = " << Rcpp_norm(GRAD) << "\n"
			<< "   NormIHessGrad = " << Rcpp_norm(iHESS * GRAD) << "\n";
		Rcpp::Rcout << "Var = "; prt_vec(iHESS.diag());
	}
	
}

// [[Rcpp::export]]
arma::mat dMrs_GRID(const arma::vec& XX,const arma::uvec& DELTA,
	const arma::vec& log_D2,const arma::vec& log_F2,const arma::vec& log_THETA,
	const arma::vec& log_ALPHA,const arma::vec& log_LAMBDA,
	const arma::vec& unc_KAPPA,const std::string& copula,
	const bool& verb = false,const int& ncores = 1){
	
	arma::uword num_alpha = log_ALPHA.n_elem, num_lambda = log_LAMBDA.n_elem,
		num_kappa = unc_KAPPA.n_elem, num_theta = log_THETA.n_elem,
		tot = num_alpha * num_lambda * num_kappa * num_theta,
		aa, ll, kk, tt, cnt = 0;
	// log_ALPHA,log_LAMBDA,unc_KAPPA,log_THETA,LL
	arma::mat DAT = arma::zeros<arma::mat>(tot,5);
	// double error_num = -999.0;
	bool verb2 = verb && ncores == 1;
	if( verb ){
		Rcpp::Rcout << "#ALPHA grid points = " << num_alpha << "\n";
		Rcpp::Rcout << "#LAMBDA grid points = " << num_lambda << "\n";
		Rcpp::Rcout << "#KAPPA grid points = " << num_kappa << "\n";
		Rcpp::Rcout << "#THETA grid points = " << num_theta << "\n";
		Rcpp::Rcout << "Num grid points = " << tot << "\n";
	}
	
	for(aa = 0; aa < num_alpha; aa++){
	for(ll = 0; ll < num_lambda; ll++){
	for(kk = 0; kk < num_kappa; kk++){
	for(tt = 0; tt < num_theta; tt++){
		
		DAT.at(cnt,0) = log_ALPHA.at(aa);
		DAT.at(cnt,1) = log_LAMBDA.at(ll);
		DAT.at(cnt,2) = unc_KAPPA.at(kk);
		DAT.at(cnt,3) = log_THETA.at(tt);
		
		cnt++;
	}}}}
	
	#ifdef _OPENMP
	# pragma omp parallel for schedule(dynamic) \
		num_threads(ncores) \
		shared(verb2,copula,tot,XX,DELTA,log_D2,log_F2,DAT)
	#endif
	for(arma::uword cnt2 = 0; cnt2 < tot; cnt2++){
		if( verb2 ){
			if( (cnt2 + 1) % 500 == 0 )
				Rcpp::Rcout << ".";
			if( (cnt2 + 1) % 5000 == 0 || (cnt2 + 1) == tot )
				Rcpp::Rcout << " " << (cnt2 + 1);
			if( (cnt2 + 1) % 10000 == 0 || (cnt2 + 1) == tot )
				Rcpp::Rcout << " out of " << tot << " done\n";
		}
		
		arma::vec PARS = DAT(cnt2,arma::span(0,3)).t();
		DAT.at(cnt2,4) = dMrs_cLL(XX,DELTA,log_D2,log_F2,PARS,copula,false);
	}
	
	// Get parameters with largest log-likelihood
	// double min_val = arma::min(DAT.col(4));
	// arma::uvec uDAT = arma::zeros<arma::uvec>(DAT.n_rows);
	// uDAT = DAT.col(4) == error_num;
	// if( min_val != error_num ){
		// DAT.col(4) = DAT.col(4) - std::abs(min_val) * arma::conv_to<arma::vec>::from(uDAT);
	// }
	
	return DAT;
}

// [[Rcpp::export]]
arma::mat dMrs_MATCH(const arma::mat& wDAT,
	const arma::mat& rDAT,const int& ncores = 1,
	const bool& verb = true){
	// wDAT columns: age, yr_diag, sex, obs_time
	// rDAT columns: Year, age, qx(hazard), sex
	
	arma::uword nn = wDAT.n_rows;
	arma::mat OUT = arma::zeros<arma::mat>(nn,2); // output log_dens_t2, log_cdf_t2
	bool verb2 = verb && ncores == 1;
	
	#ifdef _OPENMP
	# pragma omp parallel for schedule(dynamic) \
		num_threads(ncores) \
		shared(verb2,nn,wDAT,rDAT,OUT)
	#endif
	for(arma::uword ii = 0; ii < nn; ii++){
		
		if( verb2 ){
			if( (ii + 1) % 500 == 0 )
				Rcpp::Rcout << ".";
			if( (ii + 1) % 5000 == 0 || (ii + 1) == nn )
				Rcpp::Rcout << (ii + 1) << " out of " << nn << " done\n";
		}
		
		double 
			age 			= wDAT.at(ii,0),
			yr_diag 	= wDAT.at(ii,1),
			sex 			= wDAT.at(ii,2),
			obs_time 	= wDAT.at(ii,3);
		
		arma::uvec idx = arma::find( rDAT.col(1) >= age
			&& rDAT.col(0) >= yr_diag
			&& rDAT.col(3) == sex
			&& ( rDAT.col(0) - rDAT.col(1) ) == ( yr_diag - age )
			&& ( rDAT.col(0) - yr_diag ) <= obs_time
			);
			
		if( idx.n_elem == 0 ) Rcpp::stop("check this");
		
		arma::uword nr = idx.n_elem;
		
		arma::mat tmp_mat = arma::zeros<arma::mat>(nr,3); // t1, t2, qx
		arma::mat rDAT2 = rDAT.rows(idx);
		
		// calc t1
		tmp_mat.col(0) = rDAT2.col(0) - yr_diag;
		
		// calc t2
		if( nr > 1 ) 
				tmp_mat(arma::span(0,nr - 2),1) = tmp_mat(arma::span(1,nr - 1),0);
		
		// set qx
		tmp_mat.col(2) = rDAT2.col(2);
		
		double haz = tmp_mat.at(nr - 1,2), cumHaz = 0.0;
		
		if( nr > 1 ){
			cumHaz = arma::dot(tmp_mat(arma::span(0,nr - 2),2),
				tmp_mat(arma::span(0,nr - 2),1) - tmp_mat(arma::span(0,nr - 2),0)) +
				haz * (obs_time - tmp_mat.at(nr - 1,0));
		} else {
			cumHaz = haz * obs_time;
		}
		
		double log_surv_t2 = -1.0 * cumHaz;
		double log_dens_t2 = std::log(haz) + log_surv_t2;
		double log_cdf_t2 = std::log(1.0 - std::exp(log_surv_t2));
		
		if( log_cdf_t2 == -arma::datum::inf ){
			log_cdf_t2 = -1.0 * std::exp(log_surv_t2);
		}
		
		// log_cdf_t2
		OUT.at(ii,1) = log_cdf_t2;
		
		// log_dens_t2
		OUT.at(ii,0) = log_dens_t2;
		
	}
	
	return OUT;
}

