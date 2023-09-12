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
arma::vec calc_expweibull_CDF_PDF(const double& XX,
	const double& LAM,const double& ALP,const double& KAP){
	
	double log_cdf_weibull = R::pweibull(XX,ALP,LAM,1,1);
	double log_cdf_expweibull = log_cdf_weibull;
	double log_pdf_expweibull = R::dweibull(XX,ALP,LAM,1);
	
	if( KAP != 1.0 ){
		log_cdf_expweibull *= KAP;
		log_pdf_expweibull += std::log(KAP) + (KAP - 1.0) * log_cdf_weibull;
	}
	
	arma::vec out = arma::zeros<arma::vec>(2);
	out.at(0) = std::exp(log_cdf_expweibull);
	out.at(1) = std::exp(log_pdf_expweibull);
	
	return out;
	
}


// --------------------
// New Functions

// [[Rcpp::export]]
double calc_copula(const double& F1,const double& F2,
	const std::string& copula,const double& THETA){
	
	double F_T1_T2 = 0.0;
	
	// A general property of copulas
	if( F1 == 0.0 || F2 == 0.0 )
		return 0.0;
	
	// A general property of copulas
	if( F1 == 1.0 ) return F2;
	if( F2 == 1.0 ) return F1;
	
	// if( F1 == 1.0 && F2 == 1.0 )
		// return 1.0;
	
	if( copula == "Independent" ){
		
		F_T1_T2 = F1 * F2;
		
	} else if( copula == "Clayton" ){
		
		// Precision method
		arma::vec log_CDFs = { std::log(F1), std::log(F2) };
		log_CDFs *= -1.0 * THETA;
		double log_mm = arma::max(log_CDFs), log_COP;
		log_COP = -1.0 / THETA *
			( log_mm + 
			std::log( arma::sum(arma::exp(log_CDFs - log_mm)) - 
			1.0 / std::exp(log_mm)) );
		F_T1_T2 = std::exp(log_COP);
		
		// Non-precision method
		// F_T1_T2 = std::pow(F1,-THETA) + std::pow(F2,-THETA) - 1.0;
		// F_T1_T2 = std::pow(F_T1_T2,-1.0 / THETA);
		
	} else if( copula == "Gumbel" ){
		
		arma::vec log_vec = arma::zeros<arma::vec>(2);
		log_vec.at(0) = THETA * std::log(-std::log(F1));
		log_vec.at(1) = THETA * std::log(-std::log(F2));
		
		F_T1_T2 = std::exp(1.0 / THETA * Rcpp_logSumExp(log_vec));
		F_T1_T2 = std::exp(-F_T1_T2);
		
	} else {
		Rcpp::stop("Not a valid copula!");
	}
	
	if( F_T1_T2 < 0 )
		Rcpp::stop("Negative copula detected!");
	
	return F_T1_T2;
}

// [[Rcpp::export]]
double calc_copula_dens(const double& D1,const double& D2,
	const double& F1,const double& F2,
	const std::string& copula,const double& THETA,
	const double& F_T1_T2){
	
	double f_T1_T2 = 0.0, nlog_F1, nlog_F2,
		log_PDF_1, log_PDF_2;
	
	arma::vec log_vec = arma::zeros<arma::vec>(2);
	
	if( copula == "Independent" ){
		
		f_T1_T2 = D1 * F2 + D2 * F1;
	
	} else if( copula == "Clayton" ){
		
		if( F_T1_T2 == 0.0 ) return 0.0;
		
		arma::vec log_CDFs = { std::log(F1), std::log(F2) };
		log_CDFs *= -1.0 * THETA;
		double log_mm = arma::max(log_CDFs);
		log_PDF_1 = (-1.0 / THETA - 1.0) *
			( log_mm + 
			std::log( arma::sum(arma::exp(log_CDFs - log_mm)) - 
			1.0 / std::exp(log_mm)) );
		
		log_vec.at(0) = std::log(D1) - (THETA + 1.0) * std::log(F1);
		log_vec.at(1) = std::log(D2) - (THETA + 1.0) * std::log(F2);
		log_PDF_2 = Rcpp_logSumExp(log_vec);
		
		f_T1_T2 = std::exp(log_PDF_1 + log_PDF_2);
		
	} else if( copula == "Gumbel" ){
		
		if( F_T1_T2 == 0.0 ) return 0.0;
		
		// log first term
		f_T1_T2 = std::log(F_T1_T2);
		
		nlog_F1 = -std::log(F1);
		nlog_F2 = -std::log(F2);
		
		// plus log second term
		log_vec.at(0) = THETA * std::log(nlog_F1);
		log_vec.at(1) = THETA * std::log(nlog_F2);
		f_T1_T2 += (1.0/THETA - 1) * Rcpp_logSumExp(log_vec);
		
		// plus log third term
		log_vec.at(0) = (THETA - 1.0) * std::log(nlog_F1) + std::log(D1) + nlog_F1;
		log_vec.at(1) = (THETA - 1.0) * std::log(nlog_F2) + std::log(D2) + nlog_F2;
		f_T1_T2 += Rcpp_logSumExp(log_vec);
		
		// finally exponentiate
		f_T1_T2 = std::exp(f_T1_T2);
		
	} else {
		Rcpp::stop("Not a valid copula!");
	}
	
	return f_T1_T2;
}

// [[Rcpp::export]]
arma::vec calc_copula_CDF_PDF(const double& D1,const double& D2,
	const double& F1,const double& F2,
	const std::string& copula,const double& THETA){
	
	double f_T1_T2 = 0.0, 
		F_T1_T2 = calc_copula(F1,F2,copula,THETA);
	arma::vec out = arma::zeros<arma::vec>(2),
		log_P = out;
	out.at(0) = F_T1_T2;
	
	/*
	if( F_T1_T2 == 0.0 || F_T1_T2 == 1.0 ){
		out.at(1) = 0.0;
		return out;
	}
	*/
	
	f_T1_T2 = calc_copula_dens(D1,D2,F1,F2,copula,THETA,F_T1_T2);
	
	out.at(1) = f_T1_T2;
	
	return out;
}

double dMrs_LL(const arma::vec& XX,const arma::uvec& DELTA,
	const arma::vec& D2,const arma::vec& S2,const double& THETA,
	const double& ALPHA,const double& LAMBDA,const double& KAPPA,
	const std::string& copula,const bool& verb = false){
	
	arma::uword ii, NN = XX.n_elem;
	double LL = 0.0, LL_ii, D1, S1, F1, F2,
		// log_F1, log_F2, 
		D1_D2, 
		// XDL, enXDLa, 
		// U1T_U2T, U1T1, U2T1,
		F1_F2, AA, BB, error_num = -999.0;
		// KAL = KAPPA * ALPHA / LAMBDA
	arma::vec out = arma::zeros<arma::vec>(2), 
		vEW = out;
	
	for(ii = 0; ii < NN; ii++){
		
		vEW = calc_expweibull_CDF_PDF(XX.at(ii),LAMBDA,ALPHA,KAPPA);
		
		// XDL = XX.at(ii) / LAMBDA;
		// enXDLa = std::exp(-std::pow(XDL,ALPHA));
		
		// F1 = 1.0 - enXDLa;
		// if( KAPPA != 1.0 ) F1 = std::pow(F1,KAPPA);
		
		F1 = vEW.at(0);
		D1 = vEW.at(1);
		S1 = 1.0 - F1;
		
		// D1 = KAL * std::pow(XDL,ALPHA - 1.0) * enXDLa;
		// if( KAPPA != 1.0 ) D1 *= F1 / ( 1.0 - enXDLa );
		// if( F1 == 0.0 ) D1 = 0.0;
		
		F2 = 1.0 - S2.at(ii);
		
		out = calc_copula_CDF_PDF(D1,D2.at(ii),
			F1,F2,copula,THETA);
		
		if( out.has_nan() ){
			// if( verb ){
				// Rcpp::Rcout << "ii = " << ii + 1 << "; ";
				// Rcpp::Rcout << "Delta = " << DELTA.at(ii) << "; ";
				// Rcpp::Rcout << "D1 = " << D1 << "; ";
				// Rcpp::Rcout << "D2 = " << D2.at(ii) << "; ";
				// Rcpp::Rcout << "F1 = " << F1 << "; ";
				// Rcpp::Rcout << "F2 = " << F2 << "; ";
				// Rcpp::Rcout << "F1_F2 = " << out.at(0) << "; ";
				// Rcpp::Rcout << "D1_D2 = " << out.at(1) << "\n";
			// }
			return error_num;
		}
		
		F1_F2 = out.at(0);
		D1_D2 = out.at(1);
		
		if( DELTA.at(ii) == 1 ){
			AA = D1 + D2.at(ii) - D1_D2;
			if( AA <= 0.0 ) return error_num;
			LL_ii = std::log( AA );
		} else {
			BB = S1 + S2.at(ii) - 1.0 + F1_F2;
			if( BB <= 0.0 ) return error_num;
			LL_ii = std::log( BB );
		}
		
		/*
		if( verb ){
			Rcpp::Rcout << "ii = " << ii+1 << "; ";
			Rcpp::Rcout << "Delta = " << DELTA.at(ii) << "; ";
			Rcpp::Rcout << "D1 = " << D1 << "; ";
			Rcpp::Rcout << "F1 = " << F1 << "; ";
			Rcpp::Rcout << "U1T1 = " << U1T1 << "; ";
			Rcpp::Rcout << "U2T1 = " << U2T1 << "\n";
			Rcpp::Rcout << "   D2 = " << D2.at(ii) << "; ";
			Rcpp::Rcout << "F2 = " << F2 << "; ";
			Rcpp::Rcout << "U1T_U2T = " << U1T_U2T << "; ";
			Rcpp::Rcout << "F1_F2 = " << F1_F2 << "; ";
			Rcpp::Rcout << "D1_D2 = " << D1_D2 << "\n";
			Rcpp::Rcout << "   AA = " << AA << "; ";
			Rcpp::Rcout << "BB = " << BB << "; ";
			Rcpp::Rcout << "LL_ii = " << LL_ii << "; ";
			Rcpp::Rcout << "LL = " << LL << "\n";
		}
		*/
		
		LL += LL_ii;
	}
	
	return LL;
	
}

// [[Rcpp::export]]
double dMrs_cLL(const arma::vec& XX,const arma::uvec& DELTA,
	const arma::vec& D2,const arma::vec& S2,const arma::vec& PARS,
	const std::string& copula,const bool& verb = false){
	
	// For Clayton: 0 <= THETA < infinity
	// For Gumbel: 1 <= THETA < infinity
	arma::vec ePARS = arma::exp(PARS);
	double THETA = ePARS.at(3);
	// double KAPPA = 1.0 / (1.0 + std::exp(-PARS.at(2)));
	double KAPPA = ePARS.at(2);
	
	if( copula == "Gumbel" ){
		THETA += 1.0;
		/* 
		For Gumbel: theta >= 1
		For Clayton: theta >= 0
		*/
	}
	
	return dMrs_LL(XX,DELTA,D2,S2,THETA,ePARS.at(0),
		ePARS.at(1),KAPPA,copula,verb);
}

arma::vec dMrs_GRAD(const arma::vec& XX,const arma::uvec& DELTA,
	const arma::vec& D2,const arma::vec& S2,const double& THETA,
	const double& ALPHA,const double& LAMBDA,const double& KAPPA,
	const std::string& copula,const arma::vec& upPARS){
	
	arma::uword ii;
	double shift = 5e-6, old_LL, new_LL, error_num = -999.0, new_theta;
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
	
	old_LL = dMrs_LL(XX,DELTA,D2,S2,THETA,ALPHA,
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
		
		new_LL = dMrs_LL(XX,DELTA,D2,S2,new_theta,
			std::exp(new_PARS.at(0)),
			std::exp(new_PARS.at(1)),
			std::exp(new_PARS.at(2)),
			copula,false);
		if( new_LL == error_num ){
			GRAD.fill(error_num);
			return GRAD;
		}
		GRAD.at(ii) = (new_LL - old_LL) / shift;
	}
	
	return GRAD;
}

// [[Rcpp::export]]
arma::vec dMrs_cGRAD(const arma::vec& XX,const arma::uvec& DELTA,
	const arma::vec& D2,const arma::vec& S2,const arma::vec& PARS,
	const std::string& copula,const arma::vec& upPARS){
	
	arma::vec ePARS = arma::exp(PARS);
	double THETA = ePARS.at(3);
	// double KAPPA = 1.0 / (1.0 + std::exp(-PARS.at(2)));
	double KAPPA = ePARS.at(2);
	
	if( copula == "Gumbel" ){
		THETA += 1.0;
	}
	
	return dMrs_GRAD(XX,DELTA,D2,S2,THETA,ePARS.at(0),
		ePARS.at(1),KAPPA,copula,upPARS);
}

arma::mat dMrs_HESS(const arma::vec& XX,const arma::uvec& DELTA,
	const arma::vec& D2,const arma::vec& S2,const double& THETA,
	const double& ALPHA,const double& LAMBDA,const double& KAPPA,
	const std::string& copula,const arma::vec& upPARS){
	
	arma::uword ii, np = 4;
	double shift = 5e-6, error_num = -999.0, new_theta;
	arma::vec old_GRAD = dMrs_GRAD(XX,DELTA,D2,S2,
		THETA,ALPHA,LAMBDA,KAPPA,copula,upPARS),
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
		
		tmp_vec = dMrs_GRAD(XX,DELTA,D2,S2,
			new_theta,std::exp(PARS_2.at(0)),
			std::exp(PARS_2.at(1)),std::exp(PARS_2.at(2)),
			copula,upPARS);
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
	const arma::vec& D2,const arma::vec& S2,const arma::vec& PARS,
	const std::string& copula,const arma::vec& upPARS){
	
	arma::vec ePARS = arma::exp(PARS);
	double THETA = ePARS.at(3);
	// double KAPPA = 1.0 / (1.0 + std::exp(-PARS.at(2)));
	double KAPPA = ePARS.at(2);
	
	if( copula == "Gumbel" ){
		THETA += 1.0;
	}
	
	return dMrs_HESS(XX,DELTA,D2,S2,THETA,ePARS.at(0),
		ePARS.at(1),KAPPA,copula,upPARS);
}

// [[Rcpp::export]]
void dMrs_BFGS(const arma::vec& XX,const arma::uvec& DELTA,
	const arma::vec& D2,const arma::vec& S2,arma::vec& PARS,
	const std::string& copula,const arma::vec& upPARS,
	const arma::uword& max_iter = 4e3,const double& eps = 1e-7,
	const bool& verb = true){
	
	arma::uword iter = 0, jj, uu,
		reset_Bk = 0, np = PARS.n_elem;
	
	// Initialize parameters
	if( verb ) PARS.t().print("iPARS = "); // Rcpp::Rcout << "iPARS = " << PARS.t();
	arma::mat I_np = arma::eye<arma::mat>(np,np),
		inv_Bk = I_np, ISYT = I_np;
	arma::vec xk = PARS, curr_xk = arma::zeros<arma::vec>(np),
		new_xk = curr_xk, gr_k = curr_xk, gr2_k = curr_xk,
		p_k = curr_xk, s_k = curr_xk, y_k = curr_xk;
	double old_LL, new_LL, inv_norm_p_k, tmp_step, ys,
		fnscale = -1.0, curr_LL = 0.0, error_num = -999.0,
		diff_LL = 0.0;
	double orig_LL = dMrs_cLL(XX,DELTA,D2,S2,xk,copula,false);
	arma::vec orig_GRAD = dMrs_cGRAD(XX,DELTA,D2,S2,xk,copula,upPARS);
	
	if( orig_LL == error_num || arma::any(orig_GRAD == error_num) ){
		// Rcpp::Rcout << "Rerun optimization w/ new initialized parameters ...\n";
		return;
	}
	
	while(iter < max_iter){
		
		// Calculate Direction p_k
		gr_k = fnscale * dMrs_cGRAD(XX,DELTA,D2,S2,xk,copula,upPARS);
		p_k = -1.0 * inv_Bk * gr_k;
		inv_norm_p_k = 1.0 / std::max(1.0,Rcpp_norm(p_k));

		// Line search for new xk
		uu = 0;
		old_LL = fnscale * dMrs_cLL(XX,DELTA,D2,S2,xk,copula,false);
		diff_LL = std::abs(old_LL - curr_LL);
		
		if( verb ){
			Rcpp::Rcout << "Iter = " << iter + 1 << "\n";
			Rcpp::Rcout << "   LL = " << -old_LL << "\n";
			if( iter > 0 ) Rcpp::Rcout << "   abs_diff_LL = " << diff_LL << "\n";
			Rcpp::Rcout << "   PARS = " << xk.t();
			Rcpp::Rcout << "   grad = " << gr_k.t();
		}
		
		for(jj = 0; jj <= 30; jj++){
			tmp_step = inv_norm_p_k / std::pow(4,jj);
			new_xk = xk + tmp_step * p_k;
			
			new_LL = fnscale * dMrs_cLL(XX,DELTA,D2,S2,new_xk,copula,false);
			if( new_LL >= old_LL || new_LL == fnscale * error_num ) continue;
			
			gr2_k = fnscale * dMrs_cGRAD(XX,DELTA,D2,S2,new_xk,copula,upPARS);
			if( arma::any(gr2_k == fnscale * error_num) ) continue;
			
			y_k = gr2_k - gr_k;
			s_k = tmp_step * p_k;
			ys = arma::dot(y_k,s_k);
			if( ys > 0.0 ){
				ISYT = I_np - (s_k * y_k.t()) / ys;
				inv_Bk = ISYT * inv_Bk * ISYT.t() + s_k * s_k.t() / ys;
			}
			xk = new_xk;
			old_LL = new_LL;
			uu = 1;
			break;
			
		}
		
		if( uu == 0 ) { // aka no update
			if( Rcpp_norm(gr_k) > 1.0 ){
				if( verb ) printR_obj("Reset inv_Bk");
				inv_Bk = I_np;
				reset_Bk++;
			} else {
				if( verb ) printR_obj("Failed line search");
				break;
			}
		}
		
		if( reset_Bk > 5 ) break;
		
		// Check Convergence
		if( iter > 0 ){
			if( std::abs(curr_LL - old_LL) < eps &&
				// Rcpp_norm(curr_xk - xk) < eps
				Rcpp_max_abs_diff(curr_xk,xk) < eps ){
				gr_k = dMrs_cGRAD(XX,DELTA,D2,S2,xk,copula,upPARS);
				if( Rcpp_norm(gr_k) < eps ){
					break;
				}
			}
		}
		
		curr_xk = xk;
		curr_LL = old_LL;
		iter++;
	
	}
	
	old_LL = dMrs_cLL(XX,DELTA,D2,S2,xk,copula,false);
	if( old_LL > orig_LL ){ // Criteria for convergence
		PARS = xk;
	}

	if( verb ){
		Rcpp::Rcout << "####\nNum Iter = " << iter+1 << "\n";
		Rcpp::Rcout << "Params = " << xk.t();
		Rcpp::Rcout << "LL = " << old_LL << "\n";
		gr_k = dMrs_cGRAD(XX,DELTA,D2,S2,xk,copula,upPARS);
		arma::mat hh = dMrs_cHESS(XX,DELTA,D2,S2,PARS,copula,upPARS),
			ihh = arma::zeros<arma::mat>(4,4);
		arma::uvec nz = arma::find(hh.diag() != 0.0);
		ihh.submat(nz,nz) = arma::inv(-1.0 * hh.submat(nz,nz));
		Rcpp::Rcout << "GRAD = " << gr_k.t();
		Rcpp::Rcout << "Convergence Indicators: \n"
			<< "   NormGrad = " << Rcpp_norm(gr_k) << "\n"
			<< "   NormIHessGrad = " << Rcpp_norm(inv_Bk * gr_k) << "\n"
			<< "   NormIHessGrad2 = " << Rcpp_norm(ihh * gr_k) << "\n";
		Rcpp::Rcout << "Var = " << ihh.diag().t();
	}
	
}

// [[Rcpp::export]]
void dMrs_NR(const arma::vec& XX,const arma::uvec& DELTA,
	const arma::vec& D2,const arma::vec& S2,arma::vec& PARS,
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
		orig_LL = dMrs_cLL(XX,DELTA,D2,S2,xk,copula,false),
		nGRAD, old_LL, new_LL;
	arma::uvec chk = arma::zeros<arma::uvec>(np),
		idx_fin = arma::find_finite(PARS);
	
	old_LL = orig_LL;
	while( iter < max_iter ){
		GRAD = dMrs_cGRAD(XX,DELTA,D2,S2,xk,copula,upPARS);
		if( GRAD.at(0) == error_num ){
			if( verb ) Rcpp::Rcout << "Invalid pars\n";
			return;
		}
		
		HESS = dMrs_cHESS(XX,DELTA,D2,S2,xk,copula,upPARS);
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
		for(kk = 0; kk < 2; kk++){
		for(jj = 0; jj <= 40; jj++){
			
			if( kk == 0 ) new_xk = xk + pk_NR / std::pow(4.0,jj);
			if( kk == 1 ) new_xk = xk + pk_GD / std::pow(4.0,jj);
			
			new_LL = dMrs_cLL(XX,DELTA,D2,S2,new_xk,copula,false);
			if( new_LL == error_num ) continue;
			if( new_LL <= old_LL ) continue;
			
			GRAD = dMrs_cGRAD(XX,DELTA,D2,S2,new_xk,copula,upPARS);
			if( arma::any(GRAD == error_num) ) continue;
			
			HESS = dMrs_cHESS(XX,DELTA,D2,S2,new_xk,copula,upPARS);
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
	
	old_LL = dMrs_cLL(XX,DELTA,D2,S2,xk,copula,false);
	if( old_LL > orig_LL ){ // Criteria for convergence
		PARS = xk;
	}

	if( verb ){
		Rcpp::Rcout << "####\nNum Iter = " << iter+1 << "\n";
		Rcpp::Rcout << "Params = "; prt_vec(xk);
		Rcpp::Rcout << "LL = " << old_LL << "\n";
		GRAD = dMrs_cGRAD(XX,DELTA,D2,S2,xk,copula,upPARS);
		HESS = dMrs_cHESS(XX,DELTA,D2,S2,PARS,copula,upPARS);
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
	const arma::vec& D2,const arma::vec& S2,const arma::vec& log_THETA,
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
	if( verb )
		Rcpp::Rcout << "Num grid points = " << tot << "\n";
	
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
		shared(verb2,copula,tot,num_alpha,num_lambda,\
			num_kappa,num_theta,XX,DELTA,D2,S2,DAT)
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
		DAT.at(cnt2,4) = dMrs_cLL(XX,DELTA,D2,S2,PARS,copula,false);
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
arma::mat dMrs_MATCH(const arma::mat& wDAT,const arma::mat& rDAT,
	const int& ncores = 1,const bool& verb = true){
	// wDAT columns: age, yr_diag, yr_event, sex
	// rDAT columns: Year, age, qx(hazard), sex
	
	arma::uword nn = wDAT.n_rows;
	arma::mat OUT = arma::zeros<arma::mat>(nn,2); // output dens_t2, surv_t2
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
		
		double yr_diag, yr_event, age, sex;
		yr_diag = wDAT.at(ii,1);
		yr_event = wDAT.at(ii,2);
		age = wDAT.at(ii,0);
		sex = wDAT.at(ii,3);
		
		arma::uvec idx = arma::find(rDAT.col(0) >= yr_diag
			&& rDAT.col(0) <= yr_event && rDAT.col(1) == age
			&& rDAT.col(3) == sex);
		
		arma::mat tmp_mat = rDAT.rows(idx);
		arma::vec tmp_haz = tmp_mat.col(2);
		
		// Calculate survival
		OUT.at(ii,1) = std::exp(-1.0 * arma::sum(tmp_haz));
		
		// Calculate density
		OUT.at(ii,0) = tmp_haz.at(tmp_haz.n_elem - 1) / OUT.at(ii,1);
		
	}
	
	return OUT;
}





