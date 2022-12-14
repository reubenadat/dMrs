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

// --------------------
// New Functions

double dMrs_LL(const arma::vec& XX,const arma::uvec& DELTA,
	const arma::vec& D2,const arma::vec& S2,const double& THETA,
	const double& ALPHA,const double& LAMBDA,const double& KAPPA,
	const std::string& copula,const bool& show = false){
	
	arma::uword ii, NN = XX.n_elem;
	double LL = 0.0, LL_ii, D1, S1, F1, F2,
		log_F1, log_F2, D1_D2, XDL, enXDLa, U1T_U2T, U1T1, U2T1,
		F1_F2, AA, BB, error_num = -999.0,
		KAL = KAPPA * ALPHA / LAMBDA;
	
	for(ii = 0; ii < NN; ii++){
		XDL = XX.at(ii) / LAMBDA;
		enXDLa = std::exp(-std::pow(XDL,ALPHA));
		
		F1 = 1.0 - enXDLa;
		if( KAPPA != 1.0 ) F1 = std::pow(F1,KAPPA);
		
		S1 = 1.0 - F1;
		
		D1 = KAL * std::pow(XDL,ALPHA - 1.0) * enXDLa;
		if( KAPPA != 1.0 ) D1 *= F1 / ( 1.0 - enXDLa );
		
		F2 = 1.0 - S2.at(ii);
		if( F1 == 0.0 || F2 == 0.0 ) return error_num;
		
		if( copula == "Clayton" ){
			if( THETA == 0.0 ){
				F1_F2 = F1 * F2;
				D1_D2 = D1 * F2 + F1 * D2.at(ii);
			} else {
				U1T1 = std::pow(F1,THETA) * F1;
				U2T1 = std::pow(F2,THETA) * F2;
				if( U1T1 == 0.0 || U2T1 == 0.0 ) return error_num;
				U1T_U2T = std::pow(F1,-THETA) + std::pow(F2,-THETA) - 1.0;
				D1_D2 = std::pow(U1T_U2T,-1.0 / THETA - 1.0) *
					( D1 / U1T1 + D2.at(ii) / U2T1 );
				F1_F2 = std::pow(U1T_U2T,-1.0 / THETA);
			}
		} else if( copula == "Gumbel" ){
			if( THETA == 1.0 ){
				F1_F2 = F1 * F2;
				D1_D2 = D1 * F2 + F1 * D2.at(ii);
			} else {
				log_F1 = std::log(F1);
				log_F2 = std::log(F2);
				U1T_U2T = std::pow(-log_F1,THETA) + std::pow(-log_F2,THETA);
				F1_F2 = std::exp(-std::pow(U1T_U2T,1.0/THETA));
				D1_D2 = F1_F2 * std::pow(U1T_U2T,1.0/THETA - 1.0) *
					( std::pow(-log_F1,THETA - 1.0) * D1 / F1 + 
					std::pow(-log_F2,THETA - 1.0) * D2.at(ii) / F2 );
			}
		} else {
			Rcpp::stop("Not a valid copula!");
		}
		
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
		if( show ){
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
	const std::string& copula,const bool& show = false){
	
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
		ePARS.at(1),KAPPA,copula,show);
}

arma::vec dMrs_GRAD(const arma::vec& XX,const arma::uvec& DELTA,
	const arma::vec& D2,const arma::vec& S2,const double& THETA,
	const double& ALPHA,const double& LAMBDA,const double& KAPPA,
	const std::string& copula,const arma::vec& upPARS){
	
	arma::uword ii;
	double shift = 1e-5, old_LL, new_LL, error_num = -999.0;
	arma::vec GRAD = arma::zeros<arma::vec>(4),
		old_PARS = GRAD, new_PARS = GRAD;
	old_PARS.at(0) = std::log(ALPHA);
	old_PARS.at(1) = std::log(LAMBDA);
	old_PARS.at(2) = std::log(KAPPA);
	old_PARS.at(3) = std::log(THETA);
	
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
		new_LL = dMrs_LL(XX,DELTA,D2,S2,std::exp(new_PARS.at(3)),
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
	double shift = 1e-5, error_num = -999.0;
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
	PARS.at(3) = std::log(THETA);
	
	for(ii = 0; ii < np; ii++){
		if( upPARS.at(ii) == 0.0 ) continue;
		PARS_2 = PARS + shift * I_np.col(ii) % upPARS;
		tmp_vec = dMrs_GRAD(XX,DELTA,D2,S2,
			std::exp(PARS_2.at(3)),std::exp(PARS_2.at(0)),
			std::exp(PARS_2.at(1)),std::exp(PARS_2.at(2)),
			copula,upPARS) - old_GRAD;
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
	const bool& show = true){
	
	arma::uword iter = 0, jj, uu,
		reset_Bk = 0, np = PARS.n_elem;
	
	// Initialize parameters
	if( show ) PARS.t().print("iPARS = "); // Rcpp::Rcout << "iPARS = " << PARS.t();
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
		
		if( show ){
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
				if( show ) printR_obj("Reset inv_Bk");
				inv_Bk = I_np;
				reset_Bk++;
			} else {
				if( show ) printR_obj("Failed line search");
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

	if( show ){
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
Rcpp::List dMrs_GRID(const arma::vec& XX,const arma::uvec& DELTA,
	const arma::vec& D2,const arma::vec& S2,const arma::vec& log_THETA,
	const arma::vec& log_ALPHA,const arma::vec& log_LAMBDA,
	const arma::vec& unc_KAPPA,const std::string& copula,
	const bool& show = false,const int& ncores = 1){
	
	arma::uword num_alpha = log_ALPHA.n_elem, num_lambda = log_LAMBDA.n_elem,
		num_kappa = unc_KAPPA.n_elem, num_theta = log_THETA.n_elem,
		tot = num_alpha * num_lambda * num_kappa * num_theta,
		aa, ll, kk, tt, cnt = 0;
	// log_ALPHA,log_LAMBDA,unc_KAPPA,log_THETA,LL
	arma::mat DAT = arma::zeros<arma::mat>(tot,5);
	double error_num = -999.0;
	bool show2 = show && ncores == 1;
	if( show )
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
		shared(show2,copula,tot,num_alpha,num_lambda,\
			num_kappa,num_theta,XX,DELTA,D2,S2,DAT)
	#endif
	for(arma::uword cnt2 = 0; cnt2 < tot; cnt2++){
		if( show2 ){
			if( (cnt2 + 1) % 500 == 0 )
				Rcpp::Rcout << ".";
			if( (cnt2 + 1) % 5000 == 0 || (cnt2 + 1) == tot )
				Rcpp::Rcout << (cnt2 + 1) << " out of " << tot << " done\n";
		}
		
		arma::vec PARS = DAT(cnt2,arma::span(0,3)).t();
		DAT.at(cnt2,4) = dMrs_cLL(XX,DELTA,D2,S2,PARS,copula,false);
	}
	
	// Get parameters with largest log-likelihood
	double min_val = arma::min(DAT.col(4));
	arma::uvec uDAT = arma::zeros<arma::uvec>(DAT.n_rows);
	uDAT = DAT.col(4) == error_num;
	if( min_val != error_num ){
		DAT.col(4) = DAT.col(4) + min_val * arma::conv_to<arma::vec>::from(uDAT);
	}
	arma::uword max_idx = arma::index_max(DAT.col(4));
	arma::vec fPARS = arma::zeros<arma::vec>(4);
	fPARS.at(0) = DAT.at(max_idx,0);
	fPARS.at(1) = DAT.at(max_idx,1);
	fPARS.at(2) = DAT.at(max_idx,2);
	
	return Rcpp::List::create(
		Rcpp::Named("PARS",Rcpp::NumericVector(fPARS.begin(),fPARS.end())),
    Rcpp::Named("DAT",DAT));
	
}

// [[Rcpp::export]]
arma::mat dMrs_MATCH(const arma::mat& wDAT,const arma::mat& rDAT,
	const int& ncores = 1,const bool& show = true){
	// wDAT columns: age, yr_diag, yr_event, sex
	// rDAT columns: Year, age, qx(hazard), sex
	
	arma::uword nn = wDAT.n_rows;
	arma::mat OUT = arma::zeros<arma::mat>(nn,2); // output dens_t2, surv_t2
	bool show2 = show && ncores == 1;
	
	#ifdef _OPENMP
	# pragma omp parallel for schedule(dynamic) \
		num_threads(ncores) \
		shared(show2,nn,wDAT,rDAT,OUT)
	#endif
	for(arma::uword ii = 0; ii < nn; ii++){
		
		if( show2 ){
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





