// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// Rcpp_norm
double Rcpp_norm(const arma::vec& a);
RcppExport SEXP _dMrs_Rcpp_norm(SEXP aSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type a(aSEXP);
    rcpp_result_gen = Rcpp::wrap(Rcpp_norm(a));
    return rcpp_result_gen;
END_RCPP
}
// Rcpp_max_abs_diff
double Rcpp_max_abs_diff(const arma::vec& aa, const arma::vec& bb);
RcppExport SEXP _dMrs_Rcpp_max_abs_diff(SEXP aaSEXP, SEXP bbSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type aa(aaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type bb(bbSEXP);
    rcpp_result_gen = Rcpp::wrap(Rcpp_max_abs_diff(aa, bb));
    return rcpp_result_gen;
END_RCPP
}
// Rcpp_logSumExp
double Rcpp_logSumExp(const arma::vec& log_x);
RcppExport SEXP _dMrs_Rcpp_logSumExp(SEXP log_xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type log_x(log_xSEXP);
    rcpp_result_gen = Rcpp::wrap(Rcpp_logSumExp(log_x));
    return rcpp_result_gen;
END_RCPP
}
// dMrs_cLL
double dMrs_cLL(const arma::vec& XX, const arma::uvec& DELTA, const arma::vec& D2, const arma::vec& S2, const arma::vec& PARS, const std::string& copula, const bool& verb);
RcppExport SEXP _dMrs_dMrs_cLL(SEXP XXSEXP, SEXP DELTASEXP, SEXP D2SEXP, SEXP S2SEXP, SEXP PARSSEXP, SEXP copulaSEXP, SEXP verbSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type XX(XXSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type DELTA(DELTASEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type D2(D2SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type S2(S2SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type PARS(PARSSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type copula(copulaSEXP);
    Rcpp::traits::input_parameter< const bool& >::type verb(verbSEXP);
    rcpp_result_gen = Rcpp::wrap(dMrs_cLL(XX, DELTA, D2, S2, PARS, copula, verb));
    return rcpp_result_gen;
END_RCPP
}
// dMrs_cGRAD
arma::vec dMrs_cGRAD(const arma::vec& XX, const arma::uvec& DELTA, const arma::vec& D2, const arma::vec& S2, const arma::vec& PARS, const std::string& copula, const arma::vec& upPARS);
RcppExport SEXP _dMrs_dMrs_cGRAD(SEXP XXSEXP, SEXP DELTASEXP, SEXP D2SEXP, SEXP S2SEXP, SEXP PARSSEXP, SEXP copulaSEXP, SEXP upPARSSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type XX(XXSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type DELTA(DELTASEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type D2(D2SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type S2(S2SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type PARS(PARSSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type copula(copulaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type upPARS(upPARSSEXP);
    rcpp_result_gen = Rcpp::wrap(dMrs_cGRAD(XX, DELTA, D2, S2, PARS, copula, upPARS));
    return rcpp_result_gen;
END_RCPP
}
// dMrs_cHESS
arma::mat dMrs_cHESS(const arma::vec& XX, const arma::uvec& DELTA, const arma::vec& D2, const arma::vec& S2, const arma::vec& PARS, const std::string& copula, const arma::vec& upPARS);
RcppExport SEXP _dMrs_dMrs_cHESS(SEXP XXSEXP, SEXP DELTASEXP, SEXP D2SEXP, SEXP S2SEXP, SEXP PARSSEXP, SEXP copulaSEXP, SEXP upPARSSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type XX(XXSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type DELTA(DELTASEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type D2(D2SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type S2(S2SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type PARS(PARSSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type copula(copulaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type upPARS(upPARSSEXP);
    rcpp_result_gen = Rcpp::wrap(dMrs_cHESS(XX, DELTA, D2, S2, PARS, copula, upPARS));
    return rcpp_result_gen;
END_RCPP
}
// dMrs_BFGS
void dMrs_BFGS(const arma::vec& XX, const arma::uvec& DELTA, const arma::vec& D2, const arma::vec& S2, arma::vec& PARS, const std::string& copula, const arma::vec& upPARS, const arma::uword& max_iter, const double& eps, const bool& verb);
RcppExport SEXP _dMrs_dMrs_BFGS(SEXP XXSEXP, SEXP DELTASEXP, SEXP D2SEXP, SEXP S2SEXP, SEXP PARSSEXP, SEXP copulaSEXP, SEXP upPARSSEXP, SEXP max_iterSEXP, SEXP epsSEXP, SEXP verbSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type XX(XXSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type DELTA(DELTASEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type D2(D2SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type S2(S2SEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type PARS(PARSSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type copula(copulaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type upPARS(upPARSSEXP);
    Rcpp::traits::input_parameter< const arma::uword& >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< const double& >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< const bool& >::type verb(verbSEXP);
    dMrs_BFGS(XX, DELTA, D2, S2, PARS, copula, upPARS, max_iter, eps, verb);
    return R_NilValue;
END_RCPP
}
// dMrs_GRID
Rcpp::List dMrs_GRID(const arma::vec& XX, const arma::uvec& DELTA, const arma::vec& D2, const arma::vec& S2, const arma::vec& log_THETA, const arma::vec& log_ALPHA, const arma::vec& log_LAMBDA, const arma::vec& unc_KAPPA, const std::string& copula, const bool& verb, const int& ncores);
RcppExport SEXP _dMrs_dMrs_GRID(SEXP XXSEXP, SEXP DELTASEXP, SEXP D2SEXP, SEXP S2SEXP, SEXP log_THETASEXP, SEXP log_ALPHASEXP, SEXP log_LAMBDASEXP, SEXP unc_KAPPASEXP, SEXP copulaSEXP, SEXP verbSEXP, SEXP ncoresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type XX(XXSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type DELTA(DELTASEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type D2(D2SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type S2(S2SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type log_THETA(log_THETASEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type log_ALPHA(log_ALPHASEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type log_LAMBDA(log_LAMBDASEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type unc_KAPPA(unc_KAPPASEXP);
    Rcpp::traits::input_parameter< const std::string& >::type copula(copulaSEXP);
    Rcpp::traits::input_parameter< const bool& >::type verb(verbSEXP);
    Rcpp::traits::input_parameter< const int& >::type ncores(ncoresSEXP);
    rcpp_result_gen = Rcpp::wrap(dMrs_GRID(XX, DELTA, D2, S2, log_THETA, log_ALPHA, log_LAMBDA, unc_KAPPA, copula, verb, ncores));
    return rcpp_result_gen;
END_RCPP
}
// dMrs_MATCH
arma::mat dMrs_MATCH(const arma::mat& wDAT, const arma::mat& rDAT, const int& ncores, const bool& verb);
RcppExport SEXP _dMrs_dMrs_MATCH(SEXP wDATSEXP, SEXP rDATSEXP, SEXP ncoresSEXP, SEXP verbSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type wDAT(wDATSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type rDAT(rDATSEXP);
    Rcpp::traits::input_parameter< const int& >::type ncores(ncoresSEXP);
    Rcpp::traits::input_parameter< const bool& >::type verb(verbSEXP);
    rcpp_result_gen = Rcpp::wrap(dMrs_MATCH(wDAT, rDAT, ncores, verb));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_dMrs_Rcpp_norm", (DL_FUNC) &_dMrs_Rcpp_norm, 1},
    {"_dMrs_Rcpp_max_abs_diff", (DL_FUNC) &_dMrs_Rcpp_max_abs_diff, 2},
    {"_dMrs_Rcpp_logSumExp", (DL_FUNC) &_dMrs_Rcpp_logSumExp, 1},
    {"_dMrs_dMrs_cLL", (DL_FUNC) &_dMrs_dMrs_cLL, 7},
    {"_dMrs_dMrs_cGRAD", (DL_FUNC) &_dMrs_dMrs_cGRAD, 7},
    {"_dMrs_dMrs_cHESS", (DL_FUNC) &_dMrs_dMrs_cHESS, 7},
    {"_dMrs_dMrs_BFGS", (DL_FUNC) &_dMrs_dMrs_BFGS, 10},
    {"_dMrs_dMrs_GRID", (DL_FUNC) &_dMrs_dMrs_GRID, 11},
    {"_dMrs_dMrs_MATCH", (DL_FUNC) &_dMrs_dMrs_MATCH, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_dMrs(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}