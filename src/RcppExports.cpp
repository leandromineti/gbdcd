// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// rcpp_partition
Rcpp::NumericVector rcpp_partition(const arma::mat& neighbors, const Rcpp::NumericVector& centers);
RcppExport SEXP _gbdcd_rcpp_partition(SEXP neighborsSEXP, SEXP centersSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type neighbors(neighborsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type centers(centersSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_partition(neighbors, centers));
    return rcpp_result_gen;
END_RCPP
}
// RcppFreqMatrix
Rcpp::NumericMatrix RcppFreqMatrix(const Rcpp::NumericVector& partitions);
RcppExport SEXP _gbdcd_RcppFreqMatrix(SEXP partitionsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type partitions(partitionsSEXP);
    rcpp_result_gen = Rcpp::wrap(RcppFreqMatrix(partitions));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_gbdcd_rcpp_partition", (DL_FUNC) &_gbdcd_rcpp_partition, 2},
    {"_gbdcd_RcppFreqMatrix", (DL_FUNC) &_gbdcd_RcppFreqMatrix, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_gbdcd(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
