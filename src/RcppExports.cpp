// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// quant
List quant(std::string ref_path, std::vector<std::string> fastq_path);
RcppExport SEXP _CRISPRExpress_quant(SEXP ref_pathSEXP, SEXP fastq_pathSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type ref_path(ref_pathSEXP);
    Rcpp::traits::input_parameter< std::vector<std::string> >::type fastq_path(fastq_pathSEXP);
    rcpp_result_gen = Rcpp::wrap(quant(ref_path, fastq_path));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_CRISPRExpress_quant", (DL_FUNC) &_CRISPRExpress_quant, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_CRISPRExpress(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
