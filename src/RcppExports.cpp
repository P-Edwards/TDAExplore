// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// sparse_patch_top_intensities
IntegerMatrix sparse_patch_top_intensities(IntegerVector pixels_in_order, int number_of_pixels_to_sample, int data_rows, int data_cols, LogicalMatrix exclusion_disc_data);
RcppExport SEXP _TDAExplore_sparse_patch_top_intensities(SEXP pixels_in_orderSEXP, SEXP number_of_pixels_to_sampleSEXP, SEXP data_rowsSEXP, SEXP data_colsSEXP, SEXP exclusion_disc_dataSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type pixels_in_order(pixels_in_orderSEXP);
    Rcpp::traits::input_parameter< int >::type number_of_pixels_to_sample(number_of_pixels_to_sampleSEXP);
    Rcpp::traits::input_parameter< int >::type data_rows(data_rowsSEXP);
    Rcpp::traits::input_parameter< int >::type data_cols(data_colsSEXP);
    Rcpp::traits::input_parameter< LogicalMatrix >::type exclusion_disc_data(exclusion_disc_dataSEXP);
    rcpp_result_gen = Rcpp::wrap(sparse_patch_top_intensities(pixels_in_order, number_of_pixels_to_sample, data_rows, data_cols, exclusion_disc_data));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP _rcpp_module_boot_Landscape();

static const R_CallMethodDef CallEntries[] = {
    {"_TDAExplore_sparse_patch_top_intensities", (DL_FUNC) &_TDAExplore_sparse_patch_top_intensities, 5},
    {"_rcpp_module_boot_Landscape", (DL_FUNC) &_rcpp_module_boot_Landscape, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_TDAExplore(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
