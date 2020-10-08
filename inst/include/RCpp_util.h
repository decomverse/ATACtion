#ifndef RCPP_UTIL_H
#define RCPP_UTIL_H

namespace ATACtion {
    Rcpp::IntegerMatrix tabulate2dCpp(Rcpp::IntegerVector x1, int xmin, int xmax, Rcpp::IntegerVector y1, int ymin, int ymax);
}

#endif
