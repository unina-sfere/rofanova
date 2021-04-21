## usethis namespace: start
#' @importFrom Rcpp sourceCpp
#' @useDynLib rofanova, .registration = TRUE
## usethis namespace: end
NULL



#' @title Robust Functional Anlaysis of Variance
#' @details
#'
#'\tabular{ll}{
#'Package: \tab rofanova\cr
#'Type: \tab Package\cr
#'Version: \tab `r packageVersion("rofanova")` \cr
#'Date: \tab  `r Sys.Date()` \cr
#'License: \tab `r packageDescription("rofanova", fields="License")`\cr
#'}
#'
#'
#'
#'
#' @aliases {rofanova}-package
#' @author Fabio Centofanti, Antonio Lepore, Biagio Palumbo
#' @references
#' Centofanti, F., Lepore, A., & Palumbo, B. (2021).
#' Sparse and Smooth Functional Data Clustering.
#' \emph{arXiv preprint arXiv:2103.15224}.
#'
#' @seealso \code{\link{rofanova}},  \code{\link{sasfclust_cv}}
#' @examples
#' \donttest{
#' library(rofanova)
#' data_out<-simulate_data(scenario="one-way")
#' label=data_out$label
#' X_fdata<-data_out$X_fdata
#' per_list_median<-rofanova(X_fdata,label,B = 10,eff=eff,family="median")
#' pvalue_median<-per_list_median$pval
#' per_list_huber<-rofanova(X_fdata,label,B = 10,eff=eff,family="huber")
#' pvalue_huber<-per_list_huber$pval
#' per_list_bisquare<-rofanova(X_fdata,label,B = 10,eff=eff,family="bisquare")
#' pvalue_bisquare<-per_list_bisquare$pval
#' per_list_hampel<-rofanova(X_fdata,label,B = 10,eff=eff,family="hampel")
#' pvalue_hampel<-per_list_hampel$pval
#' per_list_optimal<-rofanova(X_fdata,label,B = 10,eff=eff,family="optimal")
#' pvalue_optimal<-per_list_optimal$pval
#'}
#'@import fda.usc robustbase
"_PACKAGE"



