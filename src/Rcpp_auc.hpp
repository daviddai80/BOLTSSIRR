#ifndef Rcpp_aux_hpp
#define Rcpp_aux_hpp
#include <Rcpp.h>
#include "coding.hpp"
#include "lasso_solver.hpp"
#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;

RcppExport SEXP wrap_fit(LassoFit* fit){
  Rcpp::List ret;
  ret["lambdas"] = fit -> lambdas;
  ret["beta"] = fit -> beta;
  ret["covs"] = fit -> covs;
  ret["iter"] = fit -> iter;
  ret["lambda_min"] = fit -> lambda_min;
  ret["cve"] = fit -> cve;
  ret["cvse"] = fit -> cvse;
  ret["xcenter"] = fit ->x_center;
  ret["xs"] = fit ->x_sd;
  Col<int> pairs(fit ->pairs.size() * 2);
  class InteractionSNPpairs pair;
  for(int i = 0; i < fit ->pairs.size(); i++){
    pair = fit -> pairs[i];
    pairs[2*i] = pair.j1 + 1;
    pairs[2*i + 1] = pair.j2 + 1;
  }
  ret["pairs"] = pairs;
  return ret;
}

RcppExport SEXP wrap_cvfit(LassoCVFit* fit){
  Rcpp::List ret;
  ret = wrap_fit(&(fit -> fit));
  //  ret["beta"] = fit -> fit.beta;
  //  ret["lambdas"] = fit -> fit.lambdas;
  ret["cvm"] = fit -> cvm;
  ret["cvsd"] = fit -> cvsd;
  ret["cvup"] = fit -> cvup;
  ret["cvlo"] = fit -> cvlo;
  ret["auc"] = fit -> auc;
  ret["R2"] = fit -> R2;
  ret["nzero"] = fit -> nzero;
  ret["lambda_min"] = fit -> lambda_min;
  ret["index_min"] = fit -> index_min;
  ret["lambda_1se"] = fit -> lambda_1se;
  return ret;
}



#endif /* BOLT_SSI_RR_hpp */
