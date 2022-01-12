// [[Rcpp::depends(Matrix)]]
#include <RcppArmadillo.h>
#include <R.h>
#include "Rcpp_auc.hpp"
#include <Rinternals.h>
#include "coding.hpp"
#include "lasso_solver.hpp"


using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat SSI(arma::mat x, arma::vec y){
  return lmpairwise_fast(y, x);
}


// [[Rcpp::export]]
RcppExport SEXP BOLT_SSI(arma::mat X, arma::vec y, int extra_pairs = -1, int code_num = 3,  int thread_num = 4){
  time_t rawtime;
  time(&rawtime);
  cout << "Start at " << ctime(&rawtime) << endl;
  codes::thread_num = thread_num;
  int N = X.n_rows;
  int P = X.n_cols;
  long long total = (long long)P * ((long long)P - 1) / 2;

  Col<int> yi(N,fill::zeros);
  uvec uy = find_unique(y);
  if(uy.size() > 2){
    uvec pos1 = find(y > median(y));
    yi.elem(pos1).fill(1);
  }else{
    yi = conv_to<Col<int>>::from(y);
  }
  double pvalue_thresh = codes::max_num_of_pair > total ? 0.99 : (codes::max_num_of_pair * 1.0 / total);
  double chisq_thresh = R::qchisq(1 - pvalue_thresh, 4, true, false);
  GenoDiscreteMatrix gdm(y.memptr(), yi.memptr(), X.memptr(), N, P, code_num);
  time_t t1,t2;
  time(&t1);
  vector<InteractionSNPpairs> top_pairs = getTopInteractionParis(&gdm, extra_pairs, chisq_thresh);
  time(&t2);
  cout << "Time for getting the top interaction paris : " <<difftime(t2, t1) << endl;

  int nTop = top_pairs.size();
  arma::mat result(nTop, 3);
  InteractionSNPpairs pr;
  for(int i = 0; i < nTop; i++){
    pr = top_pairs[i];
    result.at(i,0) = pr.j1 + 1;
    result.at(i,1) = pr.j2 + 1;
    result.at(i,2) = pr.InteractionMeasure;
  }
  return wrap(result);
}

// [[Rcpp::export]]
RcppExport SEXP BOLT_SSI_RR(arma::mat X, arma::vec y, int extra_pairs = -1, int code_num = 3 , int nLambda = 100, int thread_num = 1){
  time_t rawtime;
  time(&rawtime);
  codes::thread_num = thread_num;
  cout << "Start at " << ctime(&rawtime) << endl;
  Rcpp::List ret;
  int N = X.n_rows;
  int P = X.n_cols;
  Col<int> yi(N,fill::zeros);
  uvec uy = find_unique(y);
  if(uy.size() > 2){
    uvec pos1 = find(y > median(y));
    yi.elem(pos1).fill(1);
  }else{
    yi = conv_to<Col<int>>::from(y);
  }
  long long total = (long long)P * ((long long)P - 1) / 2;
  double pvalue_thresh = codes::max_num_of_pair > total ? 0.99 : (codes::max_num_of_pair * 1.0 / total);
  double chisq_thresh = R::qchisq(1 - pvalue_thresh, 4, true, false);
  GenoDiscreteMatrix gdm(y.memptr(), yi.memptr(), X.memptr(), N, P, code_num);
  vector<InteractionSNPpairs> top_pairs = getTopInteractionParis(&gdm, extra_pairs, chisq_thresh);

  if(top_pairs.size() > 0){
    class InteractionSNPpairs pr;
    cout <<"Top(at most 20) Interaction Terms for full Data!" << endl;
    for(int i = 0; i < top_pairs.size() && i < 20; i++){
      pr = top_pairs[i];
      cout<< pr.InteractionMeasure<<" " << pr.j1 + 1 << ":"  << pr.j2 + 1<< endl;
    }
  }
  mat rX = X;
  vec ry = y;
  LassoFit fit(X, y, gdm, nLambda, extra_pairs, chisq_thresh);
  double max_lam = cal_max_lambda(X, y);
  arma::vec lambdas = logspace<vec>(log10(max_lam), log10(0.0001* max_lam), nLambda);
  fit.lambdas  = lambdas;
  LassoPath(X, y, fit);
  ret = wrap_fit(&fit);
  return ret;
}

// [[Rcpp::export]]
RcppExport SEXP CV_BOLT_SSI_RR(arma::mat X, arma::vec y, int extra_pairs = -1, int code_num = 3 ,int nfold = 10, int nLambda = 100, int thread_num = 4){
  time_t rawtime;
  time(&rawtime);
  codes::thread_num = thread_num;
  cout << "Start at " << ctime(&rawtime) << endl;

  Rcpp::List ret;
  uword N = X.n_rows;
  uword P = X.n_cols;

  Col<int> yi(N,fill::zeros);
  uvec uy = find_unique(y);
  if(uy.size() > 2){
    uvec pos1 = find(y > median(y));
    yi.elem(pos1).fill(1);
  }else{
    yi = conv_to<Col<int>>::from(y);
  }
  long long total = (long long)P * ((long long)P - 1) / 2;
  double pvalue_thresh = codes::max_num_of_pair > total ? 0.99 : (codes::max_num_of_pair * 1.0 / total);
  double chisq_thresh = R::qchisq(1 - pvalue_thresh, 4, true, false);
  GenoDiscreteMatrix gdm(y.memptr(), yi.memptr(), X.memptr(), N, P, code_num);
  mat rX = X;
  vec ry = y;
  LassoFit fit(X, y, gdm, nLambda, extra_pairs, chisq_thresh);
  double max_lam = cal_max_lambda(X, y);
  arma::vec lambdas = logspace<vec>(log10(max_lam), log10(0.0001* max_lam), nLambda);
  fit.lambdas  = lambdas;
  LassoPath(X, y, fit);
  Col<uword> indices = cross_valind(N, nfold);
  mat y_hat(N, nLambda);
  cout << "In "<<nfold << " fold CV Procedure";
  for (uword i = 1; i <= nfold; i++) {
    cout << ".";
    Col<uword> train_idx = find(indices != i);
    Col<uword> test_idx = find(indices == i);
    mat trainM = rX.rows(train_idx);
    vec ytrain = ry(train_idx);
    GenoDiscreteMatrix subMatrix((int)train_idx.size(), (int)P);
    gdm.subMatrix(train_idx, subMatrix);
    subMatrix.Xc = trainM.memptr();
    subMatrix.yc = ytrain.memptr();

    LassoFit fit_i(trainM, ytrain, subMatrix,  nLambda, extra_pairs, chisq_thresh);
    fit_i.lambdas = lambdas;
    mat testM = rX.rows(test_idx);
    LassoPath(trainM, ytrain, fit_i);
    mat result_i = predictAlongPath(testM, fit_i);
    y_hat.rows(test_idx) = result_i;
  }
  time(&rawtime);
  cout << endl <<"Finish the circle at " << ctime(&rawtime) << endl;
  LassoCVFit fit0;
  int min_idx = 1;
  double min_err = INFINITY;
  vec cvm(nLambda);
  vec cvsd(nLambda);
  vec cvup(nLambda);
  vec cvlo(nLambda);
  vec auc(nLambda);
  vec R2(nLambda);
  vec div = (ry - mean(ry));
  double r2 = sum(div % div);
  for(int i = 0; i < nLambda; i++){
    vec diff = y_hat.col(i) - ry;
    R2[i] = 1 - sum(diff % diff) / r2;
    double err = sum(diff % diff) / N;
    cvm[i] = err;
    cvsd[i] = stddev(diff % diff) / sqrt(N);
    cvup[i] = cvm[i] + cvsd[i];
    cvlo[i] = cvm[i] - cvsd[i];
    auc[i] = calauc(ry, y_hat.col(i) );
    if(err < min_err){
      min_err = err;
      min_idx = i;
    }
  }

  fit0.cvm = cvm;
  fit0.cvsd = cvsd;
  fit0.cvlo = cvlo;
  fit0.cvup = cvup;
  fit0.auc = auc;
  fit0.index_min = min_idx + 1;
  fit0.lambda_min = lambdas[min_idx];
  fit0.fit = fit;
  fit0.R2 = R2;
  ret = wrap_cvfit(&fit0);
  return ret;
}
