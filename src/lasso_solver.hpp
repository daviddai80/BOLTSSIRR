//
//  lasso_solver.hpp
//  boost_new
//
//  Created by David on 2017/8/10.
//  Copyright © 2017年 David. All rights reserved.
//

#ifndef lasso_solver_hpp
#define lasso_solver_hpp

#include <stdio.h>
#include <armadillo>
#include "coding.hpp"

#define type_int 1
#define type_double 0


using namespace arma;



class LassoFit{
public:
    LassoFit();
    ~LassoFit();
    LassoFit(int P, int nLambda);
    LassoFit(mat& X, vec& y, int nLambda);
    LassoFit(mat& X, vec& y,GenoDiscreteMatrix& gdm, int nLambda, int extra_pairs, double chisq_thresh);
    LassoFit(LassoFit& fit0);
    int P;
    int nLambda;
    vec lambdas;
    sp_mat beta;
    vector<InteractionSNPpairs> pairs;
    vec covs;
    int iter;
   // int min;
    double mean_y;
    rowvec x_center;
    rowvec x_sd;
    rowvec x_var;
    double lambda_min;
    vec cve;
    vec cvse;
    vec sumq2;
};

class LassoCVFit{
public:
  LassoCVFit();
  ~LassoCVFit();
  LassoFit fit;
  vec lambda;
  vec cvm;
  vec cvsd;
  vec cvup;
  vec cvlo;
  vec nzero;
  vec auc;
  vec R2;
  double lambda_min;
  int index_min;
  double lambda_1se;

};


void lassosolver(double* X, double* y_, int N, int P, double lambda, double* xvar, double* beta0, int* mm, int* active_set, int& nin, double* beta, double* yhat_ptr, int& total_iter);
//void LassoPath(arma::mat& X, arma::vec y, arma::vec& lambdas, sp_mat& beta, vec& covs)
//void LassoPath(double* X, double* y, int N, int P, double* lamseq, int num,sp_mat& beta, vec& covs);
double shrinkage(double a, double kappa);
void LassoPath(arma::mat& X, arma::vec y, LassoFit& fit /*arma::vec& lambdas, sp_mat& beta, vec& covs */);
Mat<double> cal_means(Mat<double>& X_mat, int N);
void centering(Mat<double>& X_mat, vec& y, double& mean_y, mat& SZX, int N, int P, mat* SZX_ = nullptr, double* y_mean_ = nullptr);
vec predict(arma::mat& X,  double cov,  vec& beta);
arma::Col<uword> cross_valind(arma::uword N, arma::uword nfold);
double cal_max_lambda(arma::mat& X, arma::vec& y);
// double cal_max_lambda(arma::mat& X, arma::vec& y, vector<InteractionSNPpairs>& paris);
mat predictAlongPath(arma::mat& X,  LassoFit& fit);//sp_mat& beta, vec& covs);
void LassoPathCV(arma::mat& X, arma::vec y, LassoCVFit& fit0, int nfold = 10, int nLambda = 100);
void scale(arma::mat& X, arma::rowvec& S, arma::rowvec& M);
double calauc(arma::vec label, arma::vec pred);
// void CV_BOLT_SSI_RR(arma::mat X, arma::vec y, int extra_pairs = -1, int code_num = 3 ,int nfold = 10, int nLambda = 100);
#endif /* lasso_solver_hpp */
