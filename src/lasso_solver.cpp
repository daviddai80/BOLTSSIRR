//
//  lasso_solver.cpp
//  boost_new
//
//  Created by David on 2017/8/10.
//  Copyright © 2017年 David. All rights reserved.
//

#include "lasso_solver.hpp"

//void LassoPath(double* X, double* y_, int N, int P, double* lamseq, int num,sp_mat& beta, vec& covs)
void LassoPath(arma::mat& Xmat, arma::vec y, LassoFit& fit)
{
    int P = (int)Xmat.n_cols;
    int N = (int)Xmat.n_rows;
    Col<int> mm(P,fill::zeros);
    Col<int> active_set(P,fill::zeros);
    int nin = 0;
    double* beta0 = new double[P];
    for(int j = 0; j < P; j++){
      beta0[j] = 0;
    }
    double* beta = new double[P];
    double cov = 0;
    int num = (int)fit.lambdas.n_elem;

    rowvec x_var = sum(Xmat % Xmat, 0);

    double* yhat = new double[N];
    for(int i = 0; i < N; i++){
      yhat[i] = 0;
    }


    int total_iter = 0;
    for(int i = 0; i < num; i++){
      lassosolver(Xmat.memptr(), y.memptr(), N, P, fit.lambdas[i], x_var.memptr(), beta0, mm.memptr(), active_set.memptr(),  nin, beta,yhat, total_iter);
      double sum = 0;
      for(int j = 0; j < P; j++){
        if(fabs(beta[j]) > 1e-10 ){
          beta[j] /= fit.x_sd[j];
          sum += beta[j] * fit.x_center[j];
          fit.beta.at(j,i) = beta[j];
        }
      }
      cov = fit.mean_y - sum;
      fit.covs[i] = cov;
    }
    delete[] yhat;
    delete[] beta;
    delete[] beta0;
    fit.iter = total_iter;
    //cout << "Lasso Path Completion!" << endl;
}

void lassosolver(double* X, double* y_, int N, int P, double lambda, double* xvar, double* beta0, int* mm, int* active_set, int& nin, double* beta,double* yhat_ptr, int& total_iter)
{

  double tol = 1e-4;
  int iter_outer = 0;
  int iter_inner = 0;
  double diff_beta, diffabs;
  lambda *= N;
  while (true) {
    double dlx = 0.0;
    for (int j = 0; j < P; j ++) {
      double sum = 0;
      for(int i = 0; i < N; i++){
        sum += X[j* N + i]* (y_[i] - yhat_ptr[i]);
      }

      beta[j] = shrinkage(sum + xvar[j] * beta0[j], lambda) / xvar[j];
      diff_beta = beta0[j] - beta[j];
      diffabs  = fabs(diff_beta);
      if(diffabs > 0){
        for(int i = 0; i < N; i++){
          yhat_ptr[i] -= X[j* N + i] * diff_beta;
        }
        beta0[j] = beta[j];
        if(mm[j] == 0){
          active_set[nin] = j;
          nin += 1;
          mm[j] = nin;
        }
        dlx = fmax(dlx, diffabs);
      }
    }

    iter_outer += 1;
    if(dlx < tol){
      break;
    }

    while(true){
      dlx = 0.0;
      iter_inner += 1;
      for(int k = 0; k < nin; k++){
        int j = active_set[k];
        double sum = 0;
        for(int i = 0; i < N; i++){
          sum += X[j* N + i]* (y_[i] - yhat_ptr[i]);
        }
        beta[j] = shrinkage(sum + xvar[j] * beta0[j], lambda) / xvar[j];
        diff_beta = beta0[j] - beta[j];
        diffabs  = fabs(diff_beta);
        if(diffabs > 0){
          for(int i = 0; i < N; i++){
            yhat_ptr[i] -= X[j* N + i] * diff_beta;
          }
          beta0[j] = beta[j];
          dlx = fmax(dlx, diffabs);
        }
      }
      if(dlx < tol){
        break;
      }

    }

  }

  total_iter += iter_outer + iter_inner;

  //  cout << "Lasso solved (outer loops : inner loops ) = (" << iter_outer <<":"<<iter_inner<<")" << endl;

}


double shrinkage(double a, double kappa){
    return fmax(0, a - kappa) - fmax(0, - a - kappa);
}

double cal_max_lambda(arma::mat& X, arma::vec& y){
    int N = (int)X.n_rows;
    double  lambda = max(abs((y - mean(y)).t()  * X)) / N;
    return lambda;
}

// double cal_max_lambda(arma::mat& X, arma::vec& y, vector<InteractionSNPpairs>& pairs){
//   int N = (int)X.n_rows;
//   double  lambda = max(abs(y.t() * X)) / N;
//   for(int i = 0; i < pairs.size(); i++){
//
//   }
//   return lambda;
// }



LassoFit::LassoFit(){

};

LassoFit::~LassoFit(){

};

LassoFit::LassoFit(LassoFit& fit0){
  P = fit0.P;
  nLambda = fit0.nLambda;
  lambdas = fit0.lambdas;
  beta = fit0.beta;
  pairs = fit0.pairs;
  covs = fit0.covs;
  iter = fit0.iter;
 // min = fit0.min;
  mean_y = fit0.mean_y;
  x_center = fit0.x_center;
  x_sd = fit0.x_sd;
  x_var = fit0.x_var;
  lambda_min = fit0.lambda_min;
  cve = fit0.cve;
  cvse = fit0.cvse;
  sumq2 = fit0.sumq2;
}


LassoFit::LassoFit(int P, int nLambda){
    beta.resize(P, nLambda);
    beta.zeros();
    covs.resize(nLambda);
    covs.zeros();
}

LassoFit::LassoFit(mat& X, vec& y, GenoDiscreteMatrix& gdm, int nLambda, int extra_pairs, double chisq_thresh){
  vector<InteractionSNPpairs> top_pairs = getTopInteractionParis(&gdm, extra_pairs, chisq_thresh);
  this -> P =  (int)X.n_cols;
  if(top_pairs.size() > 0){
     int P0 = (int)X.n_cols;
     this -> P =  (int)(X.n_cols + top_pairs.size());
     X.resize(X.n_rows, P);
     class InteractionSNPpairs pr;
     for(int i = 0; i < top_pairs.size(); i++){
        pr = top_pairs[i];
        Col<double> iteract_i = X.col(pr.j1) %  X.col(pr.j2);
        X.col(P0 + i) = iteract_i;
     }
  }


  beta.resize(P, nLambda);
  beta.zeros();
  covs.resize(nLambda);
  covs.zeros();


  rowvec S(P);
  rowvec M = mean(X, 0);
  for(int i = 0; i < P; i++){
    X.col(i) -= M(i);
  }
  int N = (int)X.n_rows;
  S = sqrt(sum(X % X) / N);
  scale(X, S, M);
  double y_mean = as_scalar(mean(y));
  y -= y_mean;

  this -> mean_y = y_mean;
  this -> x_center = M;
  this -> x_sd = S;
  this -> pairs = top_pairs;
}

LassoFit::LassoFit(mat& X, vec& y, int nLambda){

  P = (int)X.n_cols;
  beta.resize(P, nLambda);
  beta.zeros();
  covs.resize(nLambda);
  covs.zeros();

  rowvec S(P);
  rowvec M = mean(X, 0);
  for(int i = 0; i < P; i++){
    X.col(i) = X.col(i) - M(i);
  }
  int N = (int)X.n_rows;
  S = sqrt(sum(X % X) / N);
  scale(X, S, M);
  double y_mean = as_scalar(mean(y));
  y -= y_mean;

  this -> mean_y = y_mean;
  this -> x_center = M;
  this -> x_sd = S;

}

LassoCVFit::LassoCVFit(){

}

LassoCVFit::~LassoCVFit(){

}
void LassoPathCV(arma::mat& X, arma::vec y, LassoCVFit& fit0, int nfold, int nLambda){
    int N = (int)X.n_rows;

     mat rX = X;
     vec ry = y;
     LassoFit fitlasso(X, y, nLambda);


     double max_lam = cal_max_lambda(X, y);
  //   cout <<"max_lambda=" << max_lam << endl;
     arma::vec lambdas = logspace<vec>(log10(max_lam), log10(0.0001* max_lam), nLambda);
     fitlasso.lambdas  = lambdas;

     LassoPath(X, y, fitlasso);
     mat y_hat(N, nLambda);

     Col<uword> indices = cross_valind(N, nfold);
     for (uword i = 1; i <= nfold; i++) {
     //  cout <<"fold i=" << i << endl;
       Col<uword> train_idx = find(indices != i);
       Col<uword> test_idx = find(indices == i);
       mat trainM = rX.rows(train_idx);
       vec ytrain = ry(train_idx);
       LassoFit fit_i(trainM, ytrain, nLambda);
       fit_i.lambdas = lambdas;
       mat testM = rX.rows(test_idx);
       LassoPath(trainM, ytrain, fit_i);
       mat result_i = predictAlongPath(testM, fit_i);
       y_hat.rows(test_idx) = result_i;
     }

     int min_idx = 1;
     double min_err = INFINITY;
     vec cvm(nLambda);
     vec auc(nLambda);
     vec cvsd(nLambda);
     vec cvup(nLambda);
     vec cvlo(nLambda);
     for(int i = 0; i < nLambda; i++){
        vec diff = y_hat.col(i) - ry;
        double err = sum(diff % diff) / N;
        cvm[i] = err;
        cvsd[i] = stddev(diff % diff) / sqrt(N);
        cvup[i] = cvm[i] + cvsd[i];
        cvlo[i] = cvm[i] - cvsd[i];
        auc[i] = calauc(ry, y_hat.col(i));
        if(err < min_err){
            min_err = err;
            min_idx = i + 1;
        }
     }

     fit0.cvm = cvm;
     fit0.cvsd = cvsd;
     fit0.cvlo = cvlo;
     fit0.cvup = cvup;
     fit0.index_min = min_idx;
     fit0.lambda_min = lambdas[min_idx];
     fit0.fit = fitlasso;
}



vec predict(arma::mat& X, double cov,  vec& beta){
    vec y_hat = cov + X * beta;
    return y_hat;
}

double calauc(arma::vec label, arma::vec pred){
  double auc = 0;
  double m = mean(label);
  vec label2 = label;
  label2(find(label >= m)).fill(1);
  label2(find(label <= m)).fill(0);
  label = label2;
  uword N = pred.size();
  uword N_pos = sum(label);
  uvec  idx = sort_index(pred,"descend");
  vec above = linspace<vec>(1, N, N) - cumsum(label(idx));
  auc = (1 - sum(above % label(idx)) / (N_pos * (N-N_pos)));
  auc = auc > 0.5?auc:(1-auc);
  return auc;
}

mat predictAlongPath(arma::mat& X,  LassoFit& fit){
  //  sp_mat& beta, vec& covs;
    int N = (int)X.n_rows;
    int nLambda = (int)fit.covs.size();
    mat result(N,nLambda);

    if(fit.pairs.size() > 0){
      int P0 = X.n_cols;
      X.resize(X.n_rows, P0 + fit.pairs.size());
      class InteractionSNPpairs pr;
      for(int i = 0; i < fit.pairs.size(); i++){
        pr = fit.pairs[i];
        Col<double> iteract_i = X.col(pr.j1) %  X.col(pr.j2);
        X.col(P0 + i) = iteract_i;
      }
    }
    for(int i = 0; i < nLambda; i++){
        double cov = fit.covs[i];
        vec beta_i(fit.beta.col(i));
      //  cout <<"beta=" << sum(beta_i) << endl;
        result.col(i) = predict(X, cov, beta_i);
    }
    return result;
}

Mat<double> cal_means(Mat<double>& X_mat, int N){
    mat sumX = sum(X_mat,0);;
    mat SZX =  sumX * 1.0 / N ;//conv_to<mat>::from(mean(X, 0)); //column means of X
    return SZX;
}

void centering(Mat<double>& X_mat, vec& y, double& mean_y, mat& SZX, int N, int P, mat* SZX_, double* y_mean_){
    SZX = (SZX_ != nullptr) ? *SZX_ : cal_means(X_mat, N);
    mean_y = (y_mean_ != nullptr) ? *y_mean_ : as_scalar(mean(y)) ;
    y -= mean_y;
    double* x_mean = SZX.memptr();
    for(int i = 0; i < P; i++){
        X_mat.col(i) -= x_mean[i];
    }
}

/**shuffle the index for cross validation*/
arma::Col<uword> cross_valind(arma::uword N, arma::uword nfold){
    arma::Col<uword> indices(N);
    arma_rng::set_seed_random();
    arma::Col<uword> vec_n = arma::shuffle(arma::linspace < arma::Col <uword> >(1, N, N));
    indices.fill(nfold);
    for(uword n = 1; n <= nfold - 1; n++){
        arma::Col<uword> in = vec_n.rows((n-1)*N/nfold,n*N/nfold - 1);
        indices.elem(in - 1 ).fill(n);
    }
    return indices;
}

void scale(arma::mat& X, arma::rowvec& S, arma::rowvec& M){
  uword P = X.n_cols;
  for(uword i = 0; i < P; i ++){
    X.col(i) = X.col(i)  / S(i);
  }
}
