//
//  coding.cpp
//  boost_new
//
//  Created by David on 2017/7/9.
//  Copyright © 2017年 David. All rights reserved.
//

#include "coding.hpp"
#include <fstream>
#include <string>
using namespace std;



unsigned char SingleCodes::wordbits[65536];
int codes::code_num = 3; // 3 for 0,1,2
int codes::thread_num = 1;
long codes::max_num_of_pair = 100000;

float* SnpCodes::marginal_distr;
int SnpCodes::P;

void SingleCodes::init(){
    // precompute the wordbits (a global variable)
    for (int i = 0; i < 65536; i++)
    {
        wordbits[i] = bitCount(i);
    }
}

int SingleCodes::bitCount(uint64_t i)
{
    i = i - ((i >> 1) & 0x5555555555555555);
    i = (i & 0x3333333333333333) + ((i >> 2) & 0x3333333333333333);
    i = (i + (i >> 4)) & 0x0f0f0f0f0f0f0f0f;
    i = i + (i >> 8);
    i = i + (i >> 16);
    i = i + (i >> 32);
    return (int)i & 0x7f;
}

SingleCodes::SingleCodes(int n){
    this -> n = n;
    code = new uint64_t[n];
    for(int i = 0; i < n; i++){
        code[i] = 0;
    }
}

SingleCodes::~SingleCodes(){
    if(code != NULL){
        delete[] code;
    }
}

SingleCodes SingleCodes::operator | (const SingleCodes& b)
{
    SingleCodes singlecode(this -> n);
    for(int i = 0; i < this -> n; i++){
        singlecode.code[i] = this -> code[i] | b.code[i];
    }
    return singlecode;
}

SingleCodes SingleCodes::operator & (const SingleCodes& b)
{
    SingleCodes singlecode(this -> n);
    for(int i = 0; i < this -> n; i++){
        singlecode.code[i] = this -> code[i] & b.code[i];
    }
    return singlecode;
}

int SingleCodes::pop_count(uint64_t i ){
    return( wordbits[i&0xFFFF] + wordbits[(i>>16)&0xFFFF] + wordbits[(i>>32)&0xFFFF] + wordbits[i>>48]);
}

int SingleCodes::count(){
    int n_count = 0;
    for(int i = 0; i < n; i++){
        n_count += pop_count(code[i]);//bitCount(code[i]); //
    }
    this -> cnt = n_count;
    return n_count;
}

void add_vector(float* v1, float* v2, int num, float* dest){
    for(int i = 0; i < num; i++){
        dest[i] = v1[i] + v2[i];
    }
}

/***************** implementation of SnpCodes *******************/
SnpCodes::SnpCodes(int n, int type, int index){
    codes = new SingleCodes * [codes::code_num];
    for(int i = 0; i < codes::code_num; i++){
        codes[i] = new SingleCodes(n);
    }
    this -> type = type;
    this -> index = index;

}

void SnpCodes::init(int P_, int code_num){
    P = P_;
    marginal_distr = new float[codes::code_num * 2 * P];
}


void SnpCodes::coding(int* col, int N){
    for(int i = 0; i < N; i++){
        codes[col[i]] -> code[i / LengthLongType] |= (mask << (i % LengthLongType));
    }

    for(int m = 0; m < codes::code_num; m++){
       marginal_distr[2 * codes::code_num * index + type * codes::code_num + m] = codes[m] -> count();
    }
}

// void SnpCodes::coding_no_zero(int* col, int N){
//   for(int i = 0; i < N; i++){
//     codes[col[i]] -> code[i / LengthLongType] |= (mask << (i % LengthLongType));
//   }
//
//   for(int m = 1; m < codes::code_num; m++){
//     marginal_distr[2 * codes::code_num * index + type * codes::code_num + m] = codes[m] -> count();
//   }
// }

SnpCodes::~SnpCodes(){
    for(int i = 0; i < codes::code_num; i++){
        delete codes[i];
    }
    delete[] codes;
}

void SnpCodes::jointDistrWith(SnpCodes* snpcodes, Mu&  GenoDistr, int code_num){
    int i1 = 0;
    int i2 = 0;
    for(i1 = 0; i1 < code_num - 1; i1 ++){
        for(i2 = 0; i2 < code_num - 1; i2 ++){
            GenoDistr.data[i1 * code_num + i2] = (*this -> codes[i1] & * snpcodes -> codes[i2]).count();
        }
    }
    i2 = code_num - 1;
    int i = 0;
    for(i1 = 0; i1 < code_num - 1; i1 ++){
        int sum = 0;
        for(i = 0; i < code_num - 1; i++){
            sum += GenoDistr.data[i1 * code_num + i];
        }
        GenoDistr.data[i1 * code_num + i2] = marginal_distr[2 * code_num * index + type * code_num + i1] - sum;
    }

    i1 = code_num - 1;
    for(i2 = 0; i2 < code_num; i2 ++ ){
        int sum = 0;
        for(i = 0; i < code_num - 1; i++){
            sum += GenoDistr.data[i * code_num + i2];
        }
        GenoDistr.data[i1 * code_num + i2] = marginal_distr[2 * code_num * snpcodes->index + type * code_num + i2] - sum;
    }

}

/***************************************************************************/
// GenoDiscreteMatrix(mat& m, int code_num){
//        this -> N = m.n_rows;
//        this -> P = m.n_cols;
//        this -> data = new int[N * P];
//        for(int i = 0; i < N * P; i++){
//
//        }
// }
GenoDiscreteMatrix::GenoDiscreteMatrix(int N, int P, bool rdm){
    this -> N = N;
    this -> P = P;
    this -> data = new int[N * P];
    this -> y = new int[N];
    this -> yc = new double[N];
  //  this -> Xc = new double[N * P];
    if(rdm){
        mat m(N, P, fill::randu);
        mat y_(N,1, fill::randu);
        for(int i = 0; i < N*P; i++){
            if(m[i] < 0.4){
                data[i] = 0;
            }else if(m[i] < 0.7){
                data[i] = 1;
            }else{
                data[i] = 2;
            }
        }

        for(int i = 0; i < N; i ++){
            if(y_[i] < 0.55){
                y[i] = 0;
            }else{
                y[i] = 1;
            }
        }
    }

}

void GenoDiscreteMatrix::sep2CCMatrix(GenoDiscreteMatrix * & matrix_case, GenoDiscreteMatrix * & matrix_ctrl){
    int code_num = codes::code_num;
    int ncase = 0;
    int nctrl = 0;
    for(int i = 0; i < N; i ++){
        ncase += y[i];
    }
    nctrl = N - ncase;
    matrix_case = new GenoDiscreteMatrix(ncase, P);
    matrix_ctrl = new GenoDiscreteMatrix(nctrl, P);

    int icase = 0;
    int ictrl = 0;
    int* numbers = new int[2 * code_num * P];
    for(int i = 0; i < 2 * code_num * P; i++){
        numbers[i] = 0;
    }
    for(int j = 0; j < P; j++)
    {
        for(int i = 0; i < N; i ++)
        {
            if(y[i] == 1){
                matrix_case -> data[icase ++] = this -> data[j * N  + i];
                numbers[code_num * P + code_num * j + (this -> data[j * N  + i])] ++;
            }
            else
            {
                matrix_ctrl -> data[ictrl ++] = this -> data[j * N  + i];
                numbers[code_num * j + (this -> data[j * N  + i])] ++;
            }
        }
    }

}

GenoDiscreteMatrix::~GenoDiscreteMatrix(){
    if(y != NULL)
        delete []y;
    if(data != NULL)
        delete []data;

//    if(yc != NULL)
//        delete []yc;

//    if(Xc != NULL)
//        delete []Xc;
//    if(data != NULL)
//        delete []data;
}
/*****************************************************************************/
void GenoBitMatrix::CalculateMarginalEntropy(double * MarginalEntropySNP, double * MarginalEntropySNP_Y){
    int i1, i2;
    double ptmp;
    int code_num = codes::code_num;
    for(i1 = 0; i1 < P; i1 ++){
        for(i2 = 0; i2 < code_num; i2 ++){
            int tmp = marginal_count[i1 * code_num + i2];
            if ( tmp > 0)
            {
                ptmp = tmp / N;
                MarginalEntropySNP[i1] += -(ptmp)*log(ptmp);
            }
        }
    }

}

void GenoBitMatrix::cal_taointMeasure(const Mu& mutmp, const Mu& mucc,double& tau, double& InteractionMeasure){
    double ptmp1, ptmp2;
    tau = 0.0;
    InteractionMeasure = 0.0;
    for (int i = 0; i < code_num; i++) {
        for (int j = 0; j < code_num; j++) {
            for (int k = 0; k < 2; k ++) { //index for C
                ptmp1 = mucc.data[k * code_num * code_num + i * code_num + j] / N;
                if(ptmp1 > 0) InteractionMeasure += ptmp1 * log(ptmp1);
                ptmp2 = mutmp.data[k * code_num * code_num + j * code_num + i] / N;
                if (ptmp2>0){
                    InteractionMeasure += -ptmp1*log(ptmp2);
                    tau += ptmp2;
                }
            }
        }
    }

    InteractionMeasure = (InteractionMeasure+log(tau)) * N *2;


}
void GenoBitMatrix::cal_taointMeasure(const Mu& pab, const Mu& pbc, const Mu& pca,const Mu& mucc, double& tau, double& InteractionMeasure){
    double ptmp1, ptmp2;
    tau = 0.0;
    InteractionMeasure = 0.0;
    for (int i = 0; i < code_num; i++) {
        for (int j = 0; j < code_num; j++) {
            for (int k = 0; k < 2; k ++) { //index for C
                ptmp1 = mucc.data[k * code_num * code_num + i * code_num + j] / N;
                if(ptmp1 > 0) InteractionMeasure += ptmp1 * log(ptmp1);
                ptmp2 = pab.data[code_num * i + j] * pbc.data[code_num* k + j] * pca.data[code_num*k + i];
                if (ptmp2>0){
                    InteractionMeasure += -ptmp1*log(ptmp2);
                    tau += ptmp2;
                }
               // cout << ptmp1 << " : " << ptmp2 <<" | " << InteractionMeasure << endl;
            }
        }
    }

    InteractionMeasure = (InteractionMeasure+log(tau)) * N * 2;
}


double GenoBitMatrix::CalInteractByThread(int thread_id){
    code_num = codes::code_num;
    Mu mucc(code_num, 2 * code_num,0);
    Mu pab(code_num, code_num,0);
    Mu pbc(code_num, 2, 0);
    Mu pca(code_num, 2, 0);
    double interatction = 0;
    InteractionSNPpairs pair0 = pairs[thread_id * 2];
    InteractionSNPpairs pair1 = pairs[thread_id * 2 + 1];

    int cj1 = pair0.j1;
    int cj2 = pair0.j2;
    interatction = CalInteractFor(cj1, cj2, mucc, pab, pbc, pca);
    InteractionSNPpairs Pr(cj1, cj2);
    Pr.InteractionMeasure = interatction;
    thread_results[thread_id].push_back(Pr);

    ofstream myfile;
    myfile.open("output_" + std::to_string(thread_id) +".log");
    int idx0 = 0;
    int idx1 = 0;
    time_t t1, t2;
    time(&t1);
    int gap = 500000;
    while (true) {
        if(cj1 == pair1.j1 && cj2 == pair1.j2){
            break;
        }else if(cj2 == P - 1){
            cj1 ++;
            cj2 = cj1 + 1;
        }else{
            cj2 ++;
        }
        interatction = CalInteractFor(cj1, cj2, mucc, pab, pbc, pca);
        if(interatction >= thresholdRecord)
        {
           idx1 ++;
           InteractionSNPpairs Pr(cj1, cj2);
           Pr.InteractionMeasure = interatction;
           thread_results[thread_id].push_back(Pr);
        }
        idx0 ++;
        if(idx0 % gap == 0){
          time(&t2);
          myfile << "Time for the last "<< gap <<" iterations is " << difftime(t2, t1) << " seconds, idx:" << idx0 <<" " << idx1  << endl;
          time(&t1);
        }


    }
    myfile.close();
    return interatction;

}



double GenoBitMatrix::CalInteractFor(int j1, int j2, Mu& mucc, Mu& pab, Mu& pbc, Mu& pca){
    //j1,j2:index for SNPs
    Mu mucase(mucc.data, code_num, code_num);
    Mu muctrl(mucc.data + code_num * code_num, code_num, code_num);
    codes_case[j1] -> jointDistrWith(codes_case[j2], mucase, code_num);
    codes_ctrl[j1] -> jointDistrWith(codes_ctrl[j2], muctrl,code_num);
    Mu mar_j2(marginal_count + j2 * code_num, code_num, 1);
    plusdivide_mu(mucase,muctrl, mar_j2, pab);

    /**to calculate PBC ***/
    Mu geno_j2(SnpCodes::marginal_distr + j2 * 2 * code_num, code_num, 2);
    divide_mu(geno_j2, numbers, pbc);

    /**to calculate PCA ***/
    Mu geno_j1(SnpCodes::marginal_distr + j1 * 2 * code_num, code_num, 2);
    Mu mar_j1(marginal_count + j1 * code_num, code_num, 1);
    divide_mu(geno_j1, mar_j1, pca);

    double tau, InteractionMeasure;
    InteractionMeasure = 0.0;
    cal_taointMeasure(pab, pbc, pca, mucc, tau, InteractionMeasure);

    return InteractionMeasure;

}

InteractionSNPpairs* GenoBitMatrix::next(){

    std::lock_guard<std::mutex> lock(_mtx);
    if(cj1 == P-2 && cj2 == P-1){
        return nullptr;
    }else if(cj2 == P-1){
        cj1 ++;
        cj2 = cj1 + 1;
    }else{
        cj2 ++;
    }
    InteractionSNPpairs* pair = new InteractionSNPpairs(cj1, cj2);
    return pair;
}

void GenoBitMatrix::recalInteractCont(vector<InteractionSNPpairs>& pairs){
    arma::mat X(this->Xc, this->N, this->P, false);
    arma::vec y(this->yc, this->N, false);
    int pair_size = (int)pairs.size();
    int Pm = 0; //indicate the max index

    for(int i = 0; i < pair_size; i++){
        if(pairs[i].j2 > Pm){
           Pm = pairs[i].j2;
        }
    }

    Pm += 1;
    std::unordered_set<int> nonzero;
    for(int i = 0; i < pair_size; i++){
       nonzero.insert(pairs[i].j1);
       nonzero.insert(pairs[i].j2);
    }

    vec x2(Pm, fill::zeros);
    vec xy(Pm,  fill::zeros);
    double y_mean = (double)mean(y);
    vec y_ = y - y_mean;
    double xi_mean = 0;

    for (const int& j: nonzero) {
       xi_mean = mean(X.col(j));
       X.col(j) -= xi_mean;
       x2.at(j) = sum(X.col(j) % X.col(j));
       xy.at(j) = sum(y_ % X.col(j));
    }



    double* mat22 = new double[4];
    double* inv_mat22 = new double[4];
    double* add3 = new double[3];
    double* inv_mat33 = new double[9];
    double* beta3 = new double[3];


    int j1 = 0;
    int j2 = 0;
    double LR = 0;
    for(int i = 0; i < pair_size; i++){
       j1 = pairs[i].j1;
       j2 = pairs[i].j2;
       LR = cal_interact_fast(X, y, this->N, j1, j2,
                              mat22, inv_mat22, add3, inv_mat33, beta3,
                              x2, xy, y_mean, y_);
       pairs[i].InteractionMeasure = LR;
    }
    delete[] mat22;
    delete[] inv_mat22;
    delete[] inv_mat33;
    delete[] beta3;
    delete[] add3;
}

void GenoBitMatrix::recalInteractDis(vector<InteractionSNPpairs>& pairs){
  ones = new float[code_num * code_num * 2];
  zeros = new float[code_num * code_num * 2];
  for (int i = 0; i < code_num * code_num * 2; i++) {
    ones[i] = 1;
    zeros[i] = 0;
  }

  double InteractionMeasure = 0.0;
  Mu mutcc(code_num, 2 * code_num,0);
  Mu mut_ca(mutcc.data, code_num, code_num);
  Mu mut_ct(mutcc.data + code_num * code_num, code_num, code_num);
  Mu mu0tcc(code_num, 2 * code_num,0);
  double muError = 0;
  Mu mu_ij(code_num, code_num, 1);//mu matrix for case
  Mu mui_ca(code_num, 1, 1);//mui matrix for case
  Mu mui_ct(code_num, 1, 1);//mui matrix for case
  Mu muj_ca(1, code_num, 1);//muj matrix for case
  Mu muj_ct(1, code_num, 1);//muj matrix for case
  Mu mucc(code_num, 2 * code_num,0);
  Mu mucase(mucc.data, code_num, code_num);
  Mu muctrl(mucc.data + code_num * code_num, code_num, code_num);

  double  tau;
  class InteractionSNPpairs pr;
  for(int i = 0; i < pairs.size(); i++)
  {
      pr = pairs[i];
      int j1 = pr.j1;
      int j2 = pr.j2;

      codes_ctrl[j1] -> jointDistrWith(codes_ctrl[j2], muctrl, code_num);
      codes_case[j1] -> jointDistrWith(codes_case[j2], mucase, code_num);

      Mu muComb = mucase + muctrl;

      //to initialize mutcc and mu0tcc
      mutcc.copy(ones, code_num * code_num * 2);
      mu0tcc.copy(zeros, code_num * code_num * 2);
      muError = (mutcc - mu0tcc).sum_abs();
      Mu muca_dim2 = mucase.marginal_dim2();
      Mu muct_dim2 = muctrl.marginal_dim2();
      Mu muca_dim1 = mucase.marginal_dim1();
      Mu muct_dim1 = muctrl.marginal_dim1();
      while(muError > 0.001){
        mu0tcc.clone(mutcc);
        mu_ij = mut_ca + mut_ct;
        update_mutmp(mut_ca, muComb, mu_ij);
        update_mutmp(mut_ct, muComb, mu_ij);

        /*** mui_ca, mui_ct **/
        mui_ca = mut_ca.marginal_dim1();
        mui_ct = mut_ct.marginal_dim1();

        update_mutmp(mut_ca, muca_dim2 , mui_ca);
        update_mutmp(mut_ct, muct_dim2 , mui_ct);
        /*** mui_ca, mui_ct **/
        muj_ca = mut_ca.marginal_dim2();
        muj_ct = mut_ct.marginal_dim2();
        update_mutmp(mut_ca, muca_dim1 , muj_ca);
        update_mutmp(mut_ct, muct_dim1 , muj_ct);

        muError = (mutcc - mu0tcc).sum_abs();//(mut_ca - mu0t_ca).sum_abs() + (mut_ct - mu0t_ct).sum_abs();

      }

      cal_taointMeasure(mutcc, mucc, tau, InteractionMeasure);
      pr.InteractionMeasure = InteractionMeasure;
  }
  delete[] ones;
  delete[] zeros;

}
void GenoBitMatrix::CalGenoJointDistr(){
    time_t timep, timee;
    time(&timep); /*获取time_t类型的当前时间*/
   // cout <<"Begining calculating the joint distributions! " << timep << endl;
    thread_num = codes::thread_num;
    std::thread* threads = new (std::nothrow) std::thread[thread_num];
    for(int thread_id = 0; thread_id < thread_num; thread_id++){
        threads[ thread_id ] = std::thread(&GenoBitMatrix::CalInteractByThread,this, thread_id);
    }

    for(uword k = 0; k < thread_num; k++){
        threads[k].join();
    }

    vector<InteractionSNPpairs> pairs(0);
    for(int i = 0; i < thread_num; i++){
      int cnt = (int)this -> thread_results[i].size();
      class InteractionSNPpairs pr;
      for(int k = 0; k < cnt; k++){
        pr = this -> thread_results[i][k];
        pairs.push_back(pr);
      }
    }

    std::sort(pairs.begin(), pairs.end(), compare_pairs);
    int max_pn = max(this->P, this -> N);
    if(max_pn < pairs.size()){
       pairs.resize(max_pn);
    }

    if(this->discrete == 1){
       this -> recalInteractDis(pairs);
    }else{
       this -> recalInteractCont(pairs);
    }

    std::sort(pairs.begin(), pairs.end(), compare_pairs);
    this -> final_results = pairs;
    time(&timee); /*获取time_t类型的当前时间*/
    delete[] threads;
    threads = nullptr;

}


GenoBitMatrix::GenoBitMatrix(GenoDiscreteMatrix& gdm){
    this -> code_num = codes::code_num;
    this -> thread_num = codes::thread_num;
    this -> N = gdm.N;
    this -> P = gdm.P;
    this -> y = new int[N];
    this -> discrete = gdm.discrete;
    memcpy(this -> y, gdm.y, sizeof(int) * N);
    ncase = 0;
    nctrl = 0;
    for(int i = 0; i < N; i++){
        ncase += y[i];
    }
    nctrl = N - ncase;
    if(numbers.dim1 * numbers.dim2 >= 2){
        numbers.data[0] = ncase;
        numbers.data[1] = nctrl;
    }
    int nlongintcase = ceil( ((double) ncase) / LengthLongType);
    int nlongintctrl = ceil( ((double) nctrl) / LengthLongType);
    codes_ctrl = new SnpCodes * [P];
    codes_case = new SnpCodes * [P];
    GenoDiscreteMatrix* matrix_case;
    GenoDiscreteMatrix* matrix_ctrl;
    gdm.sep2CCMatrix(matrix_case, matrix_ctrl);
    marginal_count = new float [P * code_num];
    float* marginal_dest_i = new float[code_num];
    SnpCodes::init(P, code_num);
    for(int i = 0; i < P; i ++){
        codes_case[i] = new SnpCodes(nlongintcase, 0, i);
        codes_ctrl[i] = new SnpCodes(nlongintctrl, 1, i);
        codes_case[i] -> coding(matrix_case -> data + (i * ncase), ncase);
        codes_ctrl[i] -> coding(matrix_ctrl -> data + (i * nctrl), nctrl);
        add_vector(SnpCodes::marginal_distr + 2 * code_num * i + code_num, SnpCodes::marginal_distr + 2 * code_num * i, code_num, marginal_dest_i);
        memcpy(marginal_count + i * code_num, marginal_dest_i, sizeof(float) * code_num);
    }
    thread_results = new vector<InteractionSNPpairs>[this -> thread_num];
    pairs = new InteractionSNPpairs[this -> thread_num * 2];
    cal_thread_bound((long long)P, this -> thread_num, pairs);
    this -> yc = gdm.yc;
    this -> Xc = gdm.Xc;
    delete[] marginal_dest_i;

}

GenoBitMatrix::GenoBitMatrix(arma::Mat<int>& amat){
    this -> N = (int)amat.n_rows;//gdm.N;
    this -> P = (int)(amat.n_cols - 1);//gdm.P;
    this -> y = new int[N];
    memcpy(this -> y, amat.memptr(), sizeof(int) * N);
    ncase = 0;
    nctrl = 0;
    for(int i = 0; i < N; i++){
        ncase += y[i];
    }
    nctrl = N - ncase;
    Col<int> ycol(this -> y,N,1);

    if(numbers.dim1 * numbers.dim2 >= 2){
        numbers.data[0] = ncase;
        numbers.data[1] = nctrl;
    }

    int nlongintcase = ceil( ((double) ncase) / LengthLongType);
    int nlongintctrl = ceil( ((double) nctrl) / LengthLongType);

    codes_ctrl = new SnpCodes * [P];
    codes_case = new SnpCodes * [P];


    arma::Mat<int>* matrix_case = nullptr;
    arma::Mat<int>* matrix_ctrl = nullptr;

    matrix_case = new arma::Mat<int>(ncase, P);
    uvec q1 = find(ycol == 1);
    for(int i = 0; i < ncase; i++){
        for(int j = 0; j < P; j++){
            matrix_case->at(i, j) = amat(q1(i),j + 1);
        }
    }

    matrix_ctrl = new arma::Mat<int>(nctrl, P);
    q1 = find(ycol == 0);
    for(int i = 0; i < nctrl; i++){
        for(int j = 0; j < P; j++){
            matrix_ctrl -> at(i, j) = amat(q1(i),j + 1);
        }
    }


    marginal_count = new float [P * code_num];

    float* marginal_dest_i = new float[code_num];
    SnpCodes::init(P, code_num);

    for(int i = 0; i < P; i ++){
        codes_case[i] = new SnpCodes(nlongintcase, 0, i);
        codes_ctrl[i] = new SnpCodes(nlongintctrl, 1, i);

        codes_case[i] -> coding(matrix_case -> memptr() + (i * ncase), ncase);
        codes_ctrl[i] -> coding(matrix_ctrl -> memptr() + (i * nctrl), nctrl);
        add_vector(SnpCodes::marginal_distr + 2 * code_num * i + code_num, SnpCodes::marginal_distr + 2 * code_num * i, code_num, marginal_dest_i);
        memcpy(marginal_count + i * code_num, marginal_dest_i, sizeof(float) * code_num);
        //    codes_ctrl[i] -> marginal_distr
        //    codes_case[i] -> marginal_distr
    }
    thread_results = new vector<InteractionSNPpairs>[thread_num];
    pairs = new InteractionSNPpairs[thread_num * 2];
    cal_thread_bound(P, thread_num, pairs);
    //    for(int i = 0; i < thread_num*2; i++){
    //        cout << pairs[i].j1 << " " <<pairs[i].j2<< endl;
    //    //    cout << pairs[i ].j1 << " " <<pairs[i].j2<< endl;
    //
    //    }
    delete matrix_case;
    delete matrix_ctrl;
    delete[] marginal_dest_i;

}

/****************************************************************************************/
Mu::Mu(){
    //default number of elements is 2,
    dim1 = 1;
    dim2 = 2;
    data = new float[dim1*dim2];
}
Mu::Mu(int d1, int d2, float data0){
    dim1 = d1;
    dim2 = d2;
    data = new float[d1 * d2];
    for(int i = 0; i < d1 * d2; i++){
        data[i] = data0;
    }
}

Mu::Mu(float* data0, int d1, int d2){
    dim1 = d1;
    dim2 = d2;
    data = data0;
}

void Mu::clean(){
    if(data != NULL){
       delete[] data;
    }
}

Mu::~Mu(){
//    if(data != NULL){
//        delete[] data;
//    }

//    if(marginal_dim1 != NULL){
//        delete[] marginal_dim1;
//    }
//    if(marginal_dim2 != NULL){
//        delete[] marginal_dim1;
//    }
}
//Mu::Mu(int* data0, int d1, int d2){
//
//}

Mu& Mu::operator=(const Mu& other){
    if (this != &other){
        this -> dim1 = other.dim1;
        this -> dim2 = other.dim2;
        delete[] this -> data;
        this -> data = new float[dim1 * dim2];
        memcpy(this -> data, other.data, dim1 * dim2 * sizeof(float));
    }
    return *this;
}

Mu::Mu(const Mu& mu){
    this -> dim1 = mu.dim1;
    this -> dim2 = mu.dim2;
    this -> data = new float[dim1 * dim2];
//    for (int i = 0; i < dim1 * dim2; i++) {
//        this -> data[i] = mu.data[i];
//    }
    memcpy(this -> data, mu.data, sizeof(float) * (mu.dim1 * mu.dim2));

}

void Mu::clone(const Mu& mu){
    if(this->dim1 == mu.dim1 && this->dim2 == mu.dim2){
        memcpy(this -> data, mu.data, sizeof(float) * (mu.dim1 * mu.dim2));
    }else{
        delete[] this -> data;
        this -> dim1 = mu.dim1;
        this -> dim2 = mu.dim2;
        this -> data = new float[dim1 * dim2];
        memcpy(this -> data, mu.data, sizeof(float) * (mu.dim1 * mu.dim2));
    }
}

void Mu::copy(float* data, int size){
    memcpy(this -> data, data, sizeof(float) * size);
}

Mu Mu::marginal_dim1()
{//update the marginal info for dim1
//    if(marginal_dim1 == NULL){
//        marginal_dim1 = new float[dim1];
//    }
    Mu mu_m1(dim1, 1, 0);
    for(int i = 0; i < dim1; i++){
        double marginal = 0;
        for (int j = 0; j < dim2; j++) {
            marginal += data[j * dim1 + i];
        }
        mu_m1.data[i] = marginal;
    }
    return mu_m1;
}

Mu Mu::marginal_dim2(){//update the marginal info for dim2
//    if(marginal_dim2 == NULL){
//        marginal_dim2 = new float[dim2];
//    }
//
    Mu mu_m2(1, dim2, 0);
    for(int j = 0; j < dim2; j++){
        double marginal = 0;
        for (int i = 0; i < dim1; i++) {
            marginal += data[j * dim1 + i];
        }
        mu_m2.data[j] = marginal;
    }
    return mu_m2;
}

Mu Mu::operator - (const Mu& other){
    if(this->dim1 != other.dim1 || this -> dim2 != other.dim2){
        cout << "inconsitence dimensions of the two matrix" << endl;
        Mu m(1,1,-1);
    }
    Mu mu_out(dim1, dim2,1);
    for (int i = 0; i < dim1 * dim2; i++) {
        mu_out.data[i] = this -> data[i] - other.data[i];
    }
    return mu_out;
}

Mu& Mu::operator /= (const Mu& other){
    if(this->dim1 == other.dim1 && this -> dim2 == other.dim2){
      //  type = 0;
        for (int i = 0; i < dim1 * dim2; i++) {
            this -> data[i] = this -> data[i] / other.data[i];
        }
    }else if( this-> dim2 == other.dim2 && other.dim1 == 1){
        for(int i = 0; i < dim1; i++){
            for(int j = 0; j < dim2; j++){
                this -> data[j * dim1 + i] /= other.data[j];
            }
        }
    }else if( this-> dim1 == other.dim1 && other.dim2 == 1){
        for(int i = 0; i < dim1; i++){
            for(int j = 0; j < dim2; j++){
                this -> data[j * dim1 + i] /= other.data[i];
            }
        }
    }else{
         cout << "inconsitence dimensions of the two matrix" << endl;
    }

    return *this;
}

ostream& operator<<(ostream& os, const Mu& mu){
    for(int i = 0; i < mu.dim1; i++){
        for (int j = 0; j < mu.dim2; j++) {
            os << mu.data[j * mu.dim1 + i] << " ";
        }
        os << endl;
    }
    os << endl;
    return os;
}

void plusdivide_mu(const Mu& mu1, const Mu& mu2, const Mu& other, Mu& res){
    if(mu1.dim1 == other.dim1 && mu1.dim2 == other.dim2){
        //  type = 0;
        for (int i = 0; i < mu1.dim1 * mu1.dim2; i++) {
            res.data[i] = (mu1.data[i] + mu2.data[i]) / other.data[i];
        }
    }else if( mu1.dim2 == other.dim2 && other.dim1 == 1){
        for(int i = 0; i < mu1.dim1; i++){
            for(int j = 0; j < mu1.dim2; j++){
                res.data[j * mu1.dim1 + i] =  (mu1.data[j * mu1.dim1 + i] + mu2.data[j * mu1.dim1 + i] ) / other.data[j];
            }
        }
    }else if( mu1.dim1 == other.dim1 && other.dim2 == 1){
        for(int i = 0; i < mu1.dim1; i++){
            for(int j = 0; j < mu1.dim2; j++){
                res.data[j * mu1.dim1 + i] = (mu1.data[j * mu1.dim1 + i] +mu2.data[j * mu1.dim1 + i]) / other.data[i];
            }
        }
    }else{
        cout << "inconsitence dimensions of the two matrix" << endl;
    }
}

void divide_mu(const Mu& mu1, const Mu& other, Mu& res){
    if(mu1.dim1 == other.dim1 && mu1.dim2 == other.dim2){
        //  type = 0;
        for (int i = 0; i < mu1.dim1 * mu1.dim2; i++) {
            res.data[i] = mu1.data[i] / other.data[i];
        }
    }else if( mu1.dim2 == other.dim2 && other.dim1 == 1){
        for(int i = 0; i < mu1.dim1; i++){
            for(int j = 0; j < mu1.dim2; j++){
                res.data[j * mu1.dim1 + i] =  mu1.data[j * mu1.dim1 + i] / other.data[j];
            }
        }
    }else if( mu1.dim1 == other.dim1 && other.dim2 == 1){
        for(int i = 0; i < mu1.dim1; i++){
            for(int j = 0; j < mu1.dim2; j++){
                res.data[j * mu1.dim1 + i] = mu1.data[j * mu1.dim1 + i] / other.data[i];
            }
        }
    }else{
        cout << "inconsitence dimensions of the two matrix" << endl;
    }
}


void update_mutmp(Mu& mutmp, Mu& multiplier, Mu& divisor){
    for (int i = 0; i < mutmp.dim1; i++) {
        for (int j = 0; j < mutmp.dim2; j++){
            if(multiplier.dim1 > 1 && multiplier.dim2 > 1){
                float v = divisor.data[j * mutmp.dim1 + i];
                if(v > 0)
                  mutmp.data[j * mutmp.dim1 + i] *= (multiplier.data[i * mutmp.dim2 + j] / v);
                else
                 mutmp.data[j * mutmp.dim1 + i] = 0;
            }else if(multiplier.dim2 == 1){
                float v = divisor.data[j];
                if(v > 0)
                    mutmp.data[j * mutmp.dim1 + i] *= (multiplier.data[j] / v);
                else
                    mutmp.data[j * mutmp.dim1 + i] = 0;

            }else if(multiplier.dim1 == 1){
                float v = divisor.data[i];
                if(v > 0)
                    mutmp.data[j * mutmp.dim1 + i] *= (multiplier.data[i] / v);
                else
                    mutmp.data[j * mutmp.dim1 + i] = 0;

            }
        }
    }

}

Mu Mu::operator / (const Mu& other){
    Mu mu_out(dim1, dim2,1);
    if(this->dim1 == other.dim1 && this -> dim2 == other.dim2){
        //  type = 0;
        for (int i = 0; i < dim1 * dim2; i++) {
            mu_out.data[i] = this -> data[i] / other.data[i];
        }
    }else if( this-> dim2 == other.dim2 && other.dim1 == 1){
        for(int i = 0; i < dim1; i++){
            for(int j = 0; j < dim2; j++){
               mu_out.data[j * dim1 + i] =  this -> data[j * dim1 + i] / other.data[j];
            }
        }
    }else if( this-> dim1 == other.dim1 && other.dim2 == 1){
        for(int i = 0; i < dim1; i++){
            for(int j = 0; j < dim2; j++){
                mu_out.data[j * dim1 + i] = this -> data[j * dim1 + i] = other.data[i];
            }
        }
    }else{
        cout << "inconsitence dimensions of the two matrix" << endl;
    }

    return mu_out;
}

Mu Mu::operator + (const Mu& other){
    if(this->dim1 != other.dim1 || this -> dim2 != other.dim2){
        cout << "inconsitence dimensions of the two matrix" << endl;
        Mu m(1,1,-1);
    }
    Mu mu_out(dim1, dim2,1);
    for (int i = 0; i < dim1 * dim2; i++) {
        mu_out.data[i] = this -> data[i] + other.data[i];
    }
    return mu_out;
}

Mu Mu::operator * (const Mu& other){
    if(this->dim1 != other.dim1 || this -> dim2 != other.dim2){
        cout << "inconsitence dimensions of the two matrix" << endl;
        Mu m(1,1,-1);
    }
    Mu mu_out(dim1, dim2,1);
    for (int i = 0; i < dim1 * dim2; i++) {
        mu_out.data[i] = this -> data[i] * other.data[i];
    }
    return mu_out;
}



double Mu::sum_abs(){
    double sum = 0;
    for(int i = 0; i < dim1 * dim2; i++){
        sum += fabs(this -> data[i]);
    }
    return sum;
}


//void Mu::transpose(){
//
//}

/****************************************************************************************/
// int InteractionSNPpairs::P;
// int InteractionSNPpairs::cj1;
// int InteractionSNPpairs::cj2;

InteractionSNPpairs::InteractionSNPpairs(int j1, int j2, double InteractionMeasure){
    this -> j1 = j1;
    this -> j2 = j2;
    this -> InteractionMeasure = InteractionMeasure;
}

InteractionSNPpairs::InteractionSNPpairs(int j1, int j2){
    this -> j1 = j1;
    this -> j2 = j2;
}

InteractionSNPpairs::InteractionSNPpairs(){
    this -> j1 = 0;
    this -> j2 = 1;
    this -> InteractionMeasure = 0;
}

InteractionSNPpairs::InteractionSNPpairs(const InteractionSNPpairs& mu){
    this -> j1 = mu.j1;
    this -> j2 = mu.j2;
    this -> InteractionMeasure = mu.InteractionMeasure;
}
//InteractionSNPpairs InteractionSNPpairs::next(){
//
//}

bool compare_pairs(const InteractionSNPpairs& a, const InteractionSNPpairs& b){
   return a.InteractionMeasure > b.InteractionMeasure;
}


void encode_matrix(const mat& m, int code_num, Mat<int>& data){
   int N = (int)m.n_rows;
   int P = (int)m.n_cols;
   Mat<int> indices(N,P);
   int* types = new int[P];
   int start = 0;
   int end = 0;
   const double* m_ptr = m.memptr();
   double value0, value1;
   for(int j = 0; j < P; j++){
        uvec ind = sort_index(m.col(j));
        for(int i = 0; i < N; i++){
          indices.at(i,j) = (int)ind[i];
        }
        start = ind.at(0);
        end = ind.at(ind.n_elem - 1);
        value0 = *(m_ptr + j * N + start);
        value1 = *(m_ptr + j * N + end);
        vec mj = m.col(j);
        // m.col(j) = mj;
        uvec  uniq = find_unique(m.col(j));
        if(uniq.size() > max_code_num){ // continuous
            //type1  to discretize
            types[j] = 0;
        }else{
          if(uniq.size() > code_num){
            // cout << "uniq.size() =" << uniq.size() << endl;
             code_num = uniq.size();
            // cout << "code_num =" << code_num << endl;

          }
          if(value1 <= code_num - 1 &&  value0 >= 0 ){
            //type2 no need to coding
            types[j] = 1;
          }else{
            //type3 to recoding
            types[j] = 2;
          }
        }


   }

   codes::code_num = code_num;
   Col<int> lines(code_num + 1);
   lines[0] = 0;
   for(int i = 1; i <= code_num; i++){
     lines[i] = i * N * 1.0 / code_num;
   }

   for(int j = 0; j < P; j++){
     if(types[j] == 1){   //type1  to discretize
        for(int i = 0; i < N; i++){
           data.at(i,j) = (int)m.at(i,j);
        }
     }else if(types[j] == 0){ //continuous, type2 no need to coding
       for(int i = 0; i < code_num; i++){
           start = lines[i];
           end = lines[i+1];
           for(int idx = start; idx < end; idx++){
               data.at(indices.at(idx, j),j) = i;
           }
       }
     }else if(types[j] == 2){  //type3 to recoding
        // uvec uniq_index = find_unique(m.col(j));
        uvec ind_sort = sort_index(m.col(j));
        vec  sort_col = sort(m.col(j));
        uvec uniq_index = find_unique(sort_col);
        for(int i = 0; i < uniq_index.size(); i++){
          int end = i < uniq_index.size() - 1 ? uniq_index[i+1] - 1 : (N - 1);
          for(int k = uniq_index[i]; k <= end; k++){
             data.at(ind_sort[k], j) = i;
          }
        }
     }
   }

   delete[] types;

}

void GenoDiscreteMatrix::subMatrix(const Col<uword>& indics, GenoDiscreteMatrix& subMatrix){
        uword subN = indics.size();
        for(uword i = 0; i < subN; i++){
            subMatrix.y[i] = this -> y[indics[i]];
            subMatrix.yc[i] = this -> yc[indics[i]];
            for(int j = 0; j < P; j++){
              subMatrix.data[j*subN + i] = this -> data[j * N + indics[i]];
            }
        }
}

GenoDiscreteMatrix::GenoDiscreteMatrix(double* yc, int* y_, double* fdata, int N, int P, int code_num){
  this -> N = N;
  this -> P = P;
  this -> discrete = 0;
  Mat<double> m(fdata, N, P, false);
  this -> data = new int[N*P];
  this -> y = new int[N];
  memcpy(y, y_, sizeof(int) * N);
  Mat<int> data_(data, N, P, false);
  uvec indics = find_unique(m);
  if(indics.size() > code_num){
     encode_matrix(m, code_num, data_);
  }else{
     this -> discrete = 1;
     for(int i = 0; i < N * P; i++){
        data[i] = (int)fdata[i];
     }
  }
  this -> yc = yc;
  this -> Xc = fdata;

}


vector<InteractionSNPpairs> getTopInteractionParis(GenoDiscreteMatrix* gdm, int extra_pairs, double chisq_thresh){
  vector<InteractionSNPpairs> pairs;
  if(extra_pairs == 0) return pairs;
  extra_pairs = (extra_pairs == -1) ? gdm -> N : extra_pairs;
  SingleCodes::init();
  GenoBitMatrix gbm(*gdm);
  gbm.thresholdRecord = chisq_thresh;
  gbm.CalGenoJointDistr();
  pairs = gbm.final_results;
  if(extra_pairs < pairs.size()){
     pairs.resize(extra_pairs);
  }
  return pairs;
}

/***************************************************************************************/
void cal_thread_bound(long long P, int thread_num, InteractionSNPpairs* pairs){
    long long* CM = new long long[P - 1];
    CM[0] = P - 1;
    for(long long i = P - 2; i >= 1; i--){
        CM[P - i - 1] = i + CM[P - i - 2];
    }

    long long total = P * (P-1) / 2;
    long long n = total / thread_num;
    // cout <<"total n " <<total <<" " << n << endl;

    n = total % thread_num == 0 ? n : (n + 1);
    // cout <<"n " << n << endl;

    long long current = 0;// current index 0,1,...,P-2
    long long start,end;
    for (int i = 0; i < thread_num; i++) {
        start = i * n + 1;
        end = start + n - 1 < total ? (start + n - 1) : total;
        while(true){
            if(start <= CM[current]){
          //      cout <<i <<"(" <<  current << ":" << P - 1 - (CM[current] - start) <<",";
            //    InteractionSNPpairs pair(current, P - 1 - (CM[current] - start));
                pairs[i*2].j1 = current;
                pairs[i*2].j2 = P - 1 - int(CM[current] - start);
                break;
            }
            current ++;
        }

        while(true){
            if(end <= CM[current]){
            //    cout << current << ":" << P - 1 - (CM[current] - end) <<")"<<endl;
                pairs[i*2 + 1].j1 = current;
                pairs[i*2 + 1].j2 = P - 1 - int(CM[current] - end);
                break;
            }
            current ++;
        }


    }

    delete[] CM;


}

void get_list(InteractionSNPpairs& pair0, InteractionSNPpairs& pair1, int P){
    int cj1 = pair0.j1;
    int cj2 = pair0.j2;
    while (true) {
        if(cj1 == pair1.j1 && cj2 == pair1.j2){
            break;
        }else if(cj2 == P - 1){
            cj1 ++;
            cj2 = cj1 + 1;
        }else{
            cj2 ++;
        }
    }

}

void inverse_2by2matrix(double* m, double* invm){
  double det = m[0]*m[3] - m[2]*m[1];
  invm[0] = m[3] / det;
  invm[1] = - m[1] / det;
  invm[2] = invm[1];//- m[2] / det;
  invm[3] = m[0] / det;
}


void inverse_3by3matrix(double* add3, double* a, double* invm){
  double s = 1 / (add3[2] - add3[0]*add3[0] * a[0] - 2*add3[0]*add3[1]*a[1] - add3[1]*add3[1]*a[3]); // for (D - B'A^{-1}B)^{-1})
  double* t = new double[2]; //for A^{-1}B
  t[0] = a[0] * add3[0] + a[1]*add3[1];
  t[1] = a[2] * add3[0] + a[3]*add3[1];

  invm[0] = a[0] + t[0] * t[0] * s;
  invm[1] = a[1] + t[0] * t[1] * s;
  invm[3] = invm[1];//a[2] + t[0] * t[1] * s;
  invm[4] = a[3] + t[1] * t[1] * s;

  invm[2] = - s * t[0];
  invm[5] = - s * t[1];

  invm[6] = invm[2];//- s * t[0];
  invm[7] = invm[5];//- s * t[1];
  invm[8] = s;

}


double cal_interact_fast(arma::mat& x, arma::vec& y, int n, int i1, int i2,
                         double* mat22, double* inv_mat22, double* add3, double* inv_mat33, double* beta3,
                         vec& x2, vec& xy, double y_mean, vec& y_){

  double LR=0;
  vec x12v = x.col(i1) % x.col(i2);
  double x12 = sum(x12v);

  mat22[0] = x2.at(i1);
  mat22[3] = x2.at(i2);
  mat22[1] = x12;
  mat22[2] = x12;

  inverse_2by2matrix(mat22, inv_mat22);

  double b1 = inv_mat22[0] * xy[i1] + inv_mat22[1] * xy[i2];
  double b2 = inv_mat22[2] * xy[i1] + inv_mat22[3] * xy[i2];
  double cov = y_mean;// - (x_mean[i1] * b1 + x_mean[i2] * b2);

  double* x_ptr = x.memptr();
  double s2_0 = 0;
  double tmp = 0;
  for(int i = 0; i < n; i++){
    tmp = y[i] - (cov + x_ptr[i1 * n + i] * b1 + x_ptr[i2 * n + i] * b2);
    s2_0 += tmp * tmp;
  }

  s2_0 /= n;

  add3[0] = sum(x.col(i1) % x12v);
  add3[1] = sum(x.col(i2) % x12v);
  add3[2] = sum(x12v % x12v);


  inverse_3by3matrix(add3, inv_mat22, inv_mat33);

  double addxy = sum(y_ % x12v);


  beta3[0] = inv_mat33[0] * xy[i1] + inv_mat33[3] * xy[i2] + inv_mat33[6] * addxy;
  beta3[1] = inv_mat33[1] * xy[i1] + inv_mat33[4] * xy[i2] + inv_mat33[7] * addxy;
  beta3[2] = inv_mat33[2] * xy[i1] + inv_mat33[5] * xy[i2] + inv_mat33[8] * addxy;

  cov = y_mean - (x12 * beta3[2] / n);

  //  arma::vec res1 = y - (cov + x.col(i1) * beta3[0] + x.col(i2) * beta3[1] + x12v * beta3[2]);
  double s2_1 = 0;

  tmp = 0;
  for(int i = 0; i < n; i++){
    tmp = y[i] - (cov + x_ptr[i1 * n + i] * beta3[0] + x_ptr[i2 * n + i] * beta3[1] + x12v[i] * beta3[2]);
    s2_1 += tmp * tmp;
  }


  s2_1 /= n;
  LR = n * log(s2_0/s2_1);

  return LR;

}

double cal_interact(arma::mat& x, arma::vec& y, int n, int i, int j){
    arma::mat x0(n,3,arma::fill::ones);
    x0.col(1) = x.col(i);
    x0.col(2) = x.col(j);

    arma::vec coef0 = arma::solve(x0,y);
    arma::vec res0 = y - x0 * coef0;
    double ip_res0 = arma::dot(res0,res0);
    double s2_0 = ip_res0 / n;              //here the estimate is mle, biased

    arma::mat x1 = join_rows(x0, x0.col(1) % x0.col(2));

    arma::vec coef1 = arma::solve(x1,y);
    arma::vec res1 = y - x1 * coef1;
    double ip_res1 = arma::dot(res1,res1);
    double s2_1 = ip_res1 / n;               //here the estimate is mle, biased

    //likelihood ratio
    double LR = n * log(s2_0/s2_1);
    return LR;

}

arma::mat lmpairwise_fast(arma::vec& y, arma::mat& x){

  int p = (int)x.n_cols;

  double y_mean = (double)mean(y);
  //  int p = (int)x.n_cols;
  vec x2(p);
  vec xy(p);
  vec y_ = y - y_mean;
  for(int i = 0; i < p; i++){
    x2.at(i) = sum(x.col(i) % x.col(i));
    xy.at(i) = sum(y_ % x.col(i));
  }

  arma::mat result(p*(p-1)/2,3);
  uword idx = 0;
  double LR = 0;
  int n = (int)x.n_rows;
  double* mat22 = new double[4];
  double* inv_mat22 = new double[4];
  double* add3 = new double[3];
  double* inv_mat33 = new double[9];
  double* beta3 = new double[3];
  for (int i = 0; i < p; i++ ) {
    for (int j = i+1; j < p; j++)
    {
      LR = cal_interact_fast(x, y, n, i, j,
                             mat22, inv_mat22, add3, inv_mat33, beta3,
                             x2, xy, y_mean, y_);
      result.at(idx,0) = LR;
      result.at(idx,1)= i;
      result.at(idx,2) = j;
      idx++;
    }
  }
  delete[] mat22;
  delete[] inv_mat22;
  delete[] inv_mat33;
  delete[] beta3;
  return result;
}
