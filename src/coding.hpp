//
//  coding.hpp
//  boost_new
//
//  Created by David on 2017/7/9.
//  Copyright © 2017年 David. All rights reserved.
//

#ifndef coding_hpp
#define coding_hpp

#include <stdio.h>
#include <armadillo>
#include <stdint.h>
#include <vector>
#include <thread>
#include <mutex>
#include <unordered_set>
//#include <boost/math/special_functions/erf.hpp>

//namespace bm = boost::math;

using namespace arma;
using namespace std;

const static uint64_t mask = 1;

const static int max_code_num = 8;

const static int LengthLongType = 64;

class codes{
public:
   static int code_num; // 3 for 0,1,2
   static int thread_num;
   static long max_num_of_pair;
};


//
// int mycodes::value = 1;

class InteractionSNPpairs;

class Mu;

double probit(double p);

class SingleCodes{
public:
    uint64_t* code;
    int n;
    int cnt; //number of 1
    static unsigned char wordbits[];
    static int bitCount(uint64_t i);
    static void init();
    SingleCodes();
    SingleCodes(int n);
    ~SingleCodes();
    SingleCodes operator | (const SingleCodes& b);
    SingleCodes operator & (const SingleCodes& b);
    static int pop_count(uint64_t i );
    int count();
};

class SnpCodes{
  public:
    SnpCodes(){

    }
    int type;//type = 0 for case, type = 1 for control
    static float* marginal_distr;
    static int P;
    static void init(int P, int code_num);
    int index;//index in the list of vairables
    int n;
    SingleCodes** codes;
    SnpCodes(int n, int type, int index);
    void coding(int* col, int N);
    ~SnpCodes();
    void jointDistrWith(SnpCodes* snpcodes, Mu&  GenoDistr, int code_num);
};



class GenoDiscreteMatrix{
public:
    int N;
    int P;
    int* data;
    int* y;
    double* yc;
    double* Xc;
    int discrete;//discrete = 1 for discrete matrix, discrete = 0 for continous matrix to be discretized
    GenoDiscreteMatrix(int N, int P, bool rdm = false);
    GenoDiscreteMatrix(double* yc, int* y_, double* fdata, int N, int P, int code_num);
    void subMatrix(const Col<uword>& indics, GenoDiscreteMatrix& submatrix);
    void sep2CCMatrix(GenoDiscreteMatrix * & matrix_case, GenoDiscreteMatrix * & matrix_ctrl);
    ~GenoDiscreteMatrix();
};


void add_vector(float* v1, float* v2, int num, float* dest);

class Mu{
public:
    int dim1;
    int dim2;
    float* data;
    Mu();
    Mu(int d1, int d2, float data0);
    Mu(float* data0, int d1, int d2);
    Mu(int* data0, int d1, int d2);
    void clean();
    Mu(const Mu& mu);
    void clone(const Mu& mu);
    void copy(float* data, int size);
    Mu& operator=(const Mu& other);
    Mu operator-(const Mu& other);
    Mu operator*(const Mu& other);
    Mu operator+(const Mu& other);
    Mu operator / (const Mu& other);
    Mu& operator/=(const Mu& other);

    ~Mu();
    double sum_abs();
    Mu marginal_dim1();//update the marginal info for dim1
    Mu marginal_dim2();//update the marginal info for dim2
    friend ostream& operator<<(ostream& os, const Mu& dt);
};

void divide_mu(const Mu& mu1, const Mu& other, Mu& res);
void plusdivide_mu(const Mu& mu1, const Mu& mu2, const Mu& other, Mu& res);
void update_mutmp(Mu& mutmp, Mu& multiplier, Mu& divisor);

class InteractionSNPpairs{
public:
    // static int P;
    // static int cj1;
    // static int cj2;
    int j1;
    int j2;
    double InteractionMeasure;
    InteractionSNPpairs();
    InteractionSNPpairs(const InteractionSNPpairs& mu);
    InteractionSNPpairs(int j1, int j2);
    InteractionSNPpairs(int j1, int j2, double InteractionMeasure);
};

bool compare_pairs(const InteractionSNPpairs& a, const InteractionSNPpairs& b);

class GenoBitMatrix{
public:
    int N;
    int P;
    double* Xc;
    double* yc;
    int discrete;
    int thread_num = 4;
    int code_num;
    float thresholdRecord = 30;
    // float pvalue_threashold = 0.5;
    vector<InteractionSNPpairs>* thread_results;
    vector<InteractionSNPpairs> final_results;
    InteractionSNPpairs* pairs;
    int* y;
    int ncase = 0;
    int nctrl = 0;
    SnpCodes** codes_ctrl;
    SnpCodes** codes_case;
    float* ones;
    float* zeros;
    std::mutex _mtx;
    Mu numbers;
    int cj1 = 0;
    int cj2 = 0;
    InteractionSNPpairs* next();
    float* marginal_count; // 3*P for case and control
    void CalculateMarginalEntropy(double * MarginalEntropySNP, double * MarginalEntropySNP_Y);
    void CalGenoJointDistr();
    double CalInteractFor(int j1, int j2, Mu& mucc, Mu& pab, Mu& pbc, Mu& pca);
    double CountNumFor(int j1, int j2, Mu& mucc, Mu& pab, Mu& pbc, Mu& pca);
    double CalInteractByThread(int thread_id);
    void cal_taointMeasure(const Mu& pab, const Mu& pbc, const Mu& pca,const Mu& mucc,double& tau, double& InteractionMeasure);
    void cal_taointMeasure(const Mu& mutmp, const Mu& mucc,double& tau, double& InteractionMeasure);
    void recalInteractDis(vector<InteractionSNPpairs>& pairs);
    void recalInteractCont(vector<InteractionSNPpairs>& pairs);
    //void cal_taointMeasure(const Mu& pab, const Mu& pbc, const Mu& pca,const Mu& mucc,double& tau, double& InteractionMeasure);
    GenoBitMatrix(GenoDiscreteMatrix& gdm);
    GenoBitMatrix(arma::Mat<int>& amat);
    ~GenoBitMatrix(){
      if(y != NULL)
          delete[] y;
    }
};


//vector<InteractionSNPpairs> getTopInteractionParis(GenoDiscreteMatrix& gdm, int N = -1, double pvalue_thresh = 0.5);


void encode_matrix(const mat& m, int code_num, Mat<int>& data);
void cal_thread_bound(long long P, int thread_num, InteractionSNPpairs* pairs);
void get_list(InteractionSNPpairs& pair0, InteractionSNPpairs& pair1, int P);
vector<InteractionSNPpairs> getTopInteractionParis(GenoDiscreteMatrix* gdm, int N, double chisq_thresh);
void calInteractCon();
void inverse_2by2matrix(double* m, double* invm);
void inverse_3by3matrix(double* add3, double* a, double* invm);
double cal_interact_fast(arma::mat& x, arma::vec& y, int n, int i1, int i2,
                         double* mat22, double* inv_mat22, double* add3, double* inv_mat33, double* beta3,
                         vec& x2, vec& xy, double y_mean, vec& y_);
//void BOLT_SSI(arma::mat X, arma::vec y, int extra_pairs = -1, double pvalue_thresh = 0,  int code_num = 3);
double cal_interact(arma::mat& x, arma::vec& y, int n, int i, int j);
arma::mat lmpairwise_fast(arma::vec& y, arma::mat& x);

//double sum_absdiff(Mu& mu1, Mu&);





//template<class T> GenotypeMatrix<T>::init();






#endif /* coding_hpp */
