#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <time.h>
#include <algorithm>
#include <NTL/LLL.h>
#include "weighted_LLL.h"

using namespace std;
using namespace NTL;

void Weight(mat_ZZ &A, vec_ZZ &weight, int dim){
    mat_ZZ W;
    W.SetDims(dim, dim);
    for(int i = 0; i < dim; i++){
        W[i][i] = weight[i];
    }
    mul(A, A, W);
}

void Weighted_LLL(mat_ZZ &A,mat_ZZ &B,vec_ZZ &weight, int dim){
  Weight(A,weight,dim);
  for(int i = 0;i < dim;i++){
    for(int j = 0;j < dim;j++) B[i][j] = A[i][j];
  }
  G_LLL_XD(B);  
}

void Weighted_BKZ(mat_ZZ &A,mat_ZZ &B,vec_ZZ &weight, int beta ,int dim){
  Weight(A,weight,dim);
  for(int i = 0;i < dim;i++){
    for(int j = 0;j < dim;j++) B[i][j] = A[i][j];
  }
  BKZ_QP1(B,0.99,beta,10,0,1);
}

