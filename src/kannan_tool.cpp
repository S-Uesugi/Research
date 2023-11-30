#include "kannan_tool.h"
#include <math.h>

using namespace std;
using namespace NTL;


ZZ ToolBox::inv_to_int(int a, ZZ q){
  ZZ inved_a;
  ZZ_p::init(q);
  ZZ_p inv_a;

  inv(inv_a,ZZ_p(a));
  inved_a = rep(inv_a);
  return inved_a;
}

vec_ZZ ToolBox::reverse_vector(vec_ZZ &b,int dim){
  vec_ZZ w;
  w = VectorCopy(b,dim);

  for(int i = 1;i < dim;i++) w[dim-i] = -b[i];
  return w;
}

vec_ZZ ToolBox::rot(vec_ZZ &b,int dim,Rot_Mode mode){
  vec_ZZ w;
  w = VectorCopy(b,dim);

  switch (mode)
  {
  case Rot_Mode::Cyclic:
    for(int i = 1;i < dim;i++) w[i] = b[i-1];
    w[0] = b[dim-1];
    break;
  
  case Rot_Mode::Power2:
    for(int i = 1;i < dim;i++) w[i] = b[i-1];
    w[0] = -b[dim-1];
    break;
  
  case Rot_Mode::Power3:
    for(int i = 1;i < dim;i++) w[i] = b[i-1];
    w[0] = -b[dim-1];
    w[dim/2] = b[dim/2-1] - b[dim-1];
    break;
  
  case Rot_Mode::Prime:
    for(int i = 1;i < dim;i++) w[i] = b[i-1]-b[dim-1];
    w[0] = -b[dim-1];
    break;

  default:
    break;
  }

  return w;
}

vec_ZZ ToolBox::inv_rot(vec_ZZ &b,int dim,Rot_Mode mode){
  vec_ZZ w;
  w = VectorCopy(b,dim);

  switch (mode)
  {
  case Rot_Mode::Cyclic:
    w[dim-1] = b[0];
    for(int i = 0;i < dim-1;i++) w[i] = b[i+1];
    break;
  
  case Rot_Mode::Power2:
    w[dim-1] = -b[0];
    for(int i = 0;i < dim-1;i++) w[i] = b[i+1];
    break;
  
  case Rot_Mode::Power3:
    for(int i = 1;i < dim;i++) w[i] = b[i+1];
    w[dim-1] = -b[0];
    w[dim/2-1] = b[dim/2-1] + b[0];
    break;
  
  case Rot_Mode::Prime:
    w[dim-1] = -b[0];
    for(int i = 0;i < dim-1;i++) w[i] = b[i+1]+b[dim-1];
    break;

  default:
    break;
  }

  return w;
}

vec_ZZ ToolBox::extension_vector(vec_ZZ &b,int dim){
  vec_ZZ w,extended_vector;
  w = VectorCopy(b,dim);
  extended_vector.SetLength(2*dim);
  extended_vector[0] = w[0];
  for(int i =1 ;i < dim;i++){
    extended_vector[i] = w[i];
    extended_vector[2*dim-i] = -w[i];
  }
  return extended_vector;
}

vec_ZZ ToolBox::extract_vector(vec_ZZ &b,int dim,ZZ q){
  vec_ZZ extracted_vector;
  extracted_vector.SetLength(dim);
  for(int i = 0;i < dim;i++) extracted_vector[i] = b[i]%q;
  return extracted_vector;
}

mat_ZZ ToolBox::ex_KannanEmbedding(mat_ZZ &lB,vec_ZZ &b1,vec_ZZ &b2,int dim,Rot_Mode mode, int rotation_cnt,double eta){
  mat_ZZ B;
  B.SetDims(2*dim+rotation_cnt,2*dim+rotation_cnt);
  for(int i = 0;i < 2*dim;i++){
      for(int j = 0;j < 2*dim;j++) B[i][j] = lB[i][j];
    }

  vec_ZZ rot_b1,rot_b2;
  rot_b1 = VectorCopy(b1,dim);
  rot_b2 = VectorCopy(b2,dim);
  for(int i = 0;i < rotation_cnt;i++){
    B[2*dim+i][2*dim+i] = eta;
    for(int j = 0;j <dim;j++){
        B[2*dim+i][j] = rot_b1[j];
        B[2*dim+i][dim+j] = rot_b2[j];
      }
    rot_b1 = rot(rot_b1,dim,mode);
    rot_b2 = rot(rot_b2,dim,mode);
  }
    return B;
}

vec_ZZ ToolBox::rot_mrs(vec_ZZ &b,int dim,Rot_Mode mode,int rotation_cnt,ZZ q){
  vec_ZZ b_hat = ToolBox::extension_vector(b,dim);
  vec_ZZ bp,bm,b12,b_bar;
  bp = VectorCopy(b_hat,2*dim);
  bm = VectorCopy(b_hat,2*dim);
  for(int i = 0;i < rotation_cnt;i++){
    bp = ToolBox::rot(bp,2*dim,mode);
    bm = ToolBox::inv_rot(bp,2*dim,mode);
  }
  if(rotation_cnt >= 1){
    b12 = bp+bm;
  }
  else{
    b12 = b_hat;
  }
  b_bar = ToolBox::extract_vector(b12,dim,q);
  return b_bar;
}

mat_ZZ ToolBox::ex_KannanEmbedding2(mat_ZZ &lB,vec_ZZ &b1,vec_ZZ &b2,int dim,Rot_Mode mode, int rotation_cnt,ZZ q,double eta){
  mat_ZZ B;
  B.SetDims(2*dim+rotation_cnt+1,2*dim+rotation_cnt+1);
  for(int i = 0;i < 2*dim;i++){
      for(int j = 0;j < 2*dim;j++) B[i][j] = lB[i][j];
    }
  
  for(int i = 0;i <dim;i++){
    B[2*dim][i] = b1[i];
    B[2*dim][dim+i] = b2[i];
  }
  B[2*dim][2*dim] = eta;

  for(int i = 1;i <= rotation_cnt;i++){
      vec_ZZ rot_b1,rot_b2;
      rot_b1 = VectorCopy(b1,dim);
      rot_b2 = VectorCopy(b2,dim);
      
      rot_b1 = ToolBox::rot_mrs(rot_b1,dim,mode,i,q);
      rot_b2 = ToolBox::rot_mrs(rot_b2,dim,mode,i,q);
      for(int j = 0;j <dim;j++){
          B[2*dim+i][j] = rot_b1[j];
          B[2*dim+i][dim+j] = rot_b2[j];
        }
      
      B[2*dim+i][2*dim+i] = eta;
  }
    return B;
}
