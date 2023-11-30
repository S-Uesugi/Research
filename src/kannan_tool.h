#ifndef KANNAN_TOOL_H_INCLUDED
#define KANNAN_TOOL_H_INCLUDED

#include <iostream>
#include <vector>
#include <NTL/mat_ZZ.h>
#include <NTL/ZZ_p.h>
#include <string>

using namespace std;
using namespace NTL;

enum class Rot_Mode{
  Cyclic,
  Power2,
  Power3,
  Prime
};

class ToolBox{
  public:
    ZZ inv_to_int(int a,ZZ q);
    vec_ZZ reverse_vector(vec_ZZ &b,int dim);
    vec_ZZ rot(vec_ZZ &b,int dim,Rot_Mode mode);
    vec_ZZ inv_rot(vec_ZZ &b,int dim,Rot_Mode mode);
    vec_ZZ extension_vector(vec_ZZ &b,int dim);
    vec_ZZ extract_vector(vec_ZZ &b,int dim,ZZ q);
    mat_ZZ ex_KannanEmbedding(mat_ZZ &lB,vec_ZZ &b1,vec_ZZ &b2,int dim,Rot_Mode mode,int rotation_cnt,double eta);
    vec_ZZ rot_mrs(vec_ZZ &b,int dim,Rot_Mode mode,int rotation_cnt,ZZ q);
    mat_ZZ ex_KannanEmbedding2(mat_ZZ &lB,vec_ZZ &b1,vec_ZZ &b2,int dim,Rot_Mode mode, int rotation_cnt,ZZ q,double eta);
};

#endif