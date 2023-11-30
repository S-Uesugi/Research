#ifndef WEIGHTED_LLL_INCLUDED
#define WEIGHTED_LLL_INCLUDED

#include <iostream>
#include <vector>
#include <NTL/mat_ZZ.h>
#include <NTL/ZZ_p.h>
#include <string>

using namespace std;
using namespace NTL;

void Weight(mat_ZZ &A, vec_ZZ &weight, int dim);
void Weighted_LLL(mat_ZZ &A,mat_ZZ &B,vec_ZZ &weight, int dim);
void Weighted_BKZ(mat_ZZ &A,mat_ZZ &B,vec_ZZ &weight, int beta ,int dim);

#endif