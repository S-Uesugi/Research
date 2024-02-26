#include <iostream>
#include <vector>
#include <algorithm>
#include <NTL/ZZX.h>
#include <NTL/matrix.h>
#include <NTL/LLL.h>
#include <fstream>
#include <sstream>
#include <string>
#include <time.h>
#include <filesystem>
#include "kannan_tool.h"
#include "weighted_LLL.h"

using namespace std;
using namespace NTL;

string MakeResultDirectory(string root_result_path,vector<string> data_list,int exp_no,int M ,int rotation){
  bool flag = true;
  int path_length = data_list[exp_no].length();
  int split_num;
  string result_path;
  if(exp_no == 2 || exp_no == 3 || exp_no == 4){
    split_num = 13;
  }
  else{
    split_num = 11;
  }
  while(flag){
    string main_path = root_result_path+data_list[exp_no].substr(split_num,path_length-split_num);
    cout << main_path << endl;
    if(filesystem::exists(main_path)){
      string sub_path = main_path+"/M"+to_string(M);
      if(filesystem::exists(sub_path)){
        string sub_sub_path = sub_path + "/rot"+to_string(rotation)+'/';
        if(filesystem::exists(sub_sub_path)){
          result_path = sub_sub_path;
          flag = false;
        }
        else{
          filesystem::create_directory(sub_sub_path);
        }
      }
      else{
        filesystem::create_directory(sub_path);
      }
    }
    else{
      filesystem::create_directory(main_path);
    }
  }
  return result_path;
}

void ImportFile(vector<int> & v,string filename){
  ifstream inputfile;
  inputfile.open(filename,ios::in);
  if(!inputfile){
    cout << "Can't open this file !!"<<endl;
    exit(1);
  }
  string val_list,val;
  while (getline(inputfile,val_list,'\n'))
  {
    stringstream vals(val_list);
    while(getline(vals,val,',')){
     v.push_back(stoi(val)); 
    }

  }
}

vector<mat_ZZ> MakeMatrixList(vector<int> & v,int dim){
  vector<mat_ZZ> matrix_list;
  mat_ZZ matrix;
  matrix.SetDims(dim,dim);
  int sample_num = v.size() / (dim*dim);
  for(int i = 0;i < sample_num;i++){
    for(int j = 0;j < dim;j++){
      for(int k= 0;k < dim;k++){
        matrix[j][k] = v[i*dim*dim + j*dim+k];
      }
    }
    matrix_list.push_back(matrix);
  }
  return matrix_list; 
}

vector<vec_ZZ> MakeVectorList(vector<int> & v,int dim){
  vector<vec_ZZ> vector_list;
  vec_ZZ V;
  V.SetLength(dim);
  int sample_num = v.size() / dim;
  for(int i = 0;i < sample_num;i++){
    for(int j = 0;j < dim;j++){
        V[j] = v[i*dim+j];
    }
    vector_list.push_back(V);
  }
  return vector_list; 
}

void CombineMatrix(mat_ZZ &A12, mat_ZZ &A1,mat_ZZ &A2,int dim){
    for(int i = 0;i < dim;i++){
        for(int j = 0;j < dim;j++){
            A12[i][j] = A1[i][j];
            A12[i][j+dim] = A2[i][j];
        }
    }
}

void CombineVector(vec_ZZ &comb,vec_ZZ &b1,vec_ZZ &b2,int dim){
  for(int i = 0;i<dim;i++){
    comb[i] = b1[i];
    comb[dim+i] = b2[i];
  }
}

void MakeQaray(mat_ZZ &q_aray, mat_ZZ &A12,int dim,ZZ q){
  for(int i = 0;i < dim;i++){
      for(int j = 0;j < 2*dim;j++) q_aray[i][j] = A12[i][j];
      }
      for(int i = 0;i < 2*dim;i++) q_aray[dim+i][i] = q;
}

void GetLLL(mat_ZZ &Basis,mat_ZZ &lB,int dim){
  if(IsZero(lB[0])){
    for(int i = 0; i < 2*dim;i++){
      for(int j = 0;j < 2*dim;j++)Basis[i][j] = lB[dim+i][j];
    }
  }
  else{
    for(int i = 0; i < 2*dim;i++){
      for(int j = 0;j < 2*dim;j++)Basis[i][j] = lB[i][j];
    }
  }
}

void Reversal(vec_ZZ &target){
      for(int i = 0;i < target.length();i++) target[i] = -1*target[i];
}

void AbsError(vec_ZZ &error, ZZ q){
  for(int i = 0;i < error.length();i++){
    if(error[i] >= q/2) error[i] = error[i] - q;
  }
}

bool CheckCode(vec_ZZ &target ,vec_ZZ & b){
  int chk = 0;
  for(int i = 0;i < target.length();i++){
    if(target[i] == b[i]) chk++;
  }
  if(chk == target.length()) return true;
  else return false;
  
}

vec_ZZ extract_target(vec_ZZ &b,int rot_cnt){
  vec_ZZ v;
  v.SetLength(b.length()-rot_cnt);
  for(int i = 0;i < b.length()-rot_cnt;i++) v[i] = b[i];
  return v;
}

bool minus_check(vec_ZZ &b,int rot_cnt,int dim){
  bool flag = false;
  for(int i = 0;i < rot_cnt;i++){
    if(b[2*dim+i] == -1) flag = true;
  }
  return flag;
}

ZZ Norm(vec_ZZ &V){
    ZZ res;
    InnerProduct(res, V, V);
    return SqrRoot(res);
}

/* argv[1] : kind of experiments, argv[2] : embedding factor, arg[3] : rotation, arg[4] : blocksize */
int main(int argc, char * argv[]){  
  if(argc != 5) exit(0);
  string path = "../data/";

  /* 0: exp1_131, 1 : exp1_256, 2 : original64, 3 : original67, 4 : original81 */
  vector<string> data_list;
  for(auto& data : filesystem::directory_iterator(path)) data_list.push_back(string(data.path()));
  sort(data_list.begin(),data_list.end());

  int exp_no = atoi(argv[1]);
  int M  =  atoi(argv[2]);
  int rotation = atoi(argv[3]);
  int beta  =  atoi(argv[4]);
  Rot_Mode mode;
  switch (exp_no)
  {
  case 0:
    mode = Rot_Mode::Prime;
    break;
  case 1:
    mode = Rot_Mode::Power2;
    break;
  case 2:
    mode = Rot_Mode::Power2;
    break;
  case 3:
    mode = Rot_Mode::Prime;
    break;
  case 4:
    mode = Rot_Mode::Power3;
    break;
  default:
    break;
  }
  
  string filename(data_list[exp_no] + "/parameters.csv");
  string filename1(data_list[exp_no] + "/Sample_A1.csv");
  string filename2(data_list[exp_no] + "/Sample_b1.csv");
  string filename3(data_list[exp_no] + "/Sample_A2.csv");
  string filename4(data_list[exp_no] + "/Sample_b2.csv");
  string filename5(data_list[exp_no] + "/error1.csv");
  string filename6(data_list[exp_no] + "/error2.csv");
  string root_result_path = "../result/search/rlwe";
  string result_path = MakeResultDirectory(root_result_path,data_list,exp_no,M,rotation);
  string rhfname(result_path+"Root_Hermite_Factor_bkz"+to_string(beta)+".csv");
  ofstream writing_result;

  vector<int> parameters,matrix_list_A1,matrix_list_A2,vector_list_b1,vector_list_b2,vector_list_e1,vector_list_e2;

  ImportFile(parameters,filename);
  ImportFile(matrix_list_A1,filename1);
  ImportFile(vector_list_b1,filename2);
  ImportFile(matrix_list_A2,filename3);
  ImportFile(vector_list_b2,filename4);
  ImportFile(vector_list_e1,filename5);
  ImportFile(vector_list_e2,filename6);


  ZZ q;
  q = parameters[0];
  int dim = parameters[1];


  vector<mat_ZZ> exp_A1 = MakeMatrixList(matrix_list_A1,dim);
  vector<vec_ZZ> exp_b1 = MakeVectorList(vector_list_b1,dim);
  vector<vec_ZZ> exp_e1 = MakeVectorList(vector_list_e1,dim);
  vector<mat_ZZ> exp_A2 = MakeMatrixList(matrix_list_A2,dim);
  vector<vec_ZZ> exp_b2 = MakeVectorList(vector_list_b2,dim);
  vector<vec_ZZ> exp_e2 = MakeVectorList(vector_list_e2,dim);

  
  ToolBox tools;
  int exp_num = 100;
  double success = 0;
  double rhf_sum = 0;
  clock_t start = clock();
  writing_result.open(rhfname,ios::out);
  for(int i = 0;i < exp_num;i++){
    vec_ZZ b12,error;
    mat_ZZ A12,q_aray,lB,B;
    b12.SetLength(2*dim);
    error.SetLength(2*dim);
    A12.SetDims(dim,2*dim);
    q_aray.SetDims(3*dim,2*dim);
    lB.SetDims(2*dim,2*dim);

    CombineMatrix(A12,exp_A1[i],exp_A2[i],dim);
    CombineVector(b12,exp_b1[i],exp_b2[i],dim);
    CombineVector(error,exp_e1[i],exp_e2[i],dim);

    
    MakeQaray(q_aray,A12,dim,q);
    long lll_qp = LLL_QP(q_aray);
    GetLLL(lB,q_aray,dim);
    B = tools.ex_KannanEmbedding(lB,exp_b1[i],exp_b2[i],dim,mode,rotation,M);
    vec_ZZ w;
    w.SetLength(2*dim+rotation);
    for(int i = 0 ; i < 2*dim+rotation;i++){
      if(i < 2*dim+1){
        w[i] = 1;
      } 
      else{
        div(w[i],Norm(B[i]) , Norm(B[2*dim]));
        if(w[i] == 0) w[i] = 1;
      }  
    }
    mat_ZZ WB;
    WB.SetDims(2*dim+rotation,2*dim+rotation);
    Weighted_BKZ(B,WB,w,beta,2*dim+rotation);
    vector<vec_ZZ> cand;
    double norm_error,norm_w,rhf;
    ZZ inner_error,inner_w;
    bool flag = true;
    int decided = 0;

    for(int j = 0;j < rotation+1;j++){
        if(minus_check(WB[j],rotation,dim)) Reversal(WB[j]);
        vec_ZZ b = extract_target(WB[j],rotation);
        cand.push_back(b);
    }

    AbsError(error,q);
    InnerProduct(inner_error,error,error);
    for(int j = 0;j < rotation;j++){
      if(CheckCode(error,cand[j])){
        success++;
        cout << "Success" << endl;
        InnerProduct(inner_w,cand[j],cand[j]);
        flag = false;
        decided = j;
      }
    }
    if(flag){
      InnerProduct(inner_w,cand[0],cand[0]);
    }
    conv(norm_error,inner_error);
    conv(norm_w,inner_w);
    rhf = pow(sqrt(norm_w)/sqrt(norm_error),double((double) 1/double(2*dim+1)));
    rhf_sum += rhf;

    cout << "error norm : "<< sqrt(norm_error) << endl;
    cout <<"Round "+ to_string(i+1) +"| RHF : " << rhf << " B[0] norm : " << sqrt(norm_w)<<endl;
    writing_result <<"Round "+ to_string(i+1)+"| RHF : " << rhf << "B[0] norm : " << sqrt(norm_w)<<endl;
    writing_result << "B[0] : " << cand[decided] << endl;

    
  }

  writing_result.close();
  clock_t end = clock();
  const double time = static_cast<double>(end-start) / CLOCKS_PER_SEC;

  cout <<"__________Result__________" << " n = " << 2*dim << endl;
  cout <<"Success :"<<success /100<<endl;
  cout <<"Time[s] :"<<time<<endl;
  cout <<"Average Time[s] :"<<time/100<<endl;
  ofstream writing_file;
  string result_filename = result_path + "result_bkz"+to_string(beta)+".csv";
  writing_file.open(result_filename,ios::out);
  writing_file << "__________Result__________"<<endl;
  writing_file << " n = " <<dim<<endl;
  writing_file << "Success :"<<success /100<<endl;
  writing_file << "Time[s] :"<<time <<endl;
  writing_file << "Average Time[s] :"<<time/100<<endl;
  writing_file << "Average RHF :"<<rhf_sum/100<<endl;
  writing_file.close();
  return 0;
}