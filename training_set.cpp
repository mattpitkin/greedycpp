#include <iostream>
#include <fstream>
#include <string.h>
#include <cmath>
#include <time.h>
#include <vector>
#include <sstream>

#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_vector_complex.h>
#include <gsl/gsl_block_complex_float.h>

#include "training_set.hpp"


// using implicit destructor... check valgrid for memory leaks

TrainingSetClass::TrainingSetClass(int dim, double * scale, int size, const char * model_name){

    strcpy(model_,model_name);
    param_scale_ = scale;
    param_dim_   = dim;
    ts_size_     = size;
    distributed_ = false; //default value is false. Sets to true if SplitTrainingSet called

    // allocate memory for parameter matrix //
    params_ = ((double **) malloc(size*sizeof(double *)));
    for(int j = 0; j < size; j++)
    {

      params_[j] = (double *)malloc(dim*sizeof(double));

      if(params_[j]==NULL){
        fprintf(stderr,"Failed to allocate memory in BuildTS\n");
        exit(1);
      }

    }

    std::cout << "Using waveform model: " << model_name << ". Training set class initialized!" << std::endl;
}

void TrainingSetClass::BuildTS(int *params_num, double *params_low, double *params_high){

    std::cout << "Building from tensor product grid" << std::endl;

    if(param_dim_ == 2){
        BuildTS_TensorProduct2D(params_num,params_low,params_high);
    }
    else{
        BuildTS_TensorProductND(params_num,params_low,params_high);
    }
}
 
void TrainingSetClass::BuildTS(const char *ts_file){

    std::cout << "Building from file input" << std::endl;

    if(param_dim_ == 2){
        BuildTS_FromFile2D(ts_file);
    }
    else{
        BuildTS_FromFileND(ts_file);
    }
}

void TrainingSetClass::BuildTS_FromFile2D(const char *ts_file)
{

    std::cout << "Reading TS points from file: " << ts_file << std::endl;

    if(param_dim_ != 2){
        fprintf(stderr,"TS from file does not yet support dimensions other than 2\n");
        exit(1);
    }

    int counter = 0;
    double p1, p2;

    FILE *data;
    data = fopen(ts_file,"r");

    if(data == NULL) {
        std::cout << "failed to open training set file" << std::endl;
        exit(1);
    }

    while(fscanf(data, "%lf %lf", &p1, &p2) != EOF)
    {
        params_[counter][0] = p1;
        params_[counter][1] = p2;

        counter = counter + 1;
    }
    fclose(data);

    std::cout << "ts size = " << counter << std::endl;

    if( counter != ts_size_){
        std::cout << "TS file size does not match expected size" << std::endl;
        exit(1);
    }

}

void TrainingSetClass::BuildTS_FromFileND(const char *ts_file)
{
// could try using fread
// cin a good solution for this? whats benefits of fscanf or other functions?
// if training set values are in strange format, will cin still work?

    std::cout << "Reading TS points from file: " << ts_file << std::endl;

    double parameter_tmp;
    int counter = 0;

    std::ifstream data(ts_file);

    if(!data.is_open()) {
        std::cout << "failed to open training set file" << std::endl;
        exit(1);
    }

    for(int i = 0; i < ts_size_;++i) { // loop over rows
        for(int j  = 0; j < param_dim_;++j) { // loop over columns

            if(data >> parameter_tmp){ 
                params_[i][j] = parameter_tmp;
                //fprintf(stdout,"%1.14e ",p1);
            }
            else{
                std::cout << "failed reading training set from file" << std::endl;
                exit(1);
            }

        }
        counter = counter + 1;
    }

    data.close();

    if( counter != ts_size_){
        std::cout << "TS file size does not match expected size" << std::endl;
        exit(1);
    }

    std::cout << "ts size = " << counter << std::endl;
}

void TrainingSetClass::SplitTrainingSet(const int size)
{

    int numprocs_worker = size - 1; // master (proc 0) will do no work
    int rows = ts_size_;
    int proc_id;

    mystart_         = (int *)malloc(numprocs_worker*sizeof(int));
    myend_           = (int *)malloc(numprocs_worker*sizeof(int));
    matrix_sub_size_ = (int *)malloc(numprocs_worker*sizeof(int));

    for(int i = 2; i <= size; i++)
    {  
        proc_id = i-2;
        mystart_[i-2] = (rows / numprocs_worker) * proc_id;
        if (rows % numprocs_worker > proc_id)
        {  
            mystart_[i-2] += proc_id;
            myend_[i-2] = mystart_[i-2] + (rows / numprocs_worker) + 1;
        }
        else
        {  
            mystart_[i-2] += rows % numprocs_worker;
            myend_[i-2]    = mystart_[i-2] + (rows / numprocs_worker);
        }
    }

    distributed_ = true;

    for(int i = 0; i < size - 1; i++){
        matrix_sub_size_[i] = myend_[i] - mystart_[i];
    }

}

int TrainingSetClass::FindRowIndxRank(const int global_row_indx)
{
    int row_rank = 0;
    bool found_rank = false;

    if(distributed_){

        while(found_rank == false)
        {  
            if( global_row_indx >= mystart_[row_rank] && global_row_indx <= myend_[row_rank] ){
                found_rank = true;
            }
            else{
                row_rank = row_rank + 1;
            }
        }
    }

    return row_rank;
}



void TrainingSetClass::uniform(const int &n, const double &a, const double &b, double *SomeArray)
{
    double factor = (b-a)/(double)(n-1);
    for(int i=0;i<n;i++)
    {
        SomeArray[i] = a + (double)i*factor;
    }
}


void TrainingSetClass::BuildTS_TensorProduct2D(const int *params_num, const double *params_low, const double *params_high)
{

    if(param_dim_ != 2){
        fprintf(stderr,"TS from file does not yet support dimensions other than 2\n");
        exit(1);
    }

    double *param_list_0, *param_list_1;
    double param_0, param_1;
    int counter = 0;

    param_list_0 = (double *)malloc(params_num[0]*sizeof(double));
    param_list_1 = (double *)malloc(params_num[1]*sizeof(double));

    if(param_list_0==NULL || param_list_1==NULL){
        fprintf(stderr,"Failed to allocate memory in BuildTS\n");
        free(param_list_0);
        free(param_list_1);
        exit(1);
    }

    // fills param_list with params_num[i] samples from params_low[i] to params_high[i] //
    uniform(params_num[0], params_low[0], params_high[0], param_list_0);
    uniform(params_num[1], params_low[1], params_high[1], param_list_1);

    // ex: 2D parameter space, param_list x param_list tensor product:
    // (params[i][0],params[i][1]) is the ith training set element
    for(int i = 0; i < params_num[0]; i++)
    {
        param_0 = param_list_0[i];

        for(int j = 0; j < params_num[1]; j++)
        {
            param_1 = param_list_1[j];

            params_[counter][0] = param_0;
            params_[counter][1] = param_1;

            counter = counter + 1;
        }
    }

    free(param_list_0);
    free(param_list_1);
}


void TrainingSetClass::BuildTS_RecursiveSetBuild(const double *params_low, const double *params_step_size, const int *params_num, const int level, double *param_vector, int &counter)
{

    if(level == param_dim_)
    {
        for(int j=0; j < param_dim_; j++){
            params_[counter][j] = param_vector[j];
        }
        counter = counter + 1;
    }
    else
    {
        for(int i=0; i < params_num[level]; i++){
            param_vector[level] = params_low[level]+((double)i)*params_step_size[level];
            BuildTS_RecursiveSetBuild(params_low,params_step_size,params_num,level+1,param_vector,counter);
        }
    }
}

void TrainingSetClass::BuildTS_TensorProductND(const int *params_num, const double *params_low, const double *params_high)
{
    double *params_step_size, *param_vector;

    param_vector     = (double *)malloc(param_dim_*sizeof(double));
    params_step_size = (double *)malloc(param_dim_*sizeof(double));

    int counter = 0; // used to fill ts as ts.params[counter][j] where j  = 1 to ts.param_dim-1

    // compute step sizes for equal spaced parameter samples //
    for (int i=0; i < param_dim_; i++){
        params_step_size[i] = ( params_high[i] - params_low[i] )/((double)(params_num[i]-1));
    }

    BuildTS_RecursiveSetBuild(params_low,params_step_size,params_num,0,param_vector,counter);

    free(param_vector);
    free(params_step_size);
}


void TrainingSetClass::WriteTrainingSet()
{

    FILE *data1;

    // TODO: should be N-Dimensional
    char filename[] = "TS_Points_New.txt";
    data1 = fopen(filename,"w");
    for(int i = 0; i < ts_size_ ; i++){   
        fprintf(data1,"%.15le %.15le\n",params_[i][0],params_[i][1]);
        //fprintf(data1,"%.15le %.15le %.15le %.15le\n",params_[i][0],params_[i][1],params_[i][2],params_[i][3]);
    }
    fclose(data1);
}


int TrainingSetClass::ts_size(){ return ts_size_; }

int TrainingSetClass::param_dim(){ return param_dim_; }

bool TrainingSetClass::distributed(){ return distributed_; }

const char * TrainingSetClass::model(){ return model_; }

const int * TrainingSetClass::mystart(){return mystart_; }

const int * TrainingSetClass::myend(){ return myend_; }

const int * TrainingSetClass::matrix_sub_size(){ return matrix_sub_size_; }

double ** TrainingSetClass::params(){ return params_; }

const double * TrainingSetClass::param_scale(){ return param_scale_; }

