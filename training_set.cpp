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


void ts_alloc(const int ts_size, const int param_dim, const char *model_name, TrainSet &ts)
{

    // set parameter space dimension and size //
    ts.param_dim = param_dim;
    ts.ts_size   = ts_size;

    // allocate memory for parameter matrix //
    double **params_tmp;
    params_tmp = ((double **) malloc(ts_size*sizeof(double *)));
    for(int j = 0; j < ts_size; j++)
    {

      params_tmp[j] = (double *)malloc(param_dim*sizeof(double));

      if(params_tmp[j]==NULL){
        fprintf(stderr,"Failed to allocate memory in BuildTS\n");
        exit(1);
      }

    }
    ts.params = params_tmp;

    ts.distributed = false; //default value is false. Sets to true if SplitTrainingSet called
    strcpy(ts.model,model_name);

    std::cout << "Using waveform model: " << ts.model << std::endl;

}

void uniform(const int &n, const double &a, const double &b, double *SomeArray)
{
    double factor = (b-a)/(double)(n-1);
    for(int i=0;i<n;i++)
    {
        SomeArray[i] = a + (double)i*factor;
    }
}

void BuildTS_TensorProduct(const int *params_num, const double *params_low, const double *params_high, TrainSet &ts)
{

    // KEEP 2D: its easier to understand and has been robustly tested. ND agrees exactly with 2d via unix diff 
    if(ts.param_dim == 2){
        BuildTS_TensorProduct2D(params_num,params_low,params_high,ts);
    }
    else{
        BuildTS_TensorProductND(params_num,params_low,params_high,ts);
    }

}

void BuildTS_TensorProduct2D(const int *params_num, const double *params_low, const double *params_high, TrainSet &ts)
{

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

            ts.params[counter][0] = param_0;
            ts.params[counter][1] = param_1;

            counter = counter + 1;
        }
    }

    free(param_list_0);
    free(param_list_1);
}


void BuildTS_FromFile(const char *ts_file, TrainSet &ts)
{

    std::cout << "Reading TS points from file: " << ts_file << std::endl;

    // TODO: use fread to make extend this limitation
    if(ts.param_dim != 2){
        fprintf(stderr,"TS from file does not yet support dimensions other than 2\n");
        exit(1);
    }

    int counter = 0;
    double p1, p2;

    FILE *data;
    data = fopen(ts_file,"r");
    while(fscanf(data, "%lf %lf", &p1, &p2) != EOF)
    {
        ts.params[counter][0] = p1;
        ts.params[counter][1] = p2;

        counter = counter + 1;
    }
    fclose(data);

    std::cout << "ts size = " << counter << std::endl;

    if( counter != ts.ts_size){
        std::cout << "TS file size does not match expected size" << std::endl;
        exit(1);
    }

}

void BuildTS_RecursiveSetBuild(const double *params_low, const double *params_step_size, const int *params_num, const int level, double *param_vector, int &counter, TrainSet &ts)
{

    if(level == ts.param_dim)
    {
        for(int j=0; j < ts.param_dim; j++){
            ts.params[counter][j] = param_vector[j];
        }
        counter = counter + 1;
    }
    else
    {
        for(int i=0; i < params_num[level]; i++){
            param_vector[level] = params_low[level]+((double)i)*params_step_size[level];
            BuildTS_RecursiveSetBuild(params_low,params_step_size,params_num,level+1,param_vector,counter,ts);
        }
    }
}

void BuildTS_TensorProductND(const int *params_num, const double *params_low, const double *params_high, TrainSet &ts)
{
    double *params_step_size, *param_vector;

    param_vector     = (double *)malloc(ts.param_dim*sizeof(double));
    params_step_size = (double *)malloc(ts.param_dim*sizeof(double));

    int counter = 0; // used to fill ts as ts.params[counter][j] where j  = 1 to ts.param_dim-1

    // compute step sizes for equal spaced parameter samples //
    for (int i=0; i < ts.param_dim; i++){
        params_step_size[i] = ( params_high[i] - params_low[i] )/((double)(params_num[i]-1));
    }

    BuildTS_RecursiveSetBuild(params_low,params_step_size,params_num,0,param_vector,counter,ts);

    free(param_vector);
    free(params_step_size);
}



void SplitTrainingSet(const int size, TrainSet &ts)
{
    int *mystart_tmp, *myend_tmp, *matrix_sub_size_tmp;

    int numprocs_worker = size - 1; // master (proc 0) will do no work
    int rows = ts.ts_size;
    int proc_id;

    mystart_tmp         = (int *)malloc(numprocs_worker*sizeof(int));
    myend_tmp           = (int *)malloc(numprocs_worker*sizeof(int));
    matrix_sub_size_tmp = (int *)malloc(numprocs_worker*sizeof(int));

    for(int i = 2; i <= size; i++)
    {  
        proc_id = i-2;
        mystart_tmp[i-2] = (rows / numprocs_worker) * proc_id;
        if (rows % numprocs_worker > proc_id)
        {  
            mystart_tmp[i-2] += proc_id;
            myend_tmp[i-2] = mystart_tmp[i-2] + (rows / numprocs_worker) + 1;
        }
        else
        {  
            mystart_tmp[i-2] += rows % numprocs_worker;
            myend_tmp[i-2] = mystart_tmp[i-2] + (rows / numprocs_worker);
        }
    }

    ts.distributed = true;

    for(int i = 0; i < size - 1; i++){
        matrix_sub_size_tmp[i] = myend_tmp[i] - mystart_tmp[i];
    }

    // -- NOTE: Free mystart and myend in main -- //
    ts.mystart         = mystart_tmp;
    ts.myend           = myend_tmp;
    ts.matrix_sub_size = matrix_sub_size_tmp;

}

void WriteTrainingSet(const TrainSet ts)
{

    FILE *data1;

    char filename[] = "TS_Points.txt";
    data1 = fopen(filename,"w");
    for(int i = 0; i < ts.ts_size ; i++){   
        fprintf(data1,"%.15le %.15le\n",ts.params[i][0],ts.params[i][1]);
    }
    fclose(data1);
}



