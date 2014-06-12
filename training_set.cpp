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

    // parameter space dimension //
    ts.param_dim = param_dim;

    ts.ts_size = ts_size;

    // allocate memory //
    double *m1_tmp, *m2_tmp;
    m1_tmp = (double *)malloc(ts.ts_size*sizeof(double));
    m2_tmp = (double *)malloc(ts.ts_size*sizeof(double));
    if(m1_tmp==NULL || m2_tmp==NULL)
    {
        std::cout << "Failed to allocate memory in BuildTS" << std::endl;
        free(m1_tmp); free(m2_tmp);
        exit(1);
    }

    ts.m2 = m2_tmp;
    ts.m1 = m1_tmp;

    // building parameter matrix //
    /*double **params_tmp;
    params_tmp = ((double **) malloc(ts_size*sizeof(double *)));
    for(int j = 0; j < ts_size; j++)
    {

      params_tmp[j]=(double *)malloc(param_dim*sizeof(double));

      if(params_tmp[j]==NULL)
      {
        std::cout << "Failed to allocate memory in BuildTS" << std::endl;
        // TODO: free the parameter matrix here
        exit(1);
      }

    }
    ts.params = params_tmp;*/

    ts.distributed = false; // by default. set to true if SplitTrainingSet is called
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

void BuildTS_tensor_product(const int &m_size, const double &m_low, const double &m_high, TrainSet &ts)
{

    double *mass_list;
    double m_i, m_j;
    int counter = 0;

    mass_list = (double *)malloc(m_size*sizeof(double));

    if(mass_list==NULL)
    {
        std::cout << "Failed to allocate memory in BuildTS" << std::endl;
        free(mass_list);
        exit(1);
    }

    uniform(m_size, m_low, m_high, mass_list); // m_size equidistant points from m_low to m_high

    // mass_list x mass_list tensor product -- (m1_temp[i],m2_temp[i]) ith training set element
    for(int i = 0; i < m_size; i++)
    {
        m_i = mass_list[i];

        for(int j = 0; j < m_size; j++)
        {
            m_j = mass_list[j];
            ts.m1[counter] = m_i;
            ts.m2[counter] = m_j;

            //ts.params[counter][0] = m_i;
            //ts.params[counter][1] = m_j;

            counter = counter + 1;
        }
    }


    /* -- NOTE: Free m1_tmp, m2_tmp in main through ts -- */
    free(mass_list);

}


void BuildTS_from_file(const char *ts_file, TrainSet &ts)
{

    std::cout << "Reading TS points from file: " << ts_file << std::endl;

    int counter = 0;

    FILE *data;
    data = fopen(ts_file,"r");
    while(fscanf(data, "%lf %lf", &ts.m1[counter],&ts.m2[counter]) != EOF){
        counter = counter + 1;
    }
    fclose(data);

    //ts.params[]

    std::cout << "ts size = " << counter << std::endl;

    if( counter != ts.ts_size){
        std::cout << "TS file size does not match expected size" << std::endl;
        exit(1);
    }

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
        fprintf(data1,"%.15le %.15le\n",ts.m1[i],ts.m2[i]);
    }
    fclose(data1);
}



