// AUTHOR :  Scott Field 
//           sfield@astro.cornell.edu
//
// DATE: Nov 10, 2014
//
// PURPOSE: mpi version of greedy (pivoted MGS) algorithm 


//-- NOTE: to add a model ONLY my_models.h needs to be modified --//


// --- DEFINE PRECOMPILER FLAG COMPILE_WITHOUT_MPI from makefile --- //

#ifndef COMPILE_WITHOUT_MPI
#include <mpi.h>
#endif

#ifdef USE_OPENMP
#include <omp.h>
#endif

#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <cmath>

#include "gsl_helper_functions.hpp"
#include "training_set.hpp"
#include "quadratures.hpp"
#include "parameters.hpp"
#include "utils.h"
#include "my_models.h"

struct errorDouble_rankInt { 
  double error; 
  int rank; 
}; 


void WriteGreedyInfo(const int dim_RB,
                     gsl_matrix_complex *RB_space, // "const": need to resize
                     gsl_matrix_complex *R_matrix, // same "const"
                     const double *app_err,
                     const int *sel_rows,
                     const TrainingSetClass &ts,
                     const char * output_dir,
                     const char *datatype)
{
  FILE *err_data, *pts_data;
  FILE *rb_data, *r_data;
  char err_filename[200];
  char pts_filename[200];
  char rb_filename[200];
  char r_filename[200];

  // directly resize matrix so only non-zero basis are output //
  // (undo this before leaving routine) //
  int RB_space_size1 = RB_space->size1; 
  int R_size1        = R_matrix->size1;
  int R_size2        = R_matrix->size2;
  RB_space->size1 = dim_RB+1;
  R_matrix->size1 = dim_RB+1;
  R_matrix->size2 = dim_RB+1; 


  strcpy(err_filename,output_dir);
  strcat(err_filename,"/ApproxErrors.txt");
  strcpy(pts_filename,output_dir);
  strcat(pts_filename,"/GreedyPoints.txt");

  //--- write errors and greedy points to text file ---//
  err_data = fopen(err_filename,"w");
  pts_data = fopen(pts_filename,"w");
  for(int i = 0; i <= dim_RB ; i++)
  {
    fprintf(err_data,"%1.15e\n",app_err[i]);
    ts.fprintf_ith(pts_data,sel_rows[i]);
    fprintf(pts_data,"\n");
  }
  fclose(err_data);
  fclose(pts_data);


  bool wrote = false;
  if(strcmp(datatype,"txt") == 0 || strcmp(datatype,"both") == 0){
    strcpy(rb_filename,output_dir);
    strcat(rb_filename,"/Basis");
    strcpy(r_filename,output_dir);
    strcat(r_filename,"/R");

    mygsl::gsl_matrix_complex_fprintf(rb_filename,RB_space);
    // TODO: valgrind reports memory errors here
    mygsl::gsl_matrix_complex_fprintf(r_filename,R_matrix);
    wrote = true;
  } 
  if(strcmp(datatype,"bin") == 0 || strcmp(datatype,"both") == 0){
    strcpy(rb_filename,output_dir);
    strcat(rb_filename,"/Basis.bin");
    strcpy(r_filename,output_dir);
    strcat(r_filename,"/R.bin");

    rb_data = fopen(rb_filename,"wb");
    gsl_matrix_complex_fwrite(rb_data,RB_space);
    fclose(rb_data);
    r_data = fopen(r_filename,"wb");
    gsl_matrix_complex_fwrite(r_data,R_matrix);
    fclose(r_data);
    wrote = true;
  }
  if(strcmp(datatype,"npy") == 0){

    strcpy(rb_filename,output_dir);
    strcat(rb_filename,"/Basis.npy");
    strcpy(r_filename,output_dir);
    strcat(r_filename,"/R.npy");
    mygsl::gsl_matrix_complex_npy_save(rb_filename,RB_space);
    mygsl::gsl_matrix_complex_npy_save(r_filename,R_matrix);
    wrote = true;
  }

  if(!wrote){
    fprintf(stderr,"file type not supported");
    exit(1);
  }

  // restore matrix to its original size //
  RB_space->size1 = RB_space_size1;
  R_matrix->size1 = R_size1;
  R_matrix->size2 = R_size2; 

}

// --- GreedyWorker and Master are removed for serial builds --- //
#ifndef COMPILE_WITHOUT_MPI
void GreedyWorker(const int rank,
                  const Parameters &params,
                  const gsl_vector_complex *wQuad,
                  const gsl_matrix_complex *A,
                  const TrainingSetClass &ts)
{

// worker routine for computing computationally intensive part of greedy //

  // -- unpack params class --//
  const int max_RB      = params.max_RB();
  const int seed_global = params.seed();
  const double tol      = params.tol();

  int dim_RB  = 1;
  int cols    = wQuad->size;
  bool continue_work = true;
  int worst_global, worst_local, worst_rank;
  double rb_inc_r, rb_inc_i, tmp;
  double *errors, *vec_real, *vec_imag;
  gsl_complex tmpc;
  gsl_vector_complex *last_rb, *row_vec;
  gsl_matrix_complex *project_coeff;

  last_rb       = gsl_vector_complex_alloc(cols);
  row_vec       = gsl_vector_complex_alloc(cols);
  project_coeff = gsl_matrix_complex_alloc(max_RB,ts.matrix_sub_size()[rank]);
  errors        = new double[ts.matrix_sub_size()[rank]];
  vec_real      = new double[cols];
  vec_imag      = new double[cols];

  errorDouble_rankInt reduce_data; 
  errorDouble_rankInt global_reduce_data;

  int *worst_workers_mpi = NULL;
  double *worst_errs_mpi = NULL;

  //bool useEuc = mygsl::IsConstantVector( wQuad );
  bool useEuc = false;

  fprintf(stdout,"Worker %i was given %i matrix elements from %i to %i\n",
          rank,
          ts.matrix_sub_size()[rank],
          ts.mystart()[rank],
          ts.myend()[rank]-1);

  // -- pass seed back to master -- //
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Bcast(&worst_rank, 1, MPI_INT,0,MPI_COMM_WORLD);
  // -- return seed to master -- //
  if( (worst_rank-1) == rank) {
      worst_local = seed_global - ts.mystart()[rank];
      gsl_matrix_complex_get_row(row_vec,A,worst_local);
      mygsl::gsl_vector_complex_parts(vec_real,vec_imag,row_vec);
      MPI_Send(vec_real,cols, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
      MPI_Send(vec_imag,cols, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
  }

  while(continue_work)
  {
    // -- wait for new basis and start next sweep -- //
    //MPI_Barrier(MPI_COMM_WORLD);

    // -- receive new rb -- //
    MPI_Bcast(vec_real, cols, MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(vec_imag, cols, MPI_DOUBLE,0,MPI_COMM_WORLD);
    mygsl::make_gsl_vector_complex_parts(vec_real,vec_imag,last_rb);

    #ifdef USE_OPENMP // due to extra allocs, avoid this code if not using omp
    #pragma omp parallel
    {

      /*#pragma omp master
      {
        std::cout<<"threads (GreedyWorker) = "<<omp_get_num_threads()<<std::endl;
      }*/

      // every variable declared here is thread private (thread-safe)
      gsl_vector_complex *ts_el_omp;
      ts_el_omp = gsl_vector_complex_alloc(cols);

      #pragma omp for
      for(int i = 0; i < ts.matrix_sub_size()[rank]; i++)
      {
        gsl_matrix_complex_get_row(ts_el_omp,A,i);
        gsl_matrix_complex_set(project_coeff,dim_RB-1,i,
                           mygsl::InnerProduct(wQuad,last_rb,ts_el_omp,useEuc));
        errors[i] = 1.0 - mygsl::SumColumn(project_coeff,i,dim_RB);
      }
      gsl_vector_complex_free(ts_el_omp);
    }
    #else

    // Compute overlaps of pieces of A with rb_new //
    for(int i = 0; i < ts.matrix_sub_size()[rank]; i++)
    {
      gsl_matrix_complex_get_row(row_vec,A,i);
      gsl_matrix_complex_set(project_coeff,dim_RB-1,i,
                         mygsl::InnerProduct(wQuad,last_rb,row_vec,useEuc));
      errors[i] = 1.0 - mygsl::SumColumn(project_coeff,i,dim_RB);
    }
    #endif

    // -- find worst error here (worst_local is worker's row index) -- //
    FindMax2(errors,ts.matrix_sub_size()[rank],tmp,worst_local);
    worst_global = ts.mystart()[rank] + worst_local; // global row index

    reduce_data.error = tmp;
    reduce_data.rank  = rank;

    MPI_Allreduce(&reduce_data,&global_reduce_data,1,
                  MPI_DOUBLE_INT,MPI_MAXLOC,MPI_COMM_WORLD);

    worst_rank       = global_reduce_data.rank+1;
    double worst_err = global_reduce_data.error;

    // -- return basis to master -- //
    if( (worst_rank-1) == rank) {
        gsl_matrix_complex_get_row(row_vec,A,worst_local);
        mygsl::gsl_vector_complex_parts(vec_real,vec_imag,row_vec);
        MPI_Send(vec_real,cols, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        MPI_Send(vec_imag,cols, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        MPI_Send(&worst_global, 1, MPI_INT, 0, 0,MPI_COMM_WORLD);
    }

    ++dim_RB;

    // -- decide if another greedy sweep will be needed -- //
    if( (dim_RB == max_RB) || worst_err < tol) {
      continue_work = false;
    }

  }

  free(worst_workers_mpi); // free NULL pointers?
  free(worst_errs_mpi);

  gsl_vector_complex_free(last_rb);
  gsl_vector_complex_free(row_vec);
  gsl_matrix_complex_free(project_coeff);
  delete [] vec_real;
  delete [] vec_imag;
  delete [] errors;
}

void GreedyMaster(const int size,
                  const gsl_vector_complex *wQuad,
                  const TrainingSetClass &ts,
                  const Parameters &params)
{
// Input: 
//          size: number of procs. if 1 serial mode assumed
//
// Output (written to file):
//          sel_rows: row index defining reduced basis. sel_rows[0] = seed
//          dim_RB: number of greedy_points

  fprintf(stdout,"Starting greedy algorithm...\n");

  // -- unpackage params class -- //
  const int max_RB = params.max_RB(); // global indexing
  const int seed   = params.seed();
  const double tol = params.tol();
  const char * output_dir = params.output_dir().c_str();
  const char * output_data_format = params.output_data_format().c_str();

  const int rows = ts.ts_size();// number of rows to approximate
  const int cols = wQuad->size; // samples (for quadrature)
  int *greedy_points;           // selectted greedy points (row selection)
  double *greedy_err;           // approximate error
  clock_t start, end;           // for algorithm timing experiments
  clock_t start1, end1;         // for algorithm timing experiments
  clock_t start_or, end_or;     // time between orthogonalizations
  clock_t start_sw, end_sw;     // greedy sweep times
  double alg_time,or_t,sw_t;    // for algorithm timing experiments
  double worst_err;             // errors in greedy sweep
  int worst_app, worst_worker, worst_rank;             // worst error stored
  double rb_inc_r, rb_inc_i;
  double *vec_real, *vec_imag;
  int *worst_workers_mpi;
  double *worst_errs_mpi;
  gsl_complex tmpc;
  int dummy_mpi_int       = -1;
  double dummy_mpi_double = -1.0;
  int dim_RB              = 1;
  double total_ortho_time = 0.0;
  double total_sweep_time = 0.0;
  bool continue_work = true;

  errorDouble_rankInt dummy_reduce_data;
  dummy_reduce_data.error = -1;
  dummy_reduce_data.rank  = -2;
  errorDouble_rankInt global_reduce_err_index_data;
  gsl_vector_complex *ortho_basis, *ru;
  gsl_matrix_complex *RB_space, *R_matrix;

  // --- this memory should be freed here --- //
  ortho_basis       = gsl_vector_complex_alloc(cols);
  vec_real          = new double [cols];
  vec_imag          = new double [cols];
  worst_workers_mpi = new int[size];
  worst_errs_mpi    = new double[size];
  RB_space          = gsl_matrix_complex_alloc(max_RB,cols); 
  R_matrix          = gsl_matrix_complex_alloc(max_RB,max_RB);
  greedy_points     = new int[max_RB];
  greedy_err        = new double[max_RB];
  ru                = gsl_vector_complex_alloc(max_RB);

  // --- initialize algorithm with seed --- //
  int seed_rank = ts.FindRowIndxRank(seed);
  fprintf(stdout,"seed index %i on proc rank %i\n",seed,seed_rank);

  // -- request seed from worker -- //
  MPI_Barrier(MPI_COMM_WORLD);
  seed_rank = seed_rank+1;
  MPI_Bcast(&seed_rank, 1, MPI_INT,0,MPI_COMM_WORLD);
  MPI_Recv(vec_real, cols, MPI_DOUBLE, seed_rank, 0,\
           MPI_COMM_WORLD,MPI_STATUS_IGNORE);
  MPI_Recv(vec_imag, cols, MPI_DOUBLE, seed_rank, 0,\
           MPI_COMM_WORLD,MPI_STATUS_IGNORE);
  mygsl::make_gsl_vector_complex_parts(vec_real,vec_imag,ortho_basis);
  gsl_matrix_complex_set_row(RB_space,0,ortho_basis);

  // -- send first basis to all workers -- // 
  mygsl::gsl_vector_complex_parts(vec_real,vec_imag,ortho_basis);
  //MPI_Barrier(MPI_COMM_WORLD);
  MPI_Bcast(vec_real, cols, MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(vec_imag, cols, MPI_DOUBLE,0,MPI_COMM_WORLD);

  GSL_SET_COMPLEX(&tmpc,1.0,0.0); // assumes normalized solutions
  gsl_matrix_complex_set(R_matrix,0,0,tmpc);
  greedy_points[0] = seed;
  dim_RB           = 1;
  greedy_err[0]    = 1.0;

  // --- Continue approximation until tolerance satisfied --- //
  start = clock();
  while(continue_work)
  {

    start1 = clock();

    // -- gather worst (local) info from workers -- //
    start_sw = clock();
    MPI_Allreduce(&dummy_reduce_data,&global_reduce_err_index_data,1,
               MPI_DOUBLE_INT,MPI_MAXLOC,MPI_COMM_WORLD);
    worst_rank = global_reduce_err_index_data.rank + 1;
    worst_err  = global_reduce_err_index_data.error;
    end_sw = clock();

    // -- receive row basis from worker proc worst_rank -- //
    MPI_Recv(vec_real, cols, MPI_DOUBLE, worst_rank, 0,\
             MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    MPI_Recv(vec_imag, cols, MPI_DOUBLE, worst_rank, 0,\
             MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    MPI_Recv(&worst_app, 1, MPI_INT, worst_rank, 0,\
             MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    mygsl::make_gsl_vector_complex_parts(vec_real,vec_imag,ortho_basis);

    // --- record worst approximated row element (basis) --- //
    greedy_points[dim_RB] = worst_app;
    greedy_err[dim_RB]    = worst_err;

    // --- orthogonalize solution/vector --- //
    start_or = clock();
    mygsl::IMGS(ru,ortho_basis,RB_space,wQuad,dim_RB); // IMGS is default
    //mygsl::MGS(ru,ortho_basis,RB_space,wQuad,dim_RB);
    end_or = clock();

    // -- decide if another greedy sweep will be needed -- //
    if( (dim_RB+1 == max_RB) || worst_err < tol){
      continue_work = false;
    }
    else {

      // -- send last orthonormal rb to work procs -- //
      mygsl::gsl_vector_complex_parts(vec_real,vec_imag,ortho_basis);
      //MPI_Barrier(MPI_COMM_WORLD);
      MPI_Bcast(vec_real, cols, MPI_DOUBLE,0,MPI_COMM_WORLD);
      MPI_Bcast(vec_imag, cols, MPI_DOUBLE,0,MPI_COMM_WORLD);

    }

    // -- update basis matrix RB_space and R matrix -- //
    gsl_matrix_complex_set_row(R_matrix,dim_RB,ru);
    gsl_matrix_complex_set_row(RB_space,dim_RB,ortho_basis);


    ++dim_RB;

    // --- update timers and output info --- //
    end1     = clock();
    or_t     = ((double) (end_or- start_or)/CLOCKS_PER_SEC);
    sw_t     = ((double) (end_sw - start_sw)/CLOCKS_PER_SEC);
    alg_time = ((double) (end1 - start1)/CLOCKS_PER_SEC);
    total_ortho_time += or_t;
    total_sweep_time += sw_t;
    fprintf(stdout,"RB %i | pivot # %i | err %1.3e | ortho time %1.3e "
                   "| sweep time %1.3e | total time %1.3e\n",\
            dim_RB,worst_app,worst_err,or_t,sw_t,alg_time);

  }
  end = clock();

  alg_time = ((double) (end - start)/CLOCKS_PER_SEC);
  fprintf(stdout,"Greedy routine took %f seconds\n",alg_time);
  fprintf(stdout,"GS routine took %f seconds\n",total_ortho_time);
  fprintf(stdout,"Greedy - GS took %f seconds\n",alg_time-total_ortho_time);
  fprintf(stdout,"Greedy sweeps took %f seconds\n",total_sweep_time);
  dim_RB = dim_RB - 1;

  // -- output relevant information -- //
  WriteGreedyInfo(dim_RB,RB_space,R_matrix,greedy_err,\
                  greedy_points,ts,output_dir,output_data_format);

  gsl_vector_complex_free(ortho_basis);
  delete [] vec_real;
  delete [] vec_imag;
  delete [] worst_workers_mpi;
  delete [] worst_errs_mpi;
  gsl_matrix_complex_free(RB_space); 
  gsl_matrix_complex_free(R_matrix);
  delete [] greedy_points;
  delete [] greedy_err;
  gsl_vector_complex_free(ru);

}
#endif // end mpi disabled code (starts with GreedyWorker)

void Greedy(const Parameters &params,
            const gsl_matrix_complex *A,
            const gsl_vector_complex *wQuad,
            const TrainingSetClass &ts)
{
// Input: 
//          A: matrix (each row is a "solution", cols are quadrature samples)
//          seed: first greedy pick
//          tol: approximation tolerance
//
// Output (written to file):
//          sel_rows: row index defining reduced basis. sel_rows[0] = seed
//          dim_RB: number of greedy_points

  fprintf(stdout,"Starting greedy algorithm in serial mode...\n");

  //bool useEuc = mygsl::IsConstantVector( wQuad );
  bool useEuc = false;

  // -- unpack parameter class here -- //
  const int seed                 = params.seed();
  const int max_RB               = params.max_RB();
  const double tol               = params.tol();
  const char *output_dir         = params.output_dir().c_str();
  const char *output_data_format = params.output_data_format().c_str();
  const int ts_size              = params.ts_size();

  const int rows = A->size1; // number of rows to approximate
  const int cols = A->size2; // samples (for quadrature)
  int *greedy_points;        // selectted greedy points (row selection)
  double *greedy_err;        // approximate error
  clock_t start, end;        // for algorithm timing experiments
  clock_t start1, end1;      // for algorithm timing experiments
  double alg_time;           // for algorithm timing experiments
  double worst_err;          // errors in greedy sweep
  int worst_app;             // worst error stored
  gsl_complex tmpc;          // worst error temp
  bool continue_work = true;


  gsl_vector_complex *ts_el, *last_rb, *ortho_basis, *ru;
  gsl_matrix_complex *RB_space, *R_matrix;
  double *errors;                    // approximation errors vs dim_RB
  gsl_matrix_complex *project_coeff; // h = coeff_i e_i


  // --- this memory should be freed here --- //
  ts_el         = gsl_vector_complex_alloc(cols);
  last_rb       = gsl_vector_complex_alloc(cols);
  ortho_basis   = gsl_vector_complex_alloc(cols);
  ru            = gsl_vector_complex_alloc(max_RB);
  errors        = new double[rows];
  project_coeff = gsl_matrix_complex_alloc(max_RB,rows);
  greedy_points = new int[max_RB];
  greedy_err    = new double[max_RB];
  RB_space      = gsl_matrix_complex_alloc(max_RB,cols); 
  R_matrix      = gsl_matrix_complex_alloc(max_RB,max_RB);

  // --- initialize algorithm with seed --- //
  gsl_matrix_complex_get_row(ts_el,A,seed);
  gsl_matrix_complex_set_row(RB_space,0,ts_el);
  GSL_SET_COMPLEX(&tmpc,1.0,0.0);
  gsl_matrix_complex_set(R_matrix,0,0,tmpc);  // assumes normalized solutions

  greedy_points[0] = seed;
  int dim_RB       = 1;
  greedy_err[0]    = 1.0;

  // --- Continue approximation until tolerance satisfied --- //
  start = clock();

  #ifdef USE_OPENMP
  double omp_start = omp_get_wtime();
  #endif
 
  while(continue_work)
  {

    start1 = clock();

    gsl_matrix_complex_get_row(last_rb,RB_space,dim_RB-1); // previous basis


    // --- Loop over training set ---//
    #ifdef USE_OPENMP // due to extra allocs, avoid this code if not using omp
    #pragma omp parallel
    {

      #pragma omp master
      {
        std::cout<<"threads (Greedy) = "<<omp_get_num_threads()<<std::endl;
      }

      // every variable declared here is thread private (thread-safe)
      gsl_vector_complex *ts_el_omp;
      ts_el_omp = gsl_vector_complex_alloc(cols);

      #pragma omp for
      for(int i = 0; i < rows; i++)
      {
        gsl_matrix_complex_get_row(ts_el_omp,A,i);
        gsl_matrix_complex_set(project_coeff,dim_RB-1,i,
                           mygsl::InnerProduct(wQuad,last_rb,ts_el_omp,useEuc));
        errors[i] = 1.0 - mygsl::SumColumn(project_coeff,i,dim_RB);
      }
      gsl_vector_complex_free(ts_el_omp);
    }
    #else
    // Compute overlaps of pieces of A with rb_new //
    for(int i = 0; i < rows; i++)
    {
      gsl_matrix_complex_get_row(ts_el,A,i);
      gsl_matrix_complex_set(project_coeff,dim_RB-1,i,
                            mygsl::InnerProduct(wQuad,last_rb,ts_el,useEuc));
      errors[i] = 1.0 - mygsl::SumColumn(project_coeff,i,dim_RB);
    }
    #endif


    // --- find worst represented ts element, add to basis --- //
    FindMax2(errors,rows,worst_err,worst_app);
    greedy_points[dim_RB] = worst_app;
    greedy_err[dim_RB]    = worst_err;

    // -- decide if another greedy sweep is needed -- //
    if( (dim_RB+1 == max_RB) || (worst_err < tol) || (ts_size == dim_RB) ){
      continue_work = false;
    }


    // --- add worst approximated solution to basis set --- //
    gsl_matrix_complex_get_row(ortho_basis,A,worst_app);
    mygsl::IMGS(ru,ortho_basis,RB_space,wQuad,dim_RB); // IMGS is default
    //mygsl::MGS(ru,ortho_basis,RB_space,wQuad,dim_RB);
    gsl_matrix_complex_set_row(RB_space,dim_RB,ortho_basis);
    gsl_matrix_complex_set_row(R_matrix,dim_RB,ru);
    dim_RB = dim_RB + 1;

    end1 = clock();
    alg_time = ((double) (end1 - start1)/CLOCKS_PER_SEC);

    fprintf(stdout,"RB %i || row selected %i || error %1.4e || time %f\n",\
            dim_RB,worst_app,worst_err,alg_time);

  }
  end = clock();
  alg_time = ((double) (end - start)/CLOCKS_PER_SEC);

  #ifdef USE_OPENMP
  double omp_end  = omp_get_wtime();
  double omp_time = omp_end - omp_start;
  #else
  double omp_time = alg_time;
  #endif
 
  fprintf(stdout,"Building approximation space took %f cpu seconds and %f wall seconds \n",alg_time,omp_time);
  dim_RB = dim_RB - 1;

  // -- output relevant information -- //
  WriteGreedyInfo(dim_RB,RB_space,R_matrix,greedy_err,\
                  greedy_points,ts,output_dir,output_data_format);

  gsl_vector_complex_free(ts_el);
  gsl_vector_complex_free(last_rb);
  gsl_vector_complex_free(ortho_basis);
  gsl_vector_complex_free(ru);
  delete [] errors;
  gsl_matrix_complex_free(project_coeff);
  delete [] greedy_points;
  delete [] greedy_err;
  gsl_matrix_complex_free(RB_space);
  gsl_matrix_complex_free(R_matrix);

}

int main (int argc, char **argv) {

  int rank     = 0;  // needed for serial mode too
  int size_mpi = 1;  // needed for serial mode too
  bool high_verbosity; // controls output based on rank

  std::string cfg_file(argv[1]);

  #ifndef COMPILE_WITHOUT_MPI
  // --- setup MPI info ---//
  MPI::Init(argc, argv);

  // # procs this job is using (size) and processor rank of this thread //
  // Ex: if "mpirun -np 2", size = 2. "CPU1" will be 0 and "CPU2" will be 1 //
  rank     = MPI::COMM_WORLD.Get_rank();
  size_mpi = MPI::COMM_WORLD.Get_size();

  char name[MPI_MAX_PROCESSOR_NAME];
  int len;
  memset(name,0,MPI_MAX_PROCESSOR_NAME);
  MPI::Get_processor_name(name,len);
  memset(name+len,0,MPI_MAX_PROCESSOR_NAME-len);

  printf("Procs = %i, My rank = %i, My name = %s\n",size_mpi,rank,name);
  MPI_Barrier(MPI_COMM_WORLD);
  #endif

  if( rank==0 ) {
    high_verbosity = true;
  } 
  else {
    high_verbosity = false;
  } 

  //----- Checking the number of Variables passed to the Executable -----//
  if (argc != 2) {
    std::cerr << "Pass only location of configuration file (ends in .cfg)" 
              << std::endl;
    exit(1);
  }
  if(high_verbosity) {
    // GIT_SHA1 is passed through compiler flags
    std::cout << "Compiled code git hash is " << GIT_SHA1 << std::endl;
    std::cout << "parameter file is: " << cfg_file << std::endl;
  }

  #ifndef COMPILE_WITHOUT_MPI
  MPI_Barrier(MPI_COMM_WORLD);
  #endif

  //--- Read input file. If there is an error, report and exit.
  //--- Parameters class contains relevant information about parameters 
  Parameters *params_from_file = new Parameters(argv,high_verbosity);
  if( rank==0 ) {
    std::cout << *params_from_file;
  }

  #ifndef COMPILE_WITHOUT_MPI
  MPI_Barrier(MPI_COMM_WORLD);
  #endif

  // Variables which need to be decarled, but not set in parameter file //
  gsl_matrix_complex *TS_gsl;
  gsl_vector_complex *wQuad;
  gsl_vector *xQuad;
  char ts_filename[100];
  char ts_pts_filename[100];
  char shell_command[200];
  int gsl_status;
  int ts_size;
  clock_t start, end;

  if(params_from_file->export_tspace() && size_mpi > 1 && high_verbosity){
      fprintf(stderr,"Training space exporting only works with 1 processor!");
      fprintf(stderr," Your training space will not be exported.\n");
  }

  // returns wQuad and xQuad for the full, not reduced, quadrature rule //
  SetupQuadratureRule(&wQuad,&xQuad,params_from_file);

  // Creating Run Directory //
  if(size_mpi == 1 || rank == 0){

    strcpy(shell_command, "mkdir -p -m700 ");
    strcat(shell_command, params_from_file->output_dir().c_str());
    int ret = system(shell_command);
    if(ret == -1) {
      std::cerr << "Could not make a run directory" << std::endl;
      exit(1);
    }

    // Copy basis building cfg file to the output folder //
    // TODO: next two tasks should be a single routine
    std::string copy_filename(params_from_file->output_dir());
    copy_filename.append("/run_settings.cfg");
    if(copy_filename.compare(cfg_file) != 0) {
      std::cout << "Input config " << cfg_file
                << " does not exist at the location " << copy_filename
                << " ... Copying now." << std::endl;

      std::ifstream src(cfg_file.c_str(),std::ios::binary);
      std::ofstream dst(copy_filename.c_str(),std::ios::binary);
      dst << src.rdbuf();
      src.close();
      dst.close();
    }
    else {
      std::cout << "Input config " << cfg_file
                << " exists at the location " << copy_filename
                << " ... No copy made." << std::endl;
    }

    // Create cfg file to be modified/used for validation studies //
    std::string validation_cfg(params_from_file->output_dir());
    validation_cfg.append("/validations_setup.cfg");
    std::ifstream src1(cfg_file.c_str(),std::ios::binary);
    std::ofstream dst1(validation_cfg.c_str(),std::ios::binary);
    dst1 << src1.rdbuf();
    src1.close();
    dst1.close();

    if(high_verbosity) {
      #ifdef USE_OPENMP
      fprintf(stdout,"openMP enabled. Max threads = %i and num procs = %i\n",
             omp_get_max_threads(), omp_get_num_procs());
      #else
      fprintf(stdout,"openMP disabled\n");
      #endif
    }
  }

  // allocates dynamic memory, fills up training set with values //
  TrainingSetClass *ptspace_class = new TrainingSetClass(params_from_file,\
                                                         size_mpi);

  // TS_gsl filled by evaluating model at ptspace_class->params_ //
  // NOTE: GSL error handler will abort if too much memory requested
  // Ex: size = 1 (serial mode) => rank=0 for 5th argument of FillTrainingSet
  #ifdef USE_OPENMP
  double omp_start = omp_get_wtime();
  #else
  start = clock();
  #endif  

  if(size_mpi == 1) {
    TS_gsl = gsl_matrix_complex_alloc(ptspace_class->ts_size(),xQuad->size);
    mymodel::FillTrainingSet(TS_gsl,xQuad,wQuad,*ptspace_class,rank);

    #ifdef USE_OPENMP
    double omp_end = omp_get_wtime();
    double ts_time = omp_end - omp_start;
    #else
    end = clock();
    double ts_time = ((double) (end - start)/CLOCKS_PER_SEC);
    #endif
    fprintf(stdout,"Filled TS in %f seconds\n",ts_time);

    Greedy(*params_from_file,TS_gsl,wQuad,*ptspace_class);
  }
  else{
    #ifndef COMPILE_WITHOUT_MPI
    if(rank != 0){
      TS_gsl = gsl_matrix_complex_alloc(ptspace_class->matrix_sub_size()[rank-1],xQuad->size);
      mymodel::FillTrainingSet(TS_gsl,xQuad,wQuad,*ptspace_class,rank-1);

      #ifdef USE_OPENMP
      double omp_end = omp_get_wtime();
      double ts_time = omp_end - omp_start;
      #else
      end = clock(); 
      double ts_time = ((double) (end - start)/CLOCKS_PER_SEC);
      #endif
      fprintf(stdout,"Rank %i filled TS in %f seconds\n",rank,ts_time);

    }

    if(rank == 0) {
      GreedyMaster(size_mpi,wQuad,*ptspace_class,*params_from_file);
    }
    else{
      GreedyWorker(rank-1,*params_from_file,wQuad,TS_gsl,*ptspace_class);
      gsl_matrix_complex_free(TS_gsl);
    }
    #else
    fprintf(stderr,"Code compiled with mpi yet size_mpi > 1...\n");
    exit(1);
    #endif
  }


  // -- output a variety of extra information if requested -- //
  if(rank == 0)
  {
    // -- output quadrature weights -- //
    FILE *outfile, *ts_pts;

    strcpy(shell_command, params_from_file->output_dir().c_str());
    strcat(shell_command,"/quad_rule.txt");

    outfile = fopen(shell_command,"w");
    for(int i = 0; i < wQuad->size ; i++) {
      fprintf(outfile,"%1.14f %1.14f\n",gsl_vector_get(xQuad,i),\
              GSL_REAL(gsl_vector_complex_get(wQuad,i)));
    }
    fclose(outfile);

    if(size_mpi == 1){
      if(params_from_file->export_tspace()){
        // -- output some waveform(s) for diagnostics -- //
        strcpy(ts_filename,params_from_file->output_dir().c_str());
        strcat(ts_filename,"/TSpace");
        mygsl::gsl_matrix_complex_fprintf(ts_filename,TS_gsl);

        // -- output these waveform(s) parameter values -- //
        strcpy(ts_pts_filename,params_from_file->output_dir().c_str());
        strcat(ts_pts_filename,"/TSpacePts.txt");
        ts_pts = fopen(ts_pts_filename,"w");
        for(int i = 0; i < ptspace_class->ts_size(); i++) {
          ptspace_class->fprintf_ith(ts_pts,i);
          fprintf(ts_pts,"\n");
        }
        fclose(ts_pts);


      }
      gsl_matrix_complex_free(TS_gsl);
    }
  }


  gsl_vector_complex_free(wQuad);
  gsl_vector_free(xQuad);

  delete ptspace_class;
  ptspace_class = NULL;

  delete params_from_file;
  params_from_file = NULL;

  #ifndef COMPILE_WITHOUT_MPI
  // Tell the MPI library to release all resources it is using
  MPI::Finalize();
  #endif
}

