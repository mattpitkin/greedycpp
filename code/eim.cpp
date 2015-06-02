#include "eim.hpp"
#include "gsl_helper_functions.hpp"
#include <assert.h>
#include <time.h>

EIM::~EIM() {

  delete[] lebesgue_;
  delete[] p_;
  delete[] rho_;
  gsl_matrix_complex_free(V_);
  gsl_matrix_complex_free(invV_);
}


EIM::EIM(const int basis_dim, const int full_dim, bool err_bounds) :
  basis_dim_(basis_dim),
  full_dim_(full_dim),
  err_bounds_(err_bounds),
  eim_size_(0)
{
  rho_      = new double[basis_dim_];
  p_        = new int[basis_dim_];
  lebesgue_ = new double[basis_dim_];
  V_        = gsl_matrix_complex_alloc(basis_dim_,basis_dim_);
  invV_     = gsl_matrix_complex_alloc(basis_dim_,basis_dim_);
}

gsl_vector_complex* 
EIM::eim_sub_vector(const gsl_vector_complex *u, const int N) {

  // cannot request eim vector larger than current approximation space
  assert(N<=eim_size_);

  // memory should be freed by the calling routine
  gsl_vector_complex *u_eim = gsl_vector_complex_alloc(N);

  for(int i=0;i<N;++i)
    gsl_vector_complex_set(u_eim,i,gsl_vector_complex_get(u,p_[i]));

  return u_eim;
}

// TODO: this routine is really 3 separate things: when building EIM points,
// use the L-version with V updated on the fly.
// when called for validation studies, use precomputed invV or LU+P to get c_eim
gsl_vector_complex* 
EIM::eim_full_vector(const gsl_vector_complex *u,
                     const gsl_matrix_complex *RB_space,
                     const int N) {

  gsl_vector_complex* u_eim = gsl_vector_complex_alloc(full_dim_);
  gsl_vector_complex_set_zero(u_eim);

  // compute the small eim subvector through evaluation (interpolation data)
  gsl_vector_complex* u_sub_eim = eim_sub_vector(u,N);

  // solve for the basis expansion coefficients c_eim
  gsl_vector_complex* c_eim = gsl_vector_complex_alloc(N);
  gsl_permutation* perm = gsl_permutation_alloc(N);
  int s;

  gsl_matrix_complex_view mv = gsl_matrix_complex_submatrix(V_,0,0,N,N);
  gsl_matrix_complex* m = gsl_matrix_complex_alloc(N,N);
  gsl_matrix_complex_memcpy(m, &mv.matrix);

  mygsl::gsl_linalg_complex_L_solve(m,u_sub_eim,c_eim);

  //gsl_linalg_complex_LU_decomp (m,perm,&s);
  //gsl_linalg_complex_LU_solve(m,perm,u_sub_eim,c_eim);


  // evaluate the interpolant on the full set of points
  gsl_complex c_i;
  gsl_vector_complex* e_i = gsl_vector_complex_alloc(full_dim_);
  for (int i=0;i<N;++i) {
    c_i = gsl_vector_complex_get(c_eim,i);
    //std::cout << "ci real " << GSL_REAL(c_i) << " ci imag " << GSL_IMAG(c_i) << std::endl;
    gsl_matrix_complex_get_row(e_i,RB_space,i);
    gsl_vector_complex_scale(e_i,gsl_vector_complex_get(c_eim,i)); // c_i*e_i
    gsl_vector_complex_add(u_eim,e_i); // add c_i*e_i term to empirical interpolant
  }

  gsl_vector_complex_free(u_sub_eim);
  gsl_vector_complex_free(e_i);
  gsl_vector_complex_free(c_eim);
  gsl_matrix_complex_free(m);
  gsl_permutation_free(perm);

  return u_eim;

}

void EIM::update_vandermonde(const gsl_matrix_complex* RB_space,
                             const int N) {

  // RB_space must be updated for residuals BEFORE calling this function.
  // See build_eim_data

  gsl_complex tmpc;
  gsl_vector_complex* r_i = gsl_vector_complex_alloc(full_dim_);

  // Add a row corresponding to the most recent point selection.
  // Notice that while each row of RB_space is a basis element e_i[j]
  // the cols of V_ are filled by e_i (or rather the residuals, see
  // build_eim_data) evaluated at the eim points
  for(int i=0;i<eim_size_;++i) {
    tmpc = gsl_matrix_complex_get(RB_space,i,p_[N-1]);
    gsl_matrix_complex_set(V_,eim_size_-1,i,tmpc);
  }

  // Add a column for the last basis "processed"
  gsl_matrix_complex_get_row(r_i,RB_space,eim_size_-1);
  gsl_vector_complex* v = eim_sub_vector(r_i,eim_size_);
  for(int i=0;i<eim_size_;++i) {
    tmpc = gsl_vector_complex_get(v,i);
    gsl_matrix_complex_set(V_,i,eim_size_-1,tmpc);
  }

  gsl_vector_complex_free(v);
  gsl_vector_complex_free(r_i);
}

void EIM::compute_ith_eim_node(const gsl_vector_complex* v, const int i)
{
  p_[i]   = mygsl::gsl_vector_complex_max_index(v);
  rho_[i] = gsl_complex_abs(gsl_vector_complex_get(v,p_[i]));
  fprintf(stdout,"i = %i, EIM index = %i\n",i,p_[i]);
} 

void EIM::replace_basis_with_residual(gsl_vector_complex* res,
                                      const int res_indx,
                                      const int basis_indx,
                                      gsl_matrix_complex *RB_space)
{
  gsl_complex tmp;
  tmp = gsl_vector_complex_get(res,res_indx);
  tmp = gsl_complex_inverse(tmp);
  gsl_vector_complex_scale(res,tmp);
  gsl_matrix_complex_set_row(RB_space,basis_indx,res);
}

void EIM::rebuild_vandermonde(const gsl_matrix_complex *RB_space)
{

  gsl_vector_complex *e_i = gsl_vector_complex_alloc(full_dim_);


  for(int i=0;i<eim_size_;++i) {


    gsl_matrix_complex_get_row(e_i,RB_space,i);
    gsl_vector_complex* e_i_eim = eim_sub_vector(e_i,eim_size_);

    // TODO: This is transpose of Priscillas code... agreement otherwise... might affect roq program
    // NOTE: when using "_row" there will be a non-trivial permutation matrix
    gsl_matrix_complex_set_col(V_,i,e_i_eim);
    gsl_vector_complex_free(e_i_eim);

  }

  gsl_vector_complex_free(e_i);

}

void EIM::compute_vandermonde_inverse() {

  // Put V into LU form and check that no permutations occurred (sanity check)
  gsl_matrix_complex *V_LU = gsl_matrix_complex_alloc(eim_size_,eim_size_);
  gsl_matrix_complex_memcpy(V_LU,V_);
  gsl_permutation *perm = gsl_permutation_alloc(eim_size_);
  int s;
  gsl_linalg_complex_LU_decomp(V_LU,perm,&s);


  for(int i=0;i<eim_size_;++i) {
    //std::cout << "i - perm[i] " << i - perm->data[i] << std::endl;
    if( (i - perm->data[i])!=0) {
      std::cerr << "Non-trivial permutation matrix...aborting" << std::endl;
      exit(1);
    }
  }

  gsl_linalg_complex_LU_invert(V_LU,perm,invV_);

  gsl_matrix_complex_free(V_LU);
  gsl_permutation_free(perm);

}

// RB_space needs to be writable -- fast EIM version will update
// its rows to be the residuals
void EIM::build_eim_data(gsl_matrix_complex* RB_space) {

  clock_t start, end;

  gsl_vector_complex *u   = gsl_vector_complex_alloc(full_dim_);

  // first iteration of the EIM algorithm
  gsl_matrix_complex_get_row(u,RB_space,0);
  compute_ith_eim_node(u,0);
  eim_size_++;
  replace_basis_with_residual(u,p_[0],0,RB_space); // fast eim -- V_ lower tri.
  update_vandermonde(RB_space,eim_size_);

  // find remaining eim points
  start = clock();
  for(int i=1;i<basis_dim_;++i) {

    gsl_matrix_complex_get_row(u,RB_space,i);

    gsl_vector_complex *res = eim_full_vector(u,RB_space,i); // assemble $I_{i-1}[u]$
    gsl_vector_complex_sub(res,u); // compute residual res = $u - I_{i-1}[u]$

    compute_ith_eim_node(res,i);
    replace_basis_with_residual(res,p_[i],i,RB_space);
    eim_size_++;
    update_vandermonde(RB_space,eim_size_);

    gsl_vector_complex_free(res);
  }

  end = clock();
  double alg_time = ((double) (end - start)/CLOCKS_PER_SEC);
  fprintf(stdout,"Finding the DEIM points took %f seconds.\n", alg_time);

  // There should be a rebuild V and invV matrix here?
  // Will need RB_space again... perhaps pass this back in
  // from the main executable 

  gsl_vector_complex_free(u);

}
