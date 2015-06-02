#ifndef EIM_hpp
#define EIM_hpp

#include <gsl/gsl_matrix_complex_double.h>
#include <gsl/gsl_vector_double.h>

// Implementation of the the empirical interpolation method (EIM) as described in
// appendix C of "Two-Step Greedy Algorithm for Reduced Order Quadratures"
// by Antil et al. which modifies the tradition EIM algorithm of Maday XXX: add citation here
//
// A notationally useful reference upon which this code is based can be found in
// "Gravitational wave parameter estimation with compressed
// likelihood evaluations" by Canizares et al. Their appendix B, however,
// uses the basis e_i as opposed to the residual r_i. This leads to 
// worst algorithmic scaling as discussed in both references.


class EIM {
  public:

    //EIM();
    EIM(const int basis_dim, const int full_dim, bool err_bounds);
    ~EIM();

    // While finding the empirical interpolation points and vandermonde matrix,
    // the matrix RB_space will be modified in place to the hold eim residuals
    // this allows N eim points to be found with O(N^3) complexity, which
    // should be compared to the naive O(N^4) complexity
    void build_eim_data(gsl_matrix_complex* RB_space);

    // Constructs an N-vector, a subvector of u, by evaluating u at the first
    // N eimpirical interpolant indicies stored in p_, u_sub = u[p_[0:N]].
    // This routine allocates memory on the heap and must be freed by the caller
    gsl_vector_complex* eim_sub_vector(const gsl_vector_complex *u, const int N);


    // build the full_dim_ sized EIM vector from the N-dimensional approximation 
    // space.
    // TODO: if generating multiple eim vectors, use a precomputed invV_ for speed
    gsl_vector_complex* eim_full_vector(const gsl_vector_complex *u,
                                        const gsl_matrix_complex *RB_space,
                                        const int N);

    // builds Vandermonde matrix V_ from basis + eim points p_
    // why "re"-build? in fast EIM the vandermonde built via 
    // calls to update_vandermonde are evaluations of the residual
    // and not the basis e_i
    void rebuild_vandermonde(const gsl_matrix_complex *RB_space);

    // compute the inverse of V_. Must be called after rebuild_vandermonde
    void compute_vandermonde_inverse();

    // accessors
    inline const int* p() const { return p_; }
    inline const gsl_matrix_complex& V() const { return *V_; }
    inline const gsl_matrix_complex& invV() const { return *invV_; }
    inline const int eim_size() const { return eim_size_; }

  private:

    // When building the EIM approximation, this routine will extend the
    // the Vandermonde matrix V_ from (i-1)-by-(i-1) to i-by-i using
    // the EIM index p_[i].
    void update_vandermonde(const gsl_matrix_complex* RB_space,
                            const int current_eim_index);

    // compute p_[i] and rho_[i] from v
    void compute_ith_eim_node(const gsl_vector_complex* v, const int i);

    // replace ith basis e_i with scaled residual res/res[res_indx]
    void replace_basis_with_residual(gsl_vector_complex* res,
                                     const int res_indx,
                                     const int basis_indx,
                                     gsl_matrix_complex *RB_space);


    double *rho_; // the maximum value of the residual EIM vector. vector of size basis_dim
    int *p_;   // EIM indexes. vector of size basis_dim
    double *lebesgue_; // error bound constant(s). vector of size basis_dim
    gsl_matrix_complex *V_;        // vandermonde matrix. matrix of size basis_dim-by-basis_dim
    gsl_matrix_complex *invV_;     // V^{-1}
    const int basis_dim_;
    const int full_dim_;
    const bool err_bounds_;
    int eim_size_; // when finding the eim points, this will count the current size
};

#endif /* eim.hpp */
