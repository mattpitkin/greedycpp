#ifndef EIM_hpp
#define EIM_hpp

#include <string>

#include <gsl/gsl_matrix_complex_double.h>
#include <gsl/gsl_vector_double.h>

// Implementation of the the empirical interpolation method (EIM) as described in
// appendix C of "Two-Step Greedy Algorithm for Reduced Order Quadratures"
// by Antil et al.
//
// A notationally useful reference upon which this code is based can be found in
// "Gravitational wave parameter estimation with compressed
// likelihood evaluations" by Canizares et al. Appendix B, however,
// uses the basis e_i as opposed to the residual r_i. This leads to a
// worst algorithmic complexity as discussed in both references.


class EIM {
  public:

    // use this constructor when *finding* eim points
    EIM(const int basis_dim, const int full_dim, bool err_bounds);

    // use this constructor when loading data needed for eim evaluations
    EIM(const int basis_dim, const int full_dim,
        gsl_matrix_complex *invV, int *indices);

    ~EIM();

    // While finding the empirical interpolation points and vandermonde matrix,
    // the matrix RB_space is modified inplace to the hold eim residuals.
    // This allows N eim points to be found with O(N^3) complexity, which
    // should be compared to the naive O(N^4) complexity.
    // However, the basis must be reloaded to compute the interpolation matrix
    void build_eim_data(gsl_matrix_complex* RB_space);

    // Constructs an N-vector, a subvector of u, by evaluating u at the first
    // N eimpirical interpolant indicies stored in p_, u_sub = u[p_[0:N]].
    // This routine allocates memory on the heap and must be freed by the caller
    gsl_vector_complex* eim_sub_vector(const gsl_vector_complex *u, const int N);


    // Compute the eim cofficients for the first N basis in RB_space
    void compute_eim_coeffs(const gsl_vector_complex *u_sub_eim,
                            const int N, gsl_vector_complex *c_eim);


    // build full_dim_ sized EIM vector from the N-dim approximation space.
    gsl_vector_complex* eim_full_vector(const gsl_vector_complex *u,
                                        const gsl_matrix_complex *RB_space,
                                        const int N);

    // re-builds the Vandermonde matrix V_ from basis + eim points p_.
    // Why "re"-build? The fast EIM fills V_ with residuals and not the 
    // basis e_i. See aforementioned refs.
    void rebuild_vandermonde(const gsl_matrix_complex *RB_space);

    // compute the inverse of V_. Must be called after rebuild_vandermonde
    void compute_vandermonde_inverse();

    // accessors
    inline const int* p() const { return p_; }
    inline const gsl_matrix_complex& V() const { return *V_; }
    inline const gsl_matrix_complex& invV() const { return *invV_; }
    inline int eim_size() const { return eim_size_; }

  private:
    // -- These three methods are only used when finding the EIM nodes -- //

    // This routine will extend the the Vandermonde matrix V_ from 
    // (i-1)-by-(i-1) to i-by-i using the EIM index p_[i]
    void update_vandermonde(const gsl_matrix_complex* RB_space,
                            const int current_eim_index);

    // Compute p_[i] and rho_[i] from v
    void compute_ith_eim_node(const gsl_vector_complex* v, const int i);

    // Replace ith basis e_i with scaled residual res/res[res_indx]
    void replace_basis_with_residual(gsl_vector_complex* res,
                                     const int res_indx,
                                     const int basis_indx,
                                     gsl_matrix_complex *RB_space);


    std::string mode_;         // whether class created in "build" or "evaluation" mode
    double *rho_;              // Maximum value of the i^th residual EIM vector
    int *p_;                   // EIM indexes. Vector of size basis_dim
    double *lebesgue_;         // error bound constant(s)
    gsl_matrix_complex *V_;    // Square Vandermonde matrix of size basis_dim_
    gsl_matrix_complex *invV_; // V^{-1}
    const int basis_dim_;      // Size of approximation space 
    const int full_dim_;       // Size of full space
    const bool err_bounds_;    // whether to compute error bounds 
    int eim_size_;             // Internal counter used by build_eim_data
};

#endif /* eim.hpp */
