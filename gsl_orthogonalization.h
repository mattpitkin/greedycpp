#ifndef GSL_ORTHOGONALIZATION_H
#define GSL_ORTHOGONALIZATION_H

void MGS(gsl_vector_complex *ru, gsl_vector_complex *ortho_basis,const gsl_matrix_complex *RB_space,const gsl_vector_complex *wQuad, const int dim_RB)
{
/*  Modified GS routine. 

   Input:  RB_space:    an existing orthonormal set of basis vectors
           ortho_basis: basis we shall orthogonalize agains RB_space
           wQuad:       quadrature weights for inner products
           dim_RB:      number of elements currently in RB_space
   Output: ortho_basis: orthonormal basis vector
           ru:           1-by-(dim_RB+1) slice of the R matrix (from QR = A)
*/


    int quad_num = RB_space->size2;
    gsl_complex L2_proj, tmp;
    gsl_vector_complex *basis;


    basis = gsl_vector_complex_alloc(quad_num);

    gsl_vector_complex_set_zero(ru); // if not done, R matrix fills up below diagonal with .5 instead of 0


    for(int i = 0; i < dim_RB; i++)
    {
        gsl_matrix_complex_get_row(basis,RB_space,i);

        /* --- ortho_basis = ortho_basis - L2_proj*basis; --- */
        L2_proj = WeightedInner(wQuad,basis,ortho_basis);
        gsl_vector_complex_set(ru,i,L2_proj);
        gsl_vector_complex_scale(basis,L2_proj); // basis <- basis*L2_proj
        gsl_vector_complex_sub(ortho_basis,basis); // ortho_basis <- ortho_basis - basis

    }

    double nrm = GetNorm_double(ortho_basis,wQuad);
    gsl_complex nrmc;
    GSL_SET_COMPLEX(&nrmc,nrm,0.0);
    gsl_vector_complex_set(ru,dim_RB,nrmc);

    NormalizeVector(ortho_basis,wQuad);

    gsl_vector_complex_free(basis);

}


void IMGS(gsl_vector_complex *ru, gsl_vector_complex *ortho_basis,const gsl_matrix_complex *RB_space,const gsl_vector_complex *wQuad, const int dim_RB)
{
/*  Iterated modified GS routine. 

   Input:  RB_space:    an existing orthonormal set of basis vectors
           ortho_basis: basis we shall orthogonalize agains RB_space
           wQuad:       quadrature weights for inner products
           dim_RB:      number of elements currently in RB_space (ortho_basis is dim_RB+1 element)
   Output: ortho_basis: orthonormal basis vector
           ru:           1-by-(dim_RB+1) slice of the R matrix (from QR = A)

   Hoffmann, "ITERATIVE ALGORITHMS FOR GRAM-SCHMIDT ORTHOGONALIZATION"

*/

    double ortho_condition = .5; // IMGS stopping condition (HARD CODED!!)

    int quad_num     = RB_space->size2;
    int r_size       = ru->size;
    double nrm_prev  = GetNorm_double(ortho_basis,wQuad);
    bool flag        = false;
    int iter         = 0;
    gsl_vector_complex *e,*r_last;
    double nrm_current;
    gsl_complex nrmc_current;

    // --- allocate memory --- //
    e = gsl_vector_complex_alloc(quad_num);
    r_last = gsl_vector_complex_alloc(r_size);

    gsl_vector_complex_memcpy(e,ortho_basis);
    NormalizeVector(e,wQuad);
    gsl_vector_complex_set_zero(ru); // if not done, R matrix fills up below diagonal with .5 instead of 0

    // TODO: code below looks to have some redundent parts -- after working cleanup

    while(!flag)
    {
        gsl_vector_complex_memcpy(ortho_basis,e);
        gsl_vector_complex_set_zero(r_last);

        MGS(r_last,ortho_basis,RB_space,wQuad,dim_RB);

        gsl_vector_complex_add(ru,r_last);
        nrmc_current = gsl_vector_complex_get(r_last,dim_RB);
        nrm_current = GSL_REAL(nrmc_current);

        gsl_vector_complex_scale(ortho_basis,nrmc_current);


        if( nrm_current/nrm_prev <= ortho_condition )
        {
            nrm_prev = nrm_current;
            iter = iter + 1;
            gsl_vector_complex_memcpy(e,ortho_basis);

        }
        else
        {
            flag = true;
        }

        nrm_current  = GetNorm_double(ortho_basis,wQuad);
        GSL_SET_COMPLEX(&nrmc_current,nrm_current,0.0);
        gsl_vector_complex_set(ru,dim_RB,nrmc_current);

        NormalizeVector(ortho_basis,wQuad);

    }

    gsl_vector_complex_free(e);
    gsl_vector_complex_free(r_last);

}

#endif


