#ifndef gsl_helper_functions_h
#define gsl_helper_functions_h

void gsl_vector_complex_parts(double *v_real, double *v_imag, const gsl_vector_complex * v_gsl)
{

    for(int i = 0; i < v_gsl->size; i++)
    {
        v_real[i] = GSL_REAL(gsl_vector_complex_get(v_gsl,i));
        v_imag[i] = GSL_IMAG(gsl_vector_complex_get(v_gsl,i));
    }
}

void gsl_vector_complex_gsl_parts(gsl_vector *v_real, gsl_vector *v_imag, const gsl_vector_complex * v_gsl)
{

    double tmp_r, tmp_i;

    for(int i = 0; i < v_gsl->size; i++)
    {
        tmp_r = GSL_REAL(gsl_vector_complex_get(v_gsl,i));
        tmp_i = GSL_IMAG(gsl_vector_complex_get(v_gsl,i));

        gsl_vector_set(v_real, i,tmp_r);
        gsl_vector_set(v_imag, i,tmp_i);
    }
}

void make_gsl_vector_complex_parts(const double *v_real, const double *v_imag,gsl_vector_complex * v_gsl)
{
    gsl_complex tmpc;

    for(int i = 0; i < v_gsl->size; i++){
        GSL_SET_COMPLEX(&tmpc,v_real[i],v_imag[i]);
        gsl_vector_complex_set(v_gsl,i,tmpc);
    }
}

void gsl_matrix_complex_fprintf_part(FILE *rb_data, const gsl_matrix_complex * m_gsl,const char *part)
{
/* output the real part of matrix m_gsl, using the first cols columns and rows rows, to file. */

// NOTE: text output so that A = R * Q is desired QR decomposition //

    const int cols = m_gsl->size2;
    const int rows = m_gsl->size1;

    gsl_vector_complex *v_gsl;
    double *vec_real, *vec_imag;
    v_gsl         = gsl_vector_complex_alloc(cols);
    vec_real      = (double *)malloc(cols*sizeof(double));
    vec_imag      = (double *)malloc(cols*sizeof(double));

    for(int ii = 0; ii < rows; ii++)
    {  
        gsl_matrix_complex_get_row(v_gsl,m_gsl, ii);
        gsl_vector_complex_parts(vec_real,vec_imag,v_gsl);

        for(int jj = 0; jj < cols; jj++)
        {  
            if(strcmp(part,"real") == 0){
                fprintf(rb_data, "%1.12e\t",vec_real[jj]);
            }
            else if(strcmp(part,"imag") == 0){
                fprintf(rb_data, "%1.12e\t",vec_imag[jj]);
            }
            else{
                fprintf(stderr,"part unknown");
                exit(1);
            }
        }

        fprintf(rb_data,"\n");
    }

    gsl_vector_complex_free(v_gsl);
    free(vec_real);
    free(vec_imag);
}

void gsl_vector_sqrt(gsl_vector_complex *out, const gsl_vector_complex *in)
{  

    gsl_complex z1;

    for(int i=0; i < in->size;i++)
    {   
         z1 = gsl_complex_sqrt(gsl_vector_complex_get(in,i));
         gsl_vector_complex_set(out,i,z1);
    }

}

#endif /* gsl_helper_functions.h */
