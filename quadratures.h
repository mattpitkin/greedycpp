#ifndef quadratures_h
#define quadratures_h


// TODO: there is an fcount in main -- should fix this //
static int fcount_pts_1(const char *infile) {
    //const char *c_str = infile.c_str();
    //FILE *fvec = fopen(c_str, "r");
    std::ifstream fvec(infile);
    // counts the number of lines in the given file not beginning with "#"
    // to get the number of points
    std::string line;
    int points = 0;
    while (std::getline(fvec, line))
    {
        if (line.compare(0, 1, "#") != 0)
            ++points;

    }
    //fclose(fvec);

    return points;
}

/* --- fill array with linear spacing --- */
void Linspace(const int &n, const double &a, const double &b, double *SomeArray)
{
    double factor = (b-a)/(double)(n-1);
    for(int i=0;i<n;i++)
    {
        SomeArray[i] = a + (double)i*factor;
    }
}

void ReimannQuad(const double a,const double b,double *xQuad,double * wQuad,const int quad_points)
{

    Linspace(quad_points, a, b, xQuad);

    for(int i = 0; i < quad_points; i++)
    {
        wQuad[i] = (b-a)/( (double) (quad_points-1) );
    }

}

void DynamicQuad(double *wQuad, const gsl_vector *xQuad, const int quad_points)
{
/* compute quadrature weights as Reimann sum-like rule with nonuniform xQuad */
    double this_f, next_f;

    this_f = gsl_vector_get(xQuad, 0);
    for(int i = 0; i < quad_points-1; i++)
    {
        next_f = gsl_vector_get(xQuad, i+1);
        wQuad[i] = next_f - this_f;
        this_f = next_f;
    }
    // set weight of last point to 0
    wQuad[quad_points-1] = 0.;
}

void MakeQuadratureRule(gsl_vector_complex *wQuad_c, gsl_vector *xQuad_c, const double a, const double b, const int quad_points,const int quad_type)
{
    double *wQuad_tmp, *xQuad_tmp;
    wQuad_tmp = new double[quad_points];
    xQuad_tmp = new double[quad_points];

    // -- Quadrature rule for inner product between rows -- //
    if(quad_type == 0) {
        gauleg(a,b,xQuad_tmp,wQuad_tmp,quad_points); // returns grid on [-1,1] from NR3
    }
    else if(quad_type == 1){
        ReimannQuad(a,b,xQuad_tmp,wQuad_tmp,quad_points);
    }
    else if(quad_type == 2){
        DynamicQuad(wQuad_tmp, xQuad_c,quad_points);
    }
    else{
        fprintf(stderr,"quadrature rule not coded\n");
        delete[] wQuad_tmp;
        delete[] xQuad_tmp;
        exit(1);
    }

    /* --- Make quad weights and points of type gsl_complex_vector --- */
    gsl_complex zM1;
    for (int i = 0; i < quad_points; i++)
    {
        GSL_SET_COMPLEX(&zM1,wQuad_tmp[i],0.0);
        gsl_vector_complex_set(wQuad_c,i,zM1);
    }

    if(quad_type != 2)
    {
        for (int i = 0; i < quad_points; i++){
            gsl_vector_set(xQuad_c,i,xQuad_tmp[i]);
        }
    }

    delete[] wQuad_tmp;
    delete[] xQuad_tmp;
}

void MakeWeightedInnerProduct(gsl_vector_complex *wQuad, FILE *weightf)
{
/* If the inner product includes a weight W(x), incorporate this by modifying the quadrature weights accordingly 
   Note: this routine is only called if 'weighted = true' from configuration (input) file. */

    // TODO: error check that input weight and wQuad are of equal length //

    int gsl_status;
    gsl_complex z;
    double a;
    gsl_vector *weight;

    // load the weight //
    fprintf(stdout, "Loading ASD (weight) ...\n");
    weight = gsl_vector_alloc(wQuad->size);
    gsl_status = gsl_vector_fscanf(weightf, weight);
    if( gsl_status == GSL_EFAILED )
    {
        fprintf(stderr, "Error reading ASD file\n");
        exit(1);
    }

    // divide by the ASD //
    for(int jj = 0; jj < wQuad->size; jj++)
    {
        z = gsl_vector_complex_get(wQuad, jj);
        a = gsl_vector_get(weight, jj);

        // TODO: should this be squared for ASD?
        z = gsl_complex_div_real(z, a);

        gsl_vector_complex_set(wQuad, jj, z);
    }

    gsl_vector_free(weight);

}

void SetupQuadratureRule(gsl_vector_complex **wQuad,gsl_vector **xQuad,const int quad_type,const bool weighted_inner,const char *cfg_file)
{
    // wQuad and xQuad are pointers to pointers, which allows memory allocation here to be passed back to main //


    double x_min;              // lower value x_min (physical domain)
    double x_max;              // upper value x_max (physical domain)
    int quad_points;
    const char *fvec_file_name;
    const char *weight_file_name;
    int gsl_status;
    bool cfg_status;

    gsl_vector_complex *wQuad_tmp;
    gsl_vector *xQuad_tmp;


    libconfig::Config cfg; // TODO: opening a new file! better way?
    try{
      cfg.readFile(cfg_file);
    }
    catch(const libconfig::FileIOException &fioex){
      std::cerr << "I/O error while reading file." << std::endl;
      exit(1);
    }
    catch(const libconfig::ParseException &pex){
      std::cerr << "Parse error " << std::endl;
      exit(1);
    }


    if (quad_type == 2) // Dynamic, need to look up frequency vector
    {
        cfg_status = cfg.lookupValue("frequency_vector_file", fvec_file_name); // specify the frequency vector to use from a file
        if (!cfg_status){
            fprintf(stderr, "frequency_vector_file not found in config file\n");
            exit(1);
        }
        quad_points = fcount_pts_1(fvec_file_name);
    }
    else
    {
        x_min = cfg.lookup("x_min");
        x_max = cfg.lookup("x_max");
        quad_points = cfg.lookup("quad_points");
    }

    if (weighted_inner){
        cfg_status = cfg.lookupValue("weight_file", weight_file_name);
        if (!cfg_status){
            fprintf(stderr, "weight_file not found in config file\n");
            exit(1);
        }
    }

    // -- allocate memory --//
    wQuad_tmp = gsl_vector_complex_alloc(quad_points);
    xQuad_tmp = gsl_vector_alloc(quad_points);

    // load frequency vector for dynamic frequency 
    if(quad_type == 2)
    {
        FILE *fvecf = fopen(fvec_file_name, "r");
        gsl_status = gsl_vector_fscanf(fvecf, xQuad_tmp);
        fclose(fvecf);
        if( gsl_status == GSL_EFAILED ){
            fprintf(stderr, "Error reading frequency vector from %s\n", fvec_file_name);
            exit(1);
        }

    }

    /* -- all procs should have a copy of the quadrature rule -- */
    // TODO: on nodes can this be shared? Whats the best way to impliment it? 
    MakeQuadratureRule(wQuad_tmp,xQuad_tmp,x_min,x_max,quad_points,quad_type);

    // Role inner product weight into wQuad //
    if(weighted_inner){
        FILE *weightf = fopen(weight_file_name, "r");
        MakeWeightedInnerProduct(wQuad_tmp, weightf);
        fclose(weightf);
    }

    *wQuad = wQuad_tmp;
    *xQuad = xQuad_tmp;

}


#endif // quadratures.h //
