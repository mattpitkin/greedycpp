#ifndef parameters_set_hpp
#define parameters_set_hpp

#include <string.h>
#include <libconfig.h++>

/* Class contains all parameters set in configuration file */

class Parameters {
  public:

    Parameters(char **,bool);
    ~Parameters();

    inline bool load_from_file() const { return load_from_file_; }
    inline bool export_tspace() const { return export_tspace_; }
    inline bool weighted_inner() const { return weighted_inner_; }
    inline bool ts_file_exists() const { return ts_file_exists_; }
    inline int param_dim() const { return param_dim_; }
    inline int seed() const { return seed_; }
    inline int max_RB() const { return max_RB_; }
    inline int ts_size() const { return ts_size_; }
    inline int quad_type() const { return quad_type_; }
    inline int quad_points() const { return quad_points_; }
    inline double tol() const { return tol_; }
    inline double x_min() const { return x_min_; }
    inline double x_max() const { return x_max_; }
    inline const int * params_num() const { return params_num_; }
    inline const double* param_scale() const { return param_scale_; }
    inline const double * params_low() const { return params_low_; }
    inline const double * params_high() const { return params_high_; }
    inline std::string output_dir() const { return output_dir_; }
    inline std::string model_name() const { return model_name_; }
    inline std::string ts_file_name() const { return ts_file_name_; }
    inline std::string output_data_format() const { return output_data_format_;}
    inline std::string quad_nodes_file() const { return quad_nodes_file_; }
    inline std::string num_weight_file() const { return num_weight_file_; }
    inline std::string quad_weight_file() const { return quad_weight_file_; }

  private:
    std::string config_file_name_; // name of cfg file. command line argument


    //**** settings which MUST be set ****//

    // load_from_file_: Load training points from file (If true, provide
    // ts_file_name. If false, set parms_low, params_high and params_num)
    bool load_from_file_;

    // quad_type_: Gaussian quadrature (0), Riemann sum (1), user-defined
    // pointset (2). If 2, provide quad_nodes_file_, num_weight_file_
    int quad_type_;

    // weighted_inner_: Whether to include a weight W(x), such that 
    //  \int W(x) f(x) g(x)  (if true, provide weight_file_name)
    bool weighted_inner_;

    // seed_: Greedy algorithm seed
    int seed_;                     // greedy algorithm seed

    // tol_: Greedy algorithm tolerance ( \| app \|^2)
    double tol_;

    // max_RB_: Estimated number of RB (reasonable upper bound)
    int max_RB_;

    // param_dim_: Number of parametric dimensions
    int param_dim_;

    // export_tspace_: If true, mpi_size must equal 1
    bool export_tspace_;

    // TODO: should default to 1 
    // param_scale_: scale each params[j][i] so that input to the model is 
    // is the value (param_scale[i] * params[j][i])
    double* param_scale_;

    // Folder to put all output files
    std::string output_dir_;

    // Format of output files (text or gsl binary supported)
    std::string output_data_format_;

    // Name of model used to select the appropriate model in FillTrainingSet
    std::string model_name_;


    //**** run settings - MAY need to be set (quadrature info) ****//

    // x_min_: Lower value x_min (physical domain). Set if quad_type != 2
    double x_min_;

    // x_min_: Upper value x_max (physical domain). Set if quad_type != 2
    double x_max_;

    // quad_points_: total number of quadrature points. Set if quad_type != 2
    int quad_points_;

    // quad_nodes_file_: File name for vector of quadrature points.
    // Set if quad_type = 2
    std::string quad_nodes_file_;

    // num_weight_file_: File name of numerical (quadrature) weights.
    // Set if quad_type = 2
    std::string num_weight_file_;

    // quad_weight_file_: File name of weights W(x)
    // Set if weighted_inner = true
    std::string quad_weight_file_;


    //**** run settings - MAY need to be set (Model info) ****//

    // Lower/upper interval of each parameter
    // Set if load_from_file = false
    double *params_low_, *params_high_;

    // Number of samplings in the interval [param_low,param_high]
    // Set if load_from_file = false.
    int *params_num_;

    // Set if load_from_file = true
    std::string ts_file_name_;


    //**** Internal data not set in the configuration file ****//

    // True if ts file loaded, otherwise false. Useful to avoid termination
    // when file not found. Useful when a basis validation study might not
    // make this available
    bool ts_file_exists_;

    int ts_size_;
};

/// Output operator for a Parameters
std::ostream& operator<<(std::ostream& os, const Parameters& p);

#endif
