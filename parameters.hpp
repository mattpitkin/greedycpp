#ifndef parameters_set_hpp
#define parameters_set_hpp

#include <string.h>
#include <libconfig.h++>

/* Class contains all parameters set in configuration file */

class Parameters {

    public:
        Parameters();
        Parameters(char **);
        ~Parameters();

        inline bool load_from_file() const { return load_from_file_; }
        inline bool export_tspace() const { return export_tspace_; }
        inline bool weighted_inner() const { return weighted_inner_; }
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
        inline std::string output_data_format() const { return output_data_format_; }
        inline std::string quad_nodes_file() const { return quad_nodes_file_; }
        inline std::string num_weight_file() const { return num_weight_file_; }
        inline std::string quad_weight_file() const { return quad_weight_file_; }


    private:

        // passed in from command line //
        std::string config_file_name_; // name of configuration file

        // settings which MUST be set //
        bool load_from_file_;          // load training points from file (if true, provide ts_file_name, if false parms low/high/num)
        int quad_type_;                // 0 = Gaussian quadrature, 1 = Riemann, 2 = user-defined pointset (if 2, provide quad_nodes_file_, num_weight_file_)
        bool weighted_inner_;          // whether to include a weight W(x): \int W(x) f(x) g(x)  (if true, provide weight_file_name)
        int seed_;                     // greedy algorithm seed
        double tol_;                   // greedy algorithm tolerance ( \| app \|^2)
        int max_RB_;                   // estimated number of RB (reasonable upper bound)
        int param_dim_;                // number of paramteric dimensions
        bool export_tspace_;           // if true, mpi_size must equal 1
        double* param_scale_;          // scale each params[j][i] so that model evaluated at param_scale[i] * params[j][i]

        // run settings - MAY need to be set //
        double x_min_;                 // lower value x_min (physical domain) --> needed if quad_type != 2
        double x_max_;                 // upper value x_max (physical domain) --> needed if quad_type != 2
        int quad_points_;              // total number of quadrature points   --> needed if quad_type != 2
        std::string quad_nodes_file_;  // file name for vector of quadrature points  --> needed if quad_type = 2
        std::string num_weight_file_;  // file name of numerical weights             --> needed in quad_type = 2
        std::string quad_weight_file_; // file name of weights --> needed if weighted_inner = true


        // run settings - MAY need to be set //
        double *params_low_, *params_high_; // this is required if load_from_file = false. lower/upper interval of each parameter 
        int *params_num_;                   // this is required if load_from_file = false. Number of samplings in the interval [param_low,param_high]
        std::string model_name_;            // mame of model -- used to select the appropriate model in FillTrainingSet
        std::string output_dir_;            // folder to put all output files
        std::string output_data_format_;    // format of output files (text or gsl binary supported)
        std::string ts_file_name_;          // this is required if load_from_file = true
        int ts_size_;

};

/// Output operator for a Parameters
std::ostream& operator<<(std::ostream& os, const Parameters& p);

#endif
