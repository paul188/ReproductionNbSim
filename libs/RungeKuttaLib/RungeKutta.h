#include "../TimeSeriesLib/TimeSeries.h"

typedef std::vector<double> (*differential_func)(double t, std::vector<double> x, std::vector<double> params);

class Integrator
{
public:
    Integrator();
    TimeSeries get_result(differential_func right_hand_side, std::vector<double> initial_values,
                          std::vector<double> params, const double t0, const double tmax, const int num_steps);
    TimeSeries get_result_adaptive(differential_func right_hand_side, std::vector<double> initial_values,
                                   std::vector<double> params, const double t0, const double tmax, const int num_steps);
    std::vector<double> cash_karp_step(differential_func right_hand_side, std::vector<double> y, const double t,
                                       double h, std::vector<double>& yerr, std::vector<double> params);
};