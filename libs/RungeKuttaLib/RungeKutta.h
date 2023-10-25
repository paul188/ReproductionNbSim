#include "../TimeSeriesLib/TimeSeries.h"
#include "../ParamsLib/params.h"

typedef std::vector<double> (*differential_func)(double t, std::vector<double> x, Parameters params);

class Integrator
{
public:
    Integrator();
    Integrator(Parameters params, std::vector<double> initial_values);
    void set_parameters(Parameters params);
    void set_initial_vals(std::vector<double> initial_vals);
    TimeSeries get_result(differential_func right_hand_side, const double t0, const double tmax, const int num_steps);

private:
    Parameters params;
    std::vector<double> initial_values;
};