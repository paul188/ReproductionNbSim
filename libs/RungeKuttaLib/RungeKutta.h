#include "../TimeSeriesLib/TimeSeries.h"

typedef std::vector<double> (*differential_func)(double t, Eigen::VectorXd x);

/*
Integrator for initial-value problems. Only allows for explicit Runge-Kutta methods to be passed.
*/
class Integrator
{
public:
    Integrator();
    Integrator(bool adaptive_step_size);
    TimeSeries get_result(differential_func right_hand_side, Eigen::VectorXd initial_values, const double t0,
                          const double tmax, const int num_steps);

private:
    bool adaptive_step_size;
};