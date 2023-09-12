#include "RungeKutta.h"

Integrator::Integrator()
{
    adaptive_step_size = false;
}

Integrator::Integrator(bool adaptive_step_size)
{
    this->adaptive_step_size = adaptive_step_size;
}

TimeSeries Integrator::get_result(differential_func right_hand_side, Eigen::VectorXd initial_values, const double t0,
                                  const double tmax, const int num_steps)
{
    //TO BE DONE
    TimeSeries result(t0, initial_values);
    double step_size        = (tmax - t0) / num_steps;
    double time             = t0;
    int step                = 0;
    size_t num_compartments = initial_values.size();
    Eigen::VectorXd k1(num_compartments), k2(num_compartments), k3(num_compartments), k4(num_compartments);
    while (time < tmax) {
        if (time + step_size > tmax) {
            step_size = tmax - time;
        }
        for (int i = 0; i < num_compartments; i++) {
            k1[i] = right_hand_side(time, result.get_value(step))[i];
            k2[i] = right_hand_side(time + step_size / 2, result.get_value(step) + (step_size / 2) * k1)[i];
            k3[i] = right_hand_side(time + step_size / 2, result.get_value(step) + (step_size / 2) * k2)[i];
            k4[i] = right_hand_side(time + step_size, result.get_value(step) + (step_size)*k3)[i];
            step++;
        }
        result.add_timepoint(time, (step_size / 6) * (k1 + 2 * k2 + 2 * k3 + k4));
    }
    return result;
}
