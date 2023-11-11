#include "RungeKutta.h"
#include <cmath>

std::vector<double> operator+(std::vector<double> a, std::vector<double> b)
{
    std::vector<double> temp;
    for (int i = 0; i < a.size(); i++) {
        temp.push_back(a[i] + b[i]);
    }
    return temp;
}

std::vector<double> operator*(std::vector<double> b, double a)
{
    std::vector<double> temp(0);
    for (int i = 0; i < b.size(); i++) {
        temp.push_back(a * b[i]);
    }
    return temp;
}

std::vector<double> operator*(double a, std::vector<double> b)
{
    std::vector<double> temp(0);
    for (int i = 0; i < b.size(); i++) {
        temp.push_back(a * b[i]);
    }
    return temp;
}

std::vector<double> operator+=(std::vector<double>& a, std::vector<double> b)
{
    for (int i = 0; i < a.size(); i++) {
        a.push_back(a[i] + b[i]);
    }
    return a;
}

Integrator::Integrator()
{
    Parameters* params = new Parameters();
    this->params       = *params;
}

void Integrator::set_initial_vals(std::vector<double> initial_vals)
{
    initial_values.resize(9);
    if (initial_vals.size() != 9) {
        std::cout << "Initial values cannot be set because sizes don't match, setting to default" << std::endl;
        initial_values[0] = 1000;
        for (int i = 1; i < initial_values.size(); i++) {
            initial_values[i] = 0;
        }
        return;
    }
    for (int i = 0; i < initial_values.size(); i++) {
        initial_values[i] = initial_vals[i];
        if (initial_vals[i] < 0) {
            std::cout << "Initial values cannot be set because value " << i << " is negative, setting value to zero"
                      << std::endl;
            initial_values[i] = 0;
        }
    }
}

void Integrator::set_parameters(Parameters params)
{
    this->params = params;
}

TimeSeries Integrator::get_result(differential_func right_hand_side, const double t0, const double tmax,
                                  const int num_timepoints)
{
    const int num_steps = num_timepoints - 1;
    TimeSeries result(t0, initial_values);
    double step_size        = (tmax - t0) / num_steps;
    double time             = t0;
    size_t step             = 1;
    size_t num_compartments = initial_values.size();
    std::vector<double> k1(num_compartments);
    std::vector<double> k2(num_compartments);
    std::vector<double> k3(num_compartments);
    std::vector<double> k4(num_compartments);
    while (time < tmax) {
        if (time + step_size > tmax) {
            step_size = tmax - time;
        }
        if (step == num_steps) {
            step_size = tmax - time;
        }
        time += step_size;
        for (int i = 0; i < num_compartments; i++) {
            k1[i] = right_hand_side(time, result.get_value(step - 1), this->params)[i];
            k2[i] = right_hand_side(time + step_size / 2, result.get_value(step - 1) + (step_size / 2) * k1, params)[i];
            k3[i] = right_hand_side(time + step_size / 2, result.get_value(step - 1) + (step_size / 2) * k2, params)[i];
            k4[i] = right_hand_side(time + step_size, result.get_value(step - 1) + (step_size)*k3, params)[i];
        }
        result.add_timepoint(time, result.get_value(step - 1) + (step_size / 6) * (k1 + (k2 * 2.0) + (k3 * 2.0) + k4));
        step++;
    }
    return result;
}