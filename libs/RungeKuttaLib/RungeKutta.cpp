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

std::vector<double> Integrator::cash_karp_step(differential_func right_hand_side, std::vector<double> y, const double t,
                                               double h, std::vector<double>& yerr, std::vector<double> params)
{
    static double a2 = 0.2, a3 = 0.3, a4 = 0.6, a5 = 1.0, a6 = 0.875, b21 = 0.2, b31 = 3.0 / 40.0, b32 = 9.0 / 40.0,
                  b41 = 0.3, b42 = -0.9, b43 = 1.2, b51 = -11.0 / 54.0, b52 = 2.5, b53 = -70.0 / 27.0,
                  b54 = 35.0 / 27.0, b61 = 1631.0 / 55296.0, b62 = 175.0 / 512.0, b63 = 575.0 / 13824.0,
                  b64 = 44275.0 / 110592.0, b65 = 253.0 / 4096.0, c1 = 37.0 / 378.0, c3 = 250.0 / 621.0,
                  c4 = 125.0 / 594.0, c6 = 512.0 / 1771.0, dc5 = -277.00 / 14336.0;
    double dc1 = c1 - 2825.0 / 27648.0, dc3 = c3 - 18575.0 / 48384.0, dc4 = c4 - 13525.0 / 55296.0, dc6 = c6 - 0.25;
    std::vector<double> ak2, ak3, ak4, ak5, ak6, ytemp, yout;
    ak2.resize(y.size());
    ak3.resize(y.size());
    ak4.resize(y.size());
    ak5.resize(y.size());
    ak6.resize(y.size());
    ytemp.resize(y.size());

    std::vector<double> dydt = right_hand_side(t, y, params);

    ytemp = y + b21 * h * dydt; //First step

    ak2   = right_hand_side(t + a2 * h, ytemp, params); //Second step
    ytemp = y + h * (b31 * dydt + b32 * ak2);

    ak3   = right_hand_side(t + a3 * h, ytemp, params); //Third step
    ytemp = y + h * (b41 * dydt + b42 * ak2 + b43 * ak3);

    ak4   = right_hand_side(t + a4 * h, ytemp, params); //Fourth step
    ytemp = h * (b51 * dydt + b52 * ak2 + b53 * ak3 + b54 * ak4);

    ak5   = right_hand_side(t + a5 * h, ytemp, params); //Fifth step
    ytemp = y + h * (b61 * dydt + b62 * ak2 + b63 * ak3 + b64 * ak4 + b65 * ak5);

    ak6   = right_hand_side(t + a6 * h, ytemp, params);
    ytemp = y + h * (c1 * dydt + c3 * ak3 + c4 * ak4 + c6 * ak6);

    yerr = h * (dc1 * dydt + dc3 * ak3 + dc4 * ak4 + dc5 * ak5 + dc6 * ak6);
    return ytemp;
}

TimeSeries Integrator::get_result(differential_func right_hand_side, std::vector<double> initial_values,
                                  std::vector<double> params, const double t0, const double tmax, const int num_steps)
{
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
        time += step_size;
        for (int i = 0; i < num_compartments; i++) {
            k1[i] = right_hand_side(time, result.get_value(step - 1), params)[i];
            k2[i] = right_hand_side(time + step_size / 2, result.get_value(step - 1) + (step_size / 2) * k1, params)[i];
            k3[i] = right_hand_side(time + step_size / 2, result.get_value(step - 1) + (step_size / 2) * k2, params)[i];
            k4[i] = right_hand_side(time + step_size, result.get_value(step - 1) + (step_size)*k3, params)[i];
        }
        result.add_timepoint(time, result.get_value(step - 1) + (step_size / 6) * (k1 + (k2 * 2.0) + (k3 * 2.0) + k4));
        step++;
    }
    return result;
}

TimeSeries Integrator::get_result_adaptive(differential_func right_hand_side, std::vector<double> initial_values,
                                           std::vector<double> params, const double t0, const double tmax,
                                           const int num_steps)
{
    TimeSeries temp_result(t0, initial_values);
    TimeSeries interpolated_result;

    double eps = 1e-10; //Local desired truncation error
    double N   = 0; //Total population
    for (int i = 0; i < initial_values.size(); i++) {
        N += initial_values[i];
    }
    double step_size     = (int)std::floor((t0 - tmax) / num_steps);
    double min_step_size = 1e-9;

    std::vector<double> err(initial_values.size());

    for (int j = 0; j < initial_values.size(); j++) {
        for (int step_counter = 1; step_counter < 10000; step_counter++) {
        }
    }

    //Linearly interpolate the results from the adaptive Integrator to obtain equidistant timepoints
    for (int i = 0; i < num_steps; i++) {
        interpolated_result.add_timepoint(i, temp_result.get_value((double)i));
    }
    return interpolated_result;
}
