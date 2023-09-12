#include "TimeSeries.h"
#include <iostream>

TimeSeries::TimeSeries()
{
    num_compartments = 0;
    num_timepoints   = 0;
    time_series_data = {};
}

TimeSeries::TimeSeries(size_t num_compartments)
    : num_compartments(num_compartments)
{
    num_timepoints   = 0;
    time_series_data = {};
}

TimeSeries::TimeSeries(double initial_time, Eigen::VectorXd initial_values)
    : num_compartments(initial_values.size())
{
    if (initial_time < 0) {
        std::cerr << "initial time cannot be negative !" << std::endl;
    }

    if (initial_values.size() == 0) {
        std::cerr << "Initial values vector must not be empty !" << std::endl;
    }

    Eigen::VectorXd temp1;
    temp1 << initial_time;
    for (size_t i = 0; i < num_compartments; i++) {
        temp1 << initial_values[i];
    }
    time_series_data.push_back(temp1);
    num_timepoints = 1;
}

bool TimeSeries::add_timepoint(double time, Eigen::VectorXd values)
{
    if (values.size() != num_compartments) {
        std::cerr << "cannot add timepoint of wrong format !" << std::endl;
    }
    Eigen::VectorXd temp1;
    temp1 << values;
    for (size_t i = 0; i < num_compartments; i++) {
        temp1 << values[i];
    }
    time_series_data.push_back(temp1);
    num_timepoints = 1;
    return true;
}

bool TimeSeries::remove_timepoint(size_t time_index)
{
    num_timepoints -= 1;
    time_series_data.erase(time_series_data.begin() + time_index);
    return true;
}

Eigen::VectorXd TimeSeries::get_value(size_t index)
{
    Eigen::VectorXd temp                 = time_series_data[index];
    temp.block(0, 0, temp.size() - 1, 0) = temp.block(1, 0, temp.size(), 0);
    temp.conservativeResize(temp.size() - 1, 1);
    return temp;
}
