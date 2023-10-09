#include "TimeSeries.h"
#include <iostream>
#include <vcruntime.h>

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

TimeSeries::TimeSeries(double initial_time, std::vector<double> initial_values)
    : num_compartments(initial_values.size())
{
    if (initial_time < 0) {
        std::cerr << "initial time cannot be negative !" << std::endl;
    }

    if (initial_values.size() == 0) {
        std::cerr << "Initial values vector must not be empty !" << std::endl;
    }

    std::vector<double> temp1;
    temp1.push_back(initial_time);
    for (size_t i = 0; i < num_compartments; i++) {
        temp1.push_back(initial_values[i]);
    }
    time_series_data.push_back(temp1);
    num_timepoints = 1;
}

TimeSeries& TimeSeries::operator=(const TimeSeries& time_series)
{
    this->num_timepoints   = time_series.num_timepoints;
    this->num_compartments = time_series.num_compartments;
    this->time_series_data.resize(num_timepoints);
    for (size_t i = 0; i < num_timepoints; i++) {
        this->time_series_data[i].resize(num_compartments + 1);
        for (size_t j = 0; j < num_compartments + 1; j++) {
            this->time_series_data[i][j] = time_series.time_series_data[i][j];
        }
    }
    return *this;
}

bool TimeSeries::add_timepoint(double time, std::vector<double> values)
{
    if (values.size() != num_compartments) {
        std::cerr << "cannot add timepoint of wrong format !" << std::endl;
    }
    std::vector<double> temp1(0);
    temp1.push_back(time);
    for (size_t i = 0; i < num_compartments; i++) {
        temp1.push_back(values[i]);
    }
    time_series_data.push_back(temp1);
    num_timepoints += 1;
    return true;
}

bool TimeSeries::remove_timepoint(size_t time_index)
{
    num_timepoints -= 1;
    time_series_data.erase(time_series_data.begin() + time_index);
    return true;
}

std::vector<double> TimeSeries::get_value(size_t index)
{
    std::vector<double> temp = time_series_data[index];
    //Now remove the first element (time) from the vector
    temp.erase(temp.begin());
    return temp;
}

std::vector<double> TimeSeries::get_times()
{
    if (num_timepoints == 0) {
        std::cerr << "No times to return, TimeSeries is empty" << std::endl;
    }
    std::vector<double> temp;
    for (int i = 0; i < num_timepoints; i++) {
        temp.push_back(time_series_data[i][0]);
    }
    return temp;
}

double linear_interpolate(double t, double t1, double t2, double y1, double y2)
{
    return y1 + (t - t1) * (y2 - y1) / (t2 - t1);
}

double TimeSeries::get_time(size_t index)
{
    return time_series_data[index][0];
}

std::vector<double> TimeSeries::get_value(double time)
{
    if (time < 0) {
        std::cerr << "Cannot get value at time < 0";
    }
    if (time == time_series_data[0][0]) {
        return time_series_data[0];
    }
    //Find the first later time
    std::vector<double> temp;
    //Find the first later time
    auto time_index_late =
        std::distance(get_times().begin(), std::lower_bound(get_times().begin(), get_times().end(), time));
    for (int i = 1; i < num_compartments + 1; i++) {
        double y1 = time_series_data[time_index_late - 1][i];
        double y2 = time_series_data[time_index_late][i];
        temp.push_back(linear_interpolate(time, get_time(time_index_late - 1), get_time(time_index_late), y1, y2));
    }
    return temp;
}

size_t TimeSeries::get_num_compartments()
{
    return num_compartments;
}

size_t TimeSeries::get_num_timepoints()
{
    return num_timepoints;
}
