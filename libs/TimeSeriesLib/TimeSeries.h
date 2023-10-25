#include <vector>
#include <iostream>

class TimeSeries
{
public:
    TimeSeries();
    TimeSeries(size_t num_compartments);
    TimeSeries(double initial_time, std::vector<double> initial_values);
    TimeSeries& operator=(const TimeSeries& time_series);
    std::vector<double> get_value(size_t index);
    //std::vector<double> get_value(double time); //Returns value with possible interpolation
    std::vector<double> get_times();
    double get_time(size_t index);
    bool add_timepoint(double time, std::vector<double> values); //add timepoint to the end
    bool remove_timepoint(size_t time_index);
    size_t get_num_compartments();
    size_t get_num_timepoints();

private:
    size_t num_compartments;
    size_t num_timepoints;
    std::vector<std::vector<double>> time_series_data;
};