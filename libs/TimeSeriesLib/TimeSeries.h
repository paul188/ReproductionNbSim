#include "include/eigen3/Eigen/Dense"
#include <vector>

class TimeSeries
{
public:
    TimeSeries();
    TimeSeries(size_t num_compartments);
    TimeSeries(double initial_time, Eigen::VectorXd initial_values);
    ~TimeSeries();
    Eigen::VectorXd get_value(size_t index);
    bool add_timepoint(double time, Eigen::VectorXd values); //add timepoint to the end
    bool remove_timepoint(size_t time_index);

private:
    size_t num_compartments;
    size_t num_timepoints;
    std::vector<Eigen::VectorXd> time_series_data;
};