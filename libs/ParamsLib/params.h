#include <string>
#include <vector>
#include <iostream>

enum class param_key
{
    R1    = 0,
    R2    = 1,
    R3    = 2,
    R4    = 3,
    R5    = 4,
    R6    = 5,
    R7    = 6,
    R8    = 7,
    R9    = 8,
    alpha = 9,
    beta  = 10,
    rho   = 11,
    theta = 12,
    delta = 13,
    d     = 14,
    count = 15
};

class Parameters
{
public:
    Parameters();
    Parameters(std::vector<double> params);
    void set_parameter(param_key key, double param_val);
    double get_parameter(param_key key);
    bool check_constraints();
    void update_params();

private:
    void standard_init();
    std::vector<double> params;
};