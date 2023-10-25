#include "params.h"

void Parameters::standard_init()
{
    params.resize((int)param_key::count);
    //Set standard initial values from literature
    params[(int)param_key::R1]    = 0.5;
    params[(int)param_key::R3]    = 1 / (4.2);
    params[(int)param_key::R2]    = 1 / (5.2 - 1 / params[(int)param_key::R3]);
    params[(int)param_key::R4]    = 1 / (14.0);
    params[(int)param_key::R5]    = 1 / (16.0);
    params[(int)param_key::R6]    = 1 / (7.0);
    params[(int)param_key::R7]    = 1 / (3.5);
    params[(int)param_key::R8]    = 1 / (16.0);
    params[(int)param_key::R9]    = 1 / (1 / (params[(int)param_key::R3]) + 0.5 * 1 / (params[(int)param_key::R4]));
    params[(int)param_key::alpha] = 0.01;
    params[(int)param_key::beta]  = 1 / 0.05;
    params[(int)param_key::rho]   = 0.10;
    params[(int)param_key::theta] = 0.15;
    params[(int)param_key::delta] = 0.77;
    params[(int)param_key::d]     = 1 / (6.5);
}

bool ::Parameters::check_constraints()
{
    for (int i = 0; i < (int)param_key::count; i++) {
        if (params[i] < 0) {
            std::cout << "At least one parameter is negative" << std::endl;
            return false;
        }

        if (9 <= i && i <= 13) {
            if (params[i] > 1) {
                std::cout << "At least one transition probability is greater than one" << std::endl;
                for (int j = 0; j < (int)param_key::count; j++) {
                    std::cout << j << ": " << params[j] << std::endl;
                }
                return false;
            }
        }
    }

    if (!(params[1] == 1 / (5.2 - 1 / (params[2])))) {
        std::cout << "R2 is not in the right relation with R3" << std::endl;
        return false;
    }

    if (!(params[8] == 1 / (1 / params[2] + (0.5 * 1 / params[3])))) {
        std::cout << "R9 is not in the right relation with R3 and R4" << std::endl;
        return false;
    }
    return true;
}

void Parameters::update_params()
{
    if (!(params[1] == 1 / (5.2 - 1 / (params[2])))) {
        params[1] = 1 / (5.2 - 1 / (params[2]));
    }
    if (!(params[8] == 1 / (1 / params[2] + (0.5 * 1 / params[3])))) {
        params[8] = 1 / (1 / params[2] + (0.5 * 1 / params[3]));
    }
}

void Parameters::set_parameter(param_key key, double param_val)
{
    params[(int)key] = param_val;
}

double Parameters::get_parameter(param_key key)
{
    return params[(int)key];
}

Parameters::Parameters()
{
    standard_init();
}

Parameters::Parameters(std::vector<double> params)
{
    if (!(params.size() == (int)param_key::count)) {
        std::cout << "Parameter size does not match size of desired parameter set, will instead be initialized with "
                     "standard values"
                  << std::endl;
        standard_init();
        return;
    }

    for (size_t i = 0; i < (size_t)param_key::count; i++) {
        this->params[i] = params[i];
    }

    if (!check_constraints()) {
        std::cout
            << "Parameters do not fulfill the mandatory constraints, will instead be initialized with standard values";
        standard_init();
    }
}