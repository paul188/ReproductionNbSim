#include "../../libs/RungeKuttaLib/RungeKutta.h"
#include <fstream>

//Set the differential equations
std::vector<double> right_hand_side(double t, std::vector<double> x, Parameters params)
{
    double R1    = params.get_parameter(param_key::R1);
    double R2    = params.get_parameter(param_key::R2);
    double R3    = params.get_parameter(param_key::R3);
    double R4    = params.get_parameter(param_key::R4);
    double R5    = params.get_parameter(param_key::R5);
    double R6    = params.get_parameter(param_key::R6);
    double R7    = params.get_parameter(param_key::R7);
    double R8    = params.get_parameter(param_key::R8);
    double R9    = params.get_parameter(param_key::R9);
    double alpha = params.get_parameter(param_key::alpha);
    double beta  = params.get_parameter(param_key::beta);
    double rho   = params.get_parameter(param_key::rho);
    double theta = params.get_parameter(param_key::theta);
    double delta = params.get_parameter(param_key::delta);
    double d     = params.get_parameter(param_key::d);
    double N0    = 0;
    for (int i = 0; i < x.size(); i++) {
        N0 += x[i] - x[7]; //don't want to count recovered from severe disease twice
    }
    std::vector<double> result(9);
    result[0] = -R1 * (x[2] + beta * x[3]) / N0 * x[0];
    result[1] = R1 * (x[2] + beta * x[3]) / N0 * x[0] - R2 * x[1];
    result[2] = R2 * x[1] - ((1 - alpha) * R3 + alpha * R9) * x[2];
    result[3] = (1 - alpha) * R3 * x[2] - ((1 - rho) * R4 + rho * R6) * x[3];
    result[4] = rho * R6 * x[3] - ((1 - theta) * R5 + theta * R7) * x[4];
    result[5] = theta * R7 * x[4] - ((1 - delta) * R8 + delta * d) * x[5];
    result[6] = alpha * R9 * x[2] + (1 - rho) * R4 * x[3] + (1 - theta) * R5 * x[4] + (1 - delta) * R8 * x[5];
    result[7] = (1 - theta) * R5 * x[4] + (1 - delta) * R8 * x[5]; //extra compartment for fitting to Italy data
    result[8] = delta * R8 * x[5];
    return result;
}

double readline_to_double(std::string filepath, unsigned int line_number)
{
    std::fstream stream;
    stream.open(filepath);
    if (stream.fail()) {
        std::cout << "File failed to open" << std::endl;
    }

    int current_line = 0;
    std::string line = "";
    while (!stream.eof()) {
        current_line++;
        getline(stream, line);
        if (current_line == line_number) {
            break;
        }
    }
    if (current_line < line_number) {
        std::cout << "line not found, file does not contain as many lines" << std::endl;
    }
    return std::stod(line);
}

std::vector<double> read_values(unsigned int region, double initial_exposed, double initial_carrier, double total_pop)
{
    std::vector<double> temp;

    std::fstream instream;

    //Read total population
    double total_population = total_pop;

    temp.push_back(total_population);

    temp.push_back(initial_exposed);
    temp.push_back(initial_carrier);

    //Read infected population
    double initial_infected = readline_to_double("C:/Users/paul1/OneDrive/Desktop/epidemiology/coding/Secihurd_Model/"
                                                 "src/cpp/data_for_simulation/initial_populations/initial_infected.txt",
                                                 region);
    temp.push_back(initial_infected);

    //Read hospitalized population
    double initial_hospitalized =
        readline_to_double("C:/Users/paul1/OneDrive/Desktop/epidemiology/coding/Secihurd_Model/src/cpp/"
                           "data_for_simulation/initial_populations/initial_hospitalized.txt",
                           region);
    temp.push_back(initial_hospitalized);

    //Read ICU data
    double initial_icu = readline_to_double("C:/Users/paul1/OneDrive/Desktop/epidemiology/coding/Secihurd_Model/src/"
                                            "cpp/data_for_simulation/initial_populations/initial_icu.txt",
                                            region);
    temp.push_back(initial_icu);

    //Read Recovered data
    double initial_recovered =
        readline_to_double("C:/Users/paul1/OneDrive/Desktop/epidemiology/coding/Secihurd_Model/src/cpp/"
                           "data_for_simulation/initial_populations/initial_recovered.txt",
                           region);
    temp.push_back(initial_recovered);

    double initial_recovered_hospital =
        readline_to_double("C:/Users/paul1/OneDrive/Desktop/epidemiology/coding/Secihurd_Model/src/cpp/"
                           "data_for_simulation/initial_populations/initial_recovered.txt",
                           region);
    temp.push_back(initial_recovered_hospital);

    //Read Dead data
    double initial_dead = readline_to_double("C:/Users/paul1/OneDrive/Desktop/epidemiology/coding/Secihurd_Model/src/"
                                             "cpp/data_for_simulation/initial_populations/initial_dead.txt",
                                             region);
    temp.push_back(initial_dead);
    return temp;
}

Parameters create_input_params(char* argv[])
{
    Parameters params;
    params.set_parameter(param_key::R1, atof(argv[4]));
    params.set_parameter(param_key::R2, 1 / (5.2 - 1 / (atof(argv[5]))));
    params.set_parameter(param_key::R3, atof(argv[5]));
    params.set_parameter(param_key::R4, atof(argv[6]));
    params.set_parameter(param_key::R5, atof(argv[7]));
    params.set_parameter(param_key::R6, atof(argv[8]));
    params.set_parameter(param_key::R7, atof(argv[9]));
    params.set_parameter(param_key::R8, atof(argv[10]));
    params.set_parameter(param_key::R9, 1 / (1 / atof(argv[5]) + 0.5 * 1 / (atof(argv[6]))));
    params.set_parameter(param_key::alpha, atof(argv[11]));
    params.set_parameter(param_key::beta, atof(argv[12]));
    params.set_parameter(param_key::rho, atof(argv[13]));
    params.set_parameter(param_key::theta, atof(argv[14]));
    params.set_parameter(param_key::delta, atof(argv[15]));
    params.set_parameter(param_key::d, atof(argv[16]));
    return params;
}

int main(int argc, char* argv[])
{
    if (argc != 18) {
        std::cout << "Wrong number of input arguments" << std::endl;
        for (int i = 0; i < argc; i++) {
            std::cout << argv[i] << std::endl;
        }
    }
    else {
        unsigned int region = atof(argv[1]);
        std::vector<double> initial_values;
        initial_values    = read_values(region, atof(argv[2]), atof(argv[3]), atof(argv[17]));
        Parameters params = create_input_params(argv);
        params.check_constraints();
        TimeSeries new_series(9);
        Integrator integrator;
        integrator.set_parameters(params);
        integrator.set_initial_vals(initial_values);
        new_series = integrator.get_result(right_hand_side, 0, 26, 104);

        for (int j = 0; j < new_series.get_num_compartments(); j++) {
            for (int i = 0; i < new_series.get_num_timepoints(); i++) {
                std::cout << new_series.get_value((double)i)[j] << std::endl;
            }
        }
        std::ofstream stream;
        std::string filename = "C:/Users/paul1/OneDrive/Desktop/epidemiology/coding/Secihurd_Model/src/python/data" +
                               std::to_string(region) + ".txt";
        stream.open(filename);
        if (stream.fail()) {
            stream.close();
            stream.open("../" + filename);
        }

        for (int j = 0; j < new_series.get_num_compartments(); j++) {
            for (int i = 0; i < new_series.get_num_timepoints(); i++) {
                stream << new_series.get_value((double)i)[(int)j] << std::endl;
            }
        }
    }
    return 0;
}