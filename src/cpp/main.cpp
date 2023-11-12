#include "../../libs/RungeKuttaLib/RungeKutta.h"
#include <fstream>
#include <direct.h> //Filesystem handling uses windows

inline std::string get_path_to_parent_folder()
{
    wchar_t* buff;
    buff = _wgetcwd(NULL, 0);
    std::wstring dir(buff);
    std::string current_working_dir(dir.begin(), dir.end());
    std::string path_to_Secihurd;
    size_t pos = 0;
    while ((pos = current_working_dir.find_last_of("\\")) != std::string::npos) {
        std::string folder = current_working_dir.substr(pos + 1, current_working_dir.size() - 1);
        if (folder != "ReproductionNbSim") {
            current_working_dir.erase(pos, current_working_dir.size() - 1);
        }
        else {
            break;
        }
    }
    return current_working_dir;
}

inline std::string get_filepath(std::string file, unsigned int region = 0)
{
    std::string path_to_Secihurd = get_path_to_parent_folder();
    if (file == "total_populations") {
        return path_to_Secihurd + "/data/initial_populations/total_populations.txt";
    }
    if (file == "initial_recovered") {
        return path_to_Secihurd + "/data/initial_populations/initial_recovered.txt";
    }
    if (file == "initial_infected") {
        return path_to_Secihurd + "/data/initial_populations/initial_infected.txt";
    }
    if (file == "initial_icu") {
        return path_to_Secihurd + "/data/initial_populations/initial_icu.txt";
    }
    if (file == "initial_hospitalized") {
        return path_to_Secihurd + "/data/initial_populations/initial_hospitalized.txt";
    }
    if (file == "initial_dead") {
        return path_to_Secihurd + "/data/initial_populations/initial_dead.txt";
    }
    if (file == "write_sim_data") {
        return path_to_Secihurd + "/data/data_simulation_runs/" + std::to_string(region) + ".txt";
    }
    return path_to_Secihurd;
}

inline double readline_to_double(std::string file, unsigned int line_number)
{
    std::fstream stream;
    stream.open(get_filepath(file));
    if (stream.fail()) {
        std::cout << "File failed to open" << std::endl;
    }

    int current_line = 0;
    std::string line = "";
    while (!stream.eof()) {
        current_line++;
        std::getline(stream, line);
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

    //Read initial values from Italian data files
    //First find out if program is called directly from exe or python fitting

    //Read infected population
    double initial_infected = readline_to_double("initial_infected", region);
    temp.push_back(initial_infected);

    //Read hospitalized population
    double initial_hospitalized = readline_to_double("initial_hospitalized", region);
    temp.push_back(initial_hospitalized);

    //Read ICU data
    double initial_icu = readline_to_double("initial_icu", region);
    temp.push_back(initial_icu);

    //Read Recovered data
    double initial_recovered = readline_to_double("initial_recovered", region);
    temp.push_back(initial_recovered);

    double initial_recovered_hospital = readline_to_double("initial_recovered", region);
    temp.push_back(initial_recovered_hospital);

    //Read Dead data
    double initial_dead = readline_to_double("initial_dead", region);
    temp.push_back(initial_dead);
    return temp;
}
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
            std::cout << "argument count: " << argc << std::endl;
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
        /*std::ofstream stream;
        std::string filename = get_filepath("write_sim_data", region);
        stream.open(filename);

        for (int j = 0; j < new_series.get_num_compartments(); j++) {
            for (int i = 0; i < new_series.get_num_timepoints(); i++) {
                stream << new_series.get_value((double)i)[(int)j] << std::endl;
            }
        }*/
    }
    return 0;
}
