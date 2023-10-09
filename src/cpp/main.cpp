#include <corecrt_math.h>
#include <iostream>
#include "../../libs/RungeKuttaLib/RungeKutta.h"
#include <fstream>
#include <string>

//Set the parameters
/*double R1    = 0.587;
double R2    = 1 / (5.2);
double R3    = 1 / (4.2);
double R4    = 1.0 / 14;
double R5    = 1.0 / 16;
double R6    = 1 / (2.5);
double R7    = 1 / (3.5);
double R8    = 1.0 / 16;
double R9    = 1 / (1 / R3 + (0.5 * 1 / R4));
double delta = 1 / (6.92);
double alpha = 0.01;
double beta  = 0.05;
double rho   = 0.35;
double theta = 0.15;
double d     = 1 / (6.5);
double N0    = 1000;*/

//Set the differential equations
std::vector<double> right_hand_side(
    double t, std::vector<double> x,
    std::vector<double>
        params) //params contains the parameters to be fitted: R1,R3,R4,R5,R6,R7,R8 and alpha, beta, rho, theta, delta, d
{
    double R1    = params[0];
    double R2    = 5.2 - 1 / params[1];
    double R3    = params[1];
    double R4    = params[2];
    double R5    = params[3];
    double R6    = params[4];
    double R7    = params[5];
    double R8    = params[6];
    double R9    = 1 / (1 / params[1] + 0.5 * 1 / params[2]);
    double alpha = params[7];
    double beta  = params[8];
    double rho   = params[9];
    double theta = params[10];
    double delta = params[11];
    double d     = params[12];
    double N0    = 0;
    for (int i = 0; i < x.size(); i++) {
        N0 += x[i];
    }
    std::vector<double> result(8);
    result[0] = -R1 * (x[2] + beta * x[3]) / N0 * x[0];
    result[1] = R1 * (x[2] + beta * x[3]) / N0 * x[0] - R2 * x[1];
    result[2] = R2 * x[1] - ((1 - alpha) * R3 + alpha * R9) * x[2];
    result[3] = (1 - alpha) * R3 * x[2] - ((1 - rho) * R4 + rho * R6) * x[3];
    result[4] = rho * R6 * x[3] - ((1 - theta) * R5 + theta * R7) * x[4];
    result[5] = theta * R7 * x[4] - ((1 - delta) * R8 + delta * d) * x[5];
    result[6] = alpha * R9 * x[2] + (1 - rho) * R4 * x[3] + (1 - theta) * R5 * x[4] + (1 - delta) * R8 * x[5];
    result[7] = delta * R8 * x[5];

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

std::vector<double> read_values(unsigned int region, double initial_exposed, double initial_carrier)
{
    std::vector<double> temp;

    std::fstream instream;

    //Read total population
    double total_population =
        readline_to_double("C:/Users/paul1/OneDrive/Desktop/epidemiology/coding/Secihurd_Model/src/cpp/"
                           "data_for_simulation/initial_populations/total_populations.txt",
                           region);
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

    //Read Dead data
    double initial_dead = readline_to_double("C:/Users/paul1/OneDrive/Desktop/epidemiology/coding/Secihurd_Model/src/"
                                             "cpp/data_for_simulation/initial_populations/initial_dead.txt",
                                             region);
    temp.push_back(initial_dead);
    std::cout << "Data read" << std::endl;
    std::cout << temp.size() << std::endl;
    return temp;
}

int main(int argc, char* argv[])
{

    unsigned int region = atof(argv[1]);

    std::vector<double> initial_values;

    initial_values = read_values(region, atof(argv[2]), atof(argv[3])); //In arguments 1,2 we have exposed and carriers
    std::vector<double> params = {atof(argv[4]),  atof(argv[5]),  atof(argv[6]),  atof(argv[7]),  atof(argv[8]),
                                  atof(argv[9]),  atof(argv[10]), atof(argv[11]), atof(argv[12]), atof(argv[13]),
                                  atof(argv[14]), atof(argv[15]), atof(argv[16])};
    TimeSeries new_series(8);
    Integrator Integrator;
    new_series = Integrator.get_result(right_hand_side, initial_values, params, 0, 26, 102);
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
        stream << std::endl;
    }
    for (int j = 0; j < new_series.get_num_compartments(); j++) {
        for (int i = 0; i < new_series.get_num_timepoints(); i++) {
            std::cout << new_series.get_value((double)i)[j] << std::endl;
        }
    }
    return 0;
}