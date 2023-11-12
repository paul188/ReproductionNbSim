# This python script fits the model to the different regions
import os
import subprocess
import lmfit
import codecs
import matplotlib.pyplot as plt
import numpy as np
import csv
import datetime
import math
import linecache
from scipy.integrate import simpson

# read in the Italian data in the region in the order: 1. Cumulative Infected, 2. Hospitalized, 3. ICU, 4. Recovered, 5. Dead
# the region is some integer from 0 to 20


def gather_italian_data(region):
    start = datetime.date(2020, 2, 24)
    end = datetime.date(2020, 3, 20)
    res_date = start

    arrayfinal = []

    for i in range(1, 22):
        temparray = []
        while (res_date <= end):
            str_day = str(res_date.day)
            if (res_date.day < 10):
                str_day = '0'+str_day
            str_month = str(res_date.month)
            if (res_date.month < 10):
                str_month = '0'+str_month
            str_year = str(res_date.year)
            str_date = str_year + str_month + str_day
            str_file = './data/dati-regioni/dpc-covid19-ita-regioni-'+str_date+'.csv'
            temparray2 = []
            with open(str_file) as csv_file:
                csv_reader = csv.reader(csv_file, delimiter=',')
                line_count = 0
                for rows in csv_reader:
                    if (line_count == i):
                        temparray2.append(rows[6])  # Hospitalized
                        temparray2.append(rows[7])  # ICU
                        temparray2.append(rows[10])  # Cumulative Infected
                        temparray2.append(rows[13])  # Recovered
                        temparray2.append(rows[14])  # deaths
                    line_count += 1
            temparray.append(temparray2)
            res_date += datetime.timedelta(days=1)
        res_date = start
        arrayfinal.append(temparray)

    index = 0
    res_date = start
    total_data_array = []

    # Store data on cumulative infected
    while (res_date <= end):
        total_data_array.append(
            float(arrayfinal[region][index][2]))
        res_date += datetime.timedelta(days=1)
        index += 1

    res_date = start
    index = 0

    # Store data on Hospitalized
    while (res_date <= end):
        total_data_array.append(float(arrayfinal[region][index][0]))
        res_date += datetime.timedelta(days=1)
        index += 1

    res_date = start
    index = 0

    # Store data on ICU patients
    while (res_date <= end):
        total_data_array.append(
            float(arrayfinal[region][index][1]))
        res_date += datetime.timedelta(days=1)
        index += 1

    res_date = start
    index = 0

    # Store data on recovered
    while (res_date <= end):
        total_data_array.append(
            float(arrayfinal[region][index][3]))
        res_date += datetime.timedelta(days=1)
        index += 1

    res_date = start
    index = 0

    # Store data on deceased
    while (res_date <= end):
        total_data_array.append(
            float(arrayfinal[region][index][4]))
        res_date += datetime.timedelta(days=1)
        index += 1
    return total_data_array


def linearize(x, xleft, xright, yleft, yright):
    return yleft+(x-xleft)*(yright-yleft)/(xright-xleft)

# Returns the integrator results of Cumulative Infected, Hospitalized, ICU, Recovered, Dead weighted by the weight parameter
# Here the region is in the range 1 to 21


def func(x, region, E0, C0, R1, R3, R4, R5, R6, R7, R8, alpha, beta, rho, theta, delta, d, total_pop):  # x should be between 0 and 75
    print("Function called")
    print(E0, C0, R1, R3, R4, R5, R6, R7, R8,
          alpha, beta, rho, theta, delta, d, total_pop)
    # Make x iterable if it is only a float/int
    if not hasattr(x, '__iter__'):
        x = np.array([x])
    result = []
    output = subprocess.check_output([r'build/Debug/model.exe', str(region), str(
        E0), str(C0), str(R1), str(R3), str(R4), str(R5), str(R6), str(R7), str(R8), str(alpha), str(beta), str(rho), str(theta), str(delta), str(d), str(total_pop)])
    output_str = codecs.decode(output)
    output_str_list = output_str.split('\r\n')
    output_str_list.pop()
    dataarray = []
    index = 0
    for word in output_str_list:
        if index in range(312, 416):  # Cumulative Infected
            dataarray.append((float(word)+float(output_str_list[index+104])+float(
                output_str_list[index+208])))  # Infected + Hospitalized + ICU
        if index in range(416, 520):  # Hospitalized
            dataarray.append(float(word))
        if index in range(520, 624):  # ICU
            dataarray.append(float(word))
        if index in range(728, 832):  # Recovered
            dataarray.append(float(word))
        if index in range(832, 936):  # Dead
            dataarray.append(float(word))
        index += 1
    for x0 in x:
        y = linearize(x0*4, math.floor(x0*4), math.floor(x0*4)+1,
                      dataarray[math.floor(x0*4)], dataarray[math.ceil(x0*4)])
        result.append(y)
    return result

# fits the parameters to the Italian data. Returns the results from the model fit
# Takes the region which we want to fit the model to (int from 0 to 20)
# Takes an lmfit.Parameters object as parameter
# Takes an array of doubles from 0 to 1 as weighting factors. The array must have length 5, corresponding to
# Cumulative Infected,...,Dead. The weight 0 should be 1 for example, if a perfect fit to the Cumulative Infected is very important
# and 0 if we don't care about how well the fit to cumulative infected is


def fit_model(region, weight):

    max_pop = float(linecache.getline(
        "./data/initial_populations/total_populations.txt", region+1))/10

    parameters = lmfit.Parameters()
    parameters.add('region', value=region+1)
    parameters.add('E0', value=60, min=1, max=1000.0)
    parameters.add('C0', value=10, min=1, max=1000.0)
    parameters.add('R1', value=0.587, min=0.01, max=1.0)
    parameters.add('R3', value=0.3125, min=0.01, max=1.0)
    parameters.add('R4', value=0.1666, min=0.01, max=1.0)
    parameters.add('R5', value=0.1, min=0.01, max=1.0)
    parameters.add('R6', value=0.2, min=0.01, max=1.0)
    parameters.add('R7', value=0.4, min=0.01, max=1.0)
    parameters.add('R8', value=0.125, min=0.01, max=1.0)
    parameters.add('alpha', value=0.09, min=0.01, max=1.0)
    parameters.add('beta', value=0.25, min=0.01, max=1.0)
    parameters.add('rho', value=0.2, min=0.01, max=1.0)
    parameters.add('theta', value=0.26, min=0.01, max=1.0)
    parameters.add('delta', value=0.77, min=0.01, max=1.0)
    parameters.add('d', value=1/(6.5), min=0.01, max=1.0)
    parameters.add('total_pop', value=3000, min=50, max=max_pop)
    parameters['region'].vary = False

    weights = np.full(26, weight[0])
    weights = np.append(weights, np.full(26, weight[1]))
    weights = np.append(weights, np.full(26, weight[2]))
    weights = np.append(weights, np.full(26, weight[3]))
    weights = np.append(weights, np.full(26, weight[4]))

    total_data_array = gather_italian_data(region)
    my_model = lmfit.Model(func)
    results = my_model.fit(total_data_array, params=parameters,
                           x=np.arange(0, 130, 1), method='leastsq', weights=weights)
    return results


def print_fit(region, weights):
    xData = np.arange(0, 130, 1)
    xData2 = np.arange(0, 26, 1)
    results = fit_model(region, weights)

    tick_labels = ['C.I.',
                   'H', 'ICU', 'Recovered', 'Dead']
    tick_locations = [13, 39, 65, 91, 117]

    result_data = func(xData, results.best_values['region'], results.best_values['E0'], results.best_values['C0'], results.best_values['R1'], results.best_values['R3'], results.best_values['R4'], results.best_values['R5'], results.best_values['R6'],
                       results.best_values['R7'], results.best_values['R8'], results.best_values['alpha'], results.best_values['beta'], results.best_values['rho'], results.best_values['theta'], results.best_values['delta'], results.best_values['d'], results.best_values['total_pop'])

    # Plot the overall fit
    plt.xticks(tick_locations, tick_labels)
    plt.plot(xData, gather_italian_data(region), color="b")
    plt.plot(xData, result_data, color="y")
    plt.savefig(r'./data/results/regional_fitting/region' +
                str(region+1)+'/leastsq_fit/fit.png', bbox_inches="tight")
    plt.clf()

    # Plot the number of Cumulative Infected
    plt.plot(xData2, gather_italian_data(region)[0:26], color="b")
    plt.plot(xData2, result_data[0:26], color="y")
    plt.savefig(r'./data/results/regional_fitting/region' +
                str(region+1)+'/plotCumulativeInfected'+str(region+1)+".png", bbox_inches="tight")
    plt.clf()

    # Plot the number of Hospitalized
    plt.plot(xData2, gather_italian_data(region)[26:52], color="b")
    plt.plot(xData2, result_data[26:52], color="y")
    plt.savefig(r'./data/results/regional_fitting/region' +
                str(region+1)+'/plotHospitalized'+str(region+1)+".png", bbox_inches="tight")
    plt.clf()

    # Plot the number of ICU patients
    plt.plot(xData2, gather_italian_data(region)[52:78], color="b")
    plt.plot(xData2, result_data[52:78], color="y")
    plt.savefig(r'./data/results/regional_fitting/region' +
                str(region+1)+'/plotICU'+str(region+1)+".png", bbox_inches="tight")
    plt.clf()

    # Plot the number of Recovered
    plt.plot(xData2, gather_italian_data(region)[78:104], color="b")
    plt.plot(xData2, result_data[78:104], color="y")
    plt.savefig(r'./data/results/regional_fitting/region' +
                str(region+1)+'/plotRecovered'+str(region+1)+'.png', bbox_inches="tight")
    plt.clf()

    # Plot the number of Dead
    plt.plot(xData2, gather_italian_data(region)[104:130], color="b")
    plt.plot(xData2, result_data[104:130], color="y")
    plt.savefig(r'./data/results/regional_fitting/region' +
                str(region+1)+'/plotDead'+str(region+1)+'.png', bbox_inches="tight")
    plt.clf()

    f = open(r'./data/results/regional_fitting/region' +
             str(region+1)+'/leastsq_fit/best_values.txt', "w+")
    f.truncate(0)
    f.write(str(results.best_values['E0']))
    f.write('\n')
    f.write(str(results.best_values['C0']))
    f.write('\n')
    f.write(str(results.best_values['R1']))
    f.write('\n')
    f.write(str(results.best_values['R3']))
    f.write('\n')
    f.write(str(results.best_values['R4']))
    f.write('\n')
    f.write(str(results.best_values['R5']))
    f.write('\n')
    f.write(str(results.best_values['R6']))
    f.write('\n')
    f.write(str(results.best_values['R7']))
    f.write('\n')
    f.write(str(results.best_values['R8']))
    f.write('\n')
    f.write(str(results.best_values['alpha']))
    f.write('\n')
    f.write(str(results.best_values['beta']))
    f.write('\n')
    f.write(str(results.best_values['rho']))
    f.write('\n')
    f.write(str(results.best_values['theta']))
    f.write('\n')
    f.write(str(results.best_values['delta']))
    f.write('\n')
    f.write(str(results.best_values['d']))
    f.write('\n')
    f.write(str(results.best_values['total_pop']))
    f.write('\n')
    f.close()
    f = open(r'./data/results/regional_fitting/region' +
             str(region+1)+'/leastsq_fit/quality.txt', "w+")
    area = simpson(gather_italian_data(region), dx=1)
    total_sqd_deviation = np.sum(np.square(results.residual))
    distance = math.sqrt(total_sqd_deviation)/(area)
    f.write(str(distance))
    f.close()


weights = np.array([1, 1, 1, 1, 1])

os.chdir(os.path.dirname(os.path.realpath(__file__)+'/../../../'))
print(os.getcwd())

for i in range(0, 21):
    print_fit(i, weights)
# Calculate the best parameter set by averaging the best parameters.
quality_array = np.empty(21)
avg_params = np.zeros((13, 21))
for i in range(0, 21):
    f = open(r'./data/results/regional_fitting/region' +
             str(i+1)+'/leastsq_fit/quality.txt', "r")
    quality = float(f.readline())
    f.close()
    quality_array[i] = quality
    f = open(r'./data/results/regional_fitting/region' +
             str(i+1)+'/leastsq_fit/best_values.txt', "r")
    index = 1
    for line in f:
        if index >= 3 and index <= 15:
            avg_params[index-3][i] = float(line)
        index = index+1

f = open(r'./data/results/regional_fitting/averages.txt', "w+")
f.write("average quality of the fit: ")
f.write("\n")
f.write(str(np.average(quality_array)))
f.write("\n")
f.write("average parameters: ")
f.write("\n")
for i in range(0, len(avg_params)):
    f.write(str(np.average(avg_params[i])))
    f.write("\n")
# f.write("average ratios of E0, C0 and N")
f.close()

# Calculate the best fits for the values from the paper to compare
paper_params = lmfit.Parameters()
paper_params.add('E0', value=60, min=1, max=1000.0)
paper_params.add('C0', value=10, min=1, max=1000.0)
paper_params.add('R1', value=0.587, min=0.01, max=1.0)
paper_params.add('R3', value=1/(4.2), min=0.01, max=1.0)
paper_params.add('R4', value=1/(14), min=0.01, max=1.0)
paper_params.add('R5', value=1/(16), min=0.01, max=1.0)
paper_params.add('R6', value=1/(2.5), min=0.01, max=1.0)
paper_params.add('R7', value=1/(3.5), min=0.01, max=1.0)
paper_params.add('R8', value=1/(16), min=0.01, max=1.0)
paper_params.add('alpha', value=0.01, min=0.01, max=1.0)
paper_params.add('beta', value=0.05, min=0.01, max=1.0)
paper_params.add('rho', value=0.35, min=0.01, max=1.0)
paper_params.add('theta', value=0.15, min=0.01, max=1.0)
paper_params.add('delta', value=0.77, min=0.01, max=1.0)
paper_params.add('d', value=1/(6.92), min=0.01, max=1.0)
paper_params['R1'].vary = False
paper_params['R3'].vary = False
paper_params['R4'].vary = False
paper_params['R5'].vary = False
paper_params['R6'].vary = False
paper_params['R7'].vary = False
paper_params['R8'].vary = False
paper_params['alpha'].vary = False
paper_params['beta'].vary = False
paper_params['rho'].vary = False
paper_params['theta'].vary = False
paper_params['delta'].vary = False
paper_params['d'].vary = False

'''
paper_quality_array = np.zeros(21)
for region in range(0, 21):
    max_pop = float(linecache.getline(
        "./src/cpp/data_for_simulation/initial_populations/total_populations.txt", region+1))/10
    paper_params.add('region', value=region+1)
    paper_params['region'].vary = False
    paper_params.add('total_pop', value=3000, min=50, max=max_pop)
    total_data_array = gather_italian_data(region)
    my_model = lmfit.Model(func)
    results = my_model.fit(total_data_array, params=paper_params,
                           x=np.arange(0, 130, 1), method='leastsq')
    average_sqd_deviation = 0
    area = simpson(gather_italian_data(region), dx=1)
    average_sqd_deviation = np.average(np.square(results.residual))
    avg_sqd_per_height = average_sqd_deviation*((130.0/area)**2)
    paper_quality_array[i] = avg_sqd_per_height

f = open(r'.data/results/regional_fitting/averages1.txt', "w+")
f.write("\n")
f.write("average quality of paper fit: ")
f.write("\n")
f.write(str(np.average(paper_quality_array)))
f.close()
'''
