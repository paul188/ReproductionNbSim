import subprocess
import lmfit
import codecs
import matplotlib.pyplot as plt
import numpy as np
import csv
import datetime
import math


def linearize(x, xleft, xright, yleft, yright):
    return yleft+(x-xleft)*(yright-yleft)/(xright-xleft)


# read in the Italian data
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
        str_file = 'C:/Users/paul1/OneDrive/Desktop/epidemiology/coding/parameter_fitting/dati-regioni/dpc-covid19-ita-regioni-'+str_date+'.csv'
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
        float(arrayfinal[0][index][2]))
    res_date += datetime.timedelta(days=1)
    index += 1

res_date = start
index = 0

# Store data on Hospitalized
while (res_date <= end):
    total_data_array.append(float(arrayfinal[0][index][0]))
    res_date += datetime.timedelta(days=1)
    index += 1

res_date = start
index = 0

# Store data on ICU patients
while (res_date <= end):
    total_data_array.append(
        float(arrayfinal[0][index][1]))
    res_date += datetime.timedelta(days=1)
    index += 1

res_date = start
index = 0

# Store data on recovered
while (res_date <= end):
    total_data_array.append(
        float(arrayfinal[0][index][3]))
    res_date += datetime.timedelta(days=1)
    index += 1

res_date = start
index = 0

# Store data on deceased
while (res_date <= end):
    total_data_array.append(
        float(arrayfinal[0][index][4]))
    res_date += datetime.timedelta(days=1)
    index += 1


def linearize(x, xleft, xright, yleft, yright):
    return yleft+(x-xleft)*(yright-yleft)/(xright-xleft)

# Define fitting function


def func(x, region, E0, C0, R1, R3, R4, R5, R6, R7, R8, alpha, beta, rho, theta, delta, d, total_pop):  # x should be between 0 and 75
    print("Function called")
    print(E0, C0, R1, R3, R4, R5, R6, R7, R8,
          alpha, beta, rho, theta, delta, d, total_pop)
    # Make x iterable if it is only a float/int
    if not hasattr(x, '__iter__'):
        x = np.array([x])
    result = []
    output = subprocess.check_output([r'C:\Users\paul1\OneDrive\Desktop\epidemiology\coding\Secihurd_Model\build\Debug\model.exe', str(region), str(
        E0), str(C0), str(R1), str(R3), str(R4), str(R5), str(R6), str(R7), str(R8), str(alpha), str(beta), str(rho), str(theta), str(delta), str(d), str(total_pop)])
    output_str = codecs.decode(output)
    output_str_list = output_str.split('\r\n')
    output_str_list.pop()
    dataarray = []
    index = 0
    for word in output_str_list:
        if index in range(312, 416):  # Cumulative Infected
            dataarray.append(float(word)+float(output_str_list[index+104])+float(
                output_str_list[index+208]))  # Infected + Hospitalized + ICU
        elif index in range(416, 832):  # Hospitalized+ICU+Recovered+Dead
            dataarray.append(float(word))
        index += 1
    for x0 in x:
        y = linearize(x0*4, math.floor(x0*4), math.floor(x0*4)+1,
                      dataarray[math.floor(x0*4)], dataarray[math.ceil(x0*4)])
        result.append(y)
    return result


parameters = lmfit.Parameters()
parameters.add('region', value=1)
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
parameters.add('total_pop', value=3000, min=500, max=20000)
parameters['region'].vary = False


xData = np.arange(0, 130, 1)
my_model = lmfit.Model(func)
results = my_model.fit(total_data_array, params=parameters,
                       x=xData, method='leastsq')
# ci = lmfit.conf_interval(results, results)
# lmfit.report_ci(ci)

tick_labels = ['Cumulative Infected',
               'Hospitalized', 'ICU', 'Recovered', 'Dead']
tick_locations = [13, 39, 65, 91, 117]

plt.xticks(tick_locations, tick_labels)

plt.plot(xData, total_data_array, color="b")
plt.plot(xData, func(xData, results.best_values['region'], results.best_values['E0'], results.best_values['C0'], results.best_values['R1'], results.best_values['R3'], results.best_values['R4'], results.best_values['R5'], results.best_values['R6'],
                     results.best_values['R7'], results.best_values['R8'], results.best_values['alpha'], results.best_values['beta'], results.best_values['rho'], results.best_values['theta'], results.best_values['delta'], results.best_values['d'], results.best_values['total_pop']), color="y")

plt.show()
