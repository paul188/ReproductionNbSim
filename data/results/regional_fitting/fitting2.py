import subprocess
import lmfit
import codecs
import matplotlib.pyplot as plt
import numpy as np
import csv
import datetime
import math

region = 1


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
                    temparray2.append(rows[6])  # Hospitalized with symptoms
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

while (res_date <= end):
    total_data_array.append(float(arrayfinal[region][index][0]))
    res_date += datetime.timedelta(days=1)
    index += 1

res_date = start
index = 0

while (res_date <= end):
    total_data_array.append(float(arrayfinal[region][index][1]))
    res_date += datetime.timedelta(days=1)
    index += 1

res_date = start
index = 0

while (res_date <= end):
    total_data_array.append(float(arrayfinal[region][index][4]))
    res_date += datetime.timedelta(days=1)
    index += 1


def linearize(x, xleft, xright, yleft, yright):
    return yleft+(x-xleft)*(yright-yleft)/(xright-xleft)

# Define fitting function


def func(x, R1, R3, R4, R5, R6, R7, R8, alpha, beta, rho, theta, delta, d, E0, C0, region):  # x should be between 0 and 75
    print("Function called")
    print(R1, R3, R4, R5, R6, R7, R8,
          alpha, beta, rho, theta, delta, d, E0, C0, region)
    # Make x iterable if it is only a float/int
    if not hasattr(x, '__iter__'):
        x = np.array([x])
    result = []
    output = subprocess.check_output([r'C:\Users\paul1\OneDrive\Desktop\epidemiology\coding\Secihurd_Model\build\Debug\model.exe', str(R1), str(R3), str(R4), str(R5), str(R6), str(R7), str(R8), str(alpha), str(beta), str(rho), str(theta), str(delta), str(d), str(
        E0), str(C0), str(region)])
    output_str = codecs.decode(output)
    output_str_list = output_str.split('\r\n')
    output_str_list.pop()
    dataarray = []
    index = 0
    for word in output_str_list:
        if index in range(412, 618) or index in range(714, 824):
            if word == '-nan(ind)':
                dataarray.append(0.0)
            else:
                dataarray.append(float(word))
        index += 1
    for x0 in x:
        # linearize(x0*4, math.floor(x0*4), math.floor(x0*4)+1,
        #          dataarray[math.floor(x0*4)], dataarray[math.ceil(x0*4)])
        y = 1
        print("length", len(dataarray))
        print("dataarray:")
        print(dataarray)
        print(math.ceil(x0*4))
        result.append(y)
    return result


parameters = lmfit.Parameters()
parameters.add('region', value=region)
parameters.add('E0', value=10, min=5, max=10000.0)
parameters.add('C0', value=5, min=5, max=10000.0)
parameters.add('R1', value=0.587, min=0.05, max=1.0)
# Lower bound to make sure R3 well-defined
parameters.add('R3', value=0.3125, min=1/(4.2), max=1.0)
parameters.add('R4', value=0.1666, min=0.05, max=1.0)
parameters.add('R5', value=0.1, min=0.05, max=1.0)
parameters.add('R6', value=0.2, min=0.05, max=1.0)
parameters.add('R7', value=0.4, min=0.05, max=1.0)
parameters.add('R8', value=0.125, min=0.05, max=1.0)
parameters.add('alpha', value=0.09, min=0.05, max=1.0)
parameters.add('beta', value=0.25, min=0.05, max=1.0)
parameters.add('rho', value=0.2, min=0.05, max=1.0)
parameters.add('theta', value=0.26, min=0.05, max=1.0)
parameters.add('delta', value=0.77, min=0.05, max=1.0)
parameters.add('d', value=0.99, min=0.05, max=1.0)
# bound to make sure R9 well defined
parameters.add('boundR9', expr='(1/R3)+0.5*(1/R4)', min=1.0)
parameters['region'].vary = False
parameters['E0'].vary = False
parameters['C0'].vary = False

xData = np.arange(0, 78, 1)
# Running this part of the code works
my_model = lmfit.Model(func)
# results = my_model.fit(total_data_array, params=parameters,
#                       x=xData, method='leastsq')
# ci = lmfit.conf_interval(results, results)
# lmfit.report_ci(ci)

func(78, 0.5, 0.238, 0.0714, 0.0625, 0.4, 0.285, 0.0625,
     0.01, 0.05, 0.35, 0.77, 0.153, 0.01, 50, 10, region)

# plt.plot(xData, total_data_array, color="y")
# plt.plot(xData, func(xData, results.best_values["R1"], results.best_values["R3"], results.best_values["R4"], results.best_values["R5"], results.best_values["R6"], results.best_values["R7"], results.best_values[
#         "R8"], results.best_values["alpha"], results.best_values["beta"], results.best_values["rho"], results.best_values["theta"], results.best_values["delta"], results.best_values["d"], 10, 5, 1), color="b")
# plt.plot(xData, func(xData, 0.587, 0.238, 0.0714, 0.0625, 0.4, 0.285,
#         0.0625, 0.01, 0.05, 0.35, 0.77, 0.153, 0.01, 50, 10, region), color="b")
# plt.plot(xData, func(xData, 0.5, 0.238, 0.0714, 0.0625, 0.4, 0.285,
#        0.0625, 0.01, 0.05, 0.35, 0.77, 0.153, 0.01, 50, 10, region), color="g")
plt.show()
