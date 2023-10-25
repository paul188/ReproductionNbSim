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

areas = np.zeros(21)
qualities = np.zeros(21)
for region in range(0,21):
    areas[region] = math.sqrt(simpson(gather_italian_data(region),dx=1))
    f = open(r'C:/Users/paul1/OneDrive/Desktop/epidemiology/coding/Secihurd_Model/src/python/regional_fitting/region' + str(region+1)+'/leastsq_fit/quality.txt',"r")
    qualities[region] = float(f.readline())
    f.close()
corr = np.corrcoef(areas, qualities)
plt.scatter(areas,qualities,color="r")
plt.xlabel("area")
plt.ylabel("fit quality metric")
plt.show()
print("Correlation: ")
print(corr)
print("\n")
print("Best fit: ")
print("region: ")
print(np.argmin(qualities)+1)
print("Worst fit: ")
print("region: ")
print(np.argmax(qualities)+1)