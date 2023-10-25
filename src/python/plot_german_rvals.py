# Plot the R-values of the RKI (Nowcasting values) from April 2020 to April 2021

import matplotlib as plt
import csv
import datetime
import math

start = datetime.date(2020, 4, 1)
end = datetime.date(2021, 4, 1)
res_date = start

str_file = "C:/Users/paul1/OneDrive/Desktop/epidemiology/data/rki_R_values_germany/Nowcast_R_aktuell.csv"
with open(str_file) as csv_file:
    
