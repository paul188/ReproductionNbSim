from distutils import text_file
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

fig, ax = plt.subplots()

region = 1

filename = "data"+str(region)+".txt"

x = list(range(0, 100))

with open(filename, "r") as text_file:
    data = text_file.readlines()[0:100]

data2 = [float(i) for i in data]  # Susceptibles

plt.yticks(np.arange(0, data2[0], data2[0]/10))

with open(filename, "r") as text_file:
    data3 = text_file.readlines()[102:202]

data4 = [float(i) for i in data3]  # Exposed

with open(filename, "r") as text_file:
    data5 = text_file.readlines()[204:304]

data6 = [float(i) for i in data5]  # Carrier

with open(filename, "r") as text_file:
    data7 = text_file.readlines()[306:406]

data8 = [float(i) for i in data7]  # Infected

with open(filename, "r") as text_file:
    data9 = text_file.readlines()[408:508]

data10 = [float(i) for i in data9]  # Severe

with open(filename, "r") as text_file:
    data11 = text_file.readlines()[510:610]

data12 = [float(i) for i in data11]  # Critical

with open(filename, "r") as text_file:
    data13 = text_file.readlines()[612:712]

data14 = [float(i) for i in data13]  # Recovered

with open(filename, "r") as text_file:
    data15 = text_file.readlines()[714:814]

data16 = [float(i) for i in data15]  # Dead

plt.plot(x, data2, "g", label="Susceptibles")
plt.plot(x, data4, "r", label="Exposed")
plt.plot(x, data6, "blue", label="Carrier")
plt.plot(x, data8, "yellow", label="Infected")
plt.plot(x, data10, "black", label="Severe")
plt.plot(x, data12, "cyan", label="ICU")
plt.plot(x, data14, "orange", label="Recovered")
plt.plot(x, data16, "brown", label="Dead")

plt.xlabel("time elapsed in days")
plt.ylabel("number of individuals")

plt.title("Plot of the SECIHURD Model data")
plt.legend()
plt.show()
