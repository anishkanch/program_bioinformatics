#!/usr/bin/env python3

#Have global temperatures changed from 1850 and now?

import sys
import csv
import matplotlib.pyplot as plt

temp_data = "global_temperature/temp_anomalies.csv"

years = []
temp = []

with open(temp_data) as f:
    sheet = csv.reader(f)
    #the actual data begins on row 5. we can ignore the first 4 rows
    for i in range(4):
        next(sheet)
    for row in sheet:
        #years are in row number 1, considered 0 
        y = int(row[0][:4])
        years.append(y)
        #temperatures are in row number 2, considered 1
        t = float(row[1])
        temp.append(t)

#experimenting with stuff i found online
    # https://matplotlib.org/stable/users/explain/customizing.html  
    # https://www.geeksforgeeks.org/python/line-chart-in-matplotlib-python/ 
plt.plot(years, temp, color="blue", linewidth=0.75, linestyle=":")  
plt.title("Global Temperature Anomalies Trend From 1850-Present")
plt.xlabel("Year")
plt.ylabel("Anomaly Temperature in Celsius")
plt.grid(True, linestyle="-", color="gray", alpha=0.3)
plt.tick_params(colors="black", width=1)
plt.savefig("plot.png")