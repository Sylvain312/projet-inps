import matplotlib.pyplot as plt
import numpy as np
import csv

Y= []

with open ("./result.csv",'r') as file:
    csvreader = csv.reader(file)
    i = 0
    for row in csvreader:
        Y.append([])
        for e in row:
            Y[i].append(float(e) * 10**80)
        i+=1
    file.close

X = np.linspace(-5, 5, 501)

plt.clf()
for i in range(len(Y)):
    plt.plot(Y[i], X, color='blue')
plt.savefig("pres/images/result_2D.png")