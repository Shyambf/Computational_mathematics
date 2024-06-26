import matplotlib.pyplot as plt
import numpy as np
import sys

m = int(sys.argv[1])
t_max =  float(sys.argv[2].replace(',', '.'))
tau = float(sys.argv[3].replace(',', '.')) 
x = np.arange(0, t_max,  tau)
try:
    for i in range(m):
        with open(f'./txt/output{i}.txt', 'r') as file:
            func1 = list()
            func2 = list()
            for line in file:
                data = line.split()
                func1.append(float(data[0].replace(',', '.')))
                func2.append(float(data[1].replace(',', '.')))
            plt.plot(x, func1, label=f'{i}_1')
            plt.plot(x, func2, label=f'{i}_2')
except BaseException:
    x = np.arange(0, t_max + tau, tau)
    for i in range(m):
        with open(f'./txt/output{i}.txt', 'r') as file:
            func1 = list()
            func2 = list()
            for line in file:
                data = line.split()
                func1.append(float(data[0].replace(',', '.')))
                func2.append(float(data[1].replace(',', '.')))
            plt.plot(x, func1, label=f'{i}_1')
            plt.plot(x, func2, label=f'{i}_2')

plt.xlabel('x')
plt.ylabel('y')
plt.title('plot')
plt.legend()

plt.show()
