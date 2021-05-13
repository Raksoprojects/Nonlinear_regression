import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import fmin

t = np.array([])
y = np.array([])

with open("lab_6_dane/data12.txt", "r") as data:
    lines = data.readlines()
    
    for l in lines:
        line = l.split(' ')
        #print(line)
        a = line[-1].split()
        #print(a)
        t = np.append(t, float(line[2]))
        y = np.append(y, float(a[0]))

def regr_nielin(arg):
    suma = 0 
    for i in range(len(t)):
        k = arg[0]
        tau = arg[1]
        zeta = arg[2]
        tau_z = arg[3]

        wn = 1/tau
        w0 = wn * np.sqrt(1 - zeta ** 2)
        alfa = np.arctan(np.sqrt(1 - zeta ** 2)/ zeta)

        g = (wn * (np.e ** (-zeta * wn * t[i])) * np.sin(w0*t[i])/np.sqrt(1 - zeta ** 2))
        h = (1 - (np.e ** (-zeta * wn * t[i])) * np.sin(w0*t[i] + alfa)/np.sqrt(1 - zeta ** 2))

        yyy = ( k * (tau_z*g + h))

        suma += (y[i] - yyy)**2
    return suma

parameters = fmin(regr_nielin, np.array([2, 1.6, 0.4, 1.5]))

print("Parametry: k, tau, zeta, tau_z")
print(parameters)

czas = np.array([])
wartosc = np.array([])

k = parameters[0]
tau = parameters[1]
zeta = parameters[2]
tau_z = parameters[3]

for i in range(len(t)): 
        wn = 1/tau
        w0 = wn * np.sqrt(1 - zeta ** 2)
        alfa = np.arctan( np.sqrt(1 - zeta ** 2)/ zeta)

        g = (wn * (np.e ** (-zeta * wn * t[i])) * np.sin(w0*t[i])/np.sqrt (1 - zeta ** 2) )
        h = (1 - (np.e ** (-zeta * wn * t[i])) * np.sin(w0*t[i] + alfa)/np.sqrt (1 - zeta ** 2) )

        yyy = ( k * (tau_z*g + h))
        czas = np.append(czas, t[i])
        wartosc = np.append(wartosc, yyy)

plt.plot(t,y)
plt.plot(czas, wartosc)
plt.show()
