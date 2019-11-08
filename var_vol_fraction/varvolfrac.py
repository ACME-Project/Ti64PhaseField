import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

files_no = 3000
var = []
for i in range(12):
	var.append([])
times = []
for i in range(0, files_no+1, 100):
	name = "iter2/var_frac_" + str(i) + ".txt"
	df = pd.read_csv(name, sep = ',', header = None)
	for j in range(12):
		var[j].append(df[1][j]/(60**3))
	times.append(i)
	
plt.plot(times, var[0], label = 'Var 1')
plt.plot(times, var[1], label = 'Var 2')
plt.plot(times, var[2], label = 'Var 3')
plt.plot(times, var[3], label = 'Var 4')
plt.plot(times, var[4], label = 'Var 5')
plt.plot(times, var[5], label = 'Var 6')
plt.plot(times, var[6], label = 'Var 7')
plt.plot(times, var[7], label = 'Var 8')
plt.plot(times, var[8], label = 'Var 9')
plt.plot(times, var[9], label = 'Var 10')
plt.plot(times, var[10], label = 'Var 11')
plt.plot(times, var[11], label = 'Var 12')
plt.legend()
plt.xlabel('Timesteps')
plt.ylabel('Volume fraction')
plt.show()

var = []
for i in range(12):
	var.append([])
for i in range(0, files_no+1, 100):
	name = "iter3/var_frac_" + str(i) + ".txt"
	df = pd.read_csv(name, sep = ',', header = None)
	for j in range(12):
		var[j].append(df[1][j]/(60**3))

	
plt.plot(times, var[0], label = 'Var 1')
plt.plot(times, var[1], label = 'Var 2')
plt.plot(times, var[2], label = 'Var 3')
plt.plot(times, var[3], label = 'Var 4')
plt.plot(times, var[4], label = 'Var 5')
plt.plot(times, var[5], label = 'Var 6')
plt.plot(times, var[6], label = 'Var 7')
plt.plot(times, var[7], label = 'Var 8')
plt.plot(times, var[8], label = 'Var 9')
plt.plot(times, var[9], label = 'Var 10')
plt.plot(times, var[10], label = 'Var 11')
plt.plot(times, var[11], label = 'Var 12')
plt.legend()
plt.xlabel('Timesteps')
plt.ylabel('Volume fraction')
plt.show()


var = []
for i in range(12):
	var.append([])
for i in range(0, files_no+1, 100):
	name = "iter4/var_frac_" + str(i) + ".txt"
	df = pd.read_csv(name, sep = ',', header = None)
	for j in range(12):
		var[j].append(df[1][j]/(60**3))

	
plt.plot(times, var[0], label = 'Var 1')
plt.plot(times, var[1], label = 'Var 2')
plt.plot(times, var[2], label = 'Var 3')
plt.plot(times, var[3], label = 'Var 4')
plt.plot(times, var[4], label = 'Var 5')
plt.plot(times, var[5], label = 'Var 6')
plt.plot(times, var[6], label = 'Var 7')
plt.plot(times, var[7], label = 'Var 8')
plt.plot(times, var[8], label = 'Var 9')
plt.plot(times, var[9], label = 'Var 10')
plt.plot(times, var[10], label = 'Var 11')
plt.plot(times, var[11], label = 'Var 12')
plt.legend()
plt.xlabel('Timesteps')
plt.ylabel('Volume fraction')
plt.show()


var = []
for i in range(12):
	var.append([])
for i in range(0, files_no+1, 100):
	name = "iter5/var_frac_" + str(i) + ".txt"
	df = pd.read_csv(name, sep = ',', header = None)
	for j in range(12):
		var[j].append(df[1][j]/(60**3))

	
plt.plot(times, var[0], label = 'Var 1')
plt.plot(times, var[1], label = 'Var 2')
plt.plot(times, var[2], label = 'Var 3')
plt.plot(times, var[3], label = 'Var 4')
plt.plot(times, var[4], label = 'Var 5')
plt.plot(times, var[5], label = 'Var 6')
plt.plot(times, var[6], label = 'Var 7')
plt.plot(times, var[7], label = 'Var 8')
plt.plot(times, var[8], label = 'Var 9')
plt.plot(times, var[9], label = 'Var 10')
plt.plot(times, var[10], label = 'Var 11')
plt.plot(times, var[11], label = 'Var 12')
plt.legend()
plt.xlabel('Timesteps')
plt.ylabel('Volume fraction')
plt.show()
