import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np

df = pd.read_csv('results.txt', delimiter = ",")
times = df["time"]
headers = list(df)
steps = len(times)

fig = plt.figure()
ax = fig.add_subplot(111)#, projection='3d')

#for header in headers[1:len(headers)]:
    #ax0.plot(times, df[header], label=header)

spaceDivisions = int((len(headers)-1)/2)
print(spaceDivisions)

# Predator population
X = []
T = []
Y = []

for s in range(1,spaceDivisions+1):
	predFirst = headers[2*s-1]
	print(predFirst)
	location = float(predFirst[6:(len(predFirst))])
	predPopFirst = df[predFirst]
	locationList = [location]*len(predPopFirst)
	X += locationList
	for t in times:
		T += [t]
	for y in predPopFirst:
		Y += [y]
	#ax.scatter(times, locationList, predPopFirst,marker='o',c='b')

x = np.reshape(X, (spaceDivisions, steps))
y = np.reshape(T, (spaceDivisions, steps))
z = np.reshape(Y, (spaceDivisions, steps))
#ax0.legend()

levels = np.linspace(min(Y), max(Y), 100)


cnt = ax.contourf(x, y, z,levels=levels,cmap='inferno')
for c in cnt.collections:
    c.set_edgecolor("face")

ax.set_xlabel("location x")
ax.set_ylabel("time t")


#plt.show()
fig.savefig("plot.pdf")
