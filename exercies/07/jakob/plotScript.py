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

spaceDivisions = int((len(headers)-1)/2)

# Predator population
X = []
T = []
YPred = []
YPrey = []

lastLocation = []
lastPrey     = []
lastPred     = []

fraction = int(spaceDivisions/10)

for s in range(1,spaceDivisions+1):
	predFirst = headers[2*s-1]
	preyFirst = headers[2*s]
	locationPred = float(predFirst[6:(len(predFirst))])
	locationPrey = float(preyFirst[6:(len(preyFirst))])

	predPopFirst = df[predFirst]
	preyPopFirst = df[preyFirst]
	locationList = [locationPred]*len(predPopFirst)

	lastLocation += [locationPred]
	lastPred     += [predPopFirst[len(predPopFirst)-1]]
	lastPrey     += [preyPopFirst[len(preyPopFirst)-1]]

	X += locationList
	for t in times:
		T += [t]
	for y in predPopFirst:
		YPred += [y]
	for y in preyPopFirst:
		YPrey += [y]

x = np.reshape(X, (spaceDivisions, steps))
y = np.reshape(T, (spaceDivisions, steps))
zPred = np.reshape(YPred, (spaceDivisions, steps))
zPrey = np.reshape(YPrey, (spaceDivisions, steps))
#ax0.legend()
print("min(Y)", min(YPred))
print("max(Y)", max(YPred))
levelsPred = np.linspace(min(YPred), max(YPred), 100)
levelsPrey = np.linspace(min(YPrey), max(YPrey), 100)



# Plotting predator population against time

figPred = plt.figure()
axPred = figPred.add_subplot(111)

cntPred = axPred.contourf(x, y, zPred,levels=levelsPred,cmap='inferno')
for c in cntPred.collections:
    c.set_edgecolor("face")

axPred.set_xlabel("location x")
axPred.set_ylabel("time t")

figPred.savefig("predatorTime.png")

# Plotting prey population against time

figPrey = plt.figure()
axPrey = figPrey.add_subplot(111)

cntPrey = axPrey.contourf(x, y, zPrey,levels=levelsPrey,cmap='inferno')
for c in cntPrey.collections:
    c.set_edgecolor("face")

axPrey.set_xlabel("location x")
axPrey.set_ylabel("time t")

figPrey.savefig("preyTime.png")


# Plotting the last iteration

figLast = plt.figure()
axLast = figLast.add_subplot(111)

axLast.plot(lastLocation, lastPred)
axLast.plot(lastLocation, lastPrey)

axLast.set_xlabel("location x")
axLast.set_ylabel("population")

figLast.savefig("lastIteration.png")
