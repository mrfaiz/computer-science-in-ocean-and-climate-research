import numpy as np
import matplotlib.pyplot as plt

path = "/home/faiz/SS_2020/Ocean/cs_in_ocean_and_climate_research_python_fortran/exercies/10/faiz"
infilename= path+'/output.txt'
outfilename= path+'/improved_euler.png'

#read labels from first line
with open(infilename) as f:
    dataLabels = f.readline().split(',')

#clean labels (whitespace/tab/newline)
for s in dataLabels:
    s.lstrip().rstrip()

#read data from row 2 onwoards
#x[i,:] data i, i=0 is time,...
x = np.loadtxt(infilename, delimiter=',', unpack=True,skiprows=1)

#plot figure
fig = plt.figure()
dataLines=[]
for ii in range(1,x.shape[0]):
    lineIi, = plt.plot(x[0,:],x[ii,:], label=dataLabels[ii])
    dataLines.append(lineIi)

plt.xlabel('time')
plt.ylabel('populations')
plt.legend(handles=dataLines, loc='best')
#plt.show()
fig.savefig(outfilename, bbox_inches='tight')

fig.savefig(outfilename)