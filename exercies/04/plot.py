import numpy as np
import matplotlib.pyplot as plt

times_file='/home/faiz/SS_2020/Ocean/exercies/04/time.txt'
prey_file='/home/faiz/SS_2020/Ocean/exercies/04/prey.txt'
preditor_file='/home/faiz/SS_2020/Ocean/exercies/04/preditor.txt'
outfilename='/home/faiz/SS_2020/Ocean/exercies/04/output.png'

#read data from row 1 onwoards

times = np.loadtxt(times_file, delimiter=',', usecols=range(1))

prey_data = np.loadtxt(prey_file, delimiter=',', usecols=range(100))
preditor_data = np.loadtxt(preditor_file, delimiter=',', usecols=range(100))

# print(len(prey_data))
# print (len(preditor_data))



# plot figure
fig = plt.figure()
# dataLines=[]
# for ii in range(1,x.shape[0]):
#     lineIi, = plt.plot(x[0,:],x[ii,:])
for i in (1,99):
    plt.plot(prey_data[:,i],times)
    plt.plot(preditor_data[:,i],times)

plt.xlabel('space x')
plt.ylabel('time')
# plt.legend(handles=dataLines, loc='best')
#plt.show()
fig.savefig("/home/faiz/SS_2020/Ocean/exercies/04/pred-prey-space.png")