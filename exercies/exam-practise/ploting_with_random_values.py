import numpy as np 
import matplotlib.pyplot as plt


x_values = np.arange(-15,16,1)
# y_values = np.zeros(len(x_values))
y_values = []

for x in x_values:
    y_values.append(x*x)

print(y_values)

fig,graph = plt.subplots()
graph.set_ylabel("square")
graph.set_xlabel("normal distribution")
graph.plot(x_values,y_values)

fig.savefig("exam-practise/square.png")
