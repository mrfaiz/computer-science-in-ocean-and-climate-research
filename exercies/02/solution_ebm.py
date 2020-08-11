from TimeIntegrator import TimeIntegrator
from model import Model
from EnergyBalancedModel import EnergyBalancedModel
from ImprovedEulerMethod import ImprovedEulerMethod
from EulerMethod import EulerMethod
import matplotlib.pyplot as plt

import numpy as np

class TestClassEBM:
    def __init__(self,T,t0,delta_t,y0,steps):
        self.steps = steps
        self.times =  np.arange(t0, T, delta_t)
        self.results = [y0]
        self.delta_t = delta_t

    def testfunction(self):
        print(self.delta_t)

    def time_loop(self, time_integrator, model_function):
        for k in range(0, self.steps-1):
            previous_index = k-1 
            yk = self.results[previous_index]
            tk = self.times[previous_index]
            result = time_integrator(tk,yk,model_function,delta_t)
            self.results.append(result)

    def getTimes(self):
        return self.times
    
    def getResult(self):
        return self.results

fig, graph = plt.subplots()
graph.set_xlabel("Time")
T = 10e6 
t0 = 0
# y0 = 100.0
steps = 1000


intial_values = np.arange(100, 700, 50)
use_imporved_euler_model = True
output_file = "euler_method.png"
timeIntgrator = EulerMethod()

if use_imporved_euler_model:
    output_file = "ebm_improved_euler.png"
    timeIntgrator = ImprovedEulerMethod()
    # RuntimeWarning: overflow encountered in double_scalars
    intial_values = np.arange(10, 50, 10)
    T = 10e4

embModel = EnergyBalancedModel()
delta_t = int((T - t0) / steps)

for y0 in intial_values:
    obj = TestClassEBM(T,t0,delta_t,y0,steps)
    obj.time_loop(timeIntgrator.integrate,embModel.model)
    data = obj.getResult()
    graph.plot(obj.getTimes(), data, label='y0 =' + str(y0))
    graph.legend()

graph.set_ylabel("Temperature")
fig.savefig(output_file)
