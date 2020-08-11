
import numpy as np
import matplotlib.pyplot as plt
from PreditorModel import PreditorModel
from PreyModel import PreyModel
from ImprovedEulerMethod import ImprovedEulerMethod
from EulerMethod import EulerMethod

class TestPreditorPrey:
    def __init__(self,t0,delta_t,steps):
        self.steps = steps
        self.times =  [t0]
        self.prey = [100]
        self.preditor = [10]
        self.delta_t = delta_t

    def testfunction(self):
        print(self.delta_t)

    def time_loop(self, time_integrator, prey_model_function,preditor_model_function):
        for k in range(0, self.steps-1):
            previous_index = k-1 
            new_time = self.times[previous_index]+ self.delta_t
            self.times.append(new_time)

            xk = self.prey[previous_index]
            yk = self.preditor[previous_index]
            prey_result = time_integrator(yk,xk,prey_model_function,self.delta_t)
            preditor_result = time_integrator(xk,yk,preditor_model_function,self.delta_t) 
            self.prey.append(prey_result)
            self.preditor.append(preditor_result)

    def getTimes(self):
        return self.times
    
    def getPrey(self):
        return self.prey
    
    def getPreditor(self):
        return self.preditor

fig, graph = plt.subplots()

### Energy Balance model


use_imporved_euler_model = False
output_file = "e02_preditor_prey_euler_method.png"
timeIntgrator = EulerMethod()

if use_imporved_euler_model:
    output_file = "e02_preditor_prey_improved_euler.png"
    timeIntgrator = ImprovedEulerMethod()
    # RuntimeWarning: overflow encountered in double_scalars


# x prey , y preditor
alpha_prey = 0.5  #birth rate
beta_prey = 0.015  # death rate
gamma_preditor = 0.5  #death rate
delta_preditor = 0.015  # birth rate
lmda_prey = 0.0
mu_preditor = 0.0

preyModel = PreyModel(alpha_prey,beta_prey,lmda_prey)
preditorModel = PreditorModel(gamma_preditor,delta_preditor,mu_preditor)

delta_time = 0.01
steps = 4500

obj = TestPreditorPrey(0.0,delta_time,steps)
obj.time_loop(timeIntgrator.integrate,preyModel.model,preditorModel.model)
times = obj.getTimes()
prey_data = obj.getPrey()
preditor_data = obj.getPreditor()
print(len(prey_data))

## Draw accross time
graph.plot(times,prey_data,label="Prey")
graph.plot(times,preditor_data,label="Preditor")
graph.legend()

graph.set_xlabel("Time")
graph.set_ylabel("Population")
fig.savefig(output_file)

## Draw Prey and Preditor
# graph.plot(prey_data,preditor_data)
# graph.set_xlabel("Prey")
# graph.set_ylabel("Preditor")
# fig.savefig("preditor_prey")



