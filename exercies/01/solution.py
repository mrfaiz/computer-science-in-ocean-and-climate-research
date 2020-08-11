import numpy as np
import matplotlib.pyplot as plt
# • the thermal coupling constant C = 9.96 × 10 6 ,
# • the emissivity  = 0.62 (corresponding to the fact that a part of the
# outgoing radiation is hold back in the atmosphere due to the green-
# house effect),
# • the Boltzmann constant σ = 5.67 × 10 −8
# • the solar constant S = 136

class Exercise1:
    def __init__(self, year, steps):
        self.c = 9.96 * (10**6)
        self.e = 0.62
        self.alpha = 5.67 * (10 ** -8)
        self.s = 136
        self.steps = steps
        self.n = steps
        self.delet_t = int(year/steps)
        self.T = year
        self.times =  np.arange(0,self.T,self.delet_t)

    def c1(self):
        return 1/(4*self.c)
    
    def c2(self):
        return (self.alpha*self.e)/self.c 
    
    # ẏ(t) = c1*S(1 − α) − c2*y(t)^4 =: f (y(t), t)
    # Energy balance model
    def energy_balance_model(self, y_k):
        return (self.c1() * (self.s * (1-self.alpha))) - (self.c2() * (y_k ** 4))

    # y k+1 = y k + ∆tf (yk , tk ),
    def simulate(self,init_temparature):
        temperature = np.zeros(int(self.n))
        temperature[0] = init_temparature
        for k in range(1, self.n):
            y_k = temperature[k-1]
            temperature[k] = y_k + self.delet_t * self.energy_balance_model(y_k)
        return temperature

    def plots_data(self,intital_temparatures):
        fig,graph = plt.subplots()
        graph.set_xlabel("Time")
        graph.set_ylabel("Temperature")   
        
        for v in intital_temparatures:
            temp = self.simulate(v)
            graph.plot(self.times,temp)
        
        fig.savefig("graph1.png")

T = 10e6
n = 100

obj = Exercise1(T,n)
intital_temparatures = np.arange(50,900,50)
obj.plots_data(intital_temparatures)