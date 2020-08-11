class EnergyBalancedModel:
    def __init__(self, delta_t, intial_value):
        self.c = 9.96 * (10**6)
        self.e = 0.62
        self.alpha = 5.67 * (10**-8)
        self.s = 136
        self.delta_t = delta_t
        self.resutls = [intial_value]

    def c1(self):
        return 1 / (4 * self.c)

    def c2(self):
        return (self.alpha * self.e) / self.c

    # ẏ(t) = c1*S(1 − α) − c2*y(t)^4 =: f (y(t), t)
    # Energy balance model
    def energy_balance_model(self, k):
        return (self.c1() * (self.s *
                             (1 - self.alpha))) - (self.c2() *
                                                   (self.resutls[k]**4))

    # y k+1 = y k + ∆tΦ(t k , y k , f, ∆t)
    def time_integrator(self, k, model_function):
        yk_plus_1 = self.resutls[k] + self.delta_t * model_function(k)
        self.resutls.append(yk_plus_1)

    def get_result(self):
        return self.resutls


# import numpy as np
# import matplotlib.pyplot as plt

# T = 10e6
# p = 100
# t0 = 0
# delta_t = int((T - t0) / p)
# intial_values = np.arange(100, 800, 100)
# times = np.arange(t0, T, delta_t)

# fig, graph = plt.subplots()
# graph.set_xlabel("Time")
# graph.set_ylabel("Temperature")

# for y0 in intial_values:
#     model = EnergyBalancedModel(delta_t, y0)
#     for k in range(0, p - 1):
#         model.time_integrator(k, model.energy_balance_model)

#     data = model.get_result()
#     graph.plot(times, data, label='y0 =' + str(y0))
#     graph.legend(loc='best')

# fig.savefig("/home/faiz/SS_2020/Ocean/exercies/02/ebm_temp.pdf")