# ẋ = x(α − βy − λx)
# ẏ = y (δx − γ − μy ),
# t ≥ 0.
# https://www.digitalbiologist.com/blog/2018/9/a-population-dynamics-model-in-five-lines-of-python

chickens = [100]
foxes = [10]

chicken_birth_rate = 0.5
chicken_death_rate = 0.015
fox_birth_rate = 0.015
fox_death_rate = 0.5

delta_time = 0.01
cycles = 4500

for t in range(0, cycles):
    updated_chickens = chickens[t] + delta_time * (
        chicken_birth_rate * chickens[t] -
        chicken_death_rate * foxes[t] * chickens[t])
    updated_foxes = foxes[t] + delta_time * (
        -fox_death_rate * foxes[t] + fox_birth_rate * foxes[t] * chickens[t])
    chickens.append(updated_chickens)
    foxes.append(updated_foxes)

import matplotlib.pyplot as plt

fig, graph = plt.subplots()
graph.set_xlabel("Time")
graph.set_ylabel("Changes")

time_points = range(cycles + 1)



graph.plot(time_points, foxes,label ='foxes')
graph.plot(time_points, chickens,label = 'chickens')
graph.legend(loc = 'best')

fig.savefig("/home/faiz/SS_2020/Ocean/exercies/preditor-prey-best-example/pred-prey-chicken-fox.pdf")