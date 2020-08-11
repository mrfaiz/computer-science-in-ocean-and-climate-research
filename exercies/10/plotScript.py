import pandas as pd
import matplotlib.pyplot as plt

fd = pd.read_csv('results.txt')
headers = list(fd)
fig, (ax0) = plt.subplots(nrows = 1)
ax0.plot(fd['time'],fd[headers[1]],label=headers[1])
ax0.plot(fd['time'],fd[headers[2]],label=headers[2])
fig.savefig("plot.png")
