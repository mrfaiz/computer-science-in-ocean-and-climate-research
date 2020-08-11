import pandas as pd
import matplotlib.pyplot as plt

fd = pd.read_csv('ensembleResults.csv')
headers = list(fd)
fig, (ax0) = plt.subplots(nrows = 1)
ax0.plot(fd[headers[0]],fd[headers[1]])
ax0.set_xlabel(headers[0])
ax0.set_ylabel(headers[1])
fig.savefig("plot.png")
