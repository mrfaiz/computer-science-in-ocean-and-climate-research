import pandas as pd
import matplotlib.pyplot as plt
import sys

df = pd.read_csv(sys.argv[1], delimiter = ",")
times = df["time"]
headers = list(df)
fig, (ax0) = plt.subplots(nrows = 1)

for header in headers[1:len(headers)]:
    ax0.plot(times, df[header], label=header)

ax0.legend()
ax0.set_xlabel("time / s")
fig.savefig("plot.pdf")
