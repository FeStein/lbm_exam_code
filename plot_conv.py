import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("error.log", names = ["time", "error"])
df.plot("time", "error")
plt.show()
#plt.savefig("convergence.png")
