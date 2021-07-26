import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.special import erf


u_0 = -0.001
nu = 3e-3
df = pd.read_csv('velo.log', index_col=0, header=None).T

dt = 5e-5

def ut_ana(y,t):
    return u_0 - u_0 * erf(y / (2 * np.sqrt(nu * t)))

(rown, coln) = df.shape

y_vec = np.linspace(0, 2.0, rown)

ax, fig = plt.subplots()

for col in df:
    print(col)
    plt.plot(df[col]/u_0, y_vec)

plt.show()

#compare to analyt solution
col = 1500000.0
curr_t = col * dt
ana = ut_ana(y_vec, curr_t) 

plt.plot(df[col]/u_0, y_vec)
plt.plot(ana/u_0, y_vec)

plt.show()
