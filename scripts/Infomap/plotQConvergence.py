"""Check how the conductance changes overall as N changes."""

import json
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

sns.set()
for i in range(5):
    with open("Qconvergence{}.json".format(i)) as flines:
        convergence = json.load(flines)
    Ns = [1, 10, 100, 1000, 10000]
    level1Modularity = [convergence[str(N)][1] for N in Ns]

    plt.plot(Ns, level1Modularity)
plt.show()