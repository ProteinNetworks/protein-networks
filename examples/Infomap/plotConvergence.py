"""Check how the conductance changes overall as N changes."""

import json
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

sns.set()
for i in range(5):
    with open("convergence{}.json".format(i + 1)) as flines:
        convergence = json.load(flines)
    means = []
    Ns = [1, 10, 100, 1000, 10000]
    for N in Ns:
        mean = [np.mean(x) for x in convergence[str(N)]]
        means.append(mean)
    meanDomainConductance = [x[1] for x in means]

    plt.plot(Ns, meanDomainConductance)
plt.show()