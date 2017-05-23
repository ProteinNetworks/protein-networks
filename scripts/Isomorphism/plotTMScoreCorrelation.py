import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv(
    "./correlationBetweenTMandMCS.dat",
    sep=" ",
    names=["pdb1", "pdb2", "MCS", "SimScore", "TM1", "TM2"])
sns.set(style="white")
sns.jointplot(x="MCS", y="TM2", data=df, space=0, kind='reg')
plt.savefig("networkvsstructurecorr.png", dpi=300)
plt.show()
