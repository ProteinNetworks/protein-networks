import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv(
    "simscorecomparison.dat",
    sep=" ",
    names=["pdb1", "pdb2", "MCS", "SimScore"])
sns.set(style="white")
sns.jointplot(x="MCS", y="SimScore", data=df, space=0)
plt.show()
