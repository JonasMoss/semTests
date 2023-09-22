import numpy as np
import seaborn as sns
import matplotlib.pylab as plt

rng = np.random.default_rng(seed=313)
x = rng.exponential(2, (10000, 10))
medians = np.median(x, axis = 1)
sns.histplot(medians, stat="density")
plt.show()
plt.clf()
rng = np.random.default_rng(seed=313)
x = rng.exponential(2, 100)
sns.ecdfplot(x)
y = np.linspace(0, np.max(x), 100)
plt.plot(y, 1 - np.exp(-y/2))
plt.show()
