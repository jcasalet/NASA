import matplotlib.pyplot as plt
import numpy as np

size=1000
#X = np.random.uniform(0, 1000, size)
X = np.arange(0, size)
sd_deviation = np.random.uniform(0, 4, size=size) + (X / 13)
error = np.random.normal(0, sd_deviation, size=size)
y_train = error + X
plt.plot(y_train, "o")
plt.show()
