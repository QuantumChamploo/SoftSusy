import matplotlib.pyplot as plt
import numpy as np

y = [(1,1,2,3,9),(1,1,2,4)]
x = [1,2]

for xe, ye in zip(x, y):
    plt.scatter([xe] * len(ye), ye)

plt.xticks([1, 2])
plt.axes().set_xticklabels(['cat1', 'cat2'])
plt.show()