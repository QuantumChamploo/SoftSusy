import matplotlib.pyplot as plt
import numpy as np

values = [{"h0":50, "H0":3300, "A0":3300, "H+-":3300},
          {"lR":1420,"t1":1410,"vT":2500, "vL":2500,"lL":2500,"T2":2500,},
         ]
x = []
y = []
n = []


for i in range(len(values)):
    for k in values[i].keys():
        x.append(i)
        y.append(values[i][k])
        n.append(k)
        #plt.scatter(i,values[i][k])
        

plt.scatter(x,y)

for i, txt in enumerate(n):
    plt.annotate(txt, (x[i], y[i]))
    
plt.xticks([0, 1])
plt.axes().set_xticklabels(['higgs', 'sleptons'])

plt.show()