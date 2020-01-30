from graphObs import *
import matplotlib.pyplot as plt

d = Graph("test.slha")

d.orgCats()

x1 = [i.x for i in d.higgs]
y1= [i.mass for i in d.higgs]

x2 = [i.x for i in d.sleptons]
y2= [i.mass for i in d.sleptons]

x3 = [i.x for i in d.squarks]
y3= [i.mass for i in d.squarks]

x4 = [i.x for i in d.gauginos]
y4= [i.mass for i in d.gauginos]

x = x1 + x2 + x3 + x4
y = y1 + y2 + y3 + y4

print (x)
plt.scatter(x,y)
plt.xlim(0,5)
plt.xticks([1,2,3,4])
plt.axes().set_xticklabels(['higgs', 'sleptons', 'gauginos', 'squarks'])

plt.show()