from graphObs import *

d = Graph("sunday.slha",[25])
#d = Graph("test.slha")


#d.plot()
#d.plot(includes=['squarks','sleptons'])
d.plotBar()

plt.ylim(-500,11000)
plt.title("This is a pretty great title")

plt.show()