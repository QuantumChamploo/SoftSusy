from graphObs import *

#d = Graph("sunday.slha",[25])
d = Graph("slha/span5Point3.slha")


d.plot()
#d.plot(includes=['squarks','sleptons'])
#d.plotBar()

plt.ylim(-500,15000)
plt.title("mMess  10^14 Lambda 10^5.77 TanBeta 20")

plt.show()