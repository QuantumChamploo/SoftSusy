from graphObs import *
import time


start = time.time()
#d = Graph("sunday.slha",[25])
d = Graph("slha/span5Point3.slha")


d.plot()
#d.plot(includes=['squarks','sleptons'])
#d.plotBar()

plt.ylim(-500,15000)
plt.title("mMess  10^14 Lambda 10^5.77 TanBeta 20")

end = time.time()
print(end - start)

plt.show()

