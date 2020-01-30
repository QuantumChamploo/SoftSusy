import matplotlib.pyplot as plt

def testfunc(x,y):
	plt.scatter(x,y)
	plt.show()

x = [1,2,3,4]
y = [1,4,9,16]

testfunc(x,y)