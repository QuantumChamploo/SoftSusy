from graphObs import *
d = Graph("test.slha")

slep = d.gauginos


print (slep)

print( "waiting...")

clustersleps = clusterFunc2(slep)
k = 0
fitcluster(clustersleps)
for i in clustersleps:
	k+=1
	print("in the" + str(k) + "th cluster")
	for j in i:
		print(j.mass)
		print(j.delta)
