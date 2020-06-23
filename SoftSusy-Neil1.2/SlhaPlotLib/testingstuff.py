from graphObs import *
d = Graph("test.slha")

slep = d.sleptons
res = [k[1] for k in slep]
hld = [abs(i-j) for i, j in zip(res[1:],res[:-1])]

print (slep)
print (hld)
print( "waiting...")

print (clusterFunc(hld,slep))