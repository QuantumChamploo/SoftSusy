from graphObs import *
d = Graph("test.slha")
for j in d.higgs:
	print (j.mass)
d.orgCats()

for i in d.higgs:
	print (i.mass)
	print (i.delta)

for i in d.sleptons:
	print (i.mass)
	print (i.delta)

for i in d.squarks:
	print (i.mass)
	print (i.delta)

for i in d.gauginos:
	print (i.mass)
	print (i.delta)