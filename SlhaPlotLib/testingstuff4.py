from graphObs import *

d = Graph("test.slha")

print("higgs")

for i in d.higgs:
	print(i.mass)

print("sleptons")

for i in d.sleptons:
	print(i.mass)

print("gauginos")

for i in d.gauginos:
	print(i.mass)

print("squarks")

for i in d.squarks:
	print(i.mass)

