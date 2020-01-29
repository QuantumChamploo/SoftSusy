import pyslha
import os

path = os.getcwd()

#todo:
# create more robust file aquitision 
# fill out self.sleptons, self.gauginos ... with cluster organization
# fully flesh out the catagories below


# Graph class wrapper around the slha file to be used
# when parsing data

cluster_thresh = 100
#pdg codes for various particles, grouped as they will be for the graph
higgs = [24, 25, 35, 36, 37]
sleptons = [1000011, 1000013, 1000015, 2000011, 2000013, 2000015, 1000012,
1000014, 1000016]
squarks = [1000001, 1000003, 1000005, 2000001, 2000003, 2000005,
1000002, 1000004, 1000006, 2000002, 2000004, 2000006]
gauginos = [1000021, 1000022, 1000023, 1000024, 1000025, 1000035,
1000037, 1000039]
def clusterFunc(sort_arr,full_arr):
	hld = []
	i = 0
	while (i < len(full_arr)):
		
		hld2 = [full_arr[i]]
		
		while ((i<len(sort_arr))):
			if(sort_arr[i]<=cluster_thresh):
				i+=1
				hld2.append(full_arr[i])
			else:
				i+=1
				hld.append(hld2)
				break
		if(i == len(sort_arr)):
			hld.append(hld2)
			break 

	return hld 


class Graph:
	def __init__(self, file):
		self.file = pyslha.read(path + "/" + file)
		self.masses = self.file.blocks['MASS'].items()
		self.higgs = [i for i in self.masses if i[0] in higgs]
		self.sleptons = [i for i in self.masses if i[0] in sleptons]
		self.squarks = [i for i in self.masses if i[0] in squarks]
		self.gauginos = [i for i in self.masses if i[0] in gauginos]

	def orgCats(self):
		hld_higgs = []
		higgs_deltas = [abs(j-i) for i,j in zip]

class Particle:
	def __init__(self,pdg,mass):
		self.pdg = pdg
		self.mass = mass
		self.delta = 0
		if(self.pdg in higgs):
			self.cat = higgs
		if(self.pdg in sleptons):
			self.cat = slepton
		if(self.pdg in squarks):
			self.cat = squark
		if(self.pdg in gauginos):
			self.cat = gaugino
		else:
			self.cat = Na
			print ("created a particle of unknown type")






