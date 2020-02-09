import pyslha
import os
import matplotlib.pyplot as plt
import numpy as np

path = os.getcwd()

#todo:
# create more robust file aquitision 
# fill out self.sleptons, self.gauginos ... with cluster organization
# fully flesh out the catagories below


# Graph class wrapper around the slha file to be used
# when parsing data

cluster_thresh = 800
shift = .1
#pdg codes for various particles, grouped as they will be for the graph
higgs = [24, 25, 35, 36, 37]
sleptons = [1000011, 1000013, 1000015, 2000011, 2000013, 2000015, 1000012,
1000014, 1000016]
squarks = [1000001, 1000003, 1000005, 2000001, 2000003, 2000005,
1000002, 1000004, 1000006, 2000002, 2000004, 2000006]
gauginos = [1000021, 1000022, 1000023, 1000024, 1000025, 1000035,
1000037, 1000039]

higgs_anno = {24:"MW", 25:r"$h^0$", 35:r"$H^0$", 36:r"$A^0$", 37:r"$H^\pm$"}
slepton_anno = {1000011:r"$\widetilde{e}_{1}$", 1000013:r"$\widetilde{e}_{2}$", 1000015:r"$\widetilde{e}_{3}$", 2000011:r"$\widetilde{e}_{4}$", 2000013:r"$\widetilde{e}_{5}$", 2000015:r"$\widetilde{e}_{6}$", 1000012:r"$\widetilde{v}_{1}$",
1000014:r"$\widetilde{v}_{2}$", 1000016:r"$\widetilde{v}_{3}$"}
squark_anno = {1000001:r"$\widetilde{d}_{1}$", 1000003:r"$\widetilde{d}_{2}$", 1000005:r"$\widetilde{d}_{3}$", 2000001:r"$\widetilde{d}_{4}$", 2000003:r"$\widetilde{d}_{5}$", 2000005:r"$\widetilde{d}_{6}$",
1000002:r"$\widetilde{u}_{1}$", 1000004:r"$\widetilde{u}_{2}$", 1000006:r"$\widetilde{u}_{3}$", 2000002:r"$\widetilde{u}_{4}$", 2000004:r"$\widetilde{u}_{5}$", 2000006:r"$\widetilde{u}_{6}$"}
gaugino_anno = {1000021:r"$\widetilde{g}$", 1000022:r"$\widetilde{X}^0_1$", 1000023:r"$\widetilde{X}^0_2$", 1000024:r"$\widetilde{X}^\pm_1$", 1000025:r"$\widetilde{X}^0_3$", 1000035:r"$\widetilde{X}^0_4$",
1000037:r"$\widetilde{X}^\pm_2$", 1000039:r"$\widetilde{gr}$"}

cat_dict = {"higgs":higgs_anno,"slepton":slepton_anno,"squark":squark_anno,"gaugino":gaugino_anno}

# initial cluster function, used just on array objects and not the class ones
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

# cleaner implementation of cluster function 
def clusterFunc2(group):
	hld = []
	i = 0
	res = [k.mass for k in group]
	diff = [abs(i-j) for i, j in zip(res[1:],res[:-1])]
	while((i < len(group))):
		hld2 = [group[i]]

		while((i < len(diff))):
			if(diff[i] <= cluster_thresh):
				i+=1
				hld2.append(group[i])
			else:
				i+=1
				hld.append(hld2)
				break
		if(i == len(diff)):

			hld.append(hld2)
			i+=1

	return hld

def clusterFunc3(group):
	hld = []
	hs = []
	hss = []
	for i in group:
		hld.append(i.mass)
	diff = np.diff(hld)
	lss = []
	ls = []
	ls.append(group[0])
	hs.append(group[0].mass)
	for i, d in enumerate(diff):
		if d < cluster_thresh:
			ls.append(group[i+1])
			hs.append(group[i+1].mass)
			if(i == len(hld) - 2):
				lss.append(ls)
				hss.append(hs)
		else:
			lss.append(ls)
			hss.append(hs)
			ls = []
			hs = []
			ls.append(group[i+1])
			hs.append(group[i+1].mass)
			if(i == len(hld) - 2):
				lss.append(ls)
				hss.append(hs)


	return lss

	
# helper function to fix points returns labels for annotation
def fitcluster(clusters):
	anno = []
	parts = []
	for cluster in clusters:
		result_string = ""
		size = len(cluster)
		start = float(-1*(size - 1)/2*shift)
		cat = cat_dict[cluster[0].cat]
		i = 0

		for part in cluster:
			if(i == (len(cluster) - 1)):
				parts.append(part)

			result_string += cat[part.pdg]



			#print(part.pdg)
			result_string += " "
			part.delta = start
			start += shift
			i+=1
		anno.append(result_string)

	return parts, anno 
def fixX(group):
	for part in group:
		part.x = part.x + part.delta


class Graph:
	def __init__(self, file):
		self.file = pyslha.read(path + "/" + file)
		self.masses = self.file.blocks['MASS'].items()
		self.higgs = [Particle(i[0],i[1]) for i in self.masses if i[0] in higgs]
		self.sleptons = [Particle(i[0],i[1]) for i in self.masses if i[0] in sleptons]
		self.squarks = [Particle(i[0],i[1]) for i in self.masses if i[0] in squarks]
		self.gauginos = [Particle(i[0],i[1]) for i in self.masses if i[0] in gauginos]

		self.sleptons.sort(key=lambda x: x.mass)
		self.squarks.sort(key=lambda x: x.mass)
		self.gauginos.sort(key=lambda x: x.mass)
		self.higgs.sort(key=lambda x: x.mass)

	def orgCats(self):



		annotations = []
		annotated_particles = []

		higgs_annos = fitcluster(clusterFunc3(self.higgs))
		slepton_annos = fitcluster(clusterFunc3(self.sleptons))
		squark_annos = fitcluster(clusterFunc3(self.squarks))
		gaugino_annos = fitcluster(clusterFunc3(self.gauginos))

		
		annotations.append(higgs_annos[1])
		annotations.append(slepton_annos[1])
		annotations.append(squark_annos[1])
		annotations.append(gaugino_annos[1])
		annotated_particles.append(higgs_annos[0])
		annotated_particles.append(slepton_annos[0])
		annotated_particles.append(squark_annos[0])
		annotated_particles.append(gaugino_annos[0])


		fixX(self.higgs)
		fixX(self.sleptons)
		fixX(self.squarks)
		fixX(self.gauginos)
		return annotations, annotated_particles


	def plot(self,includes=['sleptons','higgs','gauginos','squarks']):
		print(includes)
		annotations, annotated_particles = self.orgCats(includes)

		x1 = [i.x for i in self.higgs]
		y1 = [i.mass for i in self.higgs]

		x2 = [i.x for i in self.sleptons]
		y2 = [i.mass for i in self.sleptons]

		x3 = [i.x for i in self.squarks]
		y3 = [i.mass for i in self.squarks]

		x4 = [i.x for i in self.gauginos]
		y4 = [i.mass for i in self.gauginos]

		xs = x1 + x2 + x3 + x4
		ys = y1 + y2 + y3 + y4

		colors = [int((i)*100/len(xs)) for i in range(len(xs))]

		plt.scatter(xs,ys,c=colors, cmap='hsv')




		for x,y in zip(annotations,annotated_particles):
			for i in range(len(x)):

				plt.annotate(x[i],
						xy=(y[i].x,y[i].mass+300),
						ha="center",
						fontsize=8)

		    #label = "{:.2f}".format(y)

		    #plt.annotate(label, # this is the text
		                 #(x,y), # this is the point to label
		                 #textcoords="offset points", # how to position the text
		                 #xytext=(5,0), # distance from text to points (x,y)
		                # ha='left') # horizontal alignment can be left, right or center


		plt.xlim(0,5)

		#plt.semilogy()
		plt.grid(alpha=.5,linestyle="--")
		plt.ylabel("Mass - GeV")

		plt.xticks([1,2,3,4])
		plt.axes().set_xticklabels(['higgs', 'sleptons', 'gauginos', 'squarks'])

		#plt.show()

	def plotBar(self):
		return 


	def plotSimple(self):
		annotations, annotated_particles = self.orgCats()
		x1 = [i.x for i in self.higgs]
		y1 = [i.mass for i in self.higgs]

		x2 = [i.x for i in self.sleptons]
		y2 = [i.mass for i in self.sleptons]

		x3 = [i.x for i in self.squarks]
		y3 = [i.mass for i in self.squarks]

		x4 = [i.x for i in self.gauginos]
		y4 = [i.mass for i in self.gauginos]

		xs = x1 + x2 + x3 + x4
		ys = y1 + y2 + y3 + y4

		colors = [int((i)*100/len(xs)) for i in range(len(xs))]

		plt.scatter(xs,ys,c=colors, cmap='hsv')




		for x,y in zip(annotations,annotated_particles):
			for i in range(len(x)):

				plt.annotate(x[i],
						xy=(y[i].x,y[i].mass+300),
						ha="center",
						fontsize=8)

		    #label = "{:.2f}".format(y)

		    #plt.annotate(label, # this is the text
		                 #(x,y), # this is the point to label
		                 #textcoords="offset points", # how to position the text
		                 #xytext=(5,0), # distance from text to points (x,y)
		                # ha='left') # horizontal alignment can be left, right or center


		plt.xlim(0,5)

		#plt.semilogy()
		plt.grid(alpha=.5,linestyle="--")
		plt.ylabel("Mass - GeV")

		plt.xticks([1,2,3,4])
		plt.axes().set_xticklabels(['higgs', 'sleptons', 'gauginos', 'squarks'])

		#plt.show()

class Particle:
	def __init__(self,pdg,mass):
		self.pdg = pdg
		self.mass = mass
		self.delta = 0

		if (int(pdg) in higgs):
			self.cat = "higgs"
			self.x = 1
		elif (int(pdg) in sleptons):
			self.cat = "slepton"
			self.x = 2
		elif (int(pdg) in squarks):
			self.cat = "squark"
			self.x = 4
		elif (int(pdg) in gauginos):
			self.cat = "gaugino"
			self.x = 3
		else:
			self.cat = "Na"
			print ("created a particle of unknown type")






