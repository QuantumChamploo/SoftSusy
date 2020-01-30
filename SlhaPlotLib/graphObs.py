import pyslha
import os
import matplotlib.pyplot as plt

path = os.getcwd()

#todo:
# create more robust file aquitision 
# fill out self.sleptons, self.gauginos ... with cluster organization
# fully flesh out the catagories below


# Graph class wrapper around the slha file to be used
# when parsing data

cluster_thresh = 300
shift = .1
#pdg codes for various particles, grouped as they will be for the graph
higgs = [24, 25, 35, 36, 37]
sleptons = [1000011, 1000013, 1000015, 2000011, 2000013, 2000015, 1000012,
1000014, 1000016]
squarks = [1000001, 1000003, 1000005, 2000001, 2000003, 2000005,
1000002, 1000004, 1000006, 2000002, 2000004, 2000006]
gauginos = [1000021, 1000022, 1000023, 1000024, 1000025, 1000035,
1000037, 1000039]

higgs_anno = {24:"MW", 25:"h0", 35:"H0", 36:"A0", 37:"H+"}
slepton_anno = {1000011:"~e_1", 1000013:"~e_2", 1000015:"~e_3", 2000011:"~e_4", 2000013:"~e_5", 2000015:"~e_6", 1000012:"~nu_1",
1000014:"~nu_2", 1000016:"~nu_3"}
squark_anno = {1000001:"~d_1", 1000003:"~d_2", 1000005:"~d_3", 2000001:"~d_4", 2000003:"~d_5", 2000005:"~d_6",
1000002:"~u_1", 1000004:"~u_2", 1000006:"~u_3", 2000002:"~u_4", 2000004:"~u_5", 2000006:"~u_6"}
gaugino_anno = {1000021:"~g", 1000022:"~n1", 1000023:"~n2", 1000024:"~c1", 1000025:"~n3", 1000035:"~n4",
1000037:"~c2", 1000039:"~gr"}

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
			print(group[i].pdg)
			hld.append(hld2)
			i+=1
	return hld
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
			if(part.pdg == 1000039):
				print ("look here!!!!")

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

	def orgCats(self):
		annotations = []
		annotated_particles = []

		higgs_annos = fitcluster(clusterFunc2(self.higgs))
		slepton_annos = fitcluster(clusterFunc2(self.sleptons))
		squark_annos = fitcluster(clusterFunc2(self.squarks))
		gaugino_annos = fitcluster(clusterFunc2(self.gauginos))
		annotations.append(higgs_annos[1])
		annotations.append(slepton_annos[1])
		annotations.append(squark_annos[1])
		annotations.append(gaugino_annos[1])
		annotated_particles.append(higgs_annos[0])
		annotated_particles.append(slepton_annos[0])
		annotated_particles.append(squark_annos[0])
		annotated_particles.append(gaugino_annos[0])

		print(gaugino_annos)
		fixX(self.higgs)
		fixX(self.sleptons)
		fixX(self.squarks)
		fixX(self.gauginos)
		return annotations, annotated_particles

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

		plt.scatter(xs,ys)




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
		plt.xticks([1,2,3,4])
		plt.axes().set_xticklabels(['higgs', 'sleptons', 'gauginos', 'squarks'])

		plt.show()

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






