import itertools

""" 
	Node class contains information about the node in the network
"""

class Node(object):
	def __init__(self):
		self.name= ''
		# k_values is of the form {0:4, 1:2,...}
		self.k_values= 	{}
		self.r = 		0
		self.parentUpdateFlag= True
		self.parents= 	[]
		self.pConfigurations={}
		self.localBDeu=0
		#self.k_counts= 	kCounts
	# getters	
	def getR(self):
		return self.r
	def getName(self):
		return self.name
	def getKvalues(self):
		return self.k_values
	def getParents(self):
		return self.parents
	def getPaConfigurations(self):
		return  self.pConfigurations
	def getLocalBDeu(self):
		return self.localBDeu
	
	# setters
	def setLocalBDeu(self,bdeuScore):
		self.localBDeu= bdeuScore
	def setParents(self,pList):
		self.parents=pList
	def setName(self,name):
		self.name= name
	def setKvalues(self,kdict):
		self.k_values=dict.fromkeys(kdict)
	def setR(self,r):
		self.r= r
	def setpConfiguration(self,pConfigDict):
		self.pConfigurations=pConfigDict
			

		
		
