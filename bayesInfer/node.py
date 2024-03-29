
__author__ = "Mehmood Alam Khan"
__email__  = "malagori@kth.se"
__version__= "0.9"
__license__ = "GPL"
__credits__ = ["Mehmood Alam Khan", "Pekka Parviainen"]



""" 
	Node class contains information about the node in the network
"""

class Node(object):
	def __init__(self):
		self.name= ''
		# k_values is of the form {0:[(0,0,0),(0,0,1),..,(1,1,1)], 1:[(0,0,0),(0,0,1),..,(1,1,1)],...}
		# or simply k_values: {0:{pconfigVarValueCount}, 1: {pconfigVarValueCount}, 2: {pconfigVarValueCount}, ...}
		# keys= values of varaible, 
		# while values corresponds to the dictionary whose keys are all parent configurations and 
		# whose values are data-counts for that particular parent configuration. see details below
		# where pconfigVarValueCount is dictionary for parent's configuration:
		# keys= 	[(0,0,0),(0,0,1),...]
		# values= 	[2, 6, ....]
		
		# NOTE: keys of k_values dictionary need to be ascending order other 
		self.k_values= 	{}
		
		self.valueUpdateFlag= False
		self.r = 		0
		self.parentUpdateFlag= False
		self.parents= 	[]
		# pConfigurations: is a tuples-array whose values are different parent configuration of variable X 
		#[(0, 0, 0), (0, 0, 1), (0, 1, 0), (0, 1, 1), (0, 2, 0), (0, 2, 1), (1, 0, 0), (1, 0, 1), (1, 1, 0), (1, 1, 1), (1, 2, 0), (1, 2, 1), (0, 0, 0), (0, 0, 1), (0, 1, 0), (0, 1, 1), (0, 2, 0), (0, 2, 1), (1, 0, 0), (1, 0, 1), (1, 1, 0), (1, 1, 1), (1, 2, 0), (1, 2, 1)]
		self.pConfigurations=[]
		# dictionary for parent:variable's_value_count
		# keys= 	[(0,0,0),(0,0,1),...]
		# values= 	[2, 6, ....]
		#self.pconfigVarValueCount={} # this is extra info, need to be removed
		self.localBDeu=0.0
		self.childrenUpdateFlag=False
		self.children=[]
		# parent_k_counts is used when a node donot have any parent but it has children
		self.parent_k_counts= 	{}
	
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
	def getParentValueCount(self):
		return self.parent_k_counts
	def getChildrenUpdateFlag(self):
		return self.childrenUpdateFlag
	def getParentUpdateFlag(self):
		return self.parentUpdateFlag
	def getChildren(self):
		return self.children
	#def getsetpconfigVarValueCount(self):
	#	return self.pconfigVarValueCount
	
	# setters
	def addChild(self,child):
		self.children.append(child)
	def setParentUpdateFlag(self,flag):
		self.parentUpdateFlag= flag
	def setChildrenUpdateFlag(self, flag):
		self.childrenUpdateFlag=flag
	def setLocalBDeu(self,bdeuScore):
		self.localBDeu= bdeuScore
	def setParents(self,pList):
		self.parents=pList
	def addParent(self, singleParent):
		self.parents.append(singleParent)
	def setName(self,name):
		self.name= name
	def setKvalues(self,kdict):
		self.k_values=kdict
	def setR(self,r):
		self.r= r
	def setParentValueCounts(self, paValDict):
		self.parent_k_counts= paValDict
	def setpConfiguration(self,pConfigTupleArray):
		self.pConfigurations=pConfigTupleArray
	#def setpconfigVarValueCount(self,pConfigDict):
	#	self.pconfigVarValueCount= pConfigDict
			

		
		
