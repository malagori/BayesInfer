__author__ = "Mehmood Alam Khan"
__email__  = "malagori@kth.se"
__version__= "0.9"
__license__ = "GPL"
__credits__ = ["Mehmood Alam Khan", "Pekka Parviainen"]

import os
import numpy as np
from mlabwrap import mlab
from bayesInfer.node import Node
from bayesInfer import matlab_lib
from bayesInfer.beneWrapper import BeneWrapper

class EquivalenceClass(object):
    def __init__(self):
        self.matlabLibPath= os.path.dirname(os.path.abspath(matlab_lib.__file__))
        mlab.addpath(self.matlabLibPath)
        
    def getOptDag(self, vdFile, dataFile, outDirectory, totalVaiables):
        
        #create obj of BeneWrapper class and initialize the fields
        bwObj= BeneWrapper(vdFile, dataFile, score=1.0, outDirectory, totalVaiables) 
        
        bwObj.generateOptBnt()
        
        optDag = bwObj.readBeneBnt()
        
        optDag= np.array(optDag).astype('int')
        
        return optDag
        
    def generateCdag(self,optDag):
        """
        this function partial dag from dag represented in 2d numpy array.
        input: dag in 2d numpy array
        output: cdag in 2d numpy array
        
        """
        cDag= mlab.dag_to_cpdag(optDag)
        return cDag
    
    def generateBnt(self,dag):
        """
        this function generate populate network from dag represented in list of list.
        input: dagListofList
        output: allNodeObj (type=dictonary)
        
        """
        
        varName=0
        allNodesObj={}
        
        for i in dag:
            
            varName +=1
            cardinality=2
            parentSet= [j+1 for j in range(0,len(i)) if i[j]==1] # parent name starts from 1 not 0
            
            node= Node()
            node.setName(varName) 
            node.setR(int(cardinality)) # can update cardinality from vdFile
            node.setKvalues(dict.fromkeys(list(range(0, int(cardinality), 1))))
            node.setParents(parentSet)
            allNodesObj[varName]= node
            
        return allNodesObj
        
    def getAllDagsInPdag(self, cDag):
        """
        This function will generate all the dags in equvilance class of cDag and populate each bnet
        input: cDag
        output: dagsDict(index= int), allDagsNetworkDict(index= int)
        """
        
        dagsDict= {} # each dag is represented as list of list
        allDagsNetworkDict= {} # network related to each dag is populated 
        nDags, Dag_list = mlab.pdag_to_all_dags(cDag, nout=2)
        
        for i in range(0, nDags):
            npDagsArray= mlab.cell2mat(Dag_list[i]).astype(int)
            dagsDict[i]= npDagsArray.tolist()
            allDagsNetworkDict[i]= self.generateBnt(dagsDict[i])
        
        return dagsDict, allDagsNetworkDict
    
        
        
    
