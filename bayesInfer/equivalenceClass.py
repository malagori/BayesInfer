__author__ = "Mehmood Alam Khan"
__email__  = "malagori@kth.se"
__version__= "0.9"
__license__ = "GPL"
__credits__ = ["Mehmood Alam Khan", "Pekka Parviainen"]

import os
import numpy as np
from mlabwrap import mlab
from node import Node
import matlab_lib
from beneWrapper import BeneWrapper

class EquivalenceClass(object):
    
    def __init__(self):
        self.matlabLibPath=os.path.abspath('matlab_lib')
        #self.matlabLibPath= os.path.dirname(os.path.abspath('matlab_lib/__init__.pyc'))
        #print "os.path.abspath('matlab_lib/__init__.pyc') %s" % os.path.abspath('matlab_lib/__init__.pyc')
        #print "self.matlabLibPath: %s" % self.matlabLibPath
        #self.matlabLibPath= '/home/mehmood/glob/projects/bayesian/code/BayesInfer/bayesInfer/matlab_lib'
        #print "self.matlabLibPath: %s" % self.matlabLibPath
        mlab.addpath(self.matlabLibPath)
          
    #def __exit__(self,type, value, traceBack):
        #mlab.close()
        
        
    def getOptDag(self, vdFile, dataFile, score, outDirectory, totalVaiables, cardinality):
        
        
        #create obj of BeneWrapper class and initialize the fields
        bwObj= BeneWrapper(vdFile, dataFile, score, outDirectory, totalVaiables) 
        
        bwObj.generateOptBnt()
        
        optDag, allNodesObj = bwObj.readBeneBnt(cardinality)
        
        optDag= np.array(optDag).astype('int')
        
        return optDag, allNodesObj
        
    def generateCdag(self,optDag):
        """
        this function partial dag from dag represented in 2d numpy array.
        input: dag in 2d numpy array
        output: cdag in 2d numpy array
        
        """
   
        cDag= mlab.dag_to_cpdag(optDag)
        return cDag
    
    def generateBnt(self,dag, cardinality):
        """
        this function generate populate network from dag represented in list of list.
        input: dagListofList, cardinality for each variable [list]
        output: allNodeObj (type=dictonary)
        
        """
        
        varName=0
        allNodesObj={}
        
        print dag
        # dag's traspose to get the parent set
        myDag=map(list, zip(*dag))

        
        for i in myDag:
            
            parents=[]
            idx=1
            for j in i:
            #for j in reversed(i):
                if j == 1:
                    parents.append(idx)          
                idx+=1# parent name starts from 1 not 0
            
            node= Node()
            varName +=1
            node.setName(varName)
            node.setR(int(cardinality[varName-1])) # can update cardinality from vdFile
            node.setKvalues(dict.fromkeys(list(xrange(0, int(cardinality[varName-1]), 1))))
            node.setParents(parents)
            allNodesObj[varName]= node
            
        return allNodesObj
        
    def getAllDagsInPdag(self, cDag, cardinality):
        """
        This function will generate all the dags in equvilance class of cDag and populate each bnet
        input: cDag, cardinality for each variable [list]
        output: dagsDict(index= int), allDagsNetworkDict(index= int)
        """
  
        
        dagsDict= {} # each dag is represented as list of list
        allDagsNetworkDict= {} # network related to each dag is populated 
        signedPdag= mlab.pdag_unsigned_to_signed(cDag)
        nDags, Dag_list = mlab.pdag_to_all_dags(signedPdag, nout=2)
        
        for i in xrange(0, nDags):
            npDagsArray= mlab.cell2mat(Dag_list[i]).astype(int)
            dagsDict[i]= npDagsArray.tolist()
            allDagsNetworkDict[i]= self.generateBnt(dagsDict[i], cardinality)
        
        return dagsDict, allDagsNetworkDict
