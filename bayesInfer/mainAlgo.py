'''
Created on Sep 10, 2013

@author: malagori
'''
from __future__ import division

import time
import random as rNumber


from bayesInfer.storeRetriveSeed import RandomSeed
from bayesInfer.readDataFile import convertBeneDataFile
from bayesInfer.readDataFile import readVdFile
from bayesInfer.equivalenceClass import EquivalenceClass
from bayesInfer.node import Node
from bayesInfer.BDeuClass import BDeuClass

class MainAlgo(object):
    '''
    This class contains the main workflow of our algorithm
    '''


    def __init__(self, vdFile, dataFile, outdir, alpha, seed, steepestAsent, seedFile  ):
        '''
        Constructor
        '''
        self.vdFile          = vdFile
        self.dataFile        = dataFile
        self.outdir          = outdir
        self.alpha           = alpha
        self.seed            = seed
        self.steepestAsent   = steepestAsent
        self.seedFile        = seedFile
        
    def checkHiddenScore(self, edgeTuple, ):
        print "hi"
        
        
    def runAlgo(self):
        
        # instantiating RandomSeed object
        rs=RandomSeed()
        
        if self.seed == None and self.seedFile == None:
            self.seed= time.time()
        elif self.seed != None and self.seedFile == None:
            rs.setInitialState(self.seed)
        elif self.seedFile != None and self.seed == None:
            state= rs.getSateFromFile(self.seedFile)
            rNumber.setstate(state)
            
        # read data file
        df=convertBeneDataFile(self.dataFile)
        totalUniqueObservations= df.shape[0] # if we introduce next hidden variable, this variable would be updated
        
        # read vdFile
        variableNames, cardinality= readVdFile(self.vdFile)
        
        # create object of EquivalenceClass
        objEC= EquivalenceClass()
        # get the opt bnt from bene
        optDag= objEC.getOptDag(self.vdFile, self.dataFile, self.alpha, self.outdir, len(variableNames))
        
        # Repeat until adding a hidden variable does not increase the score
        while True:
        
            # pDag
            cDag=objEC.generateCdag(optDag)
            # generate all dags in pDag
            dagsDict, allDagsNetworkDict= objEC.getAllDagsInPdag(cDag, cardinality)
            
            # dict of dict
            cachedBDeuDict={} # keys: tuple ((parent's parent),( child's parent)) ; Values: bdeu score for the network
            edgesDict={} # keys: edge tuple (parent, child); Values: keys of Pa_C_PaPa_CPa dictionary
            
            for id, dag in dagsDict.iteritems():
                idx=1
                edges=[]
                allNodeObjects=allDagsNetworkDict[id]
                
                # instantiate CalculateBDeuClass's object 
                objCBDeu= BDeuClass(df, allNodeObjects, totalUniqueObservations)
                
                # compute initial bdeu score before adding any hidden variable
                # update the parent configurations for all variables
                # and the counts associated with the each parent configuration for each value of X
                for n in allNodeObjects:
                    objCBDeu.getUpdatedQi(allNodeObjects[n])
                    objCBDeu.populateCounts(allNodeObjects[n])
                # find the BDeu Score for the whole structure
                for n in allNodeObjects:
                    objCBDeu.nodesBDeuScore.append(objCBDeu.getBDeu(allNodeObjects[n], self.alpha))
                    
                print "BDeu Score for before adding hidden variable: %f" % sum(objCBDeu.nodesBDeuScore)
                
                for i in dag:
                    e=[edges.append((idx,j+1)) for j in range(0,len(i)) if i[j]==1]
                    idx+=1
                for i in edges:
                    # check if score for adding hidden variable at this edge is already computed in other equvilance dag
                    parentNode= Node()
                    childNode= Node()
                    
                    parentNode= objCBDeu.allNodeObjects[i[0]]
                    childNode= objCBDeu.allNodeObjects[i[1]]
                    
                    key=tuple([tuple(parentNode.getParents()), tuple(childNode.getParents())])
                    
                    
                    
                    
                    
                
        
        
        