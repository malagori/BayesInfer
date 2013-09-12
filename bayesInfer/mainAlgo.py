from __future__ import division
'''
Created on Sep 10, 2013

@author: malagori
'''


import time
import datetime
import random as rNumber


from bayesInfer.storeRetriveSeed import RandomSeed
from bayesInfer.readDataFile import convertBeneDataFile
from bayesInfer.readDataFile import readVdFile
from bayesInfer.equivalenceClass import EquivalenceClass
from bayesInfer.node import Node
from bayesInfer.bdeuClass import BDeuClass

class MainAlgo(object):
    '''
    This class contains the main workflow of our algorithm
    '''


    def __init__(self, vdFile, dataFile, outdir, alpha, seed, steepestAsent,iterations, seedFile, outPutFile  ):
        '''
        Constructor
        '''
        self.vdFile         = vdFile
        self.dataFile       = dataFile
        self.outdir         = outdir
        self.alpha          = alpha
        self.seed           = seed
        self.steepestAsent  = steepestAsent
        self.seedFile       = seedFile
        self.outputFile     = outPutFile
        self.iterations     = iterations
        self.df
        
    def printDag(self,iterations, allNodesObjects):
        with open(str("bestDag_"+iterations), 'w') as wdag:
            p=[]
            for key, node in allNodesObjects.iteritems():
                p=[]
                p=sum(node.getParents())
                if not p: # if p is empty list 
                    p.append(0)
                wdag.write("%s\n" % p)

                    
            
            
            
        
    def runAlgo(self):
        
        HIDDEN_NAME=-1
        
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
        self.df=convertBeneDataFile(self.dataFile)
        
        # read vdFile
        variableNames, cardinality= readVdFile(self.vdFile)
        
        # create object of EquivalenceClass
        objEC= EquivalenceClass()
        # get the opt bnt from bene
        optDag, allNodesObj= objEC.getOptDag(self.vdFile, self.dataFile, self.alpha, self.outdir, len(variableNames))
        
        HIDDEN_NAME= len(variableNames) +1
        
        # dict of dict
        cachedBDeuDict={} # keys: tuple ((parent's parent),( child's parent)) ; Values: bdeu score for the network
        edgesDict={} # keys: edge tuple (parent, child); Values: keys of Pa_C_PaPa_CPa dictionary
        hiddenNodesDict={}
        iterations=0
        totalPreviousBDeuScore= float("-inf")
        totalCurrentBDeuScore= float("-inf")
        currentMaxBDeu= float("-inf")
        previousMaxBDeu= float("-inf")
        
        
        with open(self.outputFile, 'w') as wf:
            
            totalUniqueObservations= self.df.shape[0]
            objCBDeu= BDeuClass(self.df, allNodesObj, totalUniqueObservations)

            # update the parent configurations for all variables
            # and the counts associated with the each parent configuration for each value of X
            for n in objCBDeu.allNodeObjects:
                objCBDeu.getUpdatedQi(objCBDeu.allNodeObjects[n])
                objCBDeu.populateCounts(objCBDeu.allNodeObjects[n])
            # find the BDeu Score for the whole structure
            for n in objCBDeu.allNodeObjects:
                objCBDeu.nodesBDeuScore.append(objCBDeu.getBDeu(objCBDeu.allNodeObjects[n], self.alpha))
            
            
            print "BDeu Score for optimal dag from Bene: %f" % sum(objCBDeu.nodesBDeuScore)
            #print initial data frame 
            self.df.to_csv('initialDF_bene'+'.csv', sep=',')
            # print the state for the random number generator
            stateOutFile= 'state_iter_'+str(iterations)+'_initialSeed_'+ str(self.seed) +'_'+self.outputFile
            rs.storeSate(stateOutFile)
            wf.write("BDeuScore for optimal dag from Bene, %f" % sum(objCBDeu.nodesBDeuScore))
            
            # Repeat until adding a hidden variable does not increase the score
            while True:
            
                
                #increment the iteration number
                iterations +=1 
            
                # pDag
                cDag=objEC.generateCdag(optDag)
                # generate all dags in pDag
                dagsDict, allDagsNetworkDict= objEC.getAllDagsInPdag(cDag, cardinality)
                           
                totalUniqueObservations= self.df.shape[0] # if we introduce next hidden variable, this variable would be updated
        
                arrayListBDeuClassObjs=[]
         
                for id, dag in dagsDict.iteritems():
                    
                    allNodeObjects=allDagsNetworkDict[id]
                    
                    # instantiate CalculateBDeuClass's object 
                    objCBDeu= BDeuClass(self.df, allNodeObjects, totalUniqueObservations)
                    
                    # compute initial bdeu score before adding any hidden variable
                    # update the parent configurations for all variables
                    # and the counts associated with the each parent configuration for each value of X
                    for n in objCBDeu.allNodeObjects:
                        objCBDeu.getUpdatedQi(objCBDeu.allNodeObjects[n])
                        objCBDeu.populateCounts(objCBDeu.allNodeObjects[n])
                    # find the BDeu Score for the whole structure
                    for n in objCBDeu.allNodeObjects:
                        objCBDeu.nodesBDeuScore.append(objCBDeu.getBDeu(allNodeObjects[n], self.alpha))
                    
                    totalPreviousBDeuScore=sum(objCBDeu.nodesBDeuScore)
                    initialBDeuScore= totalPreviousBDeuScore
                    
                    print "BDeu Score for a dag in Equivalence class before adding hidden variable: %f" % totalPreviousBDeuScore
                    
                    # find all the edges in the dag
                    idx=1
                    edges=[]
                    for i in dag:
                        e=[edges.append((idx,j+1)) for j in range(0,len(i)) if i[j]==1]
                        idx+=1
                    
                    # for each edge. replace it with hidden variable and compute the score
                    for edge in edges:
                        # check if score for adding hidden variable for this edge is already computed in other equivalence dag
                        
                        parentNode= Node()
                        childNode= Node()
                        
                        parentNode= objCBDeu.allNodeObjects[edge[0]]
                        childNode= objCBDeu.allNodeObjects[edge[1]]
                        
                        # creating key for cachedBDeuDict dictionary
                        key=tuple([tuple(parentNode.getParents()), tuple(childNode.getParents())])
                        

                        if edge in edgesDict.keys():
                            if key == edgesDict[edge]: # if true do not add hidden variable 
                                # get the score from the cachedBDeu score for this edge after being h is added to the network
                                # check if the difference is positive add the hidden variable and add the difference to the netowork bdeu score
                                # else donot add hidden variable
                                bdeuDiffScore= cachedBDeuDict[key]
                                if (totalPreviousBDeuScore + bdeuDiffScore) > totalPreviousBDeuScore:
                                    # add hidden variable to the network
                                    objCBDeu.allNodeObjects.append( hiddenNodesDict[edge] )
                                    totalPreviousBDeuScore= totalPreviousBDeuScore + bdeuDiffScore
                                    objCBDeu.dagBDeuScore= totalPreviousBDeuScore
                                    print "hidden variable %s added from cache" % hiddenNodesDict[edge].getName()
                                else:
                                    print "bdeu diff score is not greater then current BDeuScore "
                        else:
                            # copy allnodes
                            tmpAllNodesObj      = objCBDeu.allNodeObjects
                            tmpDF               = objCBDeu.df.copy()
                            tmpDagBDeuScore     = objCBDeu.dagBDeuScore
                            
                            # add hidden variable to the network
                            h=objCBDeu.addHiddenNode(HIDDEN_NAME, 2 , parentNode.getName(), childNode.getName())
                            
                            # split the dataframe counts
                            print "data frame before adding hidden variable"
                            print objCBDeu.df
                            objCBDeu.percentageHiddenCoutsSplit(h)
                            print "data frame after adding hidden variable"
                            print objCBDeu.df
                            
                            # write df to file called initialCountSplit.txt
                            #outName= self.outputFile+'_initialCountSplit_'+str((datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d-h%H-m%M-s%S')))
                            #newDF.to_csv(outName+'.csv', sep=',')
                            
                            # populate hidden value counts
                            objCBDeu.populateCounts(h)
                            
                            if self.steepestAsent == True:
                                
                                
                                rRecords=[i for i in xrange(0, objCBDeu.df.shape[0]-1)]
                                # randomly shuffle the indexes 
                                rNumber.shuffle(rRecords)
                                for i in rRecords:
                                    baseDFrame = objCBDeu.df.copy()
                                    #print "baseDFrame"
                                    #print baseDFrame
                                    # below loop is for the counts increment and decrement
                                    for flag in [True, False]:
                                        
                                        for j in rRecords:
                                                                     
                                            # perturb the counts in 
                                            objCBDeu.countPerturbation(h,j, incrementFlag=flag)
                                            
                                            nodesBDeuScore=[]
                                            # compute the BDeu score again after perturbations
                                            for n in objCBDeu.allNodeObjects:
                                                node=objCBDeu.allNodeObjects[n]
                                        
                                                if node.getParentUpdateFlag() == True or node.getChildrenUpdateFlag() == True: # if true its a child of hidden variable. so, calculate BDeu again
                                                    
                                                    # change the counts of this node according to some criteria i.e.
                                                    # for new parent configuration, assign the counts such that sum of counts of new parent
                                                    # configurations is equal to counts of old parent configuration which we split
                                                    # to get the new parent configuration. 
                                                    objCBDeu.populateCounts(node)
                                                    node.setLocalBDeu(objCBDeu.getBDeu(node, self.alpha))
                                                nodesBDeuScore.append(node.getLocalBDeu())
                                            totalCurrentBDeuScore= sum(nodesBDeuScore)
                                            
                                            if totalPreviousBDeuScore < totalCurrentBDeuScore:
                                                baseDFrame = objCBDeu.df.copy()
                                                totalPreviousBDeuScore = totalCurrentBDeuScore
                                        
                                        # copy the base dataframe to to data frame since the new base datafame is the one with high BDeuScore
                                        objCBDeu.df = baseDFrame.copy()
                            
                            if initialBDeuScore < totalCurrentBDeuScore:
                                # add hidden node to the dictionary
                                hiddenNodesDict[i]=h
                                diffBDeu= totalCurrentBDeuScore - initialBDeuScore
                                cachedBDeuDict[key]= diffBDeu
                                edgesDict[i]= key
                                # generate new name for hidden variable
                                HIDDEN_NAME += 1
                                objCBDeu.dagBDeuScore= totalCurrentBDeuScore
                                initialBDeuScore = totalCurrentBDeuScore
                                
                            else: # adding hidden variable didn't improve score, so go back to old state                              
                                objCBDeu.allNodeObjects= tmpAllNodesObj
                                objCBDeu.df= tmpDF
                                objCBDeu.dagBDeuScore=tmpDagBDeuScore
                                    
                    # store BDeu Class object
                    arrayListBDeuClassObjs.append(objCBDeu)            
                # find the Dag' with higest bdeu score and input it to find the equivalence dags for it and repeat the whole process
                currentMaxBDeu=float("-inf")
                for obj in arrayListBDeuClassObjs:
                    if currentMaxBDeu < obj.dagBDeuScore:
                        currentMaxBDeu                 = obj.dagBDeuScore
                        currentMaxAllNodesObjects      = obj.allNodeObjects
                        currentMaxDF                   = obj.df
                
                # check the looping condition
                if previousMaxBDeu < currentMaxBDeu:
                    previousMaxBDeu=currentMaxBDeu
                    #print optdag
                    self.printDag(currentMaxAllNodesObjects)
                    #print hidden counts and bdeu score for the dag with higest bdeu score in equivalance class
                    print "Iteration: %d , BDeu Score: %f" % (iterations, totalCurrentBDeuScore)
                    hValues= h.getKvalues().keys()
                    for i in xrange(0,len(hValues)-1):
                        count=currentMaxDF[currentMaxDF[h.getName()]==hValues[i]].Counts
                        for j in count:
                            wf.write(str(j)+',')
                        del count
                    wf.write(str(currentMaxBDeu) + "\n")
                    # print the state for the random number generator
                    stateOutFile= 'state_iter_'+str(iterations)+'_initialSeed_'+ str(self.seed) +'_'+self.outputFile
                    rs.storeSate(stateOutFile)
                else: 
                    break
                
                
                
                
                
                
        
        
        