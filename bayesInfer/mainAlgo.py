from __future__ import division

import time
import random as rNumber
import pandas as pd
import cProfile

from math import exp
from storeRetriveSeed import RandomSeed
from readDataFile import convertBeneDataFile
from readDataFile import readVdFile
from equivalenceClass import EquivalenceClass
from node import Node
from bdeu import BDeuClass

class MainAlgo(object):
    '''
    This class contains the main workflow of our algorithm
    '''
    def __init__(self, vdFile, dataFile, outdir, alpha, seed, steepestAsent,iterations, seedFile, outPutFile, simulatedAnealingFlag, decrementValue, mlabPath  ):
        '''
        Constructor
        '''
        self.mlabPath       = mlabPath
        self.vdFile         = vdFile
        self.dataFile       = dataFile
        self.outdir         = outdir
        self.alpha          = alpha
        self.seed           = seed
        self.steepestAsent  = steepestAsent
        self.seedFile       = seedFile
        self.outputFile     = outPutFile
        self.iterations     = iterations
        self.simAnealFlag   = simulatedAnealingFlag
        self.decrementValue = decrementValue
        self.df             = pd.DataFrame(index=None, columns=None)
        
    def printDag(self,iterations, allNodesObjects):
        '''
            this function return optDag after improvement in bdeu score either by adding hidden var or removing edges
        '''
        optDag=[]
        cardinality=[]
        with open(str("bestDag_"+iterations), 'w') as wdag:
            for key, node in allNodesObjects.iteritems():
                cardinality.append(node.getR())
                p=[0]*len(allNodesObjects)
                pa=[]
                pa=node.getParents()
                if not pa: # if p is empty list 
                    wdag.write("%d\n" % 0)
                    optDag.append(p)
                    continue
                else:
                    for i in pa:
                        p[i-1]=1
                    optDag.append(p)
                    p=p[::-1]
                    binStrParent=''.join(str(e) for e in p)
                    deciParents=int(binStrParent,2)
                    wdag.write("%s\n" % str(deciParents))
        
        optDag=map(list, zip(*optDag))
        
        return optDag, cardinality

    def computeBDeuUsingSteepestAsent(self,h ,objCBDeu, totalPreviousBDeuScore, sIndex, iterations, outFile):
        '''
        This function do the steepest asent to compute local optimum bdeu score.
        '''
        
        rRecords=[i for i in xrange(0, objCBDeu.df.shape[0]-1)]
        # randomly shuffle the indexes 
        #rNumber.shuffle(rRecords)
        temp                = rRecords[0]
        rRecords[0]         = rRecords[sIndex]
        rRecords[sIndex]    = temp
        objCBDeuBestCopy    = objCBDeu
        with open(outFile, 'w') as wf:
            for i in range(iterations):
                
                # below loop is for the counts increment and decrement
                for j in rRecords:
                    
                    for flag in [True, False]:                         
                        # perturb the counts in 
                        objCBDeu.countPerturbation(h, j, self.decrementValue, incrementFlag=flag )
                        
                        nodesBDeuScore=[]
                        # compute the BDeu score again after perturbations
                        for n in objCBDeu.allNodeObjects:
                            node=objCBDeu.allNodeObjects[n]
                            if node.getParentUpdateFlag() == True or node.getChildrenUpdateFlag() == True: # if true its a child of hidden variable. so, calculate BDeu again
                                
                                # change the counts of this node according to some criteria i.e.
                                # for new parent configuration, assign the counts such that sum of counts of new parent
                                # configurations is equal to counts of old parent configuration which we split
                                # to get the new parent configuration. 
                                objCBDeu.getUpdatedQi(node)
                                objCBDeu.populateCounts(node)
                                node.setLocalBDeu(objCBDeu.getBDeu(node, self.alpha))
                            nodesBDeuScore.append(node.getLocalBDeu())
                        
                        totalCurrentBDeuScore= sum(nodesBDeuScore)
                        
                        if totalPreviousBDeuScore <= totalCurrentBDeuScore:
                           
                            objCBDeu.dagBDeuScore   = totalCurrentBDeuScore
                            objCBDeuBestCopy        = objCBDeu
                            objCBDeuBestCopy.dagBDeuScore= totalCurrentBDeuScore
                            totalPreviousBDeuScore  = totalCurrentBDeuScore
                            
                        else:
                            objCBDeu    = objCBDeuBestCopy
                            
                        wf.write("Best bdeuscore: %f, Current bdeuscore: %f \n" % (totalPreviousBDeuScore, totalCurrentBDeuScore))
                        #print "Best bdeuscore: %f, Current bdeuscore: %f" % (totalPreviousBDeuScore, totalCurrentBDeuScore)
       
        return objCBDeuBestCopy
        
    def removeEdgesFromBnt(self, edges, previousBDeuScore, objCBDeu):
        '''
         This function removes edges from the bayesian network every time there is an improvement in the bdeu score
        '''
        # keep a copy of the objCBDeu
        objCBDeuCopy        = objCBDeu
        currentBDeuScore    = float("-inf")
        maxObjCBDeu         = objCBDeu
        
        for edge in edges:
            
            childNode= Node()
            
            childNode= objCBDeu.allNodeObjects[edge[1]]
            newParents= childNode.getParents()
            newParents.remove(edge[0]) # remove the parent from the parent set of child node
            childNode.setParents(newParents)
            objCBDeu.getUpdatedQi(childNode)
            objCBDeu.populateCounts(childNode)
            childNode.setLocalBDeu(objCBDeu.getBDeu(childNode, self.alpha))
            
            objCBDeu.allNodeObjects[edge[1]]= childNode
            
            for n in objCBDeu.allNodeObjects:
                node= Node()
                node= objCBDeu.allNodeObjects[n]
                currentBDeuScore+=node.getLocalBDeu()
                
            if previousBDeuScore <= currentBDeuScore:
                print "removing edge from bnt improves the bdeu scores"
                previousBDeuScore   = currentBDeuScore
                objCBDeu.dagBDeuScore= currentBDeuScore
                maxObjCBDeu         = objCBDeu
                maxObjCBDeu.dagBDeuScore= currentBDeuScore
            else:
                objCBDeu= objCBDeuCopy
 
        return maxObjCBDeu
                
    def temperature(self, k, kmax):
        """
            This function implements the temperature function of simulated anealing algorithm
        """
        temp= (1 - k/ float(kmax))
        return temp
    
    def probAcceptance(self, e, enew, T):
        """
            this function compute the probability of acceptance 
        """
        prob=0.0
        if (e < enew):
            prob=1.0
            #print "accept with prob = 1"
        else:
            prob= exp(-( -enew + e )/ float(T))
            #print "e : %f  enew: %f" % (e, enew)
        return prob
    
        
    def simulatedAnealing(self, objCBDeu, hiddenVar, previousScore, sIndex, iterations, outFile ):
        """
            This function implements the simulated Anealing algorithm (wiki) 
        """
        
        e               = previousScore                               # Initial state, energy.
        emax            = float('-inf') 
        ebest           = e                                     # Initial "best" solution
        k               = 1                                     # Energy evaluation count.
        kmax            = iterations
        objCBDeuBestState= objCBDeu
        objCBDeuOldState = objCBDeu
        j               = sIndex
        
        with open(outFile, 'w') as wf:
            
            while k < kmax and e > emax:                    # While time left & not good enough
                T =    self.temperature(k, kmax)              # Temperature calculation.
                
                # randomly choose hidden state zero or one
                num= rNumber.randint(0,1) 
                if num == 0:
                    flag = False
                else:
                    flag = True

                objCBDeu.countPerturbation(hiddenVar, j, self.decrementValue, incrementFlag=flag)     
                
                j=rNumber.randint(0, objCBDeu.df.shape[0]-1) # randomly select another record for next iteration
                
                nodesBDeuScore= []
                # compute the BDeu score again after perturbations
                for n in objCBDeu.allNodeObjects:
                    node=objCBDeu.allNodeObjects[n]
                    if node.getParentUpdateFlag() == True or node.getChildrenUpdateFlag() == True: # if true its a child of hidden variable. so, calculate BDeu again
                        objCBDeu.populateCounts(node)
                        node.setLocalBDeu(objCBDeu.getBDeu(node, self.alpha))
                    nodesBDeuScore.append(node.getLocalBDeu())
#                for n in objCBDeu.allNodeObjects:               # populate the counts for each node
#                    tempNode    = Node()
#                    tempNode    = objCBDeu.allNodeObjects[n]
#                    objCBDeu.getUpdatedQi(tempNode)
#                    objCBDeu.populateCounts(tempNode)
#                
#                for n in objCBDeu.allNodeObjects:
#                    tempNode    = Node()
#                    tempNode    = objCBDeu.allNodeObjects[n]
#                    tempScore   = objCBDeu.getBDeu(objCBDeu.allNodeObjects[n], self.alpha)
#                    nodesBDeuScore.append(tempScore)
                objCBDeu.dagBDeuScore= sum(nodesBDeuScore)
                
                enew = objCBDeu.dagBDeuScore                              # Compute its energy.
                #NOTE Inverse logic here using  '<' instead of '>' as in org algo
                acceptprob= self.probAcceptance(e, enew, T)
                if acceptprob < rNumber.random():# reject the current state 
                    objCBDeu= objCBDeuOldState          # go back to the old state
                else:  # accept the new state
                    objCBDeuOldState= objCBDeu
                    e               = enew
                                                        
                if enew > ebest:                              # Is this a new best?
                    objCBDeuBestState= objCBDeu
                    ebest = enew                              # Save 'new neighbour' to 'best found'.
                k = k + 1
                #print "--->iteration  %d " % k                                     # One more evaluation done
                #print "Best bdeuscore: %f and Current bdeuscore %f :" % (ebest, enew)
                wf.write("Best bdeuscore: %f, Current bdeuscore: %f, proposal bdeuscore: %f  , temp: %f, prob: %f\n" % (ebest, e, enew,T, acceptprob))  
        return objCBDeuBestState                           # Return the best solution found.

    
    def runAlgo(self):
        '''
        This function run the main algorithm
        '''
        
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
            
        
        
        # read vdFile
        variableNames, cardinality= readVdFile(self.vdFile)


        # read data file
        self.df=convertBeneDataFile(self.dataFile, len(variableNames))
        
        # create object of EquivalenceClass
        objEC= EquivalenceClass(self.mlabPath)
        # get the opt bnt from bene
        optDag, allNodesObj= objEC.getOptDag(self.vdFile, self.dataFile, self.alpha, self.outdir, len(variableNames), cardinality)
        
        #print "optDag inside class MainAlgo and function runAlgo()"
        #print optDag
        
        HIDDEN_NAME= len(variableNames) +1
        
        # dict of dict
        cachedBDeuDict={} # keys: tuple ((parent's parent),( child's parent)) ; Values: bdeu score for the network
        edgesDict={} # keys: edge tuple (parent, child); Values: keys of Pa_C_PaPa_CPa dictionary
        hiddenNodesDict={} # keys: edge tuple (parent, child): values: hidden node objects
        algoIteratios=0
        totalPreviousBDeuScore= float("-inf")
        totalCurrentBDeuScore= float("-inf")
        currentMaxBDeu= float("-inf")
        previousMaxBDeu= float("-inf")
        
        
        with open(self.outputFile, 'w') as wf:
            nodesBDeuScore=[]
            totalUniqueObservations= self.df.shape[0]
            #print "totalUniqueObservations: %d" % totalUniqueObservations
            #print "df:"
            #print self.df
            objCBDeu= BDeuClass(self.df, allNodesObj, totalUniqueObservations, variableNames)

            # update the parent configurations for all variables
            # and the counts associated with the each parent configuration for each value of X
            for n in objCBDeu.allNodeObjects:
                tmpNode= Node()
                tmpNode= objCBDeu.allNodeObjects[n]
                objCBDeu.getUpdatedQi(tmpNode)
                objCBDeu.populateCounts(tmpNode)
                
            # find the BDeu Score for the whole structure
            for n in objCBDeu.allNodeObjects:
                tmpNode= Node()
                tmpNode=objCBDeu.allNodeObjects[n]
                tmpScore= objCBDeu.getBDeu(objCBDeu.allNodeObjects[n], self.alpha)
                nodesBDeuScore.append(tmpScore)
            
            
            print "BDeu Score for optimal dag from Bene: %f" % sum(nodesBDeuScore)
            #print initial data frame 
            self.df.to_csv('initialDF_bene'+'.csv', sep=',', index=False)
            # print the state for the random number generator
            if self.seed != None:
                stateOutFile= 'state_iter_'+str(algoIteratios)+'_initialSeed_'+ str(self.seed) +'_'+self.outputFile
                rs.storeSate(stateOutFile)
            wf.write("BDeuScore for optimal dag from Bene, %f" % sum(nodesBDeuScore))
            
            # Repeat until adding a hidden variable does not increase the score
            while True:
            
                #increment the iteration number
                algoIteratios +=1 
                #print "printing optDag"
                #print optDag
                # pDag
                cDag=objEC.generateCdag(optDag)
                # generate all dags in pDag
                dagsDict, allDagsNetworkDict= objEC.getAllDagsInPdag(cDag, cardinality)
                           
                totalUniqueObservations= self.df.shape[0] # if we introduce next hidden variable, this variable would be updated
        
                arrayListBDeuClassObjs=[]
         
                for id, dag in dagsDict.iteritems():
                    
                    nodesBDeuScore=[]
                    allNodeObjects=allDagsNetworkDict[id]
                    tdf= self.df.copy()
                    print tdf
                    # instantiate CalculateBDeuClass's object 
                    objCBDeu= BDeuClass(tdf, allNodeObjects, totalUniqueObservations, variableNames)
                    
                    #print objCBDeu.df
                    
                    # compute initial bdeu score before adding any hidden variable
                    # update the parent configurations for all variables
                    # and the counts associated with the each parent configuration for each value of X
                    for n in objCBDeu.allNodeObjects:
                        objCBDeu.getUpdatedQi(objCBDeu.allNodeObjects[n])
                        objCBDeu.populateCounts(objCBDeu.allNodeObjects[n])
                    # find the BDeu Score for the whole structure
                    for n in objCBDeu.allNodeObjects:
                        nodesBDeuScore.append(objCBDeu.getBDeu(objCBDeu.allNodeObjects[n], self.alpha))
                    
                    totalPreviousBDeuScore=sum(nodesBDeuScore)
                    initialBDeuScore= totalPreviousBDeuScore
                    
                    print "BDeu Score for a dag %d in Equivalence class before adding hidden variable: %f" % (id, totalPreviousBDeuScore)
                    print dag
                    
                    # find all the edges in the dag
                    idx=1
                    edges=[]
                    for i in dag:
                        e=[edges.append((idx,j+1)) for j in xrange(0,len(i)) if i[j]==1]
                        idx+=1
                    hiddenCount= 0 # hidden variable counter
                    
                    # for each edge. replace it with hidden variable and compute the score
                    while ((not edges)!= True):
                        # check if score for adding hidden variable for this edge is already computed in other equivalence dag
                        
                        edge= edges.pop()
                        
                        parentNode= Node()
                        childNode= Node()
                        
                        parentNode= objCBDeu.allNodeObjects[edge[0]]
                        childNode= objCBDeu.allNodeObjects[edge[1]]
                        
                        # creating key for cachedBDeuDict dictionary
                        key=tuple([tuple(parentNode.getParents()), tuple(childNode.getParents())])
                        
                        print "edge: %d ---> %d" % (edge[0], edge[1]) 
                        
                        

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
                            #print "data frame before adding hidden variable"
                            #print objCBDeu.df
                            objCBDeu.percentageHiddenCoutsSplit(h)
                            #print "data frame after adding hidden variable"
                            #print objCBDeu.df
                            
                            # write df to file called initialCountSplit.txt
                            #outName= self.outputFile+'_initialCountSplit_'+str((datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d-h%H-m%M-s%S')))
                            #newDF.to_csv(outName+'.csv', sep=',')
                            
                            # populate hidden value counts
                            hiddenBDeuScore=[]
                            for n in objCBDeu.allNodeObjects:
                                node= objCBDeu.allNodeObjects[n]
                                if node.getParentUpdateFlag() == True or node.getChildrenUpdateFlag() == True:
                                    objCBDeu.populateCounts(node)
                                    node.setLocalBDeu(objCBDeu.getBDeu(node, self.alpha))
                                hiddenBDeuScore.append(objCBDeu.getBDeu(node, self.alpha))
                            #objCBDeu.populateCounts(h)
                            
                            ##for n in objCBDeu.allNodeObjects:
                            ##    hiddenBDeuScore.append(objCBDeu.getBDeu(objCBDeu.allNodeObjects[n], self.alpha))
                    
                            initialBDeuScoreAfterAddingHidden=sum(hiddenBDeuScore)
                            
                            
                            if self.steepestAsent == True:
                                print "Steepest Asent Algorithm started ...." 
                                sIndex                  = rNumber.randint(0,objCBDeu.df.shape[0]-2) 
                                output= "Iter_"+str(algoIteratios)+"_dag_"+str(id)+"_edge_"+str(edge[0])+"_"+str(edge[1])+".sa" 
                                objCBDeu                = self.computeBDeuUsingSteepestAsent(h ,objCBDeu, initialBDeuScoreAfterAddingHidden, sIndex, self.iterations, output)
                                totalCurrentBDeuScore   = objCBDeu.dagBDeuScore
                                h                       = objCBDeu.allNodeObjects[h.getName()]
                                print "Steepest Asent Algorithm finished ...." 
                                #print "BDeu Score previousBDeu: %f; CurrentBDeu: %f" % (totalPreviousBDeuScore, totalCurrentBDeuScore)
                            elif self.simAnealFlag == True:
                                print "Simulated Annealing Algorithm started ...." 
                                sIndex                  = rNumber.randint(0,objCBDeu.df.shape[0]-1)
                                output= "Iter_"+str(algoIteratios)+"_dag_"+str(id)+"_edge_"+str(edge[0])+"_"+str(edge[1])+".sim" 
                                objCBDeu                =   self.simulatedAnealing(objCBDeu, h, initialBDeuScoreAfterAddingHidden, sIndex, self.iterations,output)
                                totalCurrentBDeuScore   = objCBDeu.dagBDeuScore
                                h                       = objCBDeu.allNodeObjects[h.getName()]
                                print "Simulated Annealing Algorithm finished ...." 
                            if initialBDeuScore < totalCurrentBDeuScore:
                                # add hidden node to the dictionary
                                hiddenNodesDict[edge]=h
                                
                                hiddenCount+=1 # count the number of hidden variable added
                                print "BDeu Score for dag %d in Equivalence class after adding hidden variable %d, PreviousBDeu: %f; CurrentBDeu: %f" % (id, h.getName(),initialBDeuScore, totalCurrentBDeuScore)   
                                print objCBDeu.df
                                diffBDeu= totalCurrentBDeuScore - initialBDeuScore
                                cachedBDeuDict[key]= diffBDeu
                                edgesDict[edge]= key
                                
                                # update the variable names after adding hidden variable
                                objCBDeu.setVariableNames(h.getName())
                                
#                                # update the edges list after adding hidden variable
#                                hChildren= h.getChildren()
#                                
#                                # update edges by adding edges of hidden variable to its children
#                                edges.append((h.getName(), hChildren[0]))
#                                edges.append((h.getName(), hChildren[1]))
                                
                                # generate new name for hidden variable
                                HIDDEN_NAME += 1
                                objCBDeu.dagBDeuScore= totalCurrentBDeuScore
                                initialBDeuScore = totalCurrentBDeuScore
                                
                                # remove edges and see if we get increase in bdeu score
                                
                                objCBDeu=self.removeEdgesFromBnt(edges, totalCurrentBDeuScore, objCBDeu)
                                
                                
                            else: # adding hidden variable didn't improve score, so go back to old state                              
                                objCBDeu.allNodeObjects= tmpAllNodesObj
                                objCBDeu.df= tmpDF.copy()
                                objCBDeu.dagBDeuScore=tmpDagBDeuScore
                                print "BDeu Score for dag %d is not changed, since no hidden varialbe is added: previousBDeu: %f; CurrentBDeu: %f"    % (id,totalPreviousBDeuScore, totalCurrentBDeuScore)        
                    # store BDeu Class object
                    arrayListBDeuClassObjs.append(objCBDeu)            
                # find the Dag' with higest bdeu score and input it to find the equivalence dags for it and repeat the whole process
                currentMaxAllNodesObjects={}
                currentMaxDF= pd.DataFrame(index=None, columns=None)
                for obj in arrayListBDeuClassObjs:
                    if currentMaxBDeu < obj.dagBDeuScore:
                        currentMaxBDeu                 = obj.dagBDeuScore
                        currentMaxAllNodesObjects      = obj.allNodeObjects
                        currentMaxDF                   = obj.df.copy()
                        currentMaxVariables            = obj.variableNames
                
                # check the looping condition
                if previousMaxBDeu < currentMaxBDeu:
                    previousMaxBDeu=currentMaxBDeu
                    # update the orginal df for next iteration
                    self.df = currentMaxDF.copy()
                    # update variable set if hidden is added
                    variableNames= currentMaxVariables
                    
                    # update optdag 
                    # update cardinality 
                    #print optdag
                    optDag, cardinality = self.printDag(currentMaxAllNodesObjects)
                    #print hidden counts and bdeu score for the dag with higest bdeu score in equivalance class
                    print "Iteration: %d , BDeu Score: %f" % (algoIteratios, currentMaxBDeu)
                    hValues= h.getKvalues().keys()
                    for i in xrange(0,len(hValues)-1):
                        count=currentMaxDF[currentMaxDF[h.getName()]==hValues[i]].Counts
                        for j in count:
                            wf.write(str(j)+',')
                        del count
                    wf.write(str(currentMaxBDeu) + "\n")
                    # print the state for the random number generator
                    stateOutFile= 'state_iter_'+str(algoIteratios)+'_initialSeed_'+ str(self.seed) +'_'+self.outputFile
                    rs.storeSate(stateOutFile)
                else: 
                    break
                
                
                
                
                
                
        
        
        
