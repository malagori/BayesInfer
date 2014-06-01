from __future__ import division
__author__ = "Mehmood Alam Khan"
__email__  = "malagori@kth.se"
__version__= "0.9"
__license__ = "GPL"
__credits__ = ["Mehmood Alam Khan", "Pekka Parviainen"]


import itertools
import math
import pandas as pd
import random as rNumber
import numpy as np
from pandas import Series
import copy

from node import Node

class BDeuClass(object):
    
    def __init__(self, df, dfOriginal, allNodeObjects, totalUniqueObservations, variableNames):
        self.df                         = df
        self.dfOriginal                 = dfOriginal
        self.allNodeObjects             = copy.deepcopy(allNodeObjects)
        self.totalUniqueObservations    = totalUniqueObservations
        self.dagBDeuScore               = float("-inf")
        self.variableNames              = variableNames
        
    def setOriginalDF(self,dfOriginal):
        self.dfOriginal= copy.deepcopy(dfOriginal.copy())
    def setTotalUniqueObservations(self,totalUniqueObservations):
        self.totalUniqueObservations= totalUniqueObservations
    def setAllNodeObjects(self, allNodeObjects):
        self.allNodeObjects= copy.deepcopy(allNodeObjects)
    def setDagBDeuScore(self, BDeuScore):
        self.dagBDeuScore= BDeuScore 
    def setDF(self,df):
        self.df= df
        
    def setVariableNames(self, name):
        if name not in self.variableNames:
            self.variableNames.append(name)
      
    def getUpdatedQi(self,node):
        # total parent configurations
        """ if node is a parent then it has no parent configuration, return the function other wise compute parent configurations"""
        #print "Function: getUpdatedQi; Node Name: %s" % node.getName()
        if len(node.getParents())==0:
            return
        else:
            #dictPaConfiguration={}
            pConfig=[]
            allParentValues=[]
            for p in node.getParents():
                # get values of each parent and add it to allParentValues array
                allParentValues.append(self.allNodeObjects[p].getKvalues().keys())
            for i in itertools.product(*allParentValues):
                # pConfig is of the form {(0, 1, 1), (0, 2, 0),...,}
                pConfig.append(i)
                
            #dictPaConfiguration= dict.fromkeys(pConfig)
            # dictPaConfiguration is of the form {(0, 1, 1): None, (0, 2, 0): None,..,}
            #node.setpConfiguration(dictPaConfiguration)
            
            node.setpConfiguration(pConfig)
            self.allNodeObjects[node.getName()]=node
            
    def binaryPurterbation(self, h, firstRIndex):
        '''
        This function will swap the counts of an observation to zero or K-count for the given value of hidden variable
        '''
        if firstRIndex < self.totalUniqueObservations:
            secondRIndex  = firstRIndex + self.totalUniqueObservations
        else:
            secondRIndex  = firstRIndex - self.totalUniqueObservations 
    
        if self.df.Counts[firstRIndex] == 0:
            self.df.Counts[firstRIndex] = self.df.Counts[secondRIndex]
            self.df.Counts[secondRIndex] = 0
        else:
            self.df.Counts[secondRIndex]= self.df.Counts[firstRIndex]
            self.df.Counts[firstRIndex]= 0

    def twoRowsCountPerturbation(self, firstRowIndex, secondRowIndex,decrementValue, incrementFlag):
    
        #print "df before peruturbation"
        #print df
        # decrement the record
        if incrementFlag == False:
            
            if (self.df.Counts[firstRowIndex]- decrementValue) >= 0:
                if firstRowIndex < self.totalUniqueObservations:
                    incrementedDfIndex  = firstRowIndex + self.totalUniqueObservations
                else:
                    incrementedDfIndex  = firstRowIndex - self.totalUniqueObservations 
                self.df.Counts[incrementedDfIndex]   += decrementValue
                self.df.Counts[firstRowIndex]        -= decrementValue 
                    
                # increment the record by decrementValue
                if secondRowIndex < self.totalUniqueObservations:
                    decrementedDfIndex  = secondRowIndex  + self.totalUniqueObservations
                    dfCopyIndex         = secondRowIndex
                else:
                    decrementedDfIndex  = secondRowIndex - self.totalUniqueObservations
                    dfCopyIndex         = secondRowIndex - self.totalUniqueObservations
                #print "firstRowIndex: %d, secondRowIndex: %d" % (firstRowIndex, secondRowIndex)
                #print "(df.Counts[secondRowIndex] + decrementValue): %d <= dfCopy.Counts[dfCopyIndex]: %d and (df.Counts[decrementedDfIndex] - decrementValue) >= 0 %d" % ((df.Counts[secondRowIndex] + decrementValue), dfCopy.Counts[dfCopyIndex], (df.Counts[decrementedDfIndex] - decrementValue))
                if (self.df.Counts[secondRowIndex] + decrementValue) <= self.dfOriginal.Counts[dfCopyIndex] and (self.df.Counts[decrementedDfIndex] - decrementValue) >= 0:
                    #print "rindex: %d, incrementedDfIndex: %d, totalUniqueObservations: %d" % (rIndex, incrementedDfIndex, totalUniqueObservations)
                    self.df.Counts[decrementedDfIndex]   -= decrementValue
                    self.df.Counts[secondRowIndex]       += decrementValue
                else:
                    self.df.Counts[incrementedDfIndex]   -= decrementValue
                    self.df.Counts[firstRowIndex]        += decrementValue 
                                    
        else:
            # increment the record by decrementValue
            if firstRowIndex < self.totalUniqueObservations:
                decrementedDfIndex  = firstRowIndex  + self.totalUniqueObservations
                dfCopyIndex         = firstRowIndex
            else:
                decrementedDfIndex  = firstRowIndex - self.totalUniqueObservations
                dfCopyIndex         = firstRowIndex - self.totalUniqueObservations
            #print "firstRowIndex: %d, secondRowIndex: %d" % (firstRowIndex, secondRowIndex)
            #print "(df.Counts[firstRowIndex] + decrementValue): %d <= dfCopy.Counts[dfCopyIndex]: %d and (df.Counts[decrementedDfIndex] - decrementValue) >= 0 %d" % ((self.df.Counts[firstRowIndex] + decrementValue), self.dfCopy.Counts[dfCopyIndex], (self.df.Counts[decrementedDfIndex] - decrementValue))
#            print "(df.Counts[firstRowIndex] + decrementValue): %d " % ((self.df.Counts[firstRowIndex] + decrementValue))
#            print "<= dfCopy.Counts[dfCopyIndex]: %d " % (self.dfOriginal.Counts[dfCopyIndex])
#            print "and (df.Counts[decrementedDfIndex] - decrementValue) >= 0 %d" % (  (self.df.Counts[decrementedDfIndex] - decrementValue))
            if (self.df.Counts[firstRowIndex] + decrementValue) <= self.dfOriginal.Counts[dfCopyIndex] and (self.df.Counts[decrementedDfIndex] - decrementValue) >= 0:
                self.df.Counts[decrementedDfIndex]   -= decrementValue
                self.df.Counts[firstRowIndex]        += decrementValue
                if (self.df.Counts[secondRowIndex]- decrementValue) >= 0:
                    if secondRowIndex < self.totalUniqueObservations:
                        incrementedDfIndex  = secondRowIndex + self.totalUniqueObservations
                    else:
                        incrementedDfIndex  = secondRowIndex - self.totalUniqueObservations 
                    self.df.Counts[incrementedDfIndex]   += decrementValue
                    self.df.Counts[secondRowIndex]       -= decrementValue
                else:
                    self.df.Counts[decrementedDfIndex]   += decrementValue
                    self.df.Counts[firstRowIndex]        -= decrementValue
        #print "df after peruturbation"
        #print df

        
    def populateCounts(self, node):
        """populate the counts for this variable for different parent configurations"""
        kValueDict= node.getKvalues()
        
        # take care here if the variable have any parents or not.
        if len(node.getParents()) == 0:
            paValueDict={}
            for i in kValueDict.keys():
                paValueDict[i]= self.getDataCount(i, [], node)
            #print "paValueDict:"
            #print paValueDict
            node.setParentValueCounts(paValueDict)
            self.allNodeObjects[node.getName()]=copy.deepcopy(node)
            #print node.getParentValueCount()
        else:
            for k in kValueDict.keys():
                pConfigDict={}
                # populate counts for different values of X for each parent configuration 
                for j in node.getPaConfigurations():
                    # j is a tuple containing parent's values. we have to change it to list 
                    pConfigDict[j]=self.getDataCount(k,list(j), node)
                kValueDict[k]=pConfigDict
                
            node.setKvalues(kValueDict)
            self.allNodeObjects[node.getName()]=copy.deepcopy(node)
        
    def getDataCount(self,k, j, node):
        # remember: j here is a list not tuple
         
        # subset data according to parent nodes
        # A B D Counts
        # 0 1 1 5
        # 1 0 1 3
        # 1 1 0 2
        # 
        # check if parent configuration is zero
        if len(j)==0:
            # compute the counts for parent variable for kth value
            localDframe=self.df[self.df[node.getName()]==k].copy()
        else:
            
            # All records with var value = k
            localDframe=self.df[self.df[node.getName()]==k].copy()
            # for each parent value
            idx=0
            for pa in node.getParents():
                localDframe=localDframe[localDframe[pa]==j[idx]].copy()
                idx+=1
        # return the row count satisfiying the conditions
        return sum(localDframe.Counts)
          
    def getBDeu(self, node, alpha):
        # calculate BDeu score for one variable
        localBDeu=0.0
        # if node.valueCountFlag is true, then do the following computation
        #if node.parentUpdateFlag == True:
            # set the array containing parent configurations of a node. i.e:
            # qi is of the form [(0, 1, 1), (0, 2, 0),..,]
        #    getUpdatedQi(node)
        # you can populate the counts of node with corresponding paraent configuration if node.valueUpdateFlag == True
        #if node.valueUpdateFlag == True:
        #    populateCounts(node)
        
        # consider the case when a node has no parents
        # if len(node.pConfigurations) == 0:
        # compute localBDeu for this node. calculateLocalBDeu function will be different.
        #print "Node name: %s" % node.getName()
        #print "parent configuration %d" % len(node.pConfigurations)
        if len(node.pConfigurations) == 0 :
            #print "calculate local BDeu for parent node"
            localBDeu=self.calculateBDeuParentNode(node, alpha)
        else:
            # get node parent configurations
            qi= node.getPaConfigurations()
                
            localBDeu= self.calculateLocalBDeu(qi, node, alpha)
        # set node localBDeu score there
        node.setLocalBDeu(localBDeu)
        self.allNodeObjects[node.getName()]=node
        
        return node.localBDeu
    

    def calculateBDeuParentNode(self, node, alpha):
        """ calculate bdeu score for parent node"""
        a= math.lgamma(alpha)
        N= sum(self.df['Counts']) # total number of observations
        b= math.lgamma(alpha + N)
        c= 0.0
        z= alpha/node.getR()
        
        for k, v in node.getParentValueCount().iteritems():
            c+= (math.lgamma(z + v) - math.lgamma(z) )
        
        localBDeu = a -b + c
        
        return localBDeu
   
    def calculateLocalBDeu(self, qi, node, alpha):
        # traverse all values of qi
        # compute the following
        # z= alpha/qi
        # a = log_gamma(z)
        # b = log_gamma( z + ri_sum_k1(N_ijk) )
        # c = ri_sum_k1{ log_gamma( z.(1/ri) + N_ijk) - log_gamma(z.(1/ri) }
        # sum_qi ( a - b + c)  
        z=  alpha/len(qi)
        ri= 1/node.getR()
        zri=z*ri
        a= math.lgamma(z)
        localBDeu=0.0
        for j in qi:
            Nijk=0
            c=0.0
            b=0.0
            # iterate over different values of variable X and retrieve each dictionary containing parentConfig:K_value_count
            for k, v in node.getKvalues().iteritems():
                Nijk+=v[j]
                c += (math.lgamma(zri+ v[j]) - math.lgamma(zri))
                
            b= math.lgamma(z+ Nijk)       
            localBDeu += a - b + c
            
        return localBDeu
    
    def randomCountPerturbation(self, h):
        #print "perturb the count here"
        
        hiddenName=h.getName()
        # pick a random index
        while (True):
            rIndex=rNumber.randint(0,self.df.shape[0]-1)   # -1 because indexing starts from 0
            if self.df.Counts[rIndex] > 0:
                #print "before decrimenting index %d %d: " % (rIndex, df.Counts[rIndex])
                self.df.Counts[rIndex] -= 1 # decrement by 1 
                #print "after decrimenting index %d %d: " % (rIndex, df.Counts[rIndex])
                break
        # select the hidden value at rIndex
        valueH=self.df[hiddenName][rIndex]
        # choose other values of Hidden variable other then
        hValuesWithoutValueH=[i for i in h.getKvalues().keys() if i != valueH] # hidden value with out randomly selected valueH
        
        # choose one of the hValuesWithoutValueH values to be incremented by 1
        incrementedHvalue=rNumber.choice(hValuesWithoutValueH)
        #print "hidden value %d " % valueH
        #print "incrementedHvalue %d " % incrementedHvalue
        #print "decremented index %d" % rIndex
        #print "totalUniqueObservations %d" % totalUniqueObservations
        if incrementedHvalue <  valueH:
            #print "rIndex %d" % rIndex
            # compute the index of record in df to be incremented
            incrementedDfIndex=rIndex-((valueH-incrementedHvalue)*self.totalUniqueObservations)
            #print "incrementedDfIndex %d " % incrementedDfIndex
            #print df.Counts[incrementedDfIndex]
            self.df.Counts[incrementedDfIndex] += 1
            #print df.Counts[incrementedDfIndex]
        else:
            incrementedDfIndex = rIndex + ((incrementedHvalue- valueH)*self.totalUniqueObservations)
            self.df.Counts[incrementedDfIndex] += 1
        #print "incremented Index %d " % incrementedDfIndex
        
    def percentageHiddenCoutsSplit(self, h):
        # old dataframe was:
        # A B C Counts 
        # 0 1 1 10          
        # 0 0 1 4           
        # 0 1 0 4      
        #
        # and 
        #  the new counts for will be randomly split
        # A B C Counts H
        # 0 1 1 4      0     
        # 0 0 1 2      0     
        # 0 1 0 2      0
        # 0 1 1 3      1    
        # 0 0 1 1      1     
        # 0 1 0 1      1 
        # 0 1 1 3      2    
        # 0 0 1 1      2     
        # 0 1 0 1      2     
        #
        # for each value of hidden variable, we create a column vector storing counts.
    
        hiddenName=h.name
        hiddenColumn=Series(np.zeros(self.df.shape[0]), index=self.df.index)
        self.df[hiddenName]=hiddenColumn
        df_temp= self.df.copy()
        
        
        #df_temp.Counts=np.zeros(df_temp.shape[0]) # rows with hidden value zero is add here
        for i in h.getKvalues().keys():
            if i != 0: # rows with hidden value not equal to zero are add here
                col=[i]*df_temp.shape[0]  # fastest way to create list ;)
                copyCountList= [math.floor(i*rNumber.random()) for i in df_temp.Counts]
                df_temp[hiddenName]=col
                df_temp.Counts=copyCountList
                self.df.Counts[0:self.totalUniqueObservations]= self.df.Counts[0:self.totalUniqueObservations]-copyCountList
                self.df= self.df.append(df_temp, ignore_index=True)
        #return self.df  # delete the temporary data frame to save memory
    
    def binaryHiddenCountSplit(self, h):
        '''
        This funciton will split the counts of dataframe in a binary fashion when we introduce hidden variable
        # old dataframe was:
        # A B C Counts 
        # 0 1 1 10          
        # 0 0 1 4           
        # 0 1 0 4      
        #
        # and 
        #  the new counts for will be randomly split
        # A B C Counts H
        # 0 1 1 10     0     
        # 0 0 1 0      0     
        # 0 1 0 4      0
        # 0 1 1 0      1    
        # 0 0 1 4      1     
        # 0 1 0 0      1 
         
        '''
        hiddenName=h.name
        hiddenColumn=Series(np.zeros(self.df.shape[0]), index=self.df.index)
        self.df[hiddenName]=hiddenColumn
        df_temp= self.df.copy()
        col=[1]*df_temp.shape[0]  # fastest way to create list ;)
        df_temp[hiddenName]=col
        self.df.Counts[0:self.totalUniqueObservations]= self.df.Counts[0:self.totalUniqueObservations]-df_temp.Counts
        self.df= self.df.append(df_temp, ignore_index=True)
        
        for i in xrange(2*self.df.shape[0], 0, -1):
            idx= rNumber.randint(0,self.df.shape[0]-1)
            if idx >=  self.totalUniqueObservations:
                if self.df.Counts[idx] == 0:
                    self.df.Counts[idx] = self.df.Counts[idx-self.totalUniqueObservations]
                    self.df.Counts[idx-self.totalUniqueObservations] = 0
                else:
                    self.df.Counts[idx - self.totalUniqueObservations] = self.df.Counts[idx]
                    self.df.Counts[idx] = 0
            else:
                if self.df.Counts[idx] == 0:
                    self.df.Counts[idx] = self.df.Counts[idx + self.totalUniqueObservations]
                    self.df.Counts[idx + self.totalUniqueObservations] = 0
                else:
                    self.df.Counts[idx + self.totalUniqueObservations] = self.df.Counts[idx]
                    self.df.Counts[idx] = 0
    
    
    def addHiddenNodeToDf(self, h):
        # old dataframe was:
        # A B C Counts 
        # 0 1 1 10          
        # 0 0 1 4           
        # 0 1 0 4      
        #
        # and 
        # new dataframe would like this
        # A B C Counts H
        # 0 1 1 4      0     
        # 0 0 1 2      0     
        # 0 1 0 2      0
        # 0 1 1 3      1    
        # 0 0 1 1      1     
        # 0 1 0 1      1 
        # 0 1 1 3      2    
        # 0 0 1 1      2     
        # 0 1 0 1      2     
        #
        # for each value of hidden variable, we create a column vector storing counts.
    
        hiddenName=h.name
        hiddenColumn=Series(np.zeros(self.df.shape[0]), index=self.df.index)
        self.df[hiddenName]=hiddenColumn
        df_temp= self.df.copy()
        copyCountList= [math.floor(i/h.getR()) for i in df_temp.Counts]
        
        #df_temp.Counts=np.zeros(df_temp.shape[0]) # rows with hidden value zero is add here
        for i in h.getKvalues().keys():
            if i != 0: # rows with hidden value not equal to zero are add here
                col=[i]*df_temp.shape[0]  # fastest way to create list ;)
                df_temp[hiddenName]=col
                df_temp.Counts=copyCountList
                self.df.Counts[0:self.totalUniqueObservations]= self.df.Counts[0:self.totalUniqueObservations]-copyCountList
                self.df= self.df.append(df_temp, ignore_index=True)
        del df_temp # delete the temporary data frame to save memory
        return self.df  
    
    def addHiddenNode(self, name, cardinality, child1, child2):
        #name = raw_input("Enter Hidden Variable Name: ")
        #card = raw_input("Enter Hidden Cardinality: ")
        #cardinality= int(card)
        #child1= raw_input("Enter first Child Variable Name: ")
        #child2= raw_input("Enter first Child Variable Name: ")
        # change the structure by introducing hidden variable
        h= Node()
        h.setName(name)
        h.setR(cardinality)
        h.setKvalues(dict.fromkeys(list(xrange(0, cardinality, 1)))) 
        h.addChild(child1) # add children to hidden variable
        h.addChild(child2)
        h.setChildrenUpdateFlag(True)
        
        ######
        
        childNode1= self.allNodeObjects[child1]
        childNode2= self.allNodeObjects[child2]
        
        childNode1.setParentUpdateFlag( True) # get the children nodes and update the parentUpdateFlag
        childNode2.setParentUpdateFlag( True)
        
        childParents1= childNode1.getParents()
        childParents2= childNode2.getParents()
        
        # remove the current edge  
        #      hidden
        #    /    \
        #  a/      \
        #   \      /
        #    \--->/
        #     <---
        
        if child2 in childParents1:
            childParents1.remove(child2)
            # add hidden parent to the parentset
            childParents1.append(name)
            childParents2.append(name)
        elif child1 in childParents2:
            childParents2.remove(child1)
            # add hidden parent to the parentset
            childParents1.append(name)
            childParents2.append(name)
        else:
            childParents1.append(name)
            childParents2.append(name)
            
        childNode1.setParents(childParents1)
        childNode2.setParents(childParents2)
        
        self.allNodeObjects[child1]= childNode1
        self.allNodeObjects[child2]= childNode2
        self.allNodeObjects[h.getName()]= h

        ######     
#        self.allNodeObjects[h.getName()]= h 
#        self.allNodeObjects[child1].setParentUpdateFlag( True) # get the children nodes and update the parentUpdateFlag
#        self.allNodeObjects[child2].setParentUpdateFlag( True)
#        self.allNodeObjects[child1].addParent(name)
#        self.allNodeObjects[child2].addParent(name)
#        
        
        # compute new parent configuration set for both the children
        self.getUpdatedQi(self.allNodeObjects[child1]) 
        self.getUpdatedQi(self.allNodeObjects[child2]) 

        return h
    
    def countPerturbation(self, h, rIndex,decrementValue, incrementFlag):
        # decrement the record
        if incrementFlag == False:
            if (self.df.Counts[rIndex]- decrementValue) > 0:
                if rIndex <  self.totalUniqueObservations:
                    incrementedDfIndex  = rIndex + self.totalUniqueObservations
                    dfCopyIndex         = rIndex
                else:
                    incrementedDfIndex  = rIndex - self.totalUniqueObservations 
                    dfCopyIndex         = rIndex - self.totalUniqueObservations
                    
                #print "rindex: %d, incrementedDfIndex: %d, totalUniqueObservations: %d" % (rIndex, incrementedDfIndex, self.totalUniqueObservations)
                self.df.Counts[incrementedDfIndex] += decrementValue
                self.df.Counts[rIndex] -= decrementValue
        else:
            # increment the record by decrementValue
            if rIndex < self.totalUniqueObservations:
                decrementedDfIndex  = rIndex  + self.totalUniqueObservations
                dfCopyIndex         = rIndex
            else:
                decrementedDfIndex  = rIndex - self.totalUniqueObservations
                dfCopyIndex         = rIndex - self.totalUniqueObservations
            #print "rIndex: %d, decrementedDfIndex: %d, totalUniqueObservations: %d" % (rIndex, decrementedDfIndex, self.totalUniqueObservations)
            if (self.df.Counts[rIndex] + decrementValue) <= self.dfOriginal.Counts[dfCopyIndex] and (self.df.Counts[decrementedDfIndex] - decrementValue) >= 0:
                self.df.Counts[decrementedDfIndex]   -= decrementValue
                self.df.Counts[rIndex]               += decrementValue
                
#    def countPerturbation(self, h, rIndex, decrementedValue, incrementFlag):
#        #print "perturb the count here"
#        hiddenName=h.getName()
#        # pick a random index
#        if incrementFlag == False:
#            if self.df.Counts[rIndex] > 0:
#                #print "before decrementing index %d %d: " % (rIndex, df.Counts[rIndex])
#                self.df.Counts[rIndex] -= decrementedValue # decrement by 1 
#                #print "after decrementing index %d %d: " % (rIndex, df.Counts[rIndex])
#        
#            # select the hidden value at rIndex
#            valueH=self.df[hiddenName][rIndex]
#            # choose other values of Hidden variable other then
#            hValuesWithoutValueH=[i for i in h.getKvalues().keys() if i != valueH] # hidden value with out the privously selected valueH
#            
#            # choose one of the hValuesWithoutValueH values to be incremented by 1
#            incrementedHvalue=rNumber.choice(hValuesWithoutValueH)
#            #print "hidden value %d " % valueH
#            #print "incrementedHvalue %d " % incrementedHvalue
#            #print "decremented index %d" % rIndex
#            #print "totalUniqueObservations %d" % totalUniqueObservations
#            if incrementedHvalue <  valueH:
#                #print "rIndex %d" % rIndex
#                # compute the index of record in df to be incremented
#                incrementedDfIndex=rIndex-((valueH-incrementedHvalue)*self.totalUniqueObservations)
#                #print "incrementedDfIndex %d " % incrementedDfIndex
#                #print df.Counts[incrementedDfIndex]
#                self.df.Counts[incrementedDfIndex] += decrementedValue
#                #print df.Counts[incrementedDfIndex]
#            else:
#                incrementedDfIndex = rIndex + ((incrementedHvalue- valueH)*self.totalUniqueObservations)
#                self.df.Counts[incrementedDfIndex] += decrementedValue
#            #print "incremented Index %d " % incrementedDfIndex
#        else:    
#            # select the hidden value at rIndex
#            valueH=self.df[hiddenName][rIndex]
#            # choose other values of Hidden variable other then
#            hValuesWithoutValueH=[i for i in h.getKvalues().keys() if i != valueH] # hidden value with out the privously selected valueH
#            
#            # choose one of the hValuesWithoutValueH values to be incremented by 1
#            decrementedHvalue=rNumber.choice(hValuesWithoutValueH)
#            
#            #  loop until valid case occur other wise exit with doing nothing
#    
#            for i in xrange(0,len(hValuesWithoutValueH)):
#                if decrementedHvalue <  valueH:
#                    decrementedDfIndex=rIndex-((valueH-decrementedHvalue)*self.totalUniqueObservations)
#                    if self.df.Counts[decrementedDfIndex] == 0:
#                        noDecrement=True
#                        continue
#                    else:
#                        self.df.Counts[decrementedDfIndex] -= decrementedValue
#                        noDecrement=False
#                        break
#                else:
#                    decrementedDfIndex = rIndex + ((decrementedHvalue- valueH)*self.totalUniqueObservations)
#                    if self.df.Counts[decrementedDfIndex] == 0:
#                        noDecrement=True
#                        continue
#                    else:
#                        self.df.Counts[decrementedDfIndex] -= decrementedValue
#                        noDecrement=False
#                        break
#         
#            
#            if noDecrement == False:
#                self.df.Counts[rIndex] += decrementedValue # decrement by 1 
#            #print "after incrementing index %d %d: " % (rIndex, df.Counts[rIndex])

        