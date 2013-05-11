#!/usr/bin/env python

__author__ = "Mehmood Alam Khan"
__email__  = "malagori@kth.se"
__version__= "0.9"
__license__ = "GPL"
__credits__ = ["Mehmood Alam Khan", "Pekka Parviainen"]

import os
import sys
import itertools
import math
import random as rNumber
from __future__ import division
import numpy as np
from pandas import Series

from bayesInfer.node import Node
from bayesInfer.readGraph import readInitialStructure
from bayesInfer.readDataFile import readDataFromFile

# total parent configurations
def getUpdatedQi(node):
    """ if node is a parent then it has no parent configuration, return the function other wise compute parent configurations"""
    if len(node.getParents())==0:
        return
    else:
        dictPaConfiguration={}
        pConfig=[]
        allParentValues=[]
        for p in node.getParents():
            # get values of each parent and add it to allParentValues array
            allParentValues.append(p.getKvalues().keys())
        for i in itertools.product(*allParentValues):
            # pConfig is of the form {(0, 1, 1), (0, 2, 0),...,}
            pConfig.append(i)
            
        #dictPaConfiguration= dict.fromkeys(pConfig)
        # dictPaConfiguration is of the form {(0, 1, 1): None, (0, 2, 0): None,..,}
        #node.setpConfiguration(dictPaConfiguration)
        
        node.setpConfiguration(pConfig)
        node.valueUpdateFlag == True 
    
def populateCounts(node):
    """populate the counts for this variable for different parent configurations"""
    kValueDict= node.getKvalues()
    
    # take care here if the variable have any parents or not.
    if len(node.getParents()) == 0:
        paValueDict={}
        for i in kValueDict.keys():
            paValueDict[i]= getDataCount(i, [], node)
            
        node.setParentValueCounts(paValueDict)
    else:
        for k in kValueDict.keys():
            pConfigDict={}
            # populate counts for different values of X for each parent configuration 
            for j in node.getPaConfigurations():
                # j is a tuple containing parent's values. we have to change it to list 
                pConfigDict[j]=getDataCount(k,list(j), node)
            kValueDict[k]=pConfigDict
            
        node.setKvalues(kValueDict)
    
def getDataCount(k, j, node):
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
        localDframe=df[df[node.name]==k]
    else:
        # All records with var value = k
        localDframe=df[df[node.name]==k]
        # for each parent value
        idx=0
        for pa in node.getParents():
            localDframe=localDframe[localDframe[pa]==j[idx]]
            idx+=1
    # return the row count satisfiying the conditions
    return sum(localDframe.Counts)
    
    
# calculate BDeu score for one variable
def getBDeu(node, alpha):
    
    localBDeu=0.0
    a=0
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
    if len(node.pConfigurations) == 0:
        print "calculate local BDeu for parent node"
        localBDeu=calculateBDeuParentNode(node, alpha)
    else:
        # get node parent configurations
        qi= node.pConfigurations
            
        localBDeu= calculateLocalBDeu(qi, node, alpha)
    # set node localBDeu score there
    node.setLocalBDeu(localBDeu)
    
    return node.localBDeu

def calculateBDeuParentNode(node, alpha):
    """ calculate bdeu score for parent node"""
    a= math.lgamma(alpha)
    N= sum(df['Counts']) # total number of observations
    b= math.lgamma(alpha + N)
    c= 0.0
    z= alpha/node.getR()
    
    for k, v in node.getParentValueCount().iteritems():
        c+= (math.lgamma(z + v) - math.lgamma(z) )
    
    localBDeu = a -b + c
    
    return localBDeu
    
def calculateLocalBDeu(qi, node, alpha):
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
        for k, v in node.getKvalues(node).iteritems():
            Nijk+=v[j]
            c += (math.lgamma(zri+ v[j]) - math.lgamma(zri))
            
        b= math.lgamma(z+ Nijk)       
        localBDeu += a - b + c
        
    return localBDeu

#def genNewJandSplitEqually(j,counts, hiddenParent):
#    print "generate new J and split the counts equally"
#    oldList= [ i for i in j] # converting tuple into list
#    newList=oldList
#    parentValues= hiddenParent.getKvalues().keys()
#    for v in parentValues:
#        newList.append(v)
#        newJ[]
#        newList=oldList
#    
#def splitCounts(node):
#    Divisor= 2
#    print "split counts according to some criteria"
#    newKValueDict= {}
#    newJ=[]
#    newJFlagDone=False
#    #oldJ= newJ[0:-1] # -1 because we are considering that only one hidden variable is added in parentset
#    for key, countDict in node.getKvalues().iteritems():
#        for j in countDict.iteritems():
#            newJ=genNewJandSplitEqually(j, countDict[j], allNodeObjects[node.parents[-1]]) # send the last parent which is hidden parent
#            newKValueDict[newJ]= countDict[j]
       
def countPerturbation(h):
    print "perturb the count here"
    
    hiddenName=h.getName()
    # pick a random index
    while (True):
        rIndex=rNumber.randint(0,df.shape[0]-1)   # -1 because indexing starts from 0
        if df.Counts[rIndex] > 0:
            df.Counts[rIndex] -= 1 # decrement by 1 
    # select the hidden value at rIndex
    valueH=df[hiddenName][rIndex]
    # choose other values of Hidden variable other then
    hValuesWithoutValueH=[i for i in h.getKvalues().keys() if i != valueH] # hidden value with out randomly selected valueH
    
    # choose one of the hValuesWithoutValueH values to be incremented by 1
    incrementedHvalue=rNumber.choice(hValuesWithoutValueH)
    if incrementedHvalue <  valueH:
        # compute the index of record in df to be incremented
        incrementedDfIndex=rIndex-((valueH-incrementedHvalue)*totalInitialObservations)
        df.Counts[incrementedDfIndex] += 1
    else:
        incrementedDfIndex = rIndex + ((incrementedHvalue- valueH)*totalInitialObservations)
        df.Counts[incrementedDfIndex] += 1

    
            
def addHiddenNodeToDf(h):
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
    df[hiddenName]=Series(np.zeros(df.shape[0]), index=df.index)
    df_temp= df.copy()
    copyCountList= [math.floor(i/h.getR()) for i in df_temp.Counts]
    
    #df_temp.Counts=np.zeros(df_temp.shape[0]) # rows with hidden value zero is add here
    for i in h.getKvalues().keys():
        if i != 0: # rows with hidden value not equal to zero are add here
            col=[i]*df_temp.shape[0]  # fastest way to create list ;)
            df_temp[hiddenName]=col
            df_temp.Counts=copyCountList
            df.Counts[0:totalInitialObservations]= df.Counts[0:totalInitialObservations]-copyCountList
            df= df.append(df_temp, ignore_index=True)
    
    del df_temp  # delete the temprory data frame to save memory
    
            
def main(dataFile, structureFile):
    
    
    global df, allNodeObjects, totalInitialObservations
    df=readDataFromFile(dataFile)
    totalInitialObservations= sum(df['Counts']) # if we introduce next hidden variable, this variable would be updated
    alpha=1
    nodesBDeuScore=[]
    infile='path to file containing initial structure information'
    allNodeObjects=readInitialStructure(structureFile)
    
    # draw initial structure when you get time .. future work
    # you can update the structure here....       future work
    
    
    # update the parent configurations for all variables
    # and the counts associated with the each parent configuration for each value of X
    for n in allNodeObjects:
        getUpdatedQi(allNodeObjects[n])
        populateCounts(allNodeObjects[n])
    # find the BDeu Score for the whole structure
    for n in allNodeObjects:
        nodesBDeuScore.append(getBDeu(allNodeObjects[n], alpha))
        
    print "BDeu Score for Initial Structure: %f" % sum(nodesBDeuScore)
    
    # change the structure by introducing hidden variable
    h= Node()
    cardinality=2 # user can input this information here
    child1= 'B'
    child2= 'C'
    h.name='h1'
    h.setR(cardinality)
    h.setKvalues(dict.fromkeys(list(range(0, cardinality, 1)))) 
    h.children.append(child1) # add children to hidden variable
    h.children.append(child2) 
    h.childrenUpdateFlag= True
    allNodeObjects[child1].parentUpdateFlag= True # get the children nodes and update the parentUpdateFlag
    allNodeObjects[child2].parentUpdateFlag= True
    # compute new parent configuration set for both the children
    getUpdatedQi(allNodeObjects[child1]) 
    getUpdatedQi(allNodeObjects[child2]) 
    allNodeObjects[h.name]= h  # adding h to the structure
    
    # add hidden variable to the dataframe and  split almost counts equally:
    addHiddenNodeToDf(h)
    
    maxIter= 1000
    for iterations in xrange(0, maxIter): 
        nodesBDeuScore=[]
        countBDeuList=[]
        # perturb the counts here
        countPerturbation(h)
        # compute the BDeu score again
        for n in allNodeObjects:
            node=allNodeObjects[n]
            if node.parentUpdateFlag== True: # if true its a child of hidden variable. so, calculate BDeu again
                             
                # change the counts of this node according to some criteria i.e.
                # for new parent configuration, assign the counts such that sum of counts of new parent
                # configurations is equal to counts of old parent configuration which we split
                # to get the new parent configuration. 
                populateCounts(node)
                node.setLocalBDeu(getBDeu(node, alpha))
            nodesBDeuScore.append(node.getLocalBDeu())
        countBDeuList= h.getParentValueCount().values()
        countBDeuList.append(sum(nodesBDeuScore))
        print countBDeuList  # countBDeu: [c1, c2, c3, ..., cn, bdeu_score]
  
    
    
       
if __name__== "__main__":
    main()