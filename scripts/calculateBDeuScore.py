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
from __future__ import division

from bayesInfer.node import Node
from bayesInfer.readGraph import readInitialStructure
from bayesInfer.readDataFile import readDataFromFile

# total parent configurations
def getUpdatedQi(node):
    dictPaConfiguration={}
    pConfig=[]
    allParentValues=[]
    for p in node.getParents():
        # get values of each parent and add it to allParentValues array
        allParentValues.append(p.k_values.keys())
    for i in itertools.product(*allParentValues):
        # pConfig is of the form {(0, 1, 1), (0, 2, 0),...,}
        pConfig.append(i)
        
    #dictPaConfiguration= dict.fromkeys(pConfig)
    # dictPaConfiguration is of the form {(0, 1, 1): None, (0, 2, 0): None,..,}
    #node.setpConfiguration(dictPaConfiguration)
    
    node.setpConfiguration(pConfig)
    node.valueUpdateFlag == True
    
def populateCounts(node):
    print "populate the counts for this variable and its corresponding parent configurations"
    kValueDict= node.getKvalues()
    # take care here if the variable have any parents or not.
    #
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
    # A B D
    # 0 1 1 
    # 1 0 1
    # 1 1 0
    # remove last parent, which is hidden
    # get the counts from data 
    nodeParents=node.getParents()
    
    # All records with var value = k
    localDf=df[df[node.name]==k]
    # for each parent value
    idx=0
    for pa in nodeParents:
        localDf=localDf[localDf[pa]==j[idx]]
        idx+=1
    # return the row count satisfiying the conditions
    return len(localDf.index)
    
    
# calculate BDeu score for one variable
def getBDeu(node, alpha):
    
    bdeuScore=0.0
    a=0
    # if node.valueCountFlag is true, then do the following computation
    #if node.parentUpdateFlag == True:
        # set the array containing parent configurations of a node. i.e:
        # qi is of the form [(0, 1, 1), (0, 2, 0),..,]
    #    getUpdatedQi(node)
    # you can populate the counts of node with corresponding paraent configuration if node.valueUpdateFlag == True
    #if node.valueUpdateFlag == True:
    #    populateCounts(node)
    # get node parent configurations
    qi= node.pConfigurations
        
    localBDeu= calculateLocalBDeu(qi, node, alpha)
                
    # set node localBDeu score there
    node.setLocalBDeu(localBDeu)
    return node.localBDeu

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
    localBDeu=0.00
    for j in qi:
        Nijk=0
        c=0
        b=0
        # iterate over different values of variable X and retrive each dictionary containing parentConfig:K_value_count
        for k, v in node.getKvalues(node).iteritems():
            Nijk+=v[j]
            c += (math.lgamma(zri+ v[j]) - math.lgamma(zri))
            
        b= math.lgamma(z+ Nijk)       
        localBDeu += a - b + c
        
    return localBDeu

def main(dataFile, structureFile):
    
    global df
    df=readDataFromFile(dataFile)
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
    h.setKvalues(dict.fromkeys(list(range(0, cardinality, 1)))) 
    h.children.append(child1) # add children to hidden variable
    h.children.append(child2) 
    h.childrenUpdateFlag= True
    allNodeObjects[child1].parentUpdateFlag= True # get the children nodes and update the parentUpdateFlag
    allNodeObjects[child2].parentUpdateFlag= True
   
    # compute the BDeu score again
    for n in allNodeObjects:
        node=allNodeObjects[n]
        if node.parentUpdateFlag== True: # if true its a child of hidden variable. so, calculate BDeu again
            # compute new parent configuration set
            getUpdatedQi(node)
            orginal_k_values_counts=node.k_values # its a copy of k_values dict
            
            # change the counts of this node according to some criteria i.e.
            # for new parent configuration, assign the counts such that sum of counts of new parent
            # configurations is equal to counts of old parent configuration which we split
            # to get the new parent configuration. 
    
    
       
if __name__== "__main__":
    main()