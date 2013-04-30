#!/usr/bin/env python

import os
import sys
import itertools
import math
from __future__ import division

from bayesInfer.node import Node
from bayesInfer.readGraph import readInitialStructure

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
    
def populateCounts(node):
    print "populate the counts for this variable and its corresponding parent configurations"
    kValueDict= node.getKvalues()
    for k in kValueDict.keys():
        # populate counts for different values of X for each parent configuration 
        
# calculate BDeu score for one variable
def getBDeu(node):
    
    alpha=1
    bdeuScore=0
    a=0
    # if node.valueCountFlag is true, then do the following computation
    if node.parentUpdateFlag == True:
        # get dictionary containing parent configurations of a node with count equal to None. i.e:
        # qi is of the form {(0, 1, 1): None, (0, 2, 0): None,..,}
        getUpdatedQi(node)
        # you can populate the counts of paraent configuration here
        # read from data basically. update the node.k_values dictionary
        # populating the pConfigurations of a node. its a dictionary of parent configuration whose values is again a dictionary of var values.
        # i.e. dict(dict(values_of_variable))
        # HINT: population of data can be done in parallel in getUpdatedQi(). read HINT tag in getUpdatedQi() funciton..
        populateCounts(node)
        
    if node.valueUpdateFlag == True:
        populateCounts(node)
    # get node parent configurations
    qi= node.pConfigurations
        
    localBDeu= calculateLocalBDeu(qi, node)
                
    # set node localBDeu score there
    node.setLocalBDeu(localBDeu)
    return node.localBDeu

def calculateLocalBDeu(qi, node):
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

def main():
    nodesBDeuScore=[]
    infile='path to file containing initial structure information'
    allNodeObjects=readInitialStructure(infile)
    
    # draw initial structure when you get time
    
    # you can update the structure here.
    
    # update the parent configurations for all variables
    # and the counts associated with the each parent configuration for each value of X
    for n in allNodeObjects:
        getUpdatedQi(n)
        populateCounts(n)
    # find the BDeu Score for the whole structure
    for n in allNodeObjects:
        nodesBDeuScore.append(getBDeu(n))
        
    print "Final BDeu Score: %f" % sum(nodesBDeuScore)
      
if __name__== "__main__":
    main()