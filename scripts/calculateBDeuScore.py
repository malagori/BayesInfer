#!/usr/bin/env python

import os
import sys
import itertools

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
    dictPaConfiguration= dict.fromkeys(pConfig)
    # dictPaConfiguration is of the form {(0, 1, 1): None, (0, 2, 0): None,..,}
    node.setpConfiguration(dictPaConfiguration)

def populateCounts(node):
    print "populate the counts for this variable and its corresponding parent configurations"

# calculate BDeu score for one variable
def getBDeu(node):
    
    bdeuScore=0
    if node.parentUpdateFlag == True:
        # get dictionary containing parent configurations of a node with count equal to None. i.e:
        # qi is of the form {(0, 1, 1): None, (0, 2, 0): None,..,}
        getUpdatedQi(node)
        # you can populate the counts of paraent configuration here
        # read from data basically. update the node.k_values dictionary
        # populating the pConfigurations of a node. its a dictionary of parent configuration whose values is again a dictionary of var values.
        # i.e. dict(dict(values_of_variable))
        populateCounts(node)
            
        # get node parent configurations
        qi= node.pConfigurations
        # traverse all values of qi
        for j in qi:
            #do something here. may be one other loop to traverse r
            # 
            
            
            # set node localBDeu score there
            node.setLocalBDeu(bdeuScore)
    return node.localBDeu
        
        
def main():
    nodesBDeuScore=[]
    infile='path to file containing initial structure information'
    allNodeObjects=readInitialStructure(infile)
    
    # draw initial structure when you get time
    
    # you can update the structure here.
    
    # find the BDeu Score for the whole structure
    for n in allNodeObjects:
        nodesBDeuScore.append(getBDeu(n))

    print "hi"
      
if __name__== "__main__":
    main()