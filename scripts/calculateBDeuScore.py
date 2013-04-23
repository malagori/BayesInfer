#!/usr/bin/env python

import os
import sys
import itertools

from bayesInfer.node import Node
from bayesInfer.readGraph import readInitialStructure

# total parent configurations
def getQi(node):
    qi=1
    dictPaConfiguration={}
    pConfig=[]
    allParentValues=[]
    for p in node.getParents():
        #qi= qi*p.getR()
        # get values of each parent and add it to allParentValues array
        allParentValues.append(p.k_values.keys())
    for i in itertools.product(*allParentValues):
        # pConfig is of the form {(0, 1, 1), (0, 2, 0),...,}
        pConfig.append(i)
    dictPaConfiguration= dict.fromkeys(pConfig)
    # dictPaConfiguration is of the form {(0, 1, 1): None, (0, 2, 0): None,..,}
    return dictPaConfiguration

# calculate BDeu score for one variable
def getBDeu(node):
    
    # get dictionary containing parent configurations of a node with count equal to None. i.e:
    # qi is of the form {(0, 1, 1): None, (0, 2, 0): None,..,}
    qi= getQi(node)
    
    # you can populate the counts of paraent configuration here
    # read from data basically
    
    
    # traverse all values of qi
    for j in qi:
        #do something here. may be one other loop to traverse r
        # and 
        
        
        
        
def main():
    nodesBDeuScore=[]
    infile='path to file containing initial structure information'
    allNodeObjects=readInitialStructure(infile)
    for n in allNodeObjects:
        nodesBDeuScore.append(getBDeu(n))

    print "hi"
      
if __name__== "__main__":
    main()