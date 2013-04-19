#!/usr/bin/env python

import os
import sys
import itertools

from bayesInfer.node import Node

# total parent configurations
def getQi(node):
    qi=1
    for p in node.getParents():
        qi= qi*p.getR()
    return qi

# calculate BDeu score for one variable
def getBDeu(node):
    qi= getQi(node)
    for j in xrange(0,qi):
        
def main():
    var= Node()
    
    print "hi"
      
if __name__== "__main__":
    main()