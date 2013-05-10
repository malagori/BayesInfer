__author__ = "Mehmood Alam Khan"
__email__  = "malagori@kth.se"
__version__= "0.9"
__license__ = "GPL"
__credits__ = ["Mehmood Alam Khan", "Pekka Parviainen"]

import sys

from bayesInfer.node import Node

def readInitialStructure(infile):
    parents=[]
    var=''
    allNodesObj={}
    for line in open(infile, 'r'):
        var, cardinality, pa=line.split(':')
        if len(pa) != 0:
            parents=pa.split(',')
        node=Node()
        node.setR(cardinality)
        
        node.setKvalues(dict.fromkeys(list(range(0, cardinality, 1))))
        node.setName(var)
        node.setParents(parents)
        allNodesObj[var]=node
    return allNodesObj
