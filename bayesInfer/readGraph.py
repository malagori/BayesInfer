__author__ = "Mehmood Alam Khan"
__email__  = "malagori@kth.se"
__version__= "0.9"
__license__ = "GPL"
import sys

from bayesInfer.node import Node

def readInitialStructure(infile):
    parents=[]
    var=''
    allNodeObjects={}
    for line in open(infile, 'r'):
        var, cardinality, pa=line.split(':')
        if len(pa) != 0:
            parents=pa.split(',')
        node=Node()
        node.setR(cardinality)
        node.setKvalues(list(range(0, cardinality, 1)))
        node.setName(var)
        node.setParents(parents)
        allNodeObjects[var]=node
    return allNodeObjects
