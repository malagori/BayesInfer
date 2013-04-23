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
        node.setKvalues(dict.fromkeys(list(range(0, cardinality, 1))))
        node.setName(var)
        node.setParents(parents)
        allNodeObjects[var]=node
    return allNodeObjects
