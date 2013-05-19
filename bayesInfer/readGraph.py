__author__ = "Mehmood Alam Khan"
__email__  = "malagori@kth.se"
__version__= "0.9"
__license__ = "GPL"
__credits__ = ["Mehmood Alam Khan", "Pekka Parviainen"]

import sys
import itertools
from random import randrange
from bayesInfer.node import Node

def readInitialStructure(infile):
    parents=[]
    var=''
    allNodesObj={}
    for line in open(infile, 'r'):
        var, cardinality, pa=line.split(':')
        if len(pa) > 0:
            parents=pa.split(',')
        node=Node()
        node.setR(cardinality)
        
        node.setKvalues(dict.fromkeys(list(range(0, cardinality, 1))))
        node.setName(var)
        node.setParents(parents)
        allNodesObj[var]=node
    return allNodesObj

def generateData():
    """ generate dataset """
    wf= open('data1.txt', 'w')
    alphabits=['A','B','C','D', 'Counts']
    cardinality= [2,2,2,2]
    
    for i in xrange(0,4):
        wf.write(str(alphabits[i]+',')) # print except last column name
    wf.write(alphabits[4]) # print the last column name
    wf.write('\n')
    values=[[0, 1], [0, 1], [0, 1], [0, 1]]
    for i in itertools.product(*values):
        for j in xrange(0,len(i)):
            wf.write(str(i[j])+',')
        #wf.write(str(i[-1]))
        wf.write(str(randrange(100)))
        wf.write('\n')
    wf.close() 
        
    