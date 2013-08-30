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
    allNodesObj={}
    try:
        with open(infile) as f:
            for line in f:
                var=''
                pa=''
                parents=[]
                var, cardinality, pa=line.split(':')
                pa=pa.strip()  # removing \n characters at the end
                if len(pa) > 0:
                    parents=pa.split(',')
                node=Node()
                # cardinality is read as a string, we must convert string to integer 
                node.setR(int(cardinality))
                node.setKvalues(dict.fromkeys(list(range(0, int(cardinality), 1))))
                node.setName(var)
                node.setParents(parents)
                allNodesObj[var]=node
    except IOError:
        print "Class readGraph; Error: can\'t find structure file"
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
        
    