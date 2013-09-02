__author__ = "Mehmood Alam Khan"
__email__  = "malagori@kth.se"
__version__= "0.9"
__license__ = "GPL"
__credits__ = ["Mehmood Alam Khan", "Pekka Parviainen"]

import os
import subprocess
from bayesInfer.node import Node

class BeneWrapper(object):

    def __init__(self, vdFile, dataFile, score, outDirectory, totalVaiables):
        self.vdFile = vdFile                # variable descriptor file
        self.dataFile = dataFile            # data file
        self.score = score                  # default is 1.0 BDe score alpha parameter
        self.outDirectory = outDirectory    # path to output directory
        self.totalVaiables = totalVaiables  # number of variables

        
    def generateOptBnt(self):
        try:
            null = open("/dev/null")
            beneStdOut= os.path.join(self.outDirectory, "bene.stdout")
            subprocess.check_call([ "data2net.sh", self.vdFile, self.dataFile, self.score, self.outDirectory], stdout=beneStdOut, stderr=null)
        except IOError, e:
            print ("Class: beneWrapper; Function: generateOptBnt();  Error: " + str(e))
            
    def readBeneBnt(self):
        """
        Read bayesian network generated by bene software
        """
        # convert decimal to binary. decimalToBinary(x, n) convert x to binary and 
        # represent it in an n-bit representation
        # ['1', '2', '3']
        decimalToBinary = lambda x, n: x >= 0 and str(bin(x))[2:].zfill(n)
        
        
        
        try:
            infile= self.outDirectory+'net'
            varName=0
            allNodesObj={}
            optDag=[]
            with open(infile) as f:
                for line in f:
                    node=Node()
                    parents=[]
                    varName+=1
                    cardinality=2
                    varParentSet=list(decimalToBinary(line, self.totalVaiables))
                    optDag.append(varParentSet)
                    for i in range(0, len(varParentSet)):
                        if varParentSet[i] == '1':
                            # set parent
                            parents.append(i+1)
                            
                            
                    node.setR(int(cardinality))
                    node.setKvalues(dict.fromkeys(list(range(0, int(cardinality), 1))))
                    node.setName(varName)
                    node.setParents(parents)
                    allNodesObj[varName]=node
            return allNodesObj, optDag
            
            
        except IOError, e:
            print ("Class: beneWrapper; Function: readBeneBnt();  Error: " + str(e))
        