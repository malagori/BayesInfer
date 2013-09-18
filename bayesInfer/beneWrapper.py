__author__ = "Mehmood Alam Khan"
__email__  = "malagori@kth.se"
__version__= "0.9"
__license__ = "GPL"
__credits__ = ["Mehmood Alam Khan", "Pekka Parviainen"]

import os
import subprocess
from node import Node

class BeneWrapper(object):
    def __init__(self, vdFile, dataFile, score, outDirectory, totalVaiables):
        self.vdFile = vdFile                # variable descriptor file
        self.dataFile = dataFile            # data file
        self.score = score                  # default is 1.0 BDe score alpha parameter
        self.outDirectory = outDirectory    # path to output directory
        self.totalVaiables = totalVaiables  # number of variables

    def checkExe(self, exePath):
        return os.path.isfile(exePath) and os.access(exePath, os.X_OK)
    
    def which(self, program):
        fpath, fname = os.path.split(program)
        if fpath:
            if self.checkExe(program):
                return program
        else:
            for path in os.environ["PATH"].split(os.pathsep):
                path = path.strip('"')
                exe_file = os.path.join(path, program)
                if self.checkExe(exe_file):
                    return exe_file
        return None
    
    def generateOptBnt(self):
        try:
            null = open("/dev/null")
            beneStdOut= os.path.join(self.outDirectory, "bene.stdout")
            benePwd= self.which('data2net.sh')
            print "benePwd: %s" % benePwd
            print "vdFile %s" % self.vdFile
            print "dataFile %s" % self.dataFile
            print "outDirectory %s" % self.outDirectory
            print "score %s" % self.score
            
            
            if benePwd != None:
                subprocess.call([ str(benePwd), str(self.vdFile), str(self.dataFile), str(self.score), str(self.outDirectory)])
        except IOError, e:
            print ("Class: beneWrapper; Function: generateOptBnt();  Error: " + str(e))
            
    def readBeneBnt(self, cardinality):
        """
        Read bayesian network generated by bene software
        """
        # convert decimal to binary. decimalToBinary(x, n) convert x to binary and 
        # represent it in an n-bit representation
        # ['1', '2', '3']
        decimalToBinary = lambda x, n: x >= 0 and str(bin(x))[2:].zfill(n)
        
        
        allNodesObj={}
        optDag=[]
        try:
            infile= self.outDirectory+'/'+'net'
            print infile
            varName=0
            with open(infile) as f:
                for line in f:
                    node= Node()
                    parents=[]
                    
                    varParentSet=list(decimalToBinary(int(line), int(self.totalVaiables)))
                    p=[i for i in reversed(varParentSet)]
                    optDag.append(p)
                    idx=1
                    for i in reversed(varParentSet):
                        if i == '1':
                            parents.append(idx)
                        idx+=1
                            
                    varName+=1
                    node.setR(int(cardinality[varName-1]))
                    node.setKvalues(dict.fromkeys(list(xrange(0, int(cardinality[varName-1]), 1))))
                    node.setName(varName)
                    node.setParents(parents)
                    allNodesObj[varName]=node
            print optDag
            # taking transpose of list of list to get the required dag
            optDag=map(list, zip(*optDag))
            print optDag
            return optDag, allNodesObj
            
            
        except IOError, e:
            print ("Class: beneWrapper; Function: readBeneBnt();  Error: " + str(e))
        