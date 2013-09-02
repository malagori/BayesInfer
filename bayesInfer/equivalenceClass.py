__author__ = "Mehmood Alam Khan"
__email__  = "malagori@kth.se"
__version__= "0.9"
__license__ = "GPL"
__credits__ = ["Mehmood Alam Khan", "Pekka Parviainen"]

import os
import numpy as np
from mlabwrap import mlab
from bayesInfer import matlab_lib
from bayesInfer.beneWrapper import BeneWrapper

class EquivalenceClass(object):
    def __init__(self):
        self.matlabLibPath= os.path.dirname(os.path.abspath(matlab_lib.__file__))
        mlab.addpath(self.matlabLibPath)
        
    def getCpdag(self, vdFile, dataFile, outDirectory, totalVaiables):
        
        #create obj of BeneWrapper class and initialize the fields
        bwObj= BeneWrapper(vdFile, dataFile, score=1.0, outDirectory, totalVaiables) 
        
        bwObj.generateOptBnt()
        
        allNodesObj, optDag = bwObj.readBeneBnt()
        
        dag= np.array(optDag).astype('float')
        
        cDag= mlab.dag_to_cpdag(dag)
        
        return allNodesObj, cDag
        
    def getAllDagsInPdag(self, cDag):
        mlab.pdag_to_all_dags(cDag)
        
        
        
    
