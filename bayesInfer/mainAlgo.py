'''
Created on Sep 10, 2013

@author: malagori
'''
from __future__ import division

import time
import random as rNumber


from bayesInfer.storeRetriveSeed import RandomSeed
from bayesInfer.readDataFile import convertBeneDataFile

class MainAlgo(object):
    '''
    classdocs
    '''


    def __init__(self, vdFile, dataFile, outdir, alpha, seed, steepestAsent, seedFile  ):
        '''
        Constructor
        '''
        self.vdFile          = vdFile
        self.dataFile        = dataFile
        self.outdir          = outdir
        self.alpha           = alpha
        self.seed            = seed
        self.steepestAsent   = steepestAsent
        self.seedFile        = seedFile
        
    def runAlgo(self):
        
        # instanciate RandomSeed object
        rs=RandomSeed()
        
        if self.seed == None and self.seedFile == None:
            self.seed= time.time()
        elif self.seed != None and self.seedFile == None:
            rs.setInitialState(self.seed)
        elif self.seedFile != None and self.seed == None:
            state= rs.getSateFromFile(self.seedFile)
            rNumber.setstate(state)
            
        # read data file
        df=convertBeneDataFile(self.dataFile)
        
        # read vdFile
        
        
        
        totalUniqueObservations= df.shape[0] # if we introduce next hidden variable, this variable would be updated
        