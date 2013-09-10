#!/usr/bin/env python

__author__ = "Mehmood Alam Khan"
__email__  = "malagori@kth.se"
__version__= "0.9"
__license__ = "GPL"
__credits__ = ["Mehmood Alam Khan", "Pekka Parviainen"]

import os
import sys
import argparse

from bayesInfer.mainAlgo import MainAlgo


def runMainAlgo(vdFile, dataFile, outdir, alpha, seed, steepestAsent, seedFile ):
    
    objMainAlgo = MainAlgo(vdFile, dataFile, outdir, alpha, seed, steepestAsent, seedFile)
    objMainAlgo.runAlgo()
    
def main(argv):
    
     
    # Take care of input
    parser = argparse.ArgumentParser(description="Parse input arguments and print output.")
    parser.add_argument('-vd', metavar='vdfile' ,type=str, help='Specify path to the file containing variable information ')
    parser.add_argument('-d', metavar='dataFile',type=str, help='Specify path to the data file ')
    parser.add_argument('-od', metavar='outDir',type=str, help='Specify path to output directory ')
    parser.add_argument('-a', metavar='alpha',type=float , help='Specify alpha parameter ', default=1.0)
    parser.add_argument('-Sa', dest='SteepestAsent',action="store_true", help='Steepest Asent is used if set to True ')
    parser.add_argument('-s', metavar='initialSeed',type=int , help='Specify initial seed. if both initialSeed and loadseed option are not provided then system time will be taken as the default seed  ', default=None)
    parser.add_argument('-l', metavar='loadSeed',type=int , help='Specify path to a file containing previous state', default=None)
    parser.add_argument('-o', metavar='outfile', type=str, help='Specify the file to output the results. ', default= 'counts_bdeu_results.txt')
    args = parser.parse_args()
    
    vdFile          = args.vd
    outdir          = args.od
    dataFile        = args.d
    alpha           = args.a
    seed            = args.s
    steepestAsent   = args.SteepestAsent
    seedFile        = args.l
    
    
    print "vdfile: %s"          % vdFile
    print "outputdirectory %s"  % outdir
    print "dataFile %s"         % dataFile
    print "alpha %f"            % alpha 
    print "seed %d"             % seed
    
    runMainAlgo(vdFile, dataFile, outdir, alpha, seed, steepestAsent, seedFile)
    
    
if __name__== "__main__":
    main(sys.argv[1:])