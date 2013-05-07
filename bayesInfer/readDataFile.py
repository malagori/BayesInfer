
import pandas as pd

__author__ = "Mehmood Alam Khan"
__email__  = "malagori@kth.se"
__version__= "0.9"
__license__ = "GPL"
__credits__ = ["Mehmood Alam Khan", "Pekka Parviainen"]

def readDataFromFile(infile):
    df=pd.read_csv(infile)
   
    #       A  B  C  D Counts
    #    0  0  1  0  0 10
    #    1  1  1  1  1 5   
    #    2  1  0  0  0 3
    #
    # without hidden variable count. so for each value of hidden variable,
    # we will have separate datasets but sum of counts for different datasets should be equal to the original counts.
    # accessing specific columns: df.A
    
    return df