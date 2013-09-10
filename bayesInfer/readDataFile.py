
import pandas as pd
from collections import Counter
import numpy as np

__author__ = "Mehmood Alam Khan"
__email__  = "malagori@kth.se"
__version__= "0.9"
__license__ = "GPL"
__credits__ = ["Mehmood Alam Khan", "Pekka Parviainen"]

def readDataFrame(infile):
    """
    Read data from file containing dataframe in following format:
        A  B  C  D Counts
        0  0  1  0  0 
        1  1  1  1  1   
        2  1  0  0  0
    
    returns dataframe, number of variables
    """
    #       A  B  C  D Counts
    #    0  0  1  0  0 10
    #    1  1  1  1  1 5   
    #    2  1  0  0  0 3
    #
    # without hidden variable count. so for each value of hidden variable,
    # we will have separate datasets but sum of counts for different datasets should be equal to the original counts.
    # accessing specific columns: df.A
    try:
        df=pd.read_csv(infile)
        
             
    except IOError:
        print "Class: readDataFile; Function readDataFrame(); Error: could not read dataframe"

    return df

def convertBeneDataFile(infile):
    """
    Read data from file contain data in following format:
        0  0  1  0  0 
        1  1  1  1  1   
        2  1  0  0  0
    
    returns dataframe, number of variables
    """
    
    try:
        benDF=pd.read_csv(infile, sep='\t')
        newCols= [i for i in range(1, benDF.shape[1]+1)]
        benDF.columns= newCols
        # converting rows into tuples
        trows= [tuple(row) for row in benDF.values]
        rowCount=[(item, count) for item, count in Counter(trows).iteritems() ]
        dfDictionary= dict(rowCount)
        newCols.append('Counts')
        
        # create new data frame
        
        varColumns=np.array([key for [key,val] in dfDictionary.iteritems()])
        countCoumns=np.array([val for [key,val] in dfDictionary.iteritems()])
        
        df =pd.DataFrame(varColumns)
        df['Counts']= pd.Series(countCoumns, df.index)
        df.columns= newCols
        
    except IOError:
        print "Class: readDataFile; Function convertBeneDataFile();  Error: could not read dataframe"
         
    return df
    