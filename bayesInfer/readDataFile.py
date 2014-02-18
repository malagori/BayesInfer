
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

def readInitialHiddenConfig(infile):
    """
    Read data from file containing dataframe in following format:
        A  B  C  D Counts h
        0  0  1  0  0    0
        1  1  1  1  1    0
        2  1  0  0  0    1
    
    returns dataframe, number of variables
    """
    try:
        df=pd.read_pickle(infile)
        
             
    except IOError:
        print "Class: readDataFile; Function readDataFrame(); Error: could not read dataframe"

    return df

def convertBeneDataFile(infile, columns):
    """
    Read data from file contain data in following format:
        0  0  1  0  0 
        1  1  1  1  1   
        2  1  0  0  0
    
    returns dataframe, number of variables
    """
    dic={}
    rows=[]
    try:
        with open(infile) as f:
            for line in f:
                line=line.strip()
                tok=[]
                tok=line.split('\t')
                tok=map(int, tok)
                if tuple(tok) in dic.keys():
                    dic[tuple(tok)].append(1)
                else:
                    dic[tuple(tok)]= [1]
                    res= map(int, tok)
                    rows.append(res)  
    except IOError:
        print "Class: readDataFile; Function convertBeneDataFile();  Error: could not read dataframe"     
    counts=[]
    for i in rows:
        counts.append(sum(dic[tuple(i)]))
    # new columns
    newCols= [i for i in xrange(1, columns+1)]
    #create df
    df =pd.DataFrame(rows)
    df.columns= newCols
    df['Counts']= pd.Series(counts, df.index)
    
      
    return df
    
    
def readVdFile(infile):
    
    tokens=[]
    cardinality=[]
    variableNames=[]
    try:
        with open(infile, 'r') as f:
            for line in f:
                tokens=line.split("\t")
                variableNames.append(tokens[0])
                cardinality.append(len(tokens)-1)
    except IOError:
        print "file: readVdFile; Function readVdFile();  Error: could not read vd file"   
        
    return variableNames, cardinality
       
                
    
    
    
    
    