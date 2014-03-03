#!/usr/bin/env python
from __future__ import division
__author__ = "Mehmood Alam Khan"
__email__  = "malagori@kth.se"
__version__= "0.9"
__license__ = "GPL"
__credits__ = ["Mehmood Alam Khan", "Pekka Parviainen"]

import os
import sys
import itertools
import math
import datetime
import time
from math import exp
import argparse
import random as rNumber
import numpy as np
from pandas import Series
import pandas as pd

from bayesInfer.node import Node
from bayesInfer.readGraph import readInitialStructure
from bayesInfer.readDataFile import readDataFrame, readInitialHiddenConfig
from bayesInfer.storeRetriveSeed import RandomSeed



# total parent configurations
def getUpdatedQi(node):
    """ if node is a parent then it has no parent configuration, return the function other wise compute parent configurations"""
    if len(node.getParents())==0:
        return
    else:
        #dictPaConfiguration={}
        pConfig=[]
        allParentValues=[]
        for p in node.getParents():
            # get values of each parent and add it to allParentValues array
            allParentValues.append(allNodeObjects[p].getKvalues().keys())
        for i in itertools.product(*allParentValues):
            # pConfig is of the form {(0, 1, 1), (0, 2, 0),...,}
            pConfig.append(i)
            
        #dictPaConfiguration= dict.fromkeys(pConfig)
        # dictPaConfiguration is of the form {(0, 1, 1): None, (0, 2, 0): None,..,}
        #node.setpConfiguration(dictPaConfiguration)
        
        node.setpConfiguration(pConfig)
        allNodeObjects[node.getName()]=node
        
    
def populateCounts(node):
    """populate the counts for this variable for different parent configurations"""
    kValueDict= node.getKvalues()
    
    # take care here if the variable have any parents or not.
    if len(node.getParents()) == 0:
        paValueDict={}
        for i in kValueDict.keys():
            paValueDict[i]= getDataCount(i, [], node)
        #print "paValueDict:"
        #print paValueDict
        node.setParentValueCounts(paValueDict)
        allNodeObjects[node.getName()]=node
        #print node.getParentValueCount()
    else:
        for k in kValueDict.keys():
            pConfigDict={}
            # populate counts for different values of X for each parent configuration 
            for j in node.getPaConfigurations():
                # j is a tuple containing parent's values. we have to change it to list 
                pConfigDict[j]=getDataCount(k,list(j), node)
            kValueDict[k]=pConfigDict
            
        node.setKvalues(kValueDict)
        allNodeObjects[node.getName()]=node
    
def getDataCount(k, j, node):
    # remember: j here is a list not tuple
     
    # subset data according to parent nodes
    # A B D Counts
    # 0 1 1 5
    # 1 0 1 3
    # 1 1 0 2
    # 
    # check if parent configuration is zero
    if len(j)==0:
        # compute the counts for parent variable for kth value
        localDframe=df[df[node.getName()]==k]
    else:
        # All records with var value = k
        localDframe=df[df[node.getName()]==k]
        # for each parent value
        idx=0
        for pa in node.getParents():
            localDframe=localDframe[localDframe[pa]==j[idx]]
            idx+=1
    # return the row count satisfiying the conditions
    return sum(localDframe.Counts)
    
    
# calculate BDeu score for one variable
def getBDeu(node, alpha):
    
    localBDeu=0.0
    # if node.valueCountFlag is true, then do the following computation
    #if node.parentUpdateFlag == True:
        # set the array containing parent configurations of a node. i.e:
        # qi is of the form [(0, 1, 1), (0, 2, 0),..,]
    #    getUpdatedQi(node)
    # you can populate the counts of node with corresponding paraent configuration if node.valueUpdateFlag == True
    #if node.valueUpdateFlag == True:
    #    populateCounts(node)
    
    # consider the case when a node has no parents
    # if len(node.pConfigurations) == 0:
    # compute localBDeu for this node. calculateLocalBDeu function will be different.
    #print "Node name: %s" % node.getName()
    #print "parent configuration %d" % len(node.pConfigurations)
    if len(node.pConfigurations) == 0 :
        #print "calculate local BDeu for parent node"
        localBDeu=calculateBDeuParentNode(node, alpha)
    else:
        # get node parent configurations
        qi= node.getPaConfigurations()
            
        localBDeu= calculateLocalBDeu(qi, node, alpha)
    # set node localBDeu score there
    node.setLocalBDeu(localBDeu)
    allNodeObjects[node.getName()]=node
    
    return node.localBDeu

def calculateBDeuParentNode(node, alpha):
    """ calculate bdeu score for parent node"""
    a= math.lgamma(alpha)
    N= sum(df['Counts']) # total number of observations
    b= math.lgamma(alpha + N)
    c= 0.0
    z= alpha/node.getR()
    
    for k, v in node.getParentValueCount().iteritems():
        c+= (math.lgamma(z + v) - math.lgamma(z) )
    
    localBDeu = a -b + c
    
    return localBDeu
    
def calculateLocalBDeu(qi, node, alpha):
    # traverse all values of qi
    # compute the following
    # z= alpha/qi
    # a = log_gamma(z)
    # b = log_gamma( z + ri_sum_k1(N_ijk) )
    # c = ri_sum_k1{ log_gamma( z.(1/ri) + N_ijk) - log_gamma(z.(1/ri) }
    # sum_qi ( a - b + c)  
    z=  alpha/len(qi)
    ri= 1/node.getR()
    zri=z*ri
    a= math.lgamma(z)
    localBDeu=0.0
    for j in qi:
        Nijk=0
        c=0.0
        b=0.0
        # iterate over different values of variable X and retrieve each dictionary containing parentConfig:K_value_count
        for k, v in node.getKvalues().iteritems():
            Nijk+=v[j]
            c += (math.lgamma(zri+ v[j]) - math.lgamma(zri))
            
        b= math.lgamma(z+ Nijk)       
        localBDeu += a - b + c
        
    return localBDeu

#def genNewJandSplitEqually(j,counts, hiddenParent):
#    print "generate new J and split the counts equally"
#    oldList= [ i for i in j] # converting tuple into list
#    newList=oldList
#    parentValues= hiddenParent.getKvalues().keys()
#    for v in parentValues:
#        newList.append(v)
#        newJ[]
#        newList=oldList
#    
#def splitCounts(node):
#    Divisor= 2
#    print "split counts according to some criteria"
#    newKValueDict= {}
#    newJ=[]
#    newJFlagDone=False
#    #oldJ= newJ[0:-1] # -1 because we are considering that only one hidden variable is added in parentset
#    for key, countDict in node.getKvalues().iteritems():
#        for j in countDict.iteritems():
#            newJ=genNewJandSplitEqually(j, countDict[j], allNodeObjects[node.parents[-1]]) # send the last parent which is hidden parent
#            newKValueDict[newJ]= countDict[j]
       
def randomCountPerturbation(h):
    #print "perturb the count here"
    
    hiddenName=h.getName()
    # pick a random index
    while (True):
        rIndex=rNumber.randint(0,df.shape[0]-1)   # -1 because indexing starts from 0
        if df.Counts[rIndex] > 0:
            #print "before decrimenting index %d %d: " % (rIndex, df.Counts[rIndex])
            df.Counts[rIndex] -= 1 # decrement by 1 
            #print "after decrimenting index %d %d: " % (rIndex, df.Counts[rIndex])
            break
    # select the hidden value at rIndex
    valueH=df[hiddenName][rIndex]
    # choose other values of Hidden variable other then
    hValuesWithoutValueH=[i for i in h.getKvalues().keys() if i != valueH] # hidden value with out randomly selected valueH
    
    # choose one of the hValuesWithoutValueH values to be incremented by 1
    incrementedHvalue=rNumber.choice(hValuesWithoutValueH)
    #print "hidden value %d " % valueH
    #print "incrementedHvalue %d " % incrementedHvalue
    #print "decremented index %d" % rIndex
    #print "totalUniqueObservations %d" % totalUniqueObservations
    if incrementedHvalue <  valueH:
        #print "rIndex %d" % rIndex
        # compute the index of record in df to be incremented
        incrementedDfIndex=rIndex-((valueH-incrementedHvalue)*totalUniqueObservations)
        #print "incrementedDfIndex %d " % incrementedDfIndex
        #print df.Counts[incrementedDfIndex]
        df.Counts[incrementedDfIndex] += 1
        #print df.Counts[incrementedDfIndex]
    else:
        incrementedDfIndex = rIndex + ((incrementedHvalue- valueH)*totalUniqueObservations)
        df.Counts[incrementedDfIndex] += 1
    #print "incremented Index %d " % incrementedDfIndex
    
def percentageHiddenCoutsSplit(h,df):
    # old dataframe was:
    # A B C Counts 
    # 0 1 1 10          
    # 0 0 1 4           
    # 0 1 0 4      
    #
    # and 
    #  the new counts for will be randomly split
    # A B C Counts H
    # 0 1 1 4      0     
    # 0 0 1 2      0     
    # 0 1 0 2      0
    # 0 1 1 3      1    
    # 0 0 1 1      1     
    # 0 1 0 1      1 
    # 0 1 1 3      2    
    # 0 0 1 1      2     
    # 0 1 0 1      2     
    #
    # for each value of hidden variable, we create a column vector storing counts.

    hiddenName=h.name
    hiddenColumn=Series(np.zeros(df.shape[0]), index=df.index)
    df[hiddenName]=hiddenColumn
    df_temp= df.copy()
    
    
    #df_temp.Counts=np.zeros(df_temp.shape[0]) # rows with hidden value zero is add here
    for i in h.getKvalues().keys():
        if i != 0: # rows with hidden value not equal to zero are add here
            col=[i]*df_temp.shape[0]  # fastest way to create list ;)
            copyCountList= [math.floor(i*rNumber.random()) for i in df_temp.Counts]
            df_temp[hiddenName]=col
            df_temp.Counts=copyCountList
            df.Counts[0:totalUniqueObservations]= df.Counts[0:totalUniqueObservations]-copyCountList
            df= df.append(df_temp, ignore_index=True)
    return df  # delete the temporary data frame to save memory

def addHiddenNodeToDf(h,df):
    # old dataframe was:
    # A B C Counts 
    # 0 1 1 10          
    # 0 0 1 4           
    # 0 1 0 4      
    #
    # and 
    # new dataframe would like this
    # A B C Counts H
    # 0 1 1 4      0     
    # 0 0 1 2      0     
    # 0 1 0 2      0
    # 0 1 1 3      1    
    # 0 0 1 1      1     
    # 0 1 0 1      1 
    # 0 1 1 3      2    
    # 0 0 1 1      2     
    # 0 1 0 1      2     
    #
    # for each value of hidden variable, we create a column vector storing counts.

    hiddenName=h.name
    hiddenColumn=Series(np.zeros(df.shape[0]), index=df.index)
    df[hiddenName]=hiddenColumn
    df_temp= df.copy()
    copyCountList= [math.floor(i/h.getR()) for i in df_temp.Counts]
    
    #df_temp.Counts=np.zeros(df_temp.shape[0]) # rows with hidden value zero is add here
    for i in h.getKvalues().keys():
        if i != 0: # rows with hidden value not equal to zero are add here
            col=[i]*df_temp.shape[0]  # fastest way to create list ;)
            df_temp[hiddenName]=col
            df_temp.Counts=copyCountList
            df.Counts[0:totalUniqueObservations]= df.Counts[0:totalUniqueObservations]-copyCountList
            df= df.append(df_temp, ignore_index=True)
    del df_temp # delete the temporary data frame to save memory
    return df  

def addHiddenNode(name, cardinality, child1, child2):
    #name = raw_input("Enter Hidden Variable Name: ")
    #card = raw_input("Enter Hidden Cardinality: ")
    #cardinality= int(card)
    #child1= raw_input("Enter first Child Variable Name: ")
    #child2= raw_input("Enter first Child Variable Name: ")
    # change the structure by introducing hidden variable
    h= Node()
    h.setName(name)
    h.setR(cardinality)
    h.setKvalues(dict.fromkeys(list(range(0, cardinality, 1)))) 
    h.addChild(child1) # add children to hidden variable
    h.addChild(child2)
    h.setChildrenUpdateFlag(True)
    allNodeObjects[child1].setParentUpdateFlag( True) # get the children nodes and update the parentUpdateFlag
    allNodeObjects[child2].setParentUpdateFlag( True)
    # compute new parent configuration set for both the children
    print "allNodeObjects[child1] config before:"
    print allNodeObjects[child1].getPaConfigurations()
    getUpdatedQi(allNodeObjects[child1]) 
    print "allNodeObjects[child1] config after:"
    print allNodeObjects[child1].getPaConfigurations()
    print "allNodeObjects[child2] config before:"
    print allNodeObjects[child2].getPaConfigurations()
    getUpdatedQi(allNodeObjects[child2]) 
    print "allNodeObjects[child2] config after:"
    print allNodeObjects[child2].getPaConfigurations() 
    allNodeObjects[h.getName()]= h  # adding h to the structure
    return h
 
def countPerturbation( h, rIndex,decrementValue, incrementFlag):
    #print "df before peruturbation"
    #print df
    # decrement the record
    if incrementFlag == False:
        if (df.Counts[rIndex]- decrementValue) >= 0:
            if rIndex < totalUniqueObservations:
                incrementedDfIndex  = rIndex + totalUniqueObservations
                #dfCopyIndex         = rIndex
            else:
                incrementedDfIndex  = rIndex - totalUniqueObservations 
                #dfCopyIndex         = rIndex - totalUniqueObservations
                
            #print "rindex: %d, incrementedDfIndex: %d, totalUniqueObservations: %d" % (rIndex, incrementedDfIndex, totalUniqueObservations)
            df.Counts[incrementedDfIndex] += decrementValue
            df.Counts[rIndex] -= decrementValue
    else:
        # increment the record by decrementValue
        if rIndex < totalUniqueObservations:
            decrementedDfIndex  = rIndex  + totalUniqueObservations
            dfCopyIndex         = rIndex
        else:
            decrementedDfIndex  = rIndex - totalUniqueObservations
            dfCopyIndex         = rIndex - totalUniqueObservations
        #print "rIndex: %d, decrementedDfIndex: %d, totalUniqueObservations: %d" % (rIndex, decrementedDfIndex, totalUniqueObservations)
        if (df.Counts[rIndex] + decrementValue) <= dfCopy.Counts[dfCopyIndex] and (df.Counts[decrementedDfIndex] - decrementValue) >= 0:
            df.Counts[decrementedDfIndex]   -= decrementValue
            df.Counts[rIndex]               += decrementValue
    #print "df after perturbation"
    #print df

def twoRowsCountPerturbation( h, firstRowIndex, secondRowIndex,decrementValue, incrementFlag):

    #print "df before peruturbation"
    #print df
    # decrement the record
    if incrementFlag == False:
        
        if (df.Counts[firstRowIndex]- decrementValue) >= 0:
            if firstRowIndex < totalUniqueObservations:
                incrementedDfIndex  = firstRowIndex + totalUniqueObservations
            else:
                incrementedDfIndex  = firstRowIndex - totalUniqueObservations 
            df.Counts[incrementedDfIndex]   += decrementValue
            df.Counts[firstRowIndex]        -= decrementValue 
                
            # increment the record by decrementValue
            if secondRowIndex < totalUniqueObservations:
                decrementedDfIndex  = secondRowIndex  + totalUniqueObservations
                dfCopyIndex         = secondRowIndex
            else:
                decrementedDfIndex  = secondRowIndex - totalUniqueObservations
                dfCopyIndex         = secondRowIndex - totalUniqueObservations
            #print "firstRowIndex: %d, secondRowIndex: %d" % (firstRowIndex, secondRowIndex)
            #print "(df.Counts[secondRowIndex] + decrementValue): %d <= dfCopy.Counts[dfCopyIndex]: %d and (df.Counts[decrementedDfIndex] - decrementValue) >= 0 %d" % ((df.Counts[secondRowIndex] + decrementValue), dfCopy.Counts[dfCopyIndex], (df.Counts[decrementedDfIndex] - decrementValue))
            if (df.Counts[secondRowIndex] + decrementValue) <= dfCopy.Counts[dfCopyIndex] and (df.Counts[decrementedDfIndex] - decrementValue) >= 0:
                #print "rindex: %d, incrementedDfIndex: %d, totalUniqueObservations: %d" % (rIndex, incrementedDfIndex, totalUniqueObservations)
                df.Counts[decrementedDfIndex]   -= decrementValue
                df.Counts[secondRowIndex]       += decrementValue
            else:
                df.Counts[incrementedDfIndex]   -= decrementValue
                df.Counts[firstRowIndex]        += decrementValue 
                                
    else:
        # increment the record by decrementValue
        if firstRowIndex < totalUniqueObservations:
            decrementedDfIndex  = firstRowIndex  + totalUniqueObservations
            dfCopyIndex         = firstRowIndex
        else:
            decrementedDfIndex  = firstRowIndex - totalUniqueObservations
            dfCopyIndex         = firstRowIndex - totalUniqueObservations
        #print "firstRowIndex: %d, secondRowIndex: %d" % (firstRowIndex, secondRowIndex)
        #print "(df.Counts[firstRowIndex] + decrementValue): %d <= dfCopy.Counts[dfCopyIndex]: %d and (df.Counts[decrementedDfIndex] - decrementValue) >= 0 %d" % ((df.Counts[firstRowIndex] + decrementValue), dfCopy.Counts[dfCopyIndex], (df.Counts[decrementedDfIndex] - decrementValue))
        if (df.Counts[firstRowIndex] + decrementValue) <= dfCopy.Counts[dfCopyIndex] and (df.Counts[decrementedDfIndex] - decrementValue) >= 0:
            df.Counts[decrementedDfIndex]   -= decrementValue
            df.Counts[firstRowIndex]        += decrementValue
            if (df.Counts[secondRowIndex]- decrementValue) >= 0:
                if secondRowIndex < totalUniqueObservations:
                    incrementedDfIndex  = secondRowIndex + totalUniqueObservations
                else:
                    incrementedDfIndex  = secondRowIndex - totalUniqueObservations 
                df.Counts[incrementedDfIndex]   += decrementValue
                df.Counts[secondRowIndex]       -= decrementValue
            else:
                df.Counts[decrementedDfIndex]   += decrementValue
                df.Counts[firstRowIndex]        -= decrementValue
    #print "df after peruturbation"
    #print df
    
    
    
def countPerturbationOld( h, rIndex,decrementValue, incrementFlag):
    #print "perturb the count here"
    hiddenName=h.getName()
    # pick a random index
    if incrementFlag == False:
        noIncrement= True
        if (df.Counts[rIndex]- decrementValue) > 0: # check if record after decrementation is above zero.  
            # select the hidden value at rIndex
            valueH=df[hiddenName][rIndex]
            # choose other values of Hidden variable other then
            hValuesWithoutValueH=[i for i in h.getKvalues().keys() if i != valueH] # hidden value with out the privously selected valueH
            
            # choose one of the hValuesWithoutValueH values to be incremented by 1
            incrementedHvalue=rNumber.choice(hValuesWithoutValueH)

            for i in xrange(0,len(hValuesWithoutValueH)):
                if incrementedHvalue <  valueH:
                    # compute the index of record in df to be incremented
                    incrementedDfIndex=rIndex-((valueH-incrementedHvalue)*totalUniqueObservations)
                    # compute the index of record in dfCopy to check for constraint
                    if incrementedDfIndex <= totalUniqueObservations:
                        dfCopyIndex= incrementedDfIndex
                    else:
                        dfCopyIndex=incrementedDfIndex - totalUniqueObservations
                        
                    if (df.Counts[incrementedDfIndex] + decrementValue) > dfCopy.Counts[dfCopyIndex]:
                        noIncrement=True
                        continue
                    else:
                        df.Counts[incrementedDfIndex] += decrementValue
                        noIncrement=False
                        break
                else:
                    
                    incrementedDfIndex = rIndex + ((incrementedHvalue- valueH)*totalUniqueObservations)
                    
                    if incrementedDfIndex <= totalUniqueObservations:
                        dfCopyIndex= incrementedDfIndex
                    else:
                        dfCopyIndex=incrementedDfIndex - totalUniqueObservations
                        
                    if (df.Counts[incrementedDfIndex] + decrementValue) > dfCopy.Counts[dfCopyIndex]:
                        noIncrement=True
                        continue
                    else:
                        df.Counts[incrementedDfIndex] += decrementValue
                        noIncrement=False
                        break
                    
            if noIncrement == False:
                df.Counts[rIndex] -= decrementValue # decrement by 1
    else:    
        # compute the index of record in dfCopy to check for constraint
        if rIndex <= totalUniqueObservations:
            dfCopyIndex= rIndex
        else:
            dfCopyIndex=rIndex - totalUniqueObservations
            
        if (df.Counts[rIndex] + decrementValue) <= dfCopy.Counts[dfCopyIndex]:
            noDecrement= True
            # select the hidden value at rIndex
            valueH=df[hiddenName][rIndex]
            # choose other values of Hidden variable other then
            hValuesWithoutValueH=[i for i in h.getKvalues().keys() if i != valueH] # hidden value with out the privously selected valueH
            
            # choose one of the hValuesWithoutValueH values to be incremented by 1
            decrementedHvalue=rNumber.choice(hValuesWithoutValueH)
            
            #  loop until valid case occur other wise exit with doing nothing
    
            for i in xrange(0,len(hValuesWithoutValueH)):
                if decrementedHvalue <  valueH:
                    decrementedDfIndex=rIndex-((valueH-decrementedHvalue)*totalUniqueObservations)
                    if df.Counts[decrementedDfIndex] == 0:
                        noDecrement=True
                        continue
                    else:
                        df.Counts[decrementedDfIndex] -= decrementValue
                        noDecrement=False
                        break
                else:
                    decrementedDfIndex = rIndex + ((decrementedHvalue- valueH)*totalUniqueObservations)
                    if df.Counts[decrementedDfIndex] == 0:
                        noDecrement=True
                        continue
                    else:
                        df.Counts[decrementedDfIndex] -= decrementValue
                        noDecrement=False
                        break
         
            
            if noDecrement == False:
                df.Counts[rIndex] += decrementValue # decrement by 1 
        #print "after incrementing index %d %d: " % (rIndex, df.Counts[rIndex])

def temperature( k, kmax):
    """
        This function implements the temperature function of simulated anealing algorithm
    """
    temp= (1 - k/ float(kmax))
    return temp

def probAcceptance( e, enew, T):
    """
        this function compute the probability of acceptance 
    """
    prob=0.0
    if (e < enew):
        prob=1.0
        #print "accept with prob = 1"
    else:
        prob= exp(-( -enew + e )/ float(T))
        #print "e : %f  enew: %f" % (e, enew)
    return prob        
            
            

def main(argv):
    # global variables 
    global df, dfCopy, allNodeObjects, totalUniqueObservations
    nodesBDeuScore=[]
     
    # Take care of input
    parser = argparse.ArgumentParser(description="Parse input arguments and print output.")
    parser.add_argument('-b', metavar='initialStructureFile' ,type=str, help='Specify path to the file containing initial structure. ')
    parser.add_argument('-d', metavar='dataFile',type=str, help='Specify path to the data file ')
    parser.add_argument('-dh', metavar='iniHiddenConfigFile',type=str, help='Specify path to the file containing initial hidden counts configurations.')
    parser.add_argument('-x', metavar='hiddenName',type=str, help='Specify Name for hidden variable')
    parser.add_argument('-c', metavar='cardinality',type=int , help='Specify cardinality of hidden variable ', default=2)
    parser.add_argument('-c1', metavar='child1',type=str, help='Specify Name for first child variable')
    parser.add_argument('-c2', metavar='child2',type=str, help='Specify Name for second child variable')
    parser.add_argument('-dc', metavar='decrementValue',type=int , help='Specify the decrement value ', default=1)
    parser.add_argument('-a', metavar='alpha',type=float , help='Specify path to the data file ', default=1.0)
    parser.add_argument('-Sa', dest='SteepestAsent',action="store_true", help='Steepest Ascent is used if set to True ')
    parser.add_argument('-hx', dest='excludeHidBDeu',action="store_true", help='BDeu score for hidden variable will be excluded from total bnt score if set to True ')
    parser.add_argument('-p', dest='perturbTwoRecods',action="store_true", help='Uses two-record- perturbation function if set to True ')
    parser.add_argument('-sim', dest='SimAnnealing',action="store_true", help='Simulated annealing is used if set to True ')
    parser.add_argument('-i', metavar='iterations',type=int , help='Specify maximum number of iterations ', default=100000)
    parser.add_argument('-t', metavar='thining',type=int , help='Display BDeu Score after iterations ', default=500)
    parser.add_argument('-s', metavar='initialSeed',type=int , help='Specify initial seed. if both initialSeed and loadseed option are not provided then system time will be taken as the default seed  ', default=None)
    parser.add_argument('-l', metavar='loadSeed',type=int , help='Specify path to a file containing previous state', default=None)
    parser.add_argument('-o', metavar='outfile', type=str, help='Specify the file to output the results. ', default= 'counts_bdeu_results.txt')
    args = parser.parse_args()
    
    structureFile   = args.b
    outputFile      = args.o
    dataFile        = args.d
    hiddenConf      = args.dh
    alpha           = args.a
    maxIter         = args.i
    thining         = args.t
    name            = args.x
    cardinality     = args.c
    child1          = args.c1
    child2          = args.c2
    seed            = args.s
    steepestAsent   = args.SteepestAsent
    simAnealFlag    = args.SimAnnealing
    seedFile        = args.l
    decrementValue  = args.dc
    exHiddenBdeuFlag= args.excludeHidBDeu
    pertTowRecFlag  = args.perturbTwoRecods
    
    print "decrementValue %d" % decrementValue
    
    # instanciate RandomSeed object
    rs=RandomSeed()
    
    if seed == None and seedFile == None:
        seed= time.time()
    elif seed != None and seedFile == None:
        rs.setInitialState(seed)
    elif seedFile != None and seed == None:
        state= rs.getSateFromFile(seedFile)
        rNumber.setstate(state)
    
    print "structure: %s" % structureFile
    print "outputFile %s" % outputFile
    print "dataFile %s" % dataFile
    print "Initial Hidden Configurations %s" % hiddenConf
    print "alpha %f" % alpha
    print "maxIter %d" % maxIter 
    print "seed %d" % seed
    
    # read initial structure
    allNodeObjects=readInitialStructure(structureFile)
    
    if hiddenConf != None and dataFile != None:
        print "Error: Specify either data file or initial hidden counts configuration file, not both."
        sys.exit()
    elif hiddenConf != None and dataFile == None:
        df=readInitialHiddenConfig(hiddenConf)
        if len(df.columns)-1 == len(allNodeObjects):
            print "Error: Wrong Initial hidden configuration file"
            sys.exit()
        print df
        totalUniqueObservations= df.shape[0] / 2
    else:
        # read data file
        df=readDataFrame(dataFile)
        if len(df.columns)-1 != len(allNodeObjects):
            print "Error: Wrong input data file"
            sys.exit()
        totalUniqueObservations= df.shape[0] # if we introduce next hidden variable, this variable would be updated
    dfCopy= df.copy()
    
    # update the parent configurations for all variables
    # and the counts associated with the each parent configuration for each value of X
    for n in allNodeObjects:
        getUpdatedQi(allNodeObjects[n])
        populateCounts(allNodeObjects[n])
    # find the BDeu Score for the whole structure
    for n in allNodeObjects:
        nodesBDeuScore.append(getBDeu(allNodeObjects[n], alpha))
    
    print "BDeu Score for Initial Structure without hidden variable: %f" % sum(nodesBDeuScore)
    
    # enter information about hidden variable
    h=addHiddenNode(name, cardinality, child1, child2)
    
    # add hidden variable to the dataframe and  split almost counts equally:
    if hiddenConf == None:
        df=percentageHiddenCoutsSplit(h,df)
        
    # write df to file called initialCountSplit.txt
    #outName= outputFile+'_initialHiddenCountSplit_'+str((datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d-h%H-m%M-s%S')))
    df.to_csv(outputFile+'.initialHiddenCount', sep=',', index=False)
    # populate hidden value counts
    populateCounts(h)
    
    # open file to write the results
    wf= open(outputFile, 'w')
    
    nodesBDeuScore=[]
    # compute the BDeu score again after perturbations
    for n in allNodeObjects:
        node=allNodeObjects[n]

        if node.getParentUpdateFlag() == True or node.getChildrenUpdateFlag() == True: # if true its a child of hidden variable. so, calculate BDeu again 
            print "Node Name: %s, score before: %f" % (node.getName(), node.getLocalBDeu())
            populateCounts(node)
            node.setLocalBDeu(getBDeu(node, alpha))
            print "Node Name: %s, score after: %f" % (node.getName(), node.getLocalBDeu())
        nodesBDeuScore.append(node.getLocalBDeu())
    totalPreviousBDeuScore= sum(nodesBDeuScore)
    if hiddenConf != None:
        print "Initial BDeu Score with Hidden Varialbe: %f" % ( totalPreviousBDeuScore)
    else :
        print "Initial BDeu Score after introduction  Hidden Varialbe: %f" % ( totalPreviousBDeuScore)
    hValues=node.getKvalues().keys()
    for i in xrange(0,len(hValues)-1):
        count=df[df[h.getName()]==hValues[i]].Counts
        for j in count:
            wf.write(str(j)+',')
        del count
    wf.write(str(sum(nodesBDeuScore)) + "\n")
    stateOutFile= outputFile+'.stateSeed.'+ str(seed)
    
    if simAnealFlag == True:
        print "Simulated Anealing starts now"
        firstRowIndex                  = rNumber.randint(0,df.shape[0]-2)
        rs.storeSate(stateOutFile)
        #simulatedAnealing( h, totalPreviousBDeuScore, sIndex, maxIter, outputFile+".sim", decrementValue, alpha )
        e               = totalPreviousBDeuScore                               # Initial state, energy.
        emax            = float('-inf') 
        ebest           = e                                     # Initial "best" solution
        k               = 1                                     # Energy evaluation count.
        kmax            = maxIter
        objCBDeuBestState= allNodeObjects
        objCBDeuOldState = allNodeObjects

        #bestDf          = pd.DataFrame(index=None, columns=None)
        
        dfCurrent       = df.copy()    
        bestDf          = df.copy()
    #    dfCurrent       = df.copy()

        with open(outputFile+".sim", 'w') as wf:
            
            while k < kmax and e > emax:                    # While time left & not good enough
                #print "iterations number: %d" % (k)
                T =    temperature(k, kmax)              # Temperature calculation.
                # randomly choose hidden state zero or one
                num= rNumber.randint(0,1) 
                if num == 0:
                    flag = False
                else:
                    flag = True
                    
                if pertTowRecFlag == True:
                    secondRowIndex = rNumber.randint(0, df.shape[0]-1)
                    while(True):
                        # loop untill first row and second row are different
                        if secondRowIndex != firstRowIndex:
                            break
                        secondRowIndex = rNumber.randint(0, df.shape[0]-1)
                    
                    twoRowsCountPerturbation( h, firstRowIndex, secondRowIndex,decrementValue, flag)
                else:
                    countPerturbation(h, firstRowIndex, decrementValue, flag)     
               
                firstRowIndex=rNumber.randint(0, df.shape[0]-1) # randomly select another record for next iteration
                
                nodesBDeuScore= []
                
                for n in allNodeObjects:
                    if n == h.getName() and exHiddenBdeuFlag == True:
                        continue
                    node=allNodeObjects[n]
                    if node.getParentUpdateFlag() == True or node.getChildrenUpdateFlag() == True: # if true its a child of hidden variable. so, calculate BDeu again
                        print "Node Name: %s, score before: %f" % (node.getName(), node.getLocalBDeu())
                        populateCounts(node)
                        node.setLocalBDeu(getBDeu(node, alpha))
                        print "Node Name: %s, score after: %f" % (node.getName(), node.getLocalBDeu())
                        allNodeObjects[n]= node
                    nodesBDeuScore.append(node.getLocalBDeu())
    
                dagBDeuScore= sum(nodesBDeuScore)
                
                enew = dagBDeuScore                              # Compute its energy.
                #NOTE Inverse logic here using  '<' instead of '>' as in org algo
                acceptprob= probAcceptance(e, enew, T)
                rnum= rNumber.random()
                
                if acceptprob < rnum:# reject the current state 
                    allNodeObjects= objCBDeuOldState          # go back to the old state
                    df = dfCurrent.copy()
                    wf.write("Rejected: Best bdeuscore: %f, Current bdeuscore: %f, proposal bdeuscore: %f, coin: %d , temp: %f, prob: %f rNumber: %f\n" % (ebest, e, enew, num, T, acceptprob, rnum))
                else:  # accept the new state
                    objCBDeuOldState= allNodeObjects
                    e               = enew
                    dfCurrent       = df.copy()
                    wf.write("Accepted: Best bdeuscore: %f, Current bdeuscore: %f, proposal bdeuscore: %f, coin: %d , temp: %f, prob: %f rNumber: %f\n" % (ebest, e, enew, num, T, acceptprob, rnum))
                                                        
                if enew > ebest:                              # Is this a new best?
                    objCBDeuBestState= allNodeObjects
                    bestDf = df.copy()
                    ebest = enew                              # Save 'new neighbour' to 'best found'.
                k = k + 1
                #print "--->iteration  %d " % k                                     # One more evaluation done
                #print "Best bdeuscore: %f and Current bdeuscore %f :" % (ebest, enew)
                
            #print "Best score (%f) count configurations:" % ebest
            #print bestDf
            #timeStamp=str((datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d-h%H-m%M-s%S')))
            bestDf.to_csv(outputFile+'.bestCounts', sep=',', index=False)
            #print "Current score (%f) count configurations:" % e
            #print dfCurrent
            df.to_csv(outputFile+'.currentCounts', sep=',', index=False)
        if exHiddenBdeuFlag == True:
            print "Best BDeu Score without penalty: %f" % ( ebest)
        else:
            print "Best BDeu Score: %f" % ( ebest)
        print "Simulated Anealing Done.."
    elif steepestAsent == True:
        iterations=0
        rRecords=[i for i in xrange(0, df.shape[0]-1)]
        # randomly shuffle the indexes 
        rNumber.shuffle(rRecords)
        for i in rRecords:
            baseDFrame = df.copy()
            #print "baseDFrame"
            #print baseDFrame
            # below loop is for the counts increment and decrement
            for flag in [True, False]:
                
                for j in rRecords:
                    #increment the iteration number
                    iterations +=1                          
                    # perturb the counts in 
                    countPerturbation(h,j, incrementFlag=flag)

                    nodesBDeuScore=[]
                    # compute the BDeu score again after perturbations
                    for n in allNodeObjects:
                        node=allNodeObjects[n]
                
                        if node.getParentUpdateFlag() == True or node.getChildrenUpdateFlag() == True: # if true its a child of hidden variable. so, calculate BDeu again
                            
                            # change the counts of this node according to some criteria i.e.
                            # for new parent configuration, assign the counts such that sum of counts of new parent
                            # configurations is equal to counts of old parent configuration which we split
                            # to get the new parent configuration. 
                            populateCounts(node)
                            node.setLocalBDeu(getBDeu(node, alpha))
                        nodesBDeuScore.append(node.getLocalBDeu())
                    totalCurrentBDeuScore= sum(nodesBDeuScore)
                    #print ("Current BDeuScore %f Previous BDeuScore %f " % (totalCurrentBDeuScore, totalPreviousBDeuScore))
                    # check if current BDeuScore is greater than previous BDeuScore
                    if totalPreviousBDeuScore < totalCurrentBDeuScore:
                        baseDFrame = df.copy()
                        print "Iteration: %d , BDeu Score: %f" % (iterations, totalCurrentBDeuScore)
                        hValues=node.getKvalues().keys()
                        for i in xrange(0,len(hValues)-1):
                            count=df[df[h.getName()]==hValues[i]].Counts
                            for j in count:
                                wf.write(str(j)+',')
                            del count
                        wf.write(str(totalCurrentBDeuScore) + "\n")
                        stateOutFile= 'state_iter_'+str(iterations)+'_initialSeed_'+ str(seed) +'_'+outputFile
                        rs.storeSate(stateOutFile)
                        totalPreviousBDeuScore = totalCurrentBDeuScore
                
                # copy the base dataframe to to data frame since the new base datafame is the one with high BDeuScore
                df = baseDFrame.copy() 
    else: 
        for iterations in xrange(0, maxIter + 1 ): 
            
            # perturb the counts here
            #print "df before perturbation."
            #print df
            randomCountPerturbation(h)
            #print "df after perturbation."
            #print df
            
            nodesBDeuScore=[]
            hiddenValueCountList=[]
            # compute the BDeu score again after perturbations
            for n in allNodeObjects:
                node=allNodeObjects[n]
    
                if node.getParentUpdateFlag() == True or node.getChildrenUpdateFlag() == True: # if true its a child of hidden variable. so, calculate BDeu again
                    
                    # change the counts of this node according to some criteria i.e.
                    # for new parent configuration, assign the counts such that sum of counts of new parent
                    # configurations is equal to counts of old parent configuration which we split
                    # to get the new parent configuration. 
                    populateCounts(node)
                    node.setLocalBDeu(getBDeu(node, alpha))
                nodesBDeuScore.append(node.getLocalBDeu())
            hiddenValueCountList= allNodeObjects[h.getName()].getParentValueCount().values()
            if (iterations % thining) == 0:
                print "Iteration: %d , BDeu Score: %f" % (iterations, sum(nodesBDeuScore))
                hValues=node.getKvalues().keys()
                for i in xrange(0,len(hValues)-1):
                    count=df[df[h.getName()]==hValues[i]].Counts
                    for j in count:
                        wf.write(str(j)+',')
                    del count
                wf.write(str(sum(nodesBDeuScore)) + "\n")
                stateOutFile= 'state_iter_'+str(iterations)+'_initialSeed_'+ str(seed) +'_'+outputFile
                rs.storeSate(stateOutFile)
    wf.close()
       
if __name__== "__main__":
    main(sys.argv[1:])
