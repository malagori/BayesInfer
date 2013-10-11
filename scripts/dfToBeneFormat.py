#!/usr/bin/env python
'''
Created on Oct 3, 2013

@author: malagori
'''

import sys

def convertDfToBene(infile):
    '''
    This function convert data frame to bene's data format
    '''
    count=1
    with open(infile, 'r') as rf:
        with open(str(infile+'.bene'), 'w') as wf:
            tokens=[]
            for line in rf:
                if count != 1:
                    line.strip()
                    tokens= line.split(',')
                    record=''
                    for i in xrange(0, len(tokens)-2):
                        record+=str(tokens[i])+'\t'
                    record=record+str(tokens[-2])+'\n'

                    for i in xrange(0, int(tokens[-1].strip())):
                        wf.write(record)
                count+=1
             
def main(argv):
    
    if not argv:
        print "usage: dfToBeneFormat <path-to-data-frame-file>" 
        sys.exit()
    print argv[0]
    convertDfToBene(argv[0])
           

if __name__ == '__main__':
    main(sys.argv[1:])