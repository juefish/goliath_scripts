#!/usr/bin/env python
import sys
import os
#import set

#from set import Set

# Usage.
if len(sys.argv) < 2:
        print ""
        print "This program extracts fasta sequences from a file by either exluding or including sequence based on a provide list of sequence names"
        print "Usage: %s -i/e -list file1 -fasta file2"
        print "-list: list of sequence names"
        print "-fasta: fasta file"
        print "-i: include list contents"
	print "-e: exclude list contents"
        print ""
        sys.exit()

# Parse args.
for i in range(len(sys.argv)):
        if sys.argv[i] == "-list":
                infile1 = sys.argv[i+1]
        elif sys.argv[i] == "-fasta":
                infile2 = sys.argv[i+1]
        elif sys.argv[i] == "-i":
                switch = True
	elif sys.argv[i] == "-e":
		switch = False

# get files
fls = [infile1,infile2]
geneContigs = set([])
name = set([])
sequence = str()
proceed = False
fin2 = open(fls[1], "r")
#fin1 = open(fls[0], "r")
prev = str()
counter = 0
#print infile1

# get list of contigs of interest
if(fls[0].count('#')<1):
    fin1 = open(fls[0], "r")
    for line in fin1:
        temp=line.lstrip('>').split()
        geneContigs.add(temp[0])
    fin1.close()
else:
    geneContigs.add(fls[0].strip('"').rstrip('#'))
#print list(geneContigs)

# extract contigs from larger file inclusion
if(switch==True):
    for line in fin2:
        if(line.count('|')>0):
		name = line.lstrip('>').split('|')
	else:
		name = line.lstrip('>').split()
	#print name
        if((str(name[0]) in geneContigs) == True):
	    if len(sequence) > 0:
                print sequence
	        sequence = ""
            proceed = True
	    print ">%s" % (name[0])	
            continue
        
        elif(proceed == True)&(line[0]!='>'):
    	    sequence += line.strip()
	    continue
        else:
            proceed = False
    print sequence

# extract contigs from larger file exclusion
elif(switch==False):
    for line in fin2:
        if(line.count('|')>0):
                name = line.lstrip('>').split('|')
        else:
		name = line.lstrip('>').split()
        if(((str(name[0]) in geneContigs) == False) and (line[0]=='>')):
	    if len(sequence) > 0:
                print sequence
                sequence = ""
            proceed = True
            print ">%s" % (name)
            continue
        
        elif(proceed == True)&(line[0]!='>'):
            sequence += line.strip()
            continue
        else:
            proceed = False
    print sequence    
fin2.close()

