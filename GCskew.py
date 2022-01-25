#!/usr/bin/env python
# coding: utf-8

# In[34]:


#!/usr/bin/env python

#set up for Debugging or commandline runs
Debug=False
#import sys library/module
import sys
#import from biopython module
from Bio import SeqIO

#Usage statement:
if(len(sys.argv)<2):
    print("")
    print("This is a script to calculate GC-content and GC-skew along a sequenc")
    print("Usage: arguments.py -i <input file> -w <window size> -o <outfiles base name>")
    #Stop program after the usage statement
    sys.exit()
    
#parameter sets for test case; can use in IDE environment like Jupyter
if Debug:
    InFileName="Diaz11.fasta"
    window = 1000
    OutFileBase = "test"
#for running via commandline
else:
    #Parse my arguments
    for i in range(len(sys.argv)):
        if(sys.argv[i] == "-i"):
            InFileName = sys.argv[i+1]
        elif(sys.argv[i] == "-w"):
            window = int(sys.argv[i+1])
        elif(sys.argv[i] == "-o"):
            OutFileBase = sys.argv[i+1]
            
#write headers to the two outfiles
GCSkewOutFile = open(OutFileBase + "." + str(window) + ".GCSkew.csv", "w")
GCSkewOutFile.write("chromosome,start,stop,GC_skew" + "\n")
PercGCOutFile = open(OutFileBase + "." + str(window) + ".PercGC.csv", "w")
PercGCOutFile.write("chromosome,start,stop,PercGC" + "\n")

#loop through fasta file sequences to create output files
for seq_record in SeqIO.parse(InFileName, "fasta"):
    #calculate number of windows in sequence; cuts off last window if < window size
    intervals = int(len(seq_record.seq)/window)
    #loop through sequence to calculate values for each window and write results to output files
    for i in range(0,intervals*window,window):
        segment=seq_record.seq[i:i+window]
        NumC = segment.count('C')
        NumG = segment.count('G')
        GCSkew = (NumG-NumC)/(NumG+NumC)
        PercGC = float((NumG+NumC))/window
        OutString = "%s,%d,%d,%.2f" % (seq_record.id, i, (i+window), GCSkew)
        GCSkewOutFile.write(OutString + "\n")
        OutString = "%s,%d,%d,%.2f" % (seq_record.id, i, (i+window), PercGC)
        PercGCOutFile.write(OutString + "\n")
GCSkewOutFile.close()
PercGCOutFile.close()


# In[ ]:




