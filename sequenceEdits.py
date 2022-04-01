#!/usr/bin/env python
# coding: utf-8

# Program to facilitate manual editing of COVID-19 sequence.
# The program uses MAFFT and read base coverage to provide additional
# evidence to users to investigate evidence support the resolution of
# ambigous bases, predicted indels, and variable sites.
# Written in Python3 and require Biopython and pysam packages to be installed

#Usage statement: python COVID_seq_edits.py

#import in packages needed for script

from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from datetime import datetime
from Bio.Align.Applications import MafftCommandline
from io import StringIO
import pysam, os, subprocess
from pysam import VariantFile
import numpy as np

#Set variables for program
nucs=["A","C","G","T","N","-","R","Y","S","W","K","M","B","D","H","B"]
variants={}
Debug=False

def summarize_alignment(align):
    NumN = align[1,:].seq.count('N')
    print("Number of ambiguous bases associated with each reference nucleotide type:")
    print("A->N:", int(align.substitutions['N','A']*2))
    print("C->N:", int(align.substitutions['N','C']*2))
    print("G->N:", int(align.substitutions['N','G']*2))
    print("T->N:", int(align.substitutions['N','T']*2))
    print("Total Number of Ambiguous bases in genome assembly:", NumN)
    print("Total Length of alignment:", align.get_alignment_length())
    print("Total Number of Variable Sites:", len(variants))
    print("Total Number of Non-ambiguous Variable Sites:", len(variants)-NumN)

# Code adapted from "Base count extractor" by david.eyre@ndm.ox.ac.uk
# Method to count number of bases from reads mapped to a specific position
def Bases_At_Pos(samfile, pos, chromname, minbasequal, minmapqual):
    'Return a string of the bases at that position.'
    position = 0
    coverage = 0
    bases = ""
    for pileupcolumn in samfile.pileup(reference=chromname, start=pos-1, end=pos):
        if ((pileupcolumn.pos+1) >= pos and (pileupcolumn.pos+1) <= pos):
            position = int(pileupcolumn.pos+1)
            coverage = int(pileupcolumn.n)
            for pileupread in pileupcolumn.pileups:
                if (pileupread.is_del == 0):
                    #print(pileupread.query_position)
                    bases += pileupread.alignment.seq[pileupread.query_position]
    return position, coverage, bases

# Method to summarize pileup base coverage and make base prediction
# Key points: this methods will make a nucleotide prediction IF total
# read coverage is >100 and the percentage of a base is >85%. This is
# somewhat arbitrary, so feel free to adjust as you see fit. Also, the
# method will suggest to manually investigate ANY indel in the query
# sequence. Low coverage variable sites are identified with "???".
def base_prediction(line, variants):
    global totalBases
    global predictedNuc
    for i in range(len(line)):
        if(i>0):
            totalBases+=int(line[i])
    for i in range(len(line)):
        if(i>0):
            if(totalBases>=100):
                if((int(line[i])/float(totalBases))>=0.85 and variants[int(line[0])][1].upper()!='-'):
                    if(i==1):
                        predictedNuc = 'A'
                    elif(i==2):
                        predictedNuc = 'C'
                    elif(i==3):
                        predictedNuc = 'G'
                    elif(i==4):
                        predictedNuc = 'T'
                elif((int(line[i])/float(totalBases))<0.85 and predictedNuc==""):
                    predictedNuc = 'Investigate'

            elif(totalBases < 25):
                predictedNuc="???"
            else:
                predictedNuc = 'Investigate'

#Create file objects and import genome sequences
tempFile = open("Temp.fa", 'w')
if Debug:
    reference=SeqIO.read("MN908947.3.fasta", 'fasta')
else:
    reference=SeqIO.read(input("Reference file?"), 'fasta')
if Debug:
    query=SeqIO.read("MC01975.medaka.consensus.fasta", 'fasta')
else:
    query=SeqIO.read(input("Query file?"), 'fasta')
seq_list= [reference.upper(), query.upper()]
#Combine sequences and write to a temporary file
SeqIO.write(seq_list,tempFile, 'fasta')
tempFile.close()

#run MAFFT to generate multiple sequence alignment
mafft_exe = "/usr/bin/mafft"
inFile="Temp.fa"
mafft_cline = MafftCommandline(mafft_exe, input=inFile)
stdout, stderr = mafft_cline()
align = AlignIO.read(StringIO(stdout.upper()), "fasta")
query_seq=align[1,:]
#edit_seq=query.seq.tomutable()
subprocess.call(["rm", "Temp.fa"])
for i in range(align.get_alignment_length()):
    if("n" in align[:,i] or "-" in align[:,i] or align[:,i][0]!=align[:,i][1]):
        variants[i+1]=align[:,i]
# Summarize the variant/variable information
summary = input("Would you like to summarize data? Y/N ")
while(summary!='Y' and summary!='N'):
    summary=input("Please type either a Y or an N")
if(summary=='Y'):
    summarize_alignment(align)

align_array = np.array([list(rec) for rec in align], np.str_)

# Allow you to generate a predicted variant call based off of the multiple sequence alignment and bam/read coverage
prediction = input("Would you like to generate predicted base calls for variable sites? Y/N ")

while(prediction!='Y' and prediction!='N'):
    prediction=input("Please type either a Y or an N")
if(prediction=='Y'):

    #numeric list of positions of variable sites in reference genome, numbered from zero
    POSLIST = list(variants.keys())
    POSLIST.sort()

    #path to bam files for extraction of base counts
    BAMPATH = input('Location of folder where bam files for base count extraction: ')

    #file name for bam file, without .bam suffix, these will be recycled for output name
    file = input("File name of bam file without suffix to be used for read counting: ")

    #output path
    OUTPATH = ''

    #sets reference genome ID/name; can change if using different reference
    ref = "MN908947.3"
    #loop over each of the bam files
    inbampath = BAMPATH+file+'.bam'
    #will move to next file name if can't find bam file
    if not os.path.exists(inbampath):
        #skip this file
        print("BAM file not found")
    else:
        print("Generating pileup data for %s..." % (inbampath))
    #will create bam index file if it doesn't exist
    if not os.path.exists(inbampath+".bai"):
        pysam.index(inbampath)
    #open bam file with pySAM library
    samfile = pysam.Samfile(inbampath, 'rb')
    #open output file and write headers
    with open (OUTPATH+file+'_variant_sites.txt', 'w') as f:
        f.write('site\tA\tC\tG\tT\n')
        #for each position of interest write base counts
        for pos in POSLIST:
            hq_position, hq_coverage, hq_bases = Bases_At_Pos(samfile, int(pos), ref, 30.0, 30.0)
            hq_basecounts = [ len([x for x in hq_bases if x == 'A']), len([x for x in hq_bases if x == 'C']),                              len([x for x in hq_bases if x == 'G']), len([x for x in hq_bases if x == 'T']) ]
            if(len([i for i in hq_basecounts if i>0])>=1):
                f.write('%s\t%s\n'%(hq_position, '\t'.join([str(i) for i in hq_basecounts])))
        f.close()
print("Done with pileup...")

# Option to allow for guided editing process
edit = input("Would you like to edit as you go? Y/N ")
while(edit!="Y" and edit!='N'):
    edit = input("Please type either a Y or an N")
# Option to output data on all variable sites or just the ones that are recommended for investigation
output = input("Would you like to output all results to a file? \nOtherwise, it will just write items that need to be investigated to a file. Y/N ")
while(output!="Y" and output!='N'):
    output = input("Please type either a Y or an N")

lineNumber=0
tracker=0
prevSite=0

# examine results from pileup of reads, summarize data and make recommendations on variable sites
pile_out = open(file+"_variant_sites.txt", "r")
fout = open(file+"_output.txt", 'w')
fout.write("#Ref\tVariant\tSite\tNumReads\tA\tC\tG\tT\tPredictedCall\n")

logFile = open(file+".log", 'w')
now=datetime.now()
dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
logOut = """#Log of changes made to %s sequence.
#Date and time = %s
""" % (query.id, dt_string)
logFile.write(logOut)
header="Ref\tVariant\tSite\t#Reads\tA\tC\tG\tT\tPredictedCall"

for line in pile_out:
    if(lineNumber>0):
        totalBases=0
        line=line.split()
        predictedNuc = ""
        base_prediction(line, variants)
        OutString="%s\t%s\t%s\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%s" % (variants[int(line[0])][0].upper(),                                                                  variants[int(line[0])][1].upper(),                                                                   line[0], totalBases, int(line[1])/totalBases,                                                              int(line[2])/totalBases,                                                              int(line[3])/totalBases,                                                              int(line[4])/totalBases,                                                         predictedNuc)
        if((int(line[0])-prevSite)==1):
            tracker+=1
        else:
            tracker=0
        if(tracker>=1):
            OutString=OutString + "\nWARNING: Possible Indel; May Need Investigation"
        prevSite=int(line[0])
        start=int(line[0])
        if(edit=="Y"):
            print(header +"\n" + OutString)
        if(edit=='Y' and predictedNuc!="???"):
            confirm = input("Confirm that you would like to change the %s at position %s: Y/N " % (variants[int(line[0])][1].upper(), line[0]))
            while(confirm!="Y" and confirm!='N'):
                confirm = input("Please type either a Y or an N")
            if(confirm=='Y'):
                change = input("Would you like to change this call to the predicted %s? Y/N " % predictedNuc)
                while(change!="Y" and change!='N'):
                    change = input("Please type either a Y or an N")
                if(change=='Y'):
                    print("BEFORE:", align_array[1,start-1])
                    align_array[1,start-1]=predictedNuc
                    print("AFTER:", align_array[1,start-1])
                    logOut="%s at position %s -> %s" % (variants[int(line[0])][1].upper(),start,predictedNuc)
                    logFile.write(logOut + "\n")
                else:
                    alt = input("Is there a different call you would like to change this position to? Y/N ")
                    while(alt!="Y" and alt!='N'):
                        alt = input("Please type either a Y or an N")
                    if(alt=="Y"):
                        newNuc = input("What other bases call would you like to make at this position? A, C, G, T, or N? ")
                        while(newNuc not in nucs):
                            newNuc = input("Please choose a known base call code.")
                        print("BEFORE:", align_array[1,start-1])
                        align_array[1,start-1]=newNuc
                        print("AFTER:", align_array[1,start-1])
                        logOut="%s at position %s -> %s" % (variants[int(line[0])][1].upper(),start,newNuc)
                        logFile.write(logOut + "\n")
            elif(output=='Y'):
                fout.write(OutString + '\n')

        elif(edit=='Y' and output=='N' and (predictedNuc=="Investigate" or predictedNuc=="???")):
            fout.write(OutString + '\n')
        elif(output=='Y'):
            fout.write(OutString + '\n')
    lineNumber+=1
pile_out.close()
fout.close()


# Edit sequence "by hand" - i.e. target specific sites manually
looper = True
while(looper==True):

    byhand = input("Are there any edits you wish to make by hand? Y/N ")
    while(byhand!="Y" and byhand!='N'):
        byhand = input("Please type either a Y or an N")
    if(byhand=="Y"):
        print("YOU MUST EDIT YOUR SEQUENCE FROM THE RIGHT [i.e. highest numbered position] TO THE LEFT [i.e. lowest numbered position]")
        start=int(input("Start site on reference sequence? "))
        action=input("Action to take? C(hange),D(elete),I(nsert), or E(nd manual editing) ")
        while(action not in ['C','D', 'I','E']):
            action=input("Unrecognized option entered: Action to take? C,D, or I ")
        print("The current base call at position %s is %s" % (start, align_array[1,start-1]))
        if(action=="C"):
            new=str(input("What would you like to change this position to?"))
            align_array[1,start-1]=new
            logOut="%s at position %s -> %s" % (variants[int(line[0])][1].upper(),start,new)
            logFile.write(logOut + "\n")
        elif(action=="D"):
            indel=int(input("How many bases would you like to delete? "))
            align_array=np.concatenate([align_array[:,:(start-1)],align_array[:,(start+indel-1):]],axis=1)
            logOut="%s nucleotides were deleted from the query sequence starting at position %s" % (indel,start)
            logFile.write(logOut + "\n")
        elif(action=="I"):
            # It will insert after the position you provide
            nucl=input("What bases would you like to insert? ")
            for i in nucl:
                align_array=np.insert(align_array,start,i,axis=1)
            logOut="%s was inserted into the query sequence starting at position %s" % (nucl,start)
            logFile.write(logOut + "\n")
        else:
            continue
    if(byhand=="N"):
        looper=False
logFile.close()

#write final edited sequence to new file
edit_seq=""
for nuc in align_array[1,:]:
    edit_seq+=nuc
fasta=open(file+".fasta","w")
fasta.write(">"+query.id+"\n"+edit_seq)
fasta.close()

