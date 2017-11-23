#!/usr/bin/env python
# extraction of spacer sequences from amplicon sequencing of CRISPR arrays
import os, sys
from multiprocessing import Pool as ThreadPool

### helper functions ###
def hamming(a, b): # hamming distance between two sequences
    dist = 0
    for i, j in zip(a, b):
        if i != j:
            dist += 1
    return dist

def min_hamming(ref_list, s): # find match in list with minimum hamming dist to string
    if len(ref_list) == 1: # if there is only one item in the input list
        return hamming(ref_list,s)
    if len(ref_list) > 1:
        minimum=len(s)
        match=0
        i=0
        for bc in ref_list: # iterate through list
            d=hamming(bc,s)
            if d < minimum: # record minimum hamming distance
                minimum = d
                match = i
            i=i+1
        return minimum, match

### algorithm constants ###
dr_read1_len = 12 # bp (from 3') of direct repeat sequenced at beginning of read
spacer_len=[32,33,34] # possible acquired spacer sizes in E. coli array
dr_seq = 'GTGTTCCCCGCGCCAGCGGGGATAAACC' # direct repeat seq
s1_seq = 'GAGCACAAATATCATCGCTCAAACC' # first spacer seq (truncated due to seq. primer)
num_threads = 4 # Number of threads for later multithreading

hamming_dist = 2 # max hamming for matching direct repeat or spacer

### user input ###
try:
    folder = str(sys.argv[1]) # directory to read from
    out_dir = str(sys.argv[2]) # directory to write to
    if len(sys.argv) > 3: # optionally set a custom DR sequence
        dr_seq = str(sys.argv[3])

except IndexError:
    print "Usage:"
    print "./spacer_extraction.py [fastq_directory] [out_directory] [DR_sequence (optional)]"
    exit(1)

dr_min_match = len(dr_seq) # bp (from 5') of direct repeat required for matching
# dr_min_match = 15 #if trying to extract up to 5 spacers, because of read length limitations

if not os.path.isdir(out_dir):
    raise ValueError, "Directory " +str(out_dir)+ " not found."

### main script body ###
# iterate through files and write out extracted arrays
file_list = []
for file_name in os.listdir(folder):
    if file_name.endswith(".fastq"):
        file_list.append(file_name)

def extract_spacers(file_name):
    line=0
    header=[]
    seq=[]
    qual=[]

    # read fastq, store header, sequence, quality into lists
    for l in open(os.path.join(folder,file_name)):
        line += 1
        if line%4 == 1: header.append(l.split('@')[1].split('.')[0])
        if line%4 == 2: seq.append(l.rstrip())
        if line%4 == 0: qual.append(l.rstrip())

    # loop to extract spacer sequences from reads
    extract=[]
    for i in seq:
        read = i
        spacer_list=[] # list for extracted spacers
        exit='' # store status of read extraction

        #check beginning of read matches expected DR sequence
        if hamming(read[:dr_read1_len],dr_seq[-dr_read1_len:]) > hamming_dist:
            exit='FAIL:no_read_match_DR_seq_leader'
            extract.append(exit)
            continue
        read = read[dr_read1_len:] #strip dr sequence

        while True:
            # check if immediate sequence matches spacer 1 (unexpanded array)
            if hamming(read[:len(s1_seq)],s1_seq) < hamming_dist:
                exit='SP1_MATCH'
                break

            # check to see if there is enough of the read left to do matching
            if len(read) < (max(spacer_len)+dr_min_match):
                exit='FAIL:could_not_match_DR_not_enough_length'
                break

            # find a matching DR assuming finite spacer length
            pot_dr = [read[spacer_len[i]:spacer_len[i]+dr_min_match] for i in range(3)]
            dist, match = min_hamming(pot_dr,dr_seq[:dr_min_match])
            if dist > hamming_dist:
                exit='FAIL:DR_match_not_made'
                break

            spacer_list.append(read[:spacer_len[match]])
            read = read[spacer_len[match]+len(dr_seq):]
        # if no spacers were found, determine if unexpanded array
        if len(spacer_list) == 0:
            if exit == 'SP1_MATCH':
                extract.append('ARRAY:0')
            else:
                extract.append(exit)
        # if spacers were found, determine number of expansions
        else:
            if exit == 'SP1_MATCH':
                extract.append('ARRAY:'+str(len(spacer_list))+","+",".join(spacer_list))
            else:
                extract.append(exit+'_'+'ARRAY:'+str(len(spacer_list))+","+",".join(spacer_list))

    cnt_array=[0,0,0,0,0,0] # tracking for up to 5 expanded spacers
    extract_qc=[] # extracted reads with no errors
    for i,s in zip(extract,seq):
        if "FAIL" not in i: # no extraction errors
            if "ARRAY" in i:
                num=int(i.split('ARRAY:')[1].split(',')[0])
                if num > 0:
                    if num == 1 :
                        extract_qc.append([i.split('ARRAY:1,')[1]])
                    else:
                        extract_qc.append(i.split('ARRAY:'+str(num)+',')[1].split(','))
                cnt_array[num]+=1

    # write out identified spacer sequences to fasta file
    out_file_name = file_name.split('_')[0].split('.')[0]+'_out.fa'
    fasta_out = open(os.path.join(out_dir, out_file_name),'w')
    n=1
    sp_cnt=0
    for array in extract_qc:
        s=1
        array_len=len(array)
        for array_seq in array:
            # spacer names: >a<array_id>l<array_length>p<spacer_pos>
            fasta_out.write('>a'+str(n)+'l'+str(array_len)+'p'+str(s)+'\n')
            fasta_out.write(array_seq+'\n')
            sp_cnt+=1
            s+=1
        n+=1
    fasta_out.close()

    # write output stats: sname,reads,extracted,len0,len1,len2,len3,len4,len5
    out_text=file_name.split('_')[0].split('.')[0]+','\
        +str(len(seq))+','+str(sum(cnt_array))+','+str(cnt_array)[1:-1]+'\n'
    with open(out_dir+'/'+'extraction_stats.csv','a+') as stats:
        stats.write(out_text)

pool = ThreadPool(num_threads)
res  = pool.map(extract_spacers, file_list)
