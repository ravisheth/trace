#!/usr/bin/env python
# id uniquely mapping spacers from blast results
import os, sys
import pandas as pd
from multiprocessing import Pool as ThreadPool 

num_threads = 4 # Number of threads for later multithreading

### user input ###
try:
    folder = str(sys.argv[1]) # direct to read and write results from
except IndexError:
    print "Usage:"
    print "./unique_spacers.py [working_directory]"
    exit(1)

if not os.path.isdir(folder):
    raise ValueError, "Directory " +str(folder)+ " not found."

### main script body ###
file_list = []
for file_name in os.listdir(folder):
    if file_name.endswith("blast_out.txt"):
        file_list.append(file_name)

def get_uniq_spacers(file_name):
    col=['subject','%identity','alignment_len','mismatches','gap_opens','qstart','qend','sstart','send','evalue','bitscore']
    df=pd.read_csv(folder+'/'+file_name,names=col,index_col=0)
    uniq_array=[] # array id
    uniq_id=[] # reference aligned to
    uniq_start=[] # alignment start bp
    uniq_end=[] # alignment end bp
    current_sp_store = [] # variable to store current spacer alignments
    for i, r in df.iterrows():
        if i in [ind.name for ind in current_sp_store]: # open all alignments of spacer sequence
            current_sp_store.append(r)
        else:
            if len(current_sp_store) == 1: # if current spacer is uniquely mapping
                uniq_array.append(current_sp_store[0].name)
                uniq_id.append(current_sp_store[0]['subject'])
                uniq_start.append(current_sp_store[0]['sstart'])
                uniq_end.append(current_sp_store[0]['send'])
            current_sp_store = []
            current_sp_store.append(r)
    if len(current_sp_store) == 1: # catch last spacer that was opened
        uniq_array.append(current_sp_store[0].name)
        uniq_id.append(current_sp_store[0]['subject'])
        uniq_start.append(current_sp_store[0]['sstart'])
        uniq_end.append(current_sp_store[0]['send'])

    uniq_file_name = file_name.split('_blast_out.txt')[0]+'_uniq.txt' #write output
    uniq_out = open(folder+'/'+uniq_file_name,'w')

    # write output; sname,reference,start,stop
    for i,j,k,l in zip(uniq_array,uniq_id,uniq_start,uniq_end):
        uniq_out.write(str(i)+','+j+','+str(k)+','+str(l)+'\n')
    uniq_out.close()

pool = ThreadPool(num_threads)
res  = pool.map(get_uniq_spacers, file_list)