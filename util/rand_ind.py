#!/usr/bin/env python

# Quick script to randomly generate temporal induction profiles for given number of conditions, days, samples
# WARNING: This script isn't guaranteed to print out unique sets of sample strings

import random, sys
import numpy as np

### Helper funcs ###

# replace 0's with X's, 1's with O's (letter 0)
def replace_bits(bit_str):
    return bit_str.replace('0','X').replace('1','O')

### Take in command args ###
try:
    num_conditions = int(sys.argv[1]) # No. of conditions
    num_days = int(sys.argv[2]) # Number of propagation days
    num_samps = int(sys.argv[3]) # Number of samples to follow
except IndexError:
    print "Usage:"
    print "./rand_44_str.py <num_conditions> <num_days> <num_samps>"
    exit(1)

seed = 1357 # Seeding for replicability
random.seed(seed)

### Generate random strings based on params ###

frmt_str = "{:0" + str(num_conditions) + "b}"
samp_list = np.empty(shape=(num_samps, num_days), dtype='S'+str(num_conditions))
uniq_test_list = [''] * num_samps

for i in range(num_samps): # iterate over samples
    for j in range(num_days): # iterate over days
        tmp_rand = random.getrandbits(num_conditions)
        samp_list[i,j] = replace_bits(frmt_str.format(tmp_rand))
        uniq_test_list[i] = uniq_test_list[i] + samp_list[i,j]

uniq_index = np.unique(np.array(uniq_test_list))
if len(uniq_index) == num_samps:
    footer_line = 'all samples are unique'
else:
    footer_line = 'WARNING: Samples are not unique.'

### Print to stdout ###
print '### rows= days; cols= condition ###'
for i in range(num_samps): # iterate over samples
    for j in range(num_days): # iterate over days
        print samp_list[i,j]
    print ''
print footer_line
