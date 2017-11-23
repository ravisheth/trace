#!/usr/bin/env python

# Import libraries
import os, pickle, itertools
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style("white")
from itertools import compress
from scipy.spatial.distance import euclidean

##############################
## S.01 ## Helper functions
##############################

# function searches a str for pattern, returns list of
# all beginning indices where the pattern is found
def find_in_str(in_str, patt):
    out_index_list = []
    i = -1
    while True:
        try:
            i = in_str.index(patt, i+1)
            out_index_list.append(i)
        except ValueError:
            break
    return out_index_list # returns empty list if patt not found anywhere

# Check if variable is a float or np float
def is_float(x):
    lgcl = (type(x) == float) or (type(x) == np.float64) or \
           (type(x) == np.float32) or (type(x) == np.float16)
    return lgcl

# def dict_sum(x): # Given a dict, assume that all the entries are ints or floats; sum them
#     tot = 0
#     for k in x.keys():
#         tot += x[k]
#     return tot

# Given two dicts, merge them into a new dict as a shallow copy.
def merge_two_dicts(x, y):
    z = x.copy()
    z.update(y)
    return z

# Given an array-like of sequence origin names, output the str corresponding to their assignments
# using the provided character pair
# E.g., ['bl21de3', 'psb2k3min', 'pcas12'] -> 'XOX'
def classify_array_seq(input_array, char_pair=('X','O')):
    output=''
    for i in input_array: # genome or cas plasmid
        if i == 'bl21de3' or i == 'pcas12laci' or i == 'pcas12': # 'pcas12' in the new reference fa
            output = output+char_pair[0]
        if i == 'psb2k3min': # plasmid-derived
            output = output+char_pair[1]
    if len(output) == len(input_array):
        return output
    else:
        raise ValueError("An element of 'input_array' was not a recognized sequence origin string")

# Generate list of all sample strings of len k (e.g., 'XXXX', 'XOXO', etc. for k=4)
def make_sample_list(k, char_pair=('X','O')):
    
    frmt_str = "{:0" + str(k) + "b}"
    samples = [None] * (2**k)
    for i in range(2**k):
        samples[i] = frmt_str.format(i).replace('0',char_pair[0]).replace('1',char_pair[1])
    return samples

# For an array or state str x (i.e., x has only 2 unique chars),
# replace the old char pair with the new char pair
def replace_char_pair(x, old=('X','O'), new=('G','P')):
    for i in range(2):
        x = x.replace(old[i], new[i])
    return x
double_replace = replace_char_pair # Old name of replace_char_pair, preserved for compatibility

#############################################################################
## S.02 ## Functions for extracting arrays, plotting array vector distributions
#############################################################################

# function for parsing a spacer string (e.g., a3l2p1) into the 3 component ints
def parse_array_str(array_str):
    out_list = []
    for i in ['p','l','a']: # we read through backwards
        tmp = array_str.split(i)
        array_str = tmp[0]
        out_list.append(int(tmp[1]))
    out_list.reverse()
    return out_list # returns values [a,l,p]


# function for counting up numbers of array types from provided spacer data for given samples
# takes files of spacer strings (e.g., a3l2p1), assembles the arrays, outputs freq counts of the
# arrays for given sample(s)
#
# Returns a sample-keyed dictionary of array-length-keyed dictionaries of array-type counts as pd.Series objects
#
# Can takes some time to run for many samples
#
# e.g., a3l2p1 corresponds to position 1 in array no. 3, which has length 2
def array_counts_from_spacers(folder, # path to dir to look in
                              samples, # list of samples (expecting file name formatted as <samp>_uniq.txt), case insensitive
                              do_dedup=False, # Remove duplicate arrays
                              verbose=False, # Verbose printing only implemented for dedup right now
                              char_pair=('X','O') # How to represent the arrays (reference char, pTrig char)
                             ):
    max_array_len = 5 # maximum array len to consider; determined by sequencing scheme
    array_lens = range(1, max_array_len+1)
    freq_tables_raw = dict.fromkeys(samples)
    for title in samples:
        freq_tables_raw[title] = dict.fromkeys(array_lens)
        for i in array_lens: # construct tables
            # empty_table = dict.fromkeys(make_sample_list(i), 0)
            empty_table = pd.Series(0, index=make_sample_list(i, char_pair=char_pair),
                                    dtype=int)
            freq_tables_raw[title][i] = empty_table

        #read the uniq file for current sample, parse the line
        uniq_df = pd.read_csv(os.path.join(folder, title+'_uniq.txt'), sep=',', header=None)
        uniq_df.columns = ['array','id','start','end']
        uniq_df.sort_values('array', inplace=True)

        cur_identity = None
        id_store = []

        if do_dedup: # remove duplicates
            all_uniq_array_strs = [] # list of all unique spacer strs <spacer_id><start_pos>,<end_pos>...
            num_uniq = dict.fromkeys(array_lens, 0) # keep track of number of unique arrays counted (per len)
            num_dups = dict.fromkeys(array_lens, 0) # keep track of number of dups removed (per len)
            for row in uniq_df.as_matrix():
                array_vals = parse_array_str(row[0])
                ar_identity = array_vals[0]

                if cur_identity != ar_identity: # starting a new array block
                    cur_identity = ar_identity
                    id_store = [None] * array_vals[1]
                    cur_array_str = ''

                id_store[array_vals[2]-1] = row[1] # append new id
                cur_array_str += row[1]+str(row[2])+','+str(row[3])

                # if complete array AND unique array, classify and add to table
                if not (None in id_store):
                    if not (cur_array_str in all_uniq_array_strs):
                        all_uniq_array_strs.append(cur_array_str)
                        seq_str = classify_array_seq(id_store, char_pair=char_pair)
                        freq_tables_raw[title][array_vals[1]][seq_str] += 1
                        num_uniq[array_vals[1]] += 1
                    else: # found a complete array but it's not unique
                        num_dups[array_vals[1]] += 1
                    cur_identity = None # assign cur_identity to None to force new block
            if verbose: # verbose output
                print 'sample {}:'.format(title)
                fmt_str = '    array length {:>2} dup rate  = %{:01.1f} ({:d}/{:d})'
                for arr_l in array_lens:
                    num_total = num_dups[arr_l] + num_uniq[arr_l]
                    try:
                        perc_dup = 100*num_dups[arr_l]/float(num_total)
                    except ZeroDivisionError:
                        perc_dup = 0
                    print fmt_str.format(str(arr_l), perc_dup, num_dups[arr_l], num_total)
                num_dups_sum = sum(num_dups)
                num_uniq_sum = sum(num_uniq)
                num_total = num_dups_sum + num_uniq_sum
                try: # print the overall dup percentage
                    perc_dup = 100*num_dups_sum/float(num_total)
                except ZeroDivisionError:
                    perc_dup = 0
                fmt_str = '    {:>15} dup rate  = %{:01.1f} ({:d}/{:d})'
                print fmt_str.format('overall', perc_dup, num_dups_sum, num_total)
        else: # DON'T remove duplicates
            for row in uniq_df.as_matrix():
                array_vals = parse_array_str(row[0])
                ar_identity = array_vals[0]

                if cur_identity != ar_identity: # starting a new array block
                    cur_identity = ar_identity
                    id_store = [None] * array_vals[1]

                id_store[array_vals[2]-1] = row[1] # append new id

                if not None in id_store: # if complete array, classify and add to table
                    seq_str = classify_array_seq(id_store, char_pair=char_pair)
                    freq_tables_raw[title][array_vals[1]][seq_str] += 1
                    cur_identity = None # assign cur_identity to None to force new block
    return freq_tables_raw
compute_array_vectors = array_counts_from_spacers # Old name of array_counts_from_spacers, preserved for compatibility

# function for plotting all array vectors of a given length for each sample
def plot_array_vectors(array_vectors_raw, sim_vectors=None, array_len=None, samples=None,
                       figsize=(16,12), subplot_width=4):
    plt.figure(figsize=figsize) # prep plot
    w = 0.45 # width
    if samples is None: # Plot all samples in array_vectors if none given
        samples = array_vectors_raw.keys()
        samples.sort(); samples.reverse() # We like to sort in reverse 'XXXX', 'XXXO', etc
    else:
        samples=samples[:]

    do_plot_sim = (sim_vectors != None)

    if array_len is None:
        raise ValueError('array_len cannot be None')

    for i in range(len(samples)):
        title = samples[i]
        if i==0: # get the array type strings
            array_types = array_vectors_raw[title][array_len].index.tolist()

        # Obtain vectors of given array length
        array_type_vec = array_vectors_raw[title][array_len].values
        arr_plot_vec = array_type_vec/float(sum(array_type_vec))

        if do_plot_sim:
            sim_type_vec = sim_vectors[title][array_len].values
            sim_plot_vec = sim_type_vec/float(sum(sim_type_vec))

        subplot_row_n = int( np.ceil( len(samples)/float(subplot_width) ) )
        plt.subplot(subplot_row_n, subplot_width, i+1)

        arr_bar = plt.bar(np.arange(2**array_len),   arr_plot_vec, width=w, align='center', color='k')
        if do_plot_sim:
            sim_bar = plt.bar(np.arange(2**array_len)+w, sim_plot_vec, width=w, align='center', 
                              color='white', edgecolor='black', linewidth=1)
            plt.legend((arr_bar[0], sim_bar[0]), ('data','model'))
        plt.ylim([0,1])
        plt.xlim([ -0.5, 2**array_len])
        ticks=plt.xticks(np.arange(2**array_len)+0.5/array_len, array_types)
        plt.xticks(rotation=70)
        plt.ylabel('frequency')
        plt.title(title + ', L{}={}'.format(array_len, sum(array_type_vec))) # Total num arrays of given len
    plt.tight_layout()
    plt.show()
plot_array_vectors2 = plot_array_vectors  # Old name of plot_array_vectors, preserved for compatibility

# function for extracting all array vectors for a given sample and array length from a 
# sample-keyed dictionary of raw array counts (i.e., the output of array_counts_from_spacers)
def extract_array_vectors(array_vectors, array_len, samples=None, norm_vectors=True, out_type='np_list'):
    get_array_types = True
    if samples is None:
        samples = array_vectors.keys()
    else:
        samples = samples[:]
    # samples.sort(); samples.reverse() # We like to sort in reverse 'XXXX', 'XXXO', etc
    out_list = []
    for i in range(len(samples)):
        title = samples[i]
        if get_array_types: # get the array type strings
            array_types = array_vectors[title][array_len].index.tolist()
            get_array_types = False
        # Obtain vectors of given array length
        array_type_vec = array_vectors[title][array_len].values
        if norm_vectors: # normalize vector if needed
            array_type_vec = array_type_vec/float(sum(array_type_vec))
        out_list.append(array_type_vec)
    if out_type == 'np_list': # list of numpy arrays
        return out_list
    elif out_type == 'pd_df': # pandas dataframe
        out_df = pd.DataFrame(out_list, index=samples, columns=array_types)
        return out_df

#####################################################################################################
## S.03 ## Functions for calculating the probability distributions of arrays given a state sequence
#####################################################################################################

# function for calculating the probability of spacer incorporation patterns for a given state OFF/ON (i.e., X/O)
def calc_single_state_probs(state, # Single character representing state: X=OFF or O=ON)
                            o_exp, # Probability of expansion during ON state
                            x_exp, # Probability of expansion during OFF state
                            o_plasmid, # Probability of incorporating pTrig-derived spacer during ON state
                            x_plasmid, # Probability of incorporating pTrig-derived spacer during OFF state
                            s_prob,    # Scaling probability for multiple incorporations during state
                            n_incorp,  # Number of incorporations during given state
                            log_p=True # Whether to output log probabilities
                           ):
    gp_spacers = ['G', 'P']
    spacer_types = ['N'] + gp_spacers # types of spacers ([N]one, [G]enomic, [P]lasmid)

    # Get all unique spacer patterns with length <= n_incorp, truncate at first N spacer
    spacer_patts = [''.join(x) for x in  list(itertools.product(spacer_types, repeat=n_incorp))]
    for i, patt in enumerate(spacer_patts):
        try:
            spacer_patts[i] = patt[0:patt.index('N')+1]
        except ValueError:
            pass
    spacer_patts = list(set(spacer_patts))

    # use log since these are small probabilities
    patt_dict = dict.fromkeys(spacer_patts, float(0)) # complete expansion probs for each spacer pattern
    init_dict = dict.fromkeys(spacer_types, float(0)) # initial expansion probs
    trns_dict = dict.fromkeys(spacer_types, float(0)) # transition probs

    # Set appropriate probabilities for current state
    if state.lower() == 'o':
        p_exp     = o_exp
        p_plasmid = o_plasmid
    elif state.lower() == 'x':
        p_exp     = x_exp
        p_plasmid = x_plasmid
    else:
        print "unknown state passed"
        return None

    # Fill in initial expansion probs
    init_dict['N'] = np.log(1-p_exp)
    init_dict['G'] = np.log(p_exp) + np.log(1-p_plasmid)
    init_dict['P'] = np.log(p_exp) + np.log(p_plasmid)

    # Fill in transition expansion probs
    # Note that spacer addition ends at the first 'N' event.
    trns_dict['N'] = np.log(1 - s_prob*p_exp)
    trns_dict['G'] = np.log(s_prob) + np.log(p_exp) + np.log(1-p_plasmid)
    trns_dict['P'] = np.log(s_prob) + np.log(p_exp) + np.log(p_plasmid)

    for patt in spacer_patts: # iter over patterns
        for i, sp in enumerate(patt): # iter over spacers
            if i == 0:
                patt_dict[patt] += init_dict[sp]
            else:
                patt_dict[patt] += trns_dict[sp]
            if sp == 'N':
                break

    if log_p: # already in log p, just use identity func
        def use_func(x):
            return x
    else: # since in log p, need to exponentiate
        use_func = np.exp

    out_dict = {} # Need to remove trailing 'N's from keys
    for patt in spacer_patts:
        if patt == 'N':
            out_dict[patt] = use_func(patt_dict[patt])
        else:
            out_dict[patt.replace('N','')] = use_func(patt_dict[patt])

    return out_dict
# Return list of all possible combinations of n integers <= max_incorp 
# that sum to k
#     E.g.:
#     binom_combos(3, 2, max_incorp=1)
#     [(0, 1, 1), (1, 0, 1), (1, 1, 0)]
def binom_combos(n, k, max_incorp=1):
    x = list(itertools.product(range(max_incorp+1), repeat=n))
    out_list = []
    for i in x:
        if sum(i) == k:
            out_list.append(i)
    return out_list

# Function for calculating the probability distributions of arrays given a state sequence.
#
# Returns a sample-keyed dictionary of array-length-keyed dictionaries of array-type probabilities as pd.Series objects
def calc_state_array_probs(state_str,
                           arr_len=None,       # calculate probabilities of all crispr arrays of given len
                           o_exp=None,         # Prob of expansion during ON state
                           x_exp=None,         # Prob of expansion during OFF state
                           o_plasmid=None,     # Prob of incorporating pTrig-derived spacer during ON state;
                                               # can be list of probabilities, which is interpreted as probability for incorporation
                                               # on successive days
                           x_plasmid=None,     # Prob of incorporating pTrig-derived spacer during OFF state
                           max_incorp=1,       # Maximum number of incorporations per day
                           scaling_prob=1,  # Prob of incoprorating each additional spacer in a single day
                           char_pair=('X','O')
                          ):
    # Process array inputs
    st_len = len(state_str)

    # Set defaults
    if arr_len is None:
        arr_len = min(st_len, 4) # Default to number of states or 4, whichever is smaller
        print("arr_len not given, defaulting to "+str(arr_len))
    elif arr_len > st_len*max_incorp: # Except if can't achieve array len with given # days, max # incorp/day
        raise ValueError('arr_len must be less than or equal to (length of state_str)*(max_incorp)')

    if o_exp is None:
        o_exp = 0.0988
    if x_exp is None:
        x_exp = 0.0356
    if o_plasmid is None:
        o_plasmid_dict = { 1:0.2749, # Default values of o_plasmid value dependent
                           2:0.2457, # on the array length
                           3:0.2202,
                           4:0.1865,
                           5:0.1809 }
        o_plasmid = o_plasmid_dict[arr_len]
    if x_plasmid is None:
        x_plasmid = 0.0007

    if char_pair[0] == char_pair[1]:
        raise ValueError("char_pair must have unique elements")

    # Get the probabilities for all possible single-day incorporation patterns (depending on day state)
    prob_dict = dict.fromkeys(list(set(state_str))) # Dict with one key for each state char (usually X/O)
    for st in prob_dict.keys():
        prob_dict[st] = calc_single_state_probs(st, o_exp, x_exp, o_plasmid, x_plasmid,
                                                scaling_prob, max_incorp, log_p=True)

    arr_list = make_sample_list(arr_len, char_pair=('G','P')) # all arrays of given arr_len
    p_series = pd.Series(float(0), index=arr_list)

    for arr in arr_list: # iterate over all possible crispr arrays of given arr_len
        rev_arr = arr[::-1] # We read arrays in reverse
        incorp_list = binom_combos(st_len, arr_len, max_incorp)
        for incorp in incorp_list: # iterate over possible patterns of incorporation
            p = float(0) # use log since these are small probabilities
            cnt = 0
            for i in range(st_len): # iterate over incoproration pattern (i.e., days)
                n_incorp = incorp[i] # number of incorporated spacers on current day
                st = state_str[i] # current state

                if n_incorp == 0: # Don't incorporate
                    p += prob_dict[st]['N']
                else: # Incorporate n_incorp spacers
                    p += prob_dict[st][rev_arr[cnt:cnt+n_incorp]]
                    cnt += n_incorp
            p_series[arr] += np.exp(p)
    p_series.index = [ replace_char_pair(x, old=('G','P'), new=char_pair) for x in p_series.index.values ]
    return p_series

# Function that applies 'calc_state_array_probs' to every sample in sample_list for every array length in arr_lens
# The output is comparable to that of the function 'array_counts_from_spacers'
def calc_mult_sample_array_probs(sample_list, arr_lens, char_pair=('X','O'), o_exp=None, x_exp=None,
                                 o_plasmid=None, x_plasmid=None, max_incorp=1, scaling_prob=1):
    out_dict = dict.fromkeys(sample_list)
    for samp in sample_list:
        cur_dict = dict.fromkeys(arr_lens)
        for alen in arr_lens:
            if type(o_plasmid) == dict:
                use_o_plasmid = o_plasmid[alen]
            else:
                use_o_plasmid = o_plasmid
            cur_dict[alen] = calc_state_array_probs(samp, arr_len=alen, char_pair=char_pair, o_exp=o_exp, x_exp=x_exp,
                                                    o_plasmid=use_o_plasmid, x_plasmid=x_plasmid,
                                                    max_incorp=max_incorp, scaling_prob=scaling_prob)
        out_dict[samp] = cur_dict
    return out_dict

###################################################################################################
## S.04 ## Functions for calculating, plotting euclidean distances between data and modeled array frequencies.
###################################################################################################

# Calculate euclidean distance matrix comparing vectors of data array frequencies vs. modeled array frequencies
#
# Returns numpy matrix or pandas DataFrame of pairwise distances (depending on type of 'freq_data')
def calc_euc_distance_matrix(freq_data, freq_mod):
    if type(freq_data) == np.ndarray:
        nr = len(freq_data)
        nc = len(freq_mod)
        euc_distance_matrix=np.zeros(shape=(nr,nc))
        for i in range(len(freq_data)):
            for j in range(len(freq_mod)):
                euc_distance_matrix[i,j]=euclidean(freq_data[i],freq_mod[j])
        return euc_distance_matrix
    elif type(freq_data) == pd.DataFrame:
        nr = freq_data.shape[0]
        nc = freq_mod.shape[0]
        euc_distance_df=pd.DataFrame(np.zeros(shape=(nr,nc)),
                                     index = freq_data.index.values,
                                     columns = freq_mod.index.values)
        for i in range(nr):
            for j in range(nc):
                euc_distance_df.iloc[i,j]=euclidean(freq_data.iloc[i,:],freq_mod.iloc[j,:])
        return euc_distance_df

# Plot heatmap of euclidean distances
def plot_heatmap(euc_dist_matr, row_nms=None, col_nms=None, title_text=None, mark_matches=False, figsize=None, cmap='Oranges'):
    if figsize != None:
        plt.figure(figsize=figsize)
    if row_nms is None:
        row_nms = range(euc_dist_matr.shape[0])
    if col_nms is None:
        col_nms = range(euc_dist_matr.shape[1])
    if type(euc_dist_matr) == np.ndarray:
        to_plot=pd.DataFrame(euc_dist_matr,columns=col_nms,index=row_nms)
    elif type(euc_dist_matr) == pd.DataFrame:
        to_plot = euc_dist_matr
        to_plot.columns = col_nms
        to_plot.index = row_nms
    else:
        raise TypeError("'euc_dist_matr' must be 2d numpy array or pandas DataFrame")

    if mark_matches: # mark the squares that match the rows (i.e., correctly modeled predictions)
        rn = to_plot.index.values
        cn = to_plot.columns.values
        annot_df = pd.DataFrame(columns=col_nms,index=row_nms)
        for r in rn:
            for c in cn:
                if r == c:
                    annot_df.loc[r,c] = 'X'
                else:
                    annot_df.loc[r,c] = ''
    else:
        annot_df = False
    sns.heatmap(to_plot, square=True, annot=annot_df, fmt='',cmap=cmap)
    plt.xlabel("model")
    plt.ylabel("data")
    if title_text is not None:
        plt.title(title_text)
    plt.show()

# calculate indicator matrix of the minimum euclidean distance 
def calc_min_indicator_matrix(euc_dist_matr):
    if type(euc_dist_matr) == pd.DataFrame:
        euc_dist_matr = euc_dist_matr.as_matrix()
    dim = euc_dist_matr.shape
    euc_distance_min=np.zeros(shape=dim)
    for i in range(dim[0]):
        i_vec=np.zeros(shape=(dim[1],))
        i_vec[np.argmin(euc_dist_matr[i])]+=1
        euc_distance_min[i,:]=i_vec
    return euc_distance_min

# compare array vectors against themselves (correlation plot)
def plot_clustmap(array_vec_df, nms=None, title_text=None):
    to_plot=array_vec_df.transpose()
    if nms is not None:
        to_plot.columns = nms
    clust_plot = sns.clustermap(to_plot.corr(),square=True)
    plt.setp(clust_plot.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
    if title_text is not None:
        plt.title(title_text)
    plt.show()

# Plot an array heatmap for all samples
# This function can be run in either 'from_file' or 'from_raw_vector'
#
# 'from_file' directly accesses the files and counts up spacers without regard to whether or not 
# the complete array was captured; this mode requires providing 'samples' and 'folder' arguments at minimum
#
# 'from_raw_vector' counts up spacers from a provided dictionary of raw array vectors (generated
# from the compute_array_vectors() function); this mode ensures that only complete arrays are counted;
# this mode requires providing an 'array_vectors_raw' argument at minimum
def plot_spacer_heatmaps(array_vectors_raw=None, samples=None, folder=None, char_pair=('X','O'),
                          figsize=(12,12), subplot_width=8, method='from_file'):
    def from_file(samp, array_vectors_raw, char_pair, folder):
        #opening uniq_array and uniq_id for code below
        uniq_array=[]
        uniq_id=[]
        uniq_start=[]
        uniq_end=[]
        for i in open(os.path.join(folder, samp+'_uniq.txt')):
            uniq_array.append(i.split(',')[0])
            uniq_id.append(i.split(',')[1])
            uniq_start.append(int(i.split(',')[2]))
            uniq_end.append(int(i.split(',')[3]))

        # svec=['bl21de3','pcas12laci', 'psb2k3min']
        svec=['bl21de3','pcas12','psb2k3min']
        avec=['l1','l2','l3','l4','l5']
        pvec=['p1','p2','p3','p4','p5']
        parray=np.zeros(shape=(len(avec),len(pvec)))
        for a in range(len(avec)):
            for p in range(a+1):
                ids = [(avec[a] in i and pvec[p] in i) for i in uniq_array]
                counts = list(compress(uniq_id, ids))
                sumvec = [sum([i==source for i in counts]) for source in svec]
                try:
                    sb2k3 = sumvec[2]/float(sum(sumvec))
                except ZeroDivisionError:
                    sb2k3 = float(0)
                parray[a,p] = sb2k3
        return pd.DataFrame(parray, index=[x.replace('l','L') for x in avec], columns=pvec)
    
    def from_raw_vector(samp, array_vectors_raw, char_pair, folder):
        anum = array_vectors_raw[samp].keys(); anum.sort()
        avec = [ 'L' + str(x) for x in anum ]
        pvec = [ 'p' + str(x) for x in anum ]
        parray=np.zeros(shape=(len(avec),len(pvec)))
        for a in range(len(anum)):
            for p in range(a+1):
                is_o = [ x[p] == char_pair[1] for x in array_vectors_raw[samp][anum[a]].index.values ]
                cur_cnt = sum( array_vectors_raw[samp][anum[a]][is_o] ) # number of sb2k3 spacers seen at current position in arrays of current length
                cur_sum = sum( array_vectors_raw[samp][anum[a]] ) # sum of number of spacers seen at current position in arrays of current length

                try:
                    sb2k3 = cur_cnt/float(cur_sum)
                except ZeroDivisionError:
                    sb2k3 = float(0)
                parray[a,p] = sb2k3
        return pd.DataFrame(parray, index=avec, columns=pvec)

    plt.figure(figsize=figsize)
    if method.lower() == 'from_file':
        if (samples is None) or (folder is None):
            raise TypeError("'samples' and 'folder' args must be specified if using method 'from_file'")
        use_func = from_file
            
    elif method.lower() == 'from_raw_vector':
        if array_vectors_raw is None:
            raise TypeError("'array_vectors_raw' arg must be specified if using method 'from_raw_vector'")            
        if samples == None:
            samples = array_vectors_raw.keys()
        use_func = from_raw_vector
    
    cnt = 1
    subplot_height = np.ceil( len(samples)/float(subplot_width) )
    for samp in samples:
        to_plot = use_func(samp, array_vectors_raw, char_pair, folder)
        mask = np.zeros(to_plot.shape, dtype=bool) # mask cells above the diagonal
        mask[np.triu_indices_from(mask,k=1)] = True
        plt.subplot(subplot_height, subplot_width, cnt)
        sns.heatmap(to_plot,square=True,cmap='Reds',vmax=0.3,cbar=False, mask=mask)
        plt.title(samp)
        plt.xlabel('position')
        plt.ylabel('length')
        cnt+=1
        
    plt.tight_layout()
    plt.show()
plot_len_pos_heatmaps = plot_spacer_heatmaps # Old plot_spacer_heatmaps name, for compatibility