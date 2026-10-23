#!/usr/bin/env python3

# code by Annelies de Rooij
### adapted from Alvin X. Han (https://github.com/AMC-LAEB/Within_Host_H3vH1/blob/main/H3N2/scripts/reconstruct_haplotype_wrapper.py)

import os, bisect, sys, argparse, gzip, json, re, itertools
import pandas as pd, numpy as np, statistics as st, scipy as sp
from scipy.stats import binom
from scipy.optimize import minimize
from collections import OrderedDict
import multiprocessing as mp
from tqdm import tqdm
from Bio import SeqIO
from scipy.special import loggamma
import matplotlib.pyplot as plt
import multiprocessing as mp
from tqdm import tqdm


cwd = os.getcwd() #working dir
filedir = os.path.abspath(os.path.dirname(__file__)) #directory of this file

with open(os.path.join(filedir, "..", "data", "CDS_H3N2.json")) as json_file:
    coding_json = json.load(json_file) 

ref_seqs = {r.id:r.seq for r in SeqIO.parse(os.path.join(filedir, "..", "data", "references", "stock_consensus_cds.fasta"), "fasta")}


def map_reads_per_segment(sample_reads, segment, coding_json, bqt=30):
    mapped_reads = {}
    for read_id, reads in sample_reads.items():
        #both reads must be present
        if set(reads.keys()) != set(['R1', 'R2']):
            continue

        if segment != reads["R1"]["refname"]:
            continue

        #get info per read pair id 
        r1_positions = np.arange(reads["R1"]["start_pos"], reads["R1"]["start_pos"] + len(reads["R1"]["read_seq"]))
        r1_sequence = np.frombuffer(reads["R1"]["read_seq"].encode(), dtype="S1").copy()
        r1_bq = np.fromiter((ord(c) - 33 for c in reads["R1"]["base_qual"]), dtype=np.uint8)

        r2_positions = np.arange(reads["R2"]["start_pos"], reads["R2"]["start_pos"] + len(reads["R2"]["read_seq"]))
        r2_sequence = np.frombuffer(reads["R2"]["read_seq"].encode(), dtype="S1").copy()
        r2_bq = np.fromiter((ord(c) - 33 for c in reads["R2"]["base_qual"]), dtype=np.uint8)

        #determine overlap
        overlap, r1_idx, r2_idx = np.intersect1d(r1_positions, r2_positions, return_indices=True)

        if overlap.size > 0:

            r1_overlap_read = r1_sequence[r1_idx]
            r2_overlap_read = r2_sequence[r2_idx]
            r1_overlap_bq = r1_bq[r1_idx]
            r2_overlap_bq = r2_bq[r2_idx]

            #get positions where bases are different
            mismatch_mask = r1_overlap_read != r2_overlap_read
            mismatch_indices = np.where(mismatch_mask)[0]
        
            for i in mismatch_indices:
                r1_base = r1_overlap_read[i].decode()
                r2_base = r2_overlap_read[i].decode()

                if r1_base == "N" and re.search("[ATGC]", r2_base):
                    r1_sequence[r1_idx[i]] = r2_sequence[r2_idx[i]]

                elif re.search("[ATGC]", r1_base) and re.search("[ATGC]", r2_base):
                    if r1_overlap_bq[i] < r2_overlap_bq[i]: #better base qualitiy of read 2
                        r1_sequence[r1_idx[i]] = r2_sequence[r2_idx[i]]
                        r1_bq[r1_idx[i]] = r2_bq[r2_idx[i]]
                    elif r1_overlap_bq[i] == r2_overlap_bq[i] and r1_base != r2_base: #ambiguous
                        r1_sequence[r1_idx[i]] = b"N"
                elif r1_base == "-" and r2_base != "-":
                    r1_sequence[r1_idx[i]] = r2_sequence[r2_idx[i]]
                    r1_bq[r1_idx[i]] = r2_bq[r2_idx[i]]

            # remove overlapped section in r2
            r2_keep = np.ones(len(r2_sequence), dtype=bool)
            r2_keep[r2_idx] = False
            r2_positions = r2_positions[r2_keep]
            r2_sequence = r2_sequence[r2_keep]
            r2_bq = r2_bq[r2_keep]

        all_positions = np.concatenate([r1_positions, r2_positions])
        all_sequence = np.concatenate([r1_sequence, r2_sequence])
        all_bq = np.concatenate([r1_bq, r2_bq])
        
        all_sequence[all_bq < bqt] = b"N"
        
        positions = {pos + 1: base.decode() for pos, base in zip(all_positions, all_sequence)}
        #positions = dict(OrderedDict(sorted(positions.items())))
        
        mapped_reads[read_id] = positions

    mapped_reads = pd.DataFrame.from_dict(mapped_reads, orient="index")
    cds = list(set([i for p in coding_json[segment].keys() for i in coding_json[segment][p]["nuc_cds_pos"]]))
    mapped_reads = mapped_reads[cds]#.drop_duplicates()
    return mapped_reads


def get_partial_haplotype_df(segment_snvs, mapped_reads, ref_seq):
    #determine possible allele
    possible_alleles = {}
    for nucpos in segment_snvs.index.get_level_values(0).unique():
        possible_alleles[nucpos] = []
        #determine ref base 
        refbase = ref_seq[nucpos-1]
        #determine if refbase is accept as variant
        possible_alleles[nucpos].append(refbase)
        
        possible_alleles[nucpos].extend(segment_snvs.loc[nucpos].index.tolist())


    #tally observed haplotypes from reads
    partial_haplotype_tally = {}
    for read_id, row in mapped_reads.iterrows():
        
        #determine if read covers one or more variant positions and if read is one of the possible alleles
        positions = [nucpos for nucpos, allele in row.dropna().to_dict().items() if nucpos in possible_alleles.keys() and allele in possible_alleles[nucpos]]
        if len(positions) == 0:
            continue
        alleles = [row[nucpos] for nucpos in positions]
        

        partial_haplotype_tally.setdefault(tuple(positions), {}).setdefault(tuple(alleles), 0)
        partial_haplotype_tally[tuple(positions)][tuple(alleles)] += 1

    partial_haplotypes = []
    for positions in sorted(partial_haplotype_tally.keys()):
        for alleles, allele_count in partial_haplotype_tally[positions].items():
            partial_haplotypes.append({"positions":positions, "alleles":alleles, "count": allele_count})

    partial_haplotypes = pd.DataFrame.from_dict(partial_haplotypes)
    return partial_haplotypes

def generate_possible_haplotypes(partial_haplotypes, nuc_coverage):
    #compute variant coverage
    nucpos_allele_count = {}
    for i, row in partial_haplotypes.iterrows():
        positions = list(row["positions"])
        alleles = list(row["alleles"])

        for nucpos, allele in zip(positions, alleles):
            try:
                nucpos_allele_count[nucpos][allele] += row["count"]
            except:
                try:
                    nucpos_allele_count[nucpos][allele] = row["count"]
                except:
                    nucpos_allele_count[nucpos] = {allele:row["count"]}
        
    # genrate all possible haplotypes
    sorted_nucpos = sorted(set(nucpos_allele_count.keys()))
    #all_possible_haplotypes = list(itertools.product(*[list(nucpos_allele_count[nucpos].keys()) for nucpos in sorted_nucpos]))
    all_possible_haplotypes = itertools.product(*[list(nucpos_allele_count[nucpos].keys()) for nucpos in sorted_nucpos])


    #rank possible haplotype by expected joint probability computed from single-site frequencies
    ssjp_scores = {}
    for haplotype in all_possible_haplotypes:
        allele_freqs = []
        for i, nucpos in enumerate(sorted_nucpos):
            allele = haplotype[i]
            #get frequency 
            freq = nuc_coverage.loc[nucpos, allele]/nuc_coverage.loc[nucpos, "coverage"]
            allele_freqs.append(freq)

        ssjp_scores[haplotype] = np.product(allele_freqs)

    #rank by single-site joint proba
    haplotypes_ssjp = [k for k, v in sorted(ssjp_scores.items(), key=lambda item: item[1], reverse=True)]

    # rank haplotypes by each unique read_type
    read_types = []
    multi_mutants = False
    for positions in partial_haplotypes["positions"].unique():
        #skipping if total count is <  10
        if partial_haplotypes[partial_haplotypes["positions"]==positions]['count'].sum() < 10:
            continue
        read_types.append(positions)

        #determine if there are multimutants in reads
        if len(positions) > 1:
            for alleles in partial_haplotypes[partial_haplotypes["positions"]==positions]["alleles"].unique():
                if sum([1 if nucpos_allele_count[positions[i]][allele]/sum(nucpos_allele_count[positions[i]].values()) < 0.5 else 0 for i, allele in enumerate(alleles)]) > 1:
                    multi_mutants = True
    
    if not multi_mutants:
        #if only single mutants are found > rank by single site joint prob
        sorted_haplotypes = haplotypes_ssjp[:]
    else:
        #compute multi site joint probability 
        read_type_hierarchy = {}
        #go over read types by number of sites covered > greater first
        for positions in sorted(read_types, key=lambda nucpos_list: len(nucpos_list))[::-1]:
            added = False
            
            for prev_pos in read_type_hierarchy.keys():
                #add any read type with positions in a larger top read
                if set(positions) <= set(prev_pos):
                    read_type_hierarchy[prev_pos].append(positions)
                    added = True
            
            #if read type is unique it become the top read 
            if not added:
                read_type_hierarchy[positions] = []

        #determine the top read type that dont overlap with anyother top read types
        no_overlap_topreads = [trt1 for trt1 in read_type_hierarchy.keys() if all(trt1 == trt2 or not set(trt1) & set(trt2) for trt2 in read_type_hierarchy.keys())]
        
        haplotypes_to_remove = []
        msjp_scores = {}
        for top_positions, lower_positions_list in read_type_hierarchy.items():
            top_partial_haplotypes = partial_haplotypes[partial_haplotypes["positions"]==top_positions]
            top_alleles_count = {"".join(row["alleles"]):row["count"] for i, row in top_partial_haplotypes.iterrows()}

            #absorb lower hierarchy alleles into top hierarchy alleles
            for lower_positions in sorted(lower_positions_list, key=lambda positions: partial_haplotypes[partial_haplotypes['positions']==positions]['count'].sum(), reverse=True):
                lower_partial_haplotypes = partial_haplotypes[partial_haplotypes["positions"]==lower_positions]

                for lower_alleles in lower_partial_haplotypes["alleles"]:
                    lower_alleles_regex = ''.join([lower_alleles[lower_positions.index(top_nucpos)] if top_nucpos in lower_positions else '[A-Z]' for top_nucpos in top_positions])

                    #determine the top alleles the lower alleles could be absorded in
                    matched = []
                    for top_alleles in top_alleles_count:
                        if re.search(lower_alleles_regex, top_alleles):
                            matched.append(top_alleles)

                    lower_alleles_count = lower_partial_haplotypes[lower_partial_haplotypes['alleles']==lower_alleles]['count'].iloc[0]
                    
                    if len(matched) == 0:
                        top_alleles_count[lower_alleles_regex] = lower_alleles_count
                    elif len(matched) == 1:
                        top_alleles_count[matched[0]] += lower_alleles_count
                    else:
                        matched_total = sum([top_alleles_count[top_alleles] for top_alleles in matched])
                        for top_alleles in matched:
                            top_alleles_count[top_alleles] += lower_alleles_count*(top_alleles_count[top_alleles]/matched_total)

            #compute joint probability
            all_possible_haplotypes = itertools.product(*[list(nucpos_allele_count[nucpos].keys()) for nucpos in sorted_nucpos])
            for haplotype in all_possible_haplotypes:
                alleles = "".join([haplotype[sorted_nucpos.index(top_nucpos)] for top_nucpos in top_positions])
                score = 0
                for top_alleles, count in top_alleles_count.items():
                    if re.search(top_alleles, alleles):
                        score = count/(sum(top_alleles_count.values()))
                        break
                
                if score == 0 and top_positions in no_overlap_topreads:
                    haplotypes_to_remove.append(haplotype)
                    continue

                try:
                    msjp_scores[haplotype].append(score)
                except:
                    msjp_scores[haplotype] = [score]
                    
        msjp_scores = {haplotype:np.prod(score) for haplotype, score in msjp_scores.items()}
        haplotypes_msjp = sorted(msjp_scores, key=msjp_scores.get)[::-1]
        haplotypes_msjp += haplotypes_to_remove

        #ranked by combined single-site and multi-site joint prob
        all_possible_haplotypes = itertools.product(*[list(nucpos_allele_count[nucpos].keys()) for nucpos in sorted_nucpos])

        sorted_haplotypes = {haplotype:np.prod([haplotypes_msjp.index(haplotype), haplotypes_ssjp.index(haplotype)]) for haplotype in all_possible_haplotypes if haplotype not in haplotypes_to_remove}
        sorted_haplotypes = sorted(sorted_haplotypes, key=sorted_haplotypes.get)

    #get single site allele counts as df
    ss_allele_haplotypes = [[nucpos, allele, count]for nucpos, ac in nucpos_allele_count.items() for allele, count in ac.items()]
    ss_allele_haplotypes = pd.DataFrame.from_records(ss_allele_haplotypes, columns=["nucpos", "allele", "count"])

    #generate projection matrix
    T = []
    for haplotype in sorted_haplotypes:
        T_row = []
        for i, row in ss_allele_haplotypes.iterrows():
            if haplotype[sorted_nucpos.index(row["nucpos"])] != row["allele"]:
                T_row.append(0)
            else:
                T_row.append(1)
        T.append(T_row)
        
    T = np.array(T)

    return sorted_haplotypes, T, ss_allele_haplotypes

def initialise_possible_haplotype_frequencies(size, unknown):
    #intialise full haplotype frequencies
    haplo_freq = np.random.uniform(size=size)
    if unknown:
        # unkwown haplotype frequency (<1%)
        unknown_freq = np.random.uniform()/100
        haplo_freq = (haplo_freq/np.sum(haplo_freq))*(1-unknown_freq) #normalise
        haplo_freq = np.append(haplo_freq, 1.-np.sum(haplo_freq))
    else:
        haplo_freq = haplo_freq/np.sum(haplo_freq) #normalise

    return haplo_freq

def possible_haplo_freq_to_partial_haplo_freq(T, haplo_freq):
    partial_PI = np.dot(T.T, haplo_freq)
    partial_PI = partial_PI/partial_PI.sum() # normalise
    return partial_PI

def compute_LL(observations, partial_PI, C):
    ALPHA = partial_PI * C 
    ll = 0
    for nucpos in observations["nucpos"].unique():
        nucpos_observations = observations[observations["nucpos"]==nucpos]
        X = np.array(nucpos_observations["count"])
        alpha = ALPHA[nucpos_observations.index]
        ll += loggamma(alpha.sum()) - loggamma(X.sum()+alpha.sum()) + (loggamma(X+alpha)-loggamma(alpha)).sum()
    return ll

def sa_proposal(haplo_freq, sigma, unknown):
    if unknown:
        #propose new frequency for unknown haplotype group
        new_X_freq = haplo_freq[-1] * np.random.normal(0, 0.005)
        new_X_freq = 1e-10 if new_X_freq < 0. else new_X_freq #always positive
        new_X_freq = np.random.uniform()/100 if new_X_freq > 0.01 else new_X_freq #unknown haplot freq (<1%)

        # gaussian kernel on previous haplo frequencies
        new_freq = haplo_freq[:-1] + np.random.normal(0, sigma, len(haplo_freq[:-1]))
        new_freq[new_freq<0.] = 1e-10 # always positive
        new_freq = (new_freq/np.sum(new_freq))*(1-new_X_freq) # normalise

        # combine with new X freq
        new_freq = np.append(new_freq, new_X_freq)
    else:
        # gaussian kernel on previous haplo frequencies
        new_freq = haplo_freq + np.random.normal(0, sigma, len(haplo_freq))
        new_freq[new_freq<0.] = 1e-10 # always positive
        new_freq = new_freq/np.sum(new_freq) # normalise

    return new_freq

def sa_acceptance(new_freq, haplo_freq, LL, observations, T, C):
    #update partial PI according the new frequencies
    partial_PI = possible_haplo_freq_to_partial_haplo_freq(T, new_freq)
    #compute new LL according to updated partial PI
    new_LL = compute_LL(observations, partial_PI, C)
    
    #compute acceptance probability
    if new_LL-LL > np.log(np.random.uniform()):
        return new_freq, new_LL
    else:#no change
        return haplo_freq, LL

def optimize_haploset_frequencies(observations, haploset, unknown, T, C, iterations, searches, save_state, sigma ):
    state_df = []
    for search, seed in enumerate(np.random.choice(range(10000), searches)): 
        np.random.seed(seed) #set seed
        
        #intialize haploype frequencies > random
        haplo_freq =initialise_possible_haplotype_frequencies(len(haploset), unknown) #qi
        #convert full haplotype to partial haplotype frequencies (partial pi; qpi)
        partial_PI = possible_haplo_freq_to_partial_haplo_freq(T, haplo_freq) #qpi  = T * qi

        #compute initial LL
        LL = compute_LL(observations, partial_PI, C)
        state_df.append({"search":search+1, "state":0, "freq":haplo_freq, "LL":LL})
    
        #simulated annealing to find maximum of LL and thus optimal qi
        for j in range(iterations):
            proposed_freq = sa_proposal(haplo_freq, sigma, unknown)
            haplo_freq, LL = sa_acceptance(proposed_freq, haplo_freq, LL, observations, T, C)
            
            if j == 0 or (j+1)%save_state == 0:
                state_df.append({"search":search+1, "state":j+1, "freq":haplo_freq, "LL":LL})

    state_df = pd.DataFrame.from_dict(state_df).set_index(["search", "state"])

    return state_df

def compute_BIC(d, N, LL):
    return d*np.log(N) -2*LL

def state_plot(states_df, prefix):
    fig, ax = plt.subplots(1, 1, figsize=(8, 3.5))
    for h,s,x in states_df.set_index(["haplo length", "search", "X"]).index.unique():
        df = states_df.set_index(["haplo length", "search", "X"]).loc[(h,s,x)]
        
        ax.plot(df["state"], df["LL"], label=f"S{s}-H{h}({x})", lw=1)
    #cosmetics
    ax.set_xlabel("States")
    ax.set_ylabel("LL")
    plt.legend(loc='center left', bbox_to_anchor=(1., 0.5))
    plt.tight_layout()

    plt.savefig(prefix+".pdf", bbox_inches='tight', pad_inches=0.5)
    plt.close(fig)

def reconstruct_haplotype(partial_haplotypes, nuc_coverage,  max_iter=5000, save_state=25, n_searches=2, overdispersion=1e-3, sigma=0.05, verbose=False):
    
    #compute C factor from overdispersion
    C = (1-overdispersion)/overdispersion 
    #generate all possible hypothetical haplotypes and the
    haplotypes, T_matrix, ss_allele_haplotypes = generate_possible_haplotypes(partial_haplotypes, nuc_coverage)

    observations = ss_allele_haplotypes[["nucpos", "count"]] #counts of nucleotides at each Allele
    deltabic_abv_zero = 0
    all_states = []
    for i in range(len(haplotypes)):
        #get currect haplotype set and T matrix 
        haploset = haplotypes[:i+1]
        T = T_matrix[:i+1]

        #add additonal haplotype "X" describing the cloud of potential haplotype not described in the current haploset
        X_haplotype = False
        if not ((1-T.sum(axis=0))<=0).all():
            X_haplotype = True
            T = np.append(T,[1-T.sum(axis=0)], axis=0 )
            T[T<0] = 0

        #create state df > seperate function
        searches = 1 if i==0 else n_searches
        iterations = int(np.ceil(max_iter/2)) if i==0 else max_iter

        #optimise haploset frequencies
        state_df = optimize_haploset_frequencies(observations, haploset, X_haplotype, T, C, iterations, searches, save_state, sigma)
        state_df = state_df.reset_index()
        state_df[["haplo length", "X", "haplo set"]] = [i+1,X_haplotype,str(haploset)]
        

        # compute BIC based on best LL
        best_state = state_df.loc[state_df['LL'].idxmax()]
        best_LL = best_state["LL"]
        n = i+2 if X_haplotype else i+1
        bic = compute_BIC(n, partial_haplotypes["count"].sum(), best_LL)
        state_df["best BIC"] = bic
        if i > 0:
            delta_bic = bic-prev_bic
            if verbose:
                print ("BIC: %.3f (delta=%.3f)"%(bic, delta_bic))
            # Break if delta BIC is positive for two consecutive searches
            if delta_bic > 0:
                deltabic_abv_zero += 1
                if deltabic_abv_zero == 2:
                    break
            else:
                deltabic_abv_zero = 0
        else:
            if verbose:
                print ("BIC: %.3f"%(bic))

        prev_bic = bic
        # append to all_states
        all_states.append(state_df)

    # concat all state_df
    all_states = pd.concat(all_states, ignore_index=True)
    return all_states

def generate_ML_haplotype_freq_dist(ML_freq, ML_haploset, partial_haplotypes, nuc_coverage, unknown, N=100, overdispersion=1e-3, max_iter=5000, save_state=25, n_searches=2, sigma=0.05, verbose=False):
    #compute alpha from computed ML
    alpha = ML_freq*((1/overdispersion)-1)

    #generate N random full haplotype frequencies
    if unknown:
        pi_N = np.random.dirichlet(alpha[:-1], N)
    else:
        pi_N = np.random.dirichlet(alpha, N)

    all_positions = list(set(partial_haplotypes["positions"]))
    sorted_nucpos = sorted(set([k for v in all_positions for k in v]))

    #intialize ML frequency distributions 
    ML_haplotype_freqdist = np.zeros((N, len(ML_haploset)))

    for i, pi in enumerate(pi_N): 
        #recreate partial haplotype df
        sample_partial_haplotypes = []
        for positions in all_positions:

            #generate partial haploset
            partial_haploset = [tuple([haplotype[sorted_nucpos.index(nucpos)] for nucpos in positions]) for haplotype in ML_haploset]
            
            #multinomial sampling based on pi
            n = partial_haplotypes[partial_haplotypes["positions"]==positions]['count'].sum()
            X = np.random.multinomial(n, pi)

            #count partial haplotype
            phaplotype_count = {}
            for j, count in enumerate(X):
                if count != 0:
                    phaplotype_count[partial_haploset[j]] = phaplotype_count.get(partial_haploset[j], 0) + count

            for phap, count in phaplotype_count.items():
                sample_partial_haplotypes.append({'nucpos':positions, 'var': phap, 'count':count})
        
        sample_partial_haplotypes = pd.DataFrame.from_dict(sample_partial_haplotypes)
        all_states = reconstruct_haplotype(partial_haplotypes, nuc_coverage,  max_iter, save_state, n_searches, overdispersion, sigma, verbose)
        best_BIC = min(all_states['best BIC'])
        best_all_states = all_states[all_states['best BIC']==best_BIC]
        best_state = best_all_states.loc[best_all_states['LL'].idxmax()]
        sample_freq = best_state["freq"]
        sample_haploset = eval(best_state['haplo set'])
        

        # append sample_freq to distribution
        for j, haplotype in enumerate(ML_haploset):
            if haplotype in sample_haploset:
                ML_haplotype_freqdist[(i, j)] = sample_freq[sample_haploset.index(haplotype)]
        
    return ML_haplotype_freqdist

def run_reconstruction(segment, snvs, coverage, sample_reads, ref_seq, bqt, max_iter, save_state, n_searches, overdispersion, sigma, N, verbose):
    
    #map reads for segment 
    
    mapped_reads = map_reads_per_segment(sample_reads, segment, coding_json, bqt) 

    segment_snvs = snvs[snvs['refname']==segment].set_index(['nucpos', 'read_base'])
    segment_coverage = coverage[coverage["refname"]==segment].set_index(['nucpos'])

    #get partial haplotypes
    partial_haplotypes = get_partial_haplotype_df(segment_snvs, mapped_reads, ref_seq)

    all_states = reconstruct_haplotype(partial_haplotypes,segment_coverage, max_iter=max_iter, save_state=save_state, 
                                           n_searches=n_searches, overdispersion=overdispersion, sigma=sigma, verbose=verbose)
    best_BIC = min(all_states['best BIC'])
    best_all_states_df = all_states[all_states['best BIC']==best_BIC]
    best_state = best_all_states_df.loc[best_all_states_df['LL'].idxmax()]
    best_freq = best_state["freq"]
    best_haploset = eval(best_state['haplo set'])
    best_X_haplotype = best_state['X']
    best_LL = best_state['LL']

    ML_haplotype_freqdist = generate_ML_haplotype_freq_dist(best_freq, best_haploset, partial_haplotypes, segment_coverage, best_X_haplotype, N=N, max_iter=max_iter, 
                                                            save_state=save_state, n_searches=n_searches, sigma=sigma, overdispersion=overdispersion,
                                                            verbose=verbose)
    
    return segment, all_states, ML_haplotype_freqdist
    

def ArgumentParser():
    """Argument parser"""

    parser = argparse.ArgumentParser(prog = "haplotype reconstruction", #changing this later
        formatter_class = argparse.RawTextHelpFormatter,
        description = "Haplotype reconstruction of the LTEE samples" ) 

    parser.add_argument('-s','--snvs', required=True, action="store",type=str, help="csv file with the identified SNV the sample of interest")
    parser.add_argument('-c','--coverage', required=True, action="store",type=str, help="csv file with nucleotide coverage of the sample")
    parser.add_argument('-r','--reads', required=True, action="store", type=str, nargs="+", help="read file(s) of sample with the sorted and aligned read (must be .json.gz format > flu-tube output)")
    parser.add_argument('-o','--output', required=False, action="store", type=str, default="./output", help="output directory where results are stored (default: ./output)")
    parser.add_argument('-sn','--sample-name', required=False, action="store", type=str, default="", help="sample name if not specified sample name will be inferred from coverage file name (assuming flu-tube output)")
    parser.add_argument('-bqt','--base-quality-threshold', required=False, action='store', type=int, default=30, help="minimal Phred score for base quality (default: 30)")
    parser.add_argument('-mf', '--min-snv-freq', required=False, action="store", type=float, help="minimal SNV frequency required for haplotype reconstruction")
    parser.add_argument('-mi', '--max-iter', required=False, action='store', type=int, default=5000, help="max number of iterations for haplotype reconstruction (default:5000)")
    parser.add_argument('-ss', '--save-state', required=False, action='store', type=int, default=25, help="iteration interval to save to output (default: 25)")
    parser.add_argument('-ns', '--number-searches', required=False, action='store', type=int, default=2, help="number of searches to be performed for haplotype reconstruction (default: 2)")
    parser.add_argument('--overdispersion', required=False, action='store', type=float, default=1e-3, help="variability/noise (default: 1e-3)")
    parser.add_argument('--sigma', required=False, action='store', type=float, default=0.05, help="sigma (default: 0.05)")
    parser.add_argument('--haplotype-sample', required=False, action='store', type=int, default=100, help="number of random samples for ML freq estimation (default:100)")
    parser.add_argument('-t','--threads', required=False, action="store", type=int, default=1, help="number of threads to use (default: 1)")
    parser.add_argument('--verbose',required=False,action='store_true',help="print stuff to screen")

    
    #check if arguments were entered 
    if len(sys.argv[1:])<1: 
        parser.print_help()
        sys.exit(-1)
    else:
        return parser.parse_args()
    
def main():
    args = ArgumentParser()
    bqt = args.base_quality_threshold
    N = args.haplotype_sample
    max_iter= args.max_iter
    save_state= args.save_state
    n_searches= args.number_searches
    overdispersion= args.overdispersion
    sigma= args.sigma
    threads = args.threads


    #check if depth file and allele count files exist
    snv_file = os.path.join(cwd, args.snvs)
    if not os.path.isfile(snv_file):
        sys.stderr.write(F"Error: cannot find {snv_file}.\n Check if file is properly specified")
        sys.exit(-1)

    coverage_file = os.path.join(cwd, args.coverage)
    if not os.path.isfile(coverage_file):
        sys.stderr.write(F"Error: cannot find {coverage_file}.\n Check if file is properly specified")
        sys.exit(-1)

    read_files = [os.path.join(cwd, r) for r in args.reads]
    for f in read_files:
        if not os.path.isfile(f):
            sys.stderr.write(F"Error: cannot find {f}.\n Check if file is properly specified")
            sys.exit(-1)

    #get output directory 
    output_dir = os.path.join(cwd, args.output)
    if not os.path.exists(output_dir):
        try:
            os.mkdir(output_dir)
        except:
            sys.stderr.write(f"Error: cannot create directory {output_dir}.\nPlease try creating the output directory manually.")
            sys.exit(-1)

    #determine sample name 
    if len(args.sample_name) > 1:
        sample= args.sample_name
    else:
        sample= coverage_file.split("nuc_coverage_")[-1].rstrip(".csv")

    snvs = pd.read_csv(snv_file, keep_default_na=False, na_values=[""])
    if args.min_snv_freq:
        snvs = snvs[snvs["nuc_freq"]>=args.min_snv_freq]



    coverage = pd.read_csv(coverage_file, keep_default_na=False, na_values=[""])

    if len(snvs) == 0:
        sys.stdout.write(f"No SNVs found in SNV file")
        sys.exit(-1)

    sample_reads = {}
    for f in read_files:
        with gzip.open(f) as json_file:
            sample_reads.update(json.load(json_file))

    with mp.Pool(processes=threads) as p:
        for result in tqdm([p.apply_async(run_reconstruction, args=(segment,  snvs, coverage, sample_reads, ref_seqs[segment], bqt, max_iter, 
                                                                    save_state, n_searches, overdispersion, sigma, N, args.verbose)) for segment in snvs["refname"].unique()]):
            
            segment, all_states, ML_haplotype_freqdist = result.get()

            all_states.to_csv(os.path.join(output_dir, f"reconstructed_states_{sample}_{segment}.csv"), index=False)
            np.savetxt(os.path.join(output_dir, f'freq_dist_{sample}_{segment}.csv'), ML_haplotype_freqdist, delimiter=",")
            



if __name__ == "__main__":
    main()
