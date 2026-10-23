#!/usr/bin/env python3

#Code by Annelies J.H. de Rooij (Jan 2026)
#partially adapted from code by Alvin X. Han (https://github.com/AMC-LAEB/Within_Host_H3vH1/blob/main/simulations/WHDEL.py)

import argparse, sys, os
import numpy as np, pandas as pd
from numba import jit
from scipy.integrate import odeint
from scipy.sparse import csr_matrix, vstack
from scipy.stats import poisson, multinomial, binom
import multiprocessing as mp
from tqdm import tqdm
import datetime

#np.random.seed(40038)

def ArgumentParser():
    """Argument parser"""

    parser = argparse.ArgumentParser(prog = "LTEE drift in-silico evolution", 
        formatter_class = argparse.RawTextHelpFormatter,
        description = "simulated the LTEE on binary genomic to generate drift like frequency trajectories " ) 

    parser.add_argument('-o', '--output', required=False, action='store', type=str, default="simulation.csv", help="path to out file where simulation results should be written to (default:./simulation.csv)" )
    parser.add_argument('-n','--n-sims', required=False, action="store",type=int, default=1, help="number of simulations to perform")
    parser.add_argument('-t', '--threads', required=False, action='store', type=int, default=1, help="number of threads")
    parser.add_argument('-v', '--verbose', required=False, action='store_true',help="print stuff to screen")
    
    
    ### GROWTH component argument
    passage_args = parser.add_argument_group("passage/growth arguments", "paramters specific to within passage viral growth")
    #### time-based arguments
    passage_args.add_argument('-p', '--passages', required=False, action='store', type=int, default=10, help="number of passage to simulated (default:10)")
    passage_args.add_argument('--generation-time', required=False, action="store", type=int, default=6, help="generation time in hours (default: 6h)")
    passage_args.add_argument('--hours-per-passage', required=False, action='store', type=int, default=72, help="hours per passage (default:72h)")
    #### virion arguments
    passage_args.add_argument('--input-virions', required=False, action='store', type=int, default=1000, help="number of input virions for each passage (default: 1000)")
    passage_args.add_argument('--successful-infection-frac', required=False, action='store', type=float, default=1., help="fractions of input virions that successfully infect each passage (default: 1.)")
    
    passage_args.add_argument('--max-virions', required=False, action='store', type=int, default=10e8, help="max number of virions present in well (viral growth upper limit) (default:10e8)")
    #### target-cell argument
    passage_args.add_argument('--target-cells', required=False, action='store', type=int, default=2e5, help="number of (target-)cells in HAE well (default:2e5)")
    #### model arguments
    passage_args.add_argument('--per-capita-decay-rate', required=False, action='store', type=int, default=2, help="per capita decay rate, d (default: 2 per generation)")
    passage_args.add_argument("--within-host-basic-reproduction-number", required=False, action='store', type=int, default=10,
                              help="within host basic reproduction number, R0 (default: 10)")
    passage_args.add_argument('--virion-production-number', required=False, action='store', type=int, default=100, 
                              help="number of virions produced by a cell before it dies/becomes unproductive (default:100)")
    
    ### GENOME arguments
    genome_args = parser.add_argument_group("genome arguments", "paramters specific to genome that will evolve under drift")
    genome_args.add_argument('--genome-length', required=False, action='store', type=int, default=200, help="length of the binary genome (default:200)")
    genome_args.add_argument('-mu', '--mutation-rate', required=False, action='store', type=float, default=3e-6, help="mutation rate (default: 3e-6)")
    
    ### SAMPLING components arguments
    sampling_args = parser.add_argument_group("sampling arguments", "parameters specific for sampling from one passage to the next")
    sampling_args.add_argument('--apical-release-fraction', required=False, action='store', type=float, default=0.5, 
                               help="fraction of total virus population that is release apically and can thus be sampled (default:0.5)")
    sampling_args.add_argument('--n-apical-washes', required=False, action='store', type=int, default=2, help="number of apical washes (default:2)")
    sampling_args.add_argument('--apical-wash-vol', required=False, action='store', type=float, default=0.1, help="volume of each individual apical wash perfomed in mL (default:0.1 mL(100 uL))")
    sampling_args.add_argument('--wash-recovery-frac', required=False, action='store', type=float, default=1, 
                               help="fraction of total virus recovered from the apical using the apical washes (default:1.0)")
    sampling_args.add_argument('--infection-vol', required=False, action='store', type=float, default=0.05, 
                               help="volume taken from the apical wash pool to dilute and infect the next passage in mL (default: 0.05 mL (50 uL))")
    
    
    #### SEQUENCING components arguments
    sequencing_args = parser.add_argument_group("sequencing arguments", "parameters specifc for sampling for sequencing")
    sequencing_args.add_argument('--sequence-every-nth-passage', required=False, action='store', type=int, default=10, help="sample for sequencing every 10th (default) passage")
    sequencing_args.add_argument('--n-sequence-samples', required=False, action='store', type=int, default=2, help="number of technical replicate for sequencing (default:2)")
    sequencing_args.add_argument('--sequence-sample-vol', required=False, action='store', type=float, default=0.025, help="volume of apical wash takes for each technical replicate in mL (default:0.025 mL (25 uL))")
   
    
    #check if arguments were entered 
    if len(sys.argv)<1: 
        parser.print_help()
        sys.exit(-1)
    else:
        return parser.parse_args()


@jit(nopython=True)
def resample_genotype_dist(V, genotype_prop,size=None):
    genotype_dist = np.zeros(len(genotype_prop), dtype=np.int32)
    remaining_V = int(V)
    remaining_p = 1.0
    
    for i in range(len(genotype_prop) - 1):
        if remaining_V <= 0 or remaining_p <= 0:
            break
        # Current probability relative to what is left
        p_current = genotype_prop[i] / remaining_p
        # Clip p to ensure floating point errors don't exceed 1.0
        p_current = min(max(p_current, 0.0), 1.0)
        
        
        # Rewrote binomial robust function within this function
        Np = remaining_V * p_current #expected number of success
        #poisson approxiamation when many trails and rare success > faster and numerically stable
        if remaining_V > 20 and Np < 5:
            count = np.random.poisson(Np, size)
        else:
            count = np.random.binomial(remaining_V, p_current, size)
        
        genotype_dist[i] = count
        
        remaining_V -= count
        remaining_p -= genotype_prop[i]
        
        # If remaining probability is negligible, stop early
        if remaining_p <= 1e-10:
            break
            
    # Sum what remains to the s
    # usually no
    if remaining_V > 0:
        genotype_dist[-1] += remaining_V
        
    return genotype_dist
    
def select_random_N_individuals_by_genotype(N, N_dist):

    #generate initial random array 
    result = np.random.randint(low=0, high=N_dist)
    #ensure initial array isn't all zero 
    b = 0 #breaking while loop after 10 tries
    while result.sum() == 0 and b < 10:
        result = np.random.randint(low=0, high=N_dist)
        b += 1
    
    #need to ajust to N to get result
    if result.sum() > 0:
        result = np.rint(N*(result/result.sum())).astype(np.int32)
       
    #result sum should the size of N 
    diff = result.sum() - N

    while abs(diff) != 0:
        
        #determine adjustment to be made so that results matches the size of N
        if abs(diff) < 5: #if erorr is small just fix it in one step 
            random_adjust = diff
        else:
            random_adjust = np.random.randint(low=diff, high=0) if diff < 0 else np.random.randint(low=0, high=diff)
        
        #chose 100 or 90% of the total results random result entries (by index)
        n_to_select = 100 if len(result) >= 100 else np.int8(0.9*len(result))
        result_idx = np.random.choice(np.arange(len(result)),n_to_select, replace=False)
        np.random.shuffle(result_idx)
        
        #add the adjustments
        for idx in result_idx:
            value = result[idx]
            #ensure adjustment is betweel N_dist and 0
            if value - random_adjust < N_dist[idx] and value - random_adjust >= 0:
                #add random ajdustment and break
                result[idx] = value - random_adjust
                break
        
        diff = result.sum() - N
        
    return result

def remove_extinct_genotypes(genotypes, genotype_dist):
    surviving_genotypes = np.where(genotype_dist>0)[0]
    if len(surviving_genotypes) < len(genotype_dist):
        genotype_dist = genotype_dist[surviving_genotypes]
        genotypes = genotypes[surviving_genotypes]
        
    return genotype_dist, genotypes

@jit(nopython=True) 
def build_mutated_genotypes(genotype_indptr, genotype_indices, mutator_indices, n_mutators, L, mV, total_entries):
    
    #create rows tracker and col tracker for sparse matrix based on the calculated maximum non-zero entry
    rows_out = np.empty(total_entries, dtype=np.int32)
    cols_out = np.empty(total_entries, dtype=np.int32)
    
    idx = 0 #index tracker for row and columns trackers 
    mutant_counter = -1 #counter for when mutant switches in row tracker
    for i, n_mutants in enumerate(n_mutators):
        g= mutator_indices[i]
        
        mutated_loci = genotype_indices[genotype_indptr[g]:genotype_indptr[g+1]]
        
        #for each mutated randomly select a locus to mutate
        for _ in range(n_mutants):
            mutant_counter +=1 #new mutant 
            locus_to_mutate = np.random.randint(0, L)
            for locus in mutated_loci:
                if locus != locus_to_mutate: #if back mutation > 0 so not writting it to the non-zero tracker
                    rows_out[idx] = mutant_counter
                    cols_out[idx] = locus
                    idx += 1
            #add new mutation if no back mutation
            if locus_to_mutate not in mutated_loci:
                rows_out[idx] = mutant_counter
                cols_out[idx] = locus_to_mutate
                idx += 1

    return rows_out[:idx], cols_out[:idx]

def generate_mutated_genotypes(genotype_dist, genotypes, mV, L):
    #randomly select individuals from the genotyps to mutate
    g_mutate = select_random_N_individuals_by_genotype(mV, genotype_dist)
    
    
    #remove mutated individuals from original genotype
    genotype_dist -= g_mutate
    
    #get index of mutators and number of mutators per genotype
    mutator_mask  = g_mutate > 0
    mutator_indices = np.where(mutator_mask)[0]
    n_mutators =  g_mutate[mutator_mask]
    
    #determine the number of mutations per genotype 
    mutations_per_genotype = genotypes.indptr[1:] - genotypes.indptr[:-1]
    #and how many each mutated genotype will have (assuming a forward mutation)
    #summing this to create a overview of all non-zero positions in the spars matrix
    max_mutations_per_genotype = np.sum((mutations_per_genotype[mutator_mask] + 1) * n_mutators)
    
    # build the CSR matrix for the new mutant genotypes
    new_rows, new_cols= build_mutated_genotypes(genotypes.indptr, genotypes.indices, mutator_indices, n_mutators, L, mV, max_mutations_per_genotype)
    mutant_genotypes = csr_matrix( (np.ones(len(new_rows), dtype=np.int8), (new_rows, new_cols)), 
        shape=(mV, genotypes.shape[1])
    )

    #update genotypes and genotype_dist with the mutated genotypes 
    genotypes = vstack([genotypes, mutant_genotypes])
    genotype_dist = np.hstack([genotype_dist, np.ones(mV, dtype=np.int32)])

    return genotype_dist, genotypes

@jit(nopython=True) 
def find_unique_genotypes(indptr, indices, genotype_dist):
    #also collapsing , non-unique genotypes so that they will become a unique genotype
    unique_genotype_mask = np.ones(genotype_dist.shape[0], dtype=np.bool_)
    
    last_unique_index = 0
    for i in np.arange(1, genotype_dist.shape[0]):
        #see if the length of the previous genotype match this ones
        if (indptr[last_unique_index+1] - indptr[last_unique_index]) != (indptr[i+1]-indptr[i]):
            last_unique_index = i
            continue #genotype don't match in length
        
        genotype_match = True
        #check if the actual genotype match HAS TO BE SORTED
        for j in range(indptr[i+1] - indptr[i]):
            if indices[indptr[last_unique_index] + j] != indices[indptr[i] + j]:
                genotype_match = False
                break
            
        if genotype_match:
            #add duplicate genotype counts to original genotype counts in dist
            genotype_dist[last_unique_index] += genotype_dist[i]
            #reset old one 
            genotype_dist[i] = 0
            unique_genotype_mask[i] = False 
        else:
            last_unique_index = i

    return unique_genotype_mask, genotype_dist

def combine_identical_genotypes(genotype_dist, genotypes):
    #sorting genotypes lexicographically > will sort identical genotype right next to one another
    
    #genotypes is a sparse matrix so getting both indices and indptr
    genotypes.sort_indices() #sorting first
    indices = genotypes.indices #columns of non-zero elements in the sparse matrix
    indptr = genotypes.indptr #rows of non-zero elements in the sparse matrix
    
    #for each genotype grab the mutated loci 
    genotype_mloci = [tuple(indices[indptr[g]:indptr[g+1]]) for g in np.arange(genotype_dist.size)]
    
    #sorting genotypes lexicographically
    sort_order = np.array(sorted(range(len(genotype_mloci)), key=lambda k: genotype_mloci[k]))
    #reorder according to sort 
    genotype_dist = genotype_dist[sort_order]
    genotypes = genotypes[sort_order]
   
    #again getting both indices and indptr
    indices = genotypes.indices
    indptr = genotypes.indptr
    
    keep_mask, genotype_dist = find_unique_genotypes(indptr, indices, genotype_dist)
    
    genotype_dist, genotypes = genotype_dist[keep_mask], genotypes[keep_mask] 
    return genotype_dist, genotypes

@jit(nopython=True) 
def within_passage_ODE(x, t, beta, d, r):
    C, V= x
    dC = -(beta*C*V)
    dV = (beta*C*V*r)-(d*V) 
    return dC, dV

def simulate_growth_during_passage(genotype_dist, genotypes, Cmax, Vmax, Vlim, Tpassage_vec, mu, L, beta, d, r):
    #get number of virions over time in the passage
    C_vec, V_vec = odeint(within_passage_ODE, (Cmax, Vmax), Tpassage_vec, args=(beta, d, r)).transpose()
    
    #expand monoclonal population, first replication cycle
    V = V_vec[1]
    genotype_prop = genotype_dist / np.sum(genotype_dist) #current genotype proportions
    genotype_dist = resample_genotype_dist(V, genotype_prop) #resample based on current proportion
    genotype_dist, genotypes = combine_identical_genotypes(genotype_dist, genotypes)
    #need to remove extinct genotypes
    genotype_dist, genotypes = remove_extinct_genotypes(genotypes, genotype_dist)
  
    
    for t in Tpassage_vec[1:-1]:
        #randomly select mV number of individuals to mutate 
        mV = np.random.poisson(V*mu*L)
        #print (mV)
        if mV > 0:
            genotype_dist, genotypes = generate_mutated_genotypes(genotype_dist, genotypes, mV, L)
            #need to combine identical genotypes
            genotype_dist, genotypes = combine_identical_genotypes(genotype_dist, genotypes)

        #resampling for next time point
        V = V_vec[int(t+1)] if V_vec[int(t+1)] <= Vlim else Vlim #upper limit on population size
        genotype_prop = genotype_dist / np.sum(genotype_dist) #current genotype proportions
        genotype_dist = resample_genotype_dist(V, genotype_prop) #resample based on current proportion

        #need to remove extinct genotypes
        genotype_dist, genotypes = remove_extinct_genotypes(genotypes, genotype_dist)
    
    return genotype_dist, genotypes

def get_allele_counts(genotype_dist, genotypes, loci):
    
    #initialize array count
    counts = np.zeros(loci, dtype=np.int32)
    
    #genotypes is a sparse matrix so getting both indices and indptr
    indices = genotypes.indices #columns of non-zero elements in the sparse matrix
    indptr = genotypes.indptr #rows of non-zero elements in the sparse matrix
    
    #for each genotype grab the mutated loci
    genotype_mloci = [tuple(indices[indptr[g]:indptr[g+1]]) for g in np.arange(genotype_dist.size)]
    
    #count the observed number of loci
    for genotype_id, loci in enumerate(genotype_mloci):
        for locus in loci:
            counts[locus] += genotype_dist[genotype_id]

    return counts

def simulation(args, sim_id, verbose):
    #sim_id to track the simulation
    
    ###### PARAMETERS ###########
    n_passages = args.passages 
    generation_time = args.generation_time #default 6 Baccam, et al 2006; McCrone, et al. bioRxiv 2020
    hours_per_passage = args.hours_per_passage #default 72 hours
    #generations within a passage
    Tpassage = int(hours_per_passage)/generation_time # nu
    Tpassage_vec = np.arange(Tpassage+1) 

    Vinput = args.input_virions #default 10e3
    Vmax = Vinput* args.successful_infection_frac #assuming all virion successfully infect cells
    Vlim = args.max_virions #default 10e8 hard upper limit to viral growth
    
    Cmax = args.target_cells #default 2e5
    d = args.per_capita_decay_rate*(generation_time/24)  #viral clearance rate per generation, default 2 per generation (Han, et al. 2021)
    R0 = args.within_host_basic_reproduction_number #default 10, within-host basic reproduction number (Morris, et al. 2020)
    r = args.virion_production_number #default 10, average number of virions produced by a cell (Morris, et al. 2020)
    beta = (R0*d)/(Cmax*r)
    
    L = args.genome_length #length of binary genome 
    mu = args.mutation_rate #default 3e-6 copied from Han, et al. 2021
    
    apical_release_frac = args.apical_release_fraction #fraction of virions release in the apical wash, default 0.5 assumption
    n_apical_washes = args.n_apical_washes #default 2
    apical_wash_volume_ml = args.apical_wash_vol #default 0.1 > 100 uL
    wash_recovery_frac = args.wash_recovery_frac #default 1.0 assumption
    sampling_volume_ml = args.infection_vol #default 0.05 > 50 uL 
    #determine the fraction for passaging 
    passage_frac = sampling_volume_ml/(n_apical_washes*apical_wash_volume_ml)
    
    sample_every_nth_passage_seq = args.sequence_every_nth_passage #default every 10th 
    n_sequence_samples = args.n_sequence_samples #default sequencing in duplicate
    sequence_sample_volume_ml = args.sequence_sample_vol
    #sequence sample gets taken from the remainder of the apical washes, after passage sample has been taken
    #and then gets split into two
    sequencing_frac = (n_sequence_samples*sequence_sample_volume_ml)/(n_apical_washes*apical_wash_volume_ml)
    remaining_wash_volume = (n_apical_washes*apical_wash_volume_ml)*(1-sequencing_frac)
    passage_frac = sampling_volume_ml/remaining_wash_volume
    
    
    ###### SIMULATION ###########
    # initialise monoclonal population 
    genotype_dist = np.zeros(1) + Vmax # distribution of genotypes
    # we start the experiment with a pure population > all genotype are the same 
    genotypes = csr_matrix((1, L), dtype=np.int8) # genotype boolean array (if < 0 = all WT)
    
    allele_counts = [] #track allele count for each passage that we end up 'sequencing'
    for p in np.arange(n_passages):
        if verbose:
            print (f"growing passage {p+1}, simulation {sim_id+1}")
            
        #let the virus grow
        new_genotype_dist, new_genotypes = simulate_growth_during_passage(genotype_dist, genotypes, Cmax, Vmax, Vlim, Tpassage_vec, mu, L, beta, d, r)
    
        #determine the number of virions after growth 
        currV = new_genotype_dist.sum()

        
        #first checking after sampling from apical and taking a part of the wash sample there is still enough virus present in 
        #the wash sample > see below the while loop for more detailed explaination of the steps
        npV = np.int32(np.int32(np.int32(currV * apical_release_frac)* wash_recovery_frac)*passage_frac)
        #sample used to infect the next passage must have N virions of all least  Vmax or greater
        while npV < Vinput:
            if verbose:
                print (f"target viral growth not matched in simulation {sim_id} passage {p}")
            new_genotype_dist, new_genotypes = simulate_growth_during_passage(genotype_dist, genotypes, Cmax, Vmax, Vlim, Tpassage_vec, mu, L, beta, d, r)
            currV = new_genotype_dist.sum()
            npV = np.int32(np.int32(np.int32(currV * apical_release_frac)* wash_recovery_frac)*passage_frac)

        #### SAMPLING ########

        #determine number of V virions released apically
        apicalV = np.int32(currV * apical_release_frac)
        #resampling genotypes present in apical wash based on genotype proportion 
        genotype_prop = new_genotype_dist / np.sum(new_genotype_dist)
        apical_genotype_dist = resample_genotype_dist(apicalV, genotype_prop)
        #removing extinct genotypes > genotype that got lost during sampling for the passage
        apical_genotype_dist, apical_genotypes = remove_extinct_genotypes(new_genotypes, apical_genotype_dist)
        
        #determine number of virions recovered in the apical washed 
        washV = np.int32(apicalV* wash_recovery_frac)
        if washV != apicalV:
            #once again resampling based on proportion if number of virions is not the same 
            genotype_prop = apical_genotype_dist / np.sum(apical_genotype_dist)
            wash_genotype_dist = resample_genotype_dist(washV, genotype_prop)
            wash_genotype_dist, wash_genotypes = remove_extinct_genotypes(apical_genotypes, wash_genotype_dist)
        else:
            wash_genotype_dist, wash_genotypes = apical_genotype_dist, apical_genotypes
            

        #### SAMPLING FOR SEQUENCING ########
        #every passage sequencing sample gets taken, even if the passage its self is not sequenced 
        seqV = np.int32(washV*sequencing_frac)
            
        #get genotype distirbution in total sequencing sample 
        genotype_prop = wash_genotype_dist / np.sum(wash_genotype_dist)
        seq_genotype_dist = resample_genotype_dist(seqV, genotype_prop)
        
        #remove sequencing sample from the wash genotype dist > genotypes in sequencing sample cannot be passaged anymore
        remaining_genotype_dist = wash_genotype_dist - seq_genotype_dist
        remainingV = washV - seqV
        #removing extinct genotypes > genotype that got sampled to the next passage by change and now are not present in the remaining sample anymore
        remaining_genotype_dist, remaining_genotypes = remove_extinct_genotypes(wash_genotypes, remaining_genotype_dist)
        seq_genotype_dist, seq_genotypes = remove_extinct_genotypes(wash_genotypes, seq_genotype_dist)
        
        #if passage needs to be sequenced performing 
        if (p+1) % sample_every_nth_passage_seq == 0:
            if verbose:
                print (f"sequencing passage {p+1}, simulation {sim_id+1}")

            #for each technical replicate get a sample and determine mutation frequencies, as well as count and depths
            to_sample_V = seqV #safety check 
            techrep_prop = 1/n_sequence_samples 
            for techrep in np.arange(n_sequence_samples):
                techrepV = np.int32(seqV*techrep_prop)
                
                #sample technical replicate
                genotype_prop = seq_genotype_dist / np.sum(seq_genotype_dist)
                techrep_genotype_dist = resample_genotype_dist(techrepV, genotype_prop)
            
                #update seq genotype dist by removing the genotypes that have just been sampled
                new_seq_genotype_dist = seq_genotype_dist -techrep_genotype_dist
                
                #remove extinct genotype for the techincal replicate
                techrep_genotype_dist, techrep_genotypes = remove_extinct_genotypes(seq_genotypes, techrep_genotype_dist)
                
                #finishing updating the seq sampling distrubtion by determining the genotypes that remain
                seq_genotype_dist, seq_genotypes = remove_extinct_genotypes(seq_genotypes, new_seq_genotype_dist)
                
                #get allele counts and depths
                counts = get_allele_counts(techrep_genotype_dist, techrep_genotypes, L)
                allele_counts.append([sim_id, p+1, techrep+1, np.sum(techrep_genotype_dist)]+list(counts))
                
                #safety check ensuring we don't sample more than that we have
                to_sample_V -= techrepV
                if to_sample_V <= 0 :
                    #if all virus are sampled breaking
                    break

        #### SAMPLING FOR NEXT PASSAGE ########
        #from the remaining wash volume we take the sample volume and dilute it to 
        nextpV = remainingV * passage_frac
        genotype_prop = remaining_genotype_dist / np.sum(remaining_genotype_dist)
        nextp_genotype_dist = resample_genotype_dist(nextpV, genotype_prop)
        nextp_genotype_dist, nextp_genotypes = remove_extinct_genotypes(remaining_genotypes, nextp_genotype_dist)
        
        #from to the next p sample we dilute to vinput
        genotype_prop = nextp_genotype_dist / np.sum(nextp_genotype_dist)
        genotype_dist = resample_genotype_dist(Vinput, genotype_prop)
        genotype_dist, genotypes = remove_extinct_genotypes(nextp_genotypes, genotype_dist)
        
        #need to get successfull infection proportion
        genotype_prop = genotype_dist / np.sum(genotype_dist)
        genotype_dist = resample_genotype_dist(Vmax, genotype_prop)
        genotype_dist, genotypes = remove_extinct_genotypes(genotypes, genotype_dist)
        
    return allele_counts

def main():
    args = ArgumentParser()
    verbose = args.verbose
    
    time_start = datetime.datetime.now()
    
    #number of threads
    threads = args.threads 
    if verbose:
        print (f"number of threads: {threads}")

        
    #get output file
    outfile = os.path.join(os.getcwd(), args.output)
    #ensure output file can be made 
    try:
        os.system(f'touch {outfile}')
    except:
        sys.stderr.write(f"Error: cannot create {outfile}. Please ensure path to desired output is correctly specified")
        sys.exit(-1)

    #get number of simulation to perform
    n_sims= args.n_sims
    allele_counts= []
    with mp.Pool(processes=threads) as p:
        for result in tqdm([p.apply_async(simulation, args=(args, i, verbose)) for i in range(n_sims)]):
            counts = result.get()
            allele_counts.extend(counts)

    L = args.genome_length
    allele_counts = pd.DataFrame.from_records(allele_counts, columns=["simulation", "passage", "replicate", "depth"]+list(range(L)))  
    allele_counts.to_csv(outfile, index=False)
    
    time_end = datetime.datetime.now()
    print (f"simulations started at {time_start} and ended at {time_end}. total time: {time_end-time_start}")
     
if __name__ == "__main__":
    main()