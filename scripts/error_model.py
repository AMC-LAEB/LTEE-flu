#!/usr/bin/env python3


import os, bisect, sys, argparse
import pandas as pd, numpy as np, statistics as st, scipy as sp
from scipy.stats import binom
from scipy.optimize import minimize
from collections import OrderedDict
import multiprocessing as mp
from tqdm import tqdm

cwd = os.getcwd() #working dir
filedir = os.path.abspath(os.path.dirname(__file__)) #directory of this file


def simulate_trajectory(A, D, timepoints):
    "Simulate trajectory under error model as define by good, et al. 2017 "
    #get the emperical frequencies for all snvs in the population
    F = A/D
    #for each snv calculate the average frequency in the population
    avg_f = np.sum(A) / np.sum(D)
    #calculate scale factor
    #c_factor = np.maximum(1, ((1/2)*avg_f)) #based on good,et al. text
    c_factor = 0.5/avg_f if avg_f > 0.5 else 1 #based on good, et al code 
   
    #renormalize frequencies
    cF = F * c_factor
    #scramble temporal information via random permutation
    sigma_hat = np.random.permutation(len(timepoints)) #indices of the random permutation of the timepoints
    #define error probabilities 
    P_hat = cF[sigma_hat]
    #allele count under the error model
    A_hat = np.random.binomial(D, P_hat)

    return A_hat

def calculate_autocorrelation(A, D, timepoints):
    "calculate the modified autocorrelation "
    #calculate the average frequency in the population WEIGHTED BY Dmt
    avg_f = ((np.sum(A) / np.sum(D ))* np.sum(D))/np.sum(D)
    #correct for the polarization constraint; cap avg frequency at 0.5 
    avg_f = np.minimum(avg_f, 1/2)

    autocovariance = 0 #numeration of the modified autocorrelation formula
    covariance_weight = 0 
    #calculate auto covariance 
    for i in range(1, len(timepoints)):
        f = A[i] / D[i]
        fprev = A[i-1] / D[i-1]
        weight = np.sqrt(D[i]+D[i-1])
        #get difference between f and f_macron
        read_difference = abs(A[i]-D[i]*avg_f)
        prev_read_difference =  abs(A[i-1]-D[i-1]*avg_f)
        if read_difference >= 1 and prev_read_difference >=1: #heaviside function
            autocovariance += (f-avg_f)*(fprev-avg_f)*weight
        covariance_weight += weight
    
    #calculate autocorrelation
    autocorrelation = np.sqrt(max(autocovariance/avg_f/avg_f, 0))

    return autocorrelation

def calculate_relaxation_time(A,D,timepoints, bootstrap, relaxation_cutoff=2):
    #for each snv calculate the average frequency in the population WEIGHTED BY Dmt
    avg_f = ((np.sum(A) / np.sum(D))* np.sum(D))/np.sum(D)

    #determine number of timepoints and of zero timepoints
    n_timepoints, n_zero_timepoints = len(timepoints), 0
    for i in range(len(timepoints)):
        if A[i] < 1:
            n_zero_timepoints += 1

    # T = 0 if 
    if n_zero_timepoints > 0.3*n_timepoints:
        if bootstrap:
            return len(timepoints)
        else:
            return -1
    
    #calculate relaxation time T
    running_avg_f, running_avg_weight, t = 0, 0, -1
    for i in range(1, len(timepoints)):
        f = A[i]/D[i]
        running_avg_f += f
        running_avg_weight += 1

        #skipping over the first n time points
        if i>relaxation_cutoff and (running_avg_f/running_avg_weight > avg_f*0.6):
            t = i
            break

    if t < 0:
        t = len(timepoints)

    return t

def calculate_max_run(A, D, timepoints):
    avg_f = np.sum(A) / np.sum(D)
    
    avg_f = np.minimum(avg_f, 1/2) #correct for the polarization constraint; cap avg frequency at 0.5 
    
    #determine number of timepoints and of zero timepoints
    n_timepoints, n_zero_timepoints = len(timepoints), 0
    for i in range(len(timepoints)):
        if A[i] < 1:
            n_zero_timepoints += 1

    if n_zero_timepoints > 0.3*n_timepoints:
        threshold_f = avg_f/(1+np.exp((n_zero_timepoints-0.3*n_timepoints)/5))
    else:
        threshold_f = avg_f/(1+np.exp(-(0.3*n_timepoints-n_zero_timepoints)/5))

    #calculate the max run (I)
    current_run, current_run_length , max_run = 0,0,0
    for i in range(len(timepoints)):
        #Amt[i,j] must be large that the threshold f * Depth
        if A[i] - (D[i]*threshold_f) > 1:
            current_run += (A[i]/D[i])-threshold_f
            current_run_length += 1
        else:
            #end the current run
            if current_run > max_run and current_run_length > 1: #run length should be atleast two time points
                max_run = current_run
            current_run = 0
            current_run_length = 0
        if max_run == 0 and current_run_length > 1:
            max_run = current_run

    return max_run

def calculate_pvalue_from_list(stat_list, query_stat):
    stat_list = sorted(stat_list)
    n_greater = len(stat_list) - bisect.bisect_left(stat_list, query_stat)
    return (n_greater+1)/(len(stat_list)+1)

def clip_pvalue(pvalue, pvalue_threshold=(0.05)**(1/3)):
    return 1 if pvalue > pvalue_threshold else pvalue

def calculate_composite_statistic(pvalue_list):
    composite_statistic = clip_pvalue(pvalue_list[0])
    for pvalue in pvalue_list[1:]:
        composite_statistic *= clip_pvalue(pvalue)#*np.log10(1/pvalue)
    return composite_statistic

def calculate_final_pvalue(A, D, timepoints, bootstraps=10000, relaxation_cutoff=5, i=None, max_extra_bootstraps=1000000):
    observed_autocorrelation = calculate_autocorrelation(A, D, timepoints)
    observed_max_run = calculate_max_run(A, D, timepoints)
    observed_relaxation_time = calculate_relaxation_time(A, D, timepoints, False, relaxation_cutoff=relaxation_cutoff)

    #simulate large number of draws
    ac_greater, mr_greater, rt_greater = 0,0,0

    bootstraps_autocorrelation = []
    bootstraps_max_run = []
    bootstraps_relaxation_time = []
    #calculate observation pvalue
    for draw in range(bootstraps):
        A_hat = simulate_trajectory(A, D, timepoints)
        
        bootstrap_autocorrelation = calculate_autocorrelation(A_hat,D, timepoints)
        bootstrap_max_run = calculate_max_run(A_hat, D, timepoints)
        bootstrap_relaxation_time = calculate_relaxation_time(A_hat, D, timepoints, True, relaxation_cutoff=relaxation_cutoff)
        
        bootstraps_autocorrelation.append(bootstrap_autocorrelation)
        bootstraps_max_run.append(bootstrap_max_run)
        bootstraps_relaxation_time.append(bootstrap_relaxation_time)

        #determine if error is greater than observation
        if bootstrap_autocorrelation > observed_autocorrelation:
            ac_greater += 1
        if bootstrap_max_run > observed_max_run:
            mr_greater += 1
        if bootstrap_relaxation_time > observed_relaxation_time:
            rt_greater += 1

    #calculate single statistic p-value numerically #psuedo count to avoid zero divisions
    observed_autocorrelation_pvalue = (ac_greater+1)/(bootstraps+1)
    observed_max_run_pvalue = (mr_greater+1)/(bootstraps+1)
    observed_relaxation_time_pvalue = (rt_greater+1)/(bootstraps+1)

    #calculate combined pvalue
    observed_combined_pvalue = calculate_composite_statistic([observed_autocorrelation_pvalue, observed_max_run_pvalue, observed_relaxation_time_pvalue])
    for draw in range(bootstraps):
        A_hat = simulate_trajectory(A, D, timepoints)
        
        bootstrap_autocorrelation = calculate_autocorrelation(A_hat,D, timepoints)
        bootstrap_max_run = calculate_max_run(A_hat, D, timepoints)
        bootstrap_relaxation_time = calculate_relaxation_time(A_hat, D, timepoints, True, relaxation_cutoff=relaxation_cutoff)
        
        bootstraps_autocorrelation.append(bootstrap_autocorrelation)
        bootstraps_max_run.append(bootstrap_max_run)
        bootstraps_relaxation_time.append(bootstrap_relaxation_time)

    observed_autocorrelation_pvalue = calculate_pvalue_from_list(bootstraps_autocorrelation, observed_autocorrelation)
    observed_max_run_pvalue = calculate_pvalue_from_list(bootstraps_max_run, observed_max_run)
    observed_relaxation_time_pvalue = calculate_pvalue_from_list(bootstraps_relaxation_time, observed_relaxation_time)
    #calculate combined pvalue
    observed_combined_pvalue = calculate_composite_statistic([observed_autocorrelation_pvalue, observed_max_run_pvalue, observed_relaxation_time_pvalue])

    n_greater = 0
    current_draw = 0
    for draw in range(bootstraps):
        current_draw += 1
        A_hat = simulate_trajectory(A, D, timepoints)
       
        
        bootstrap_autocorrelation = calculate_autocorrelation(A_hat,D, timepoints)
        bootstrap_max_run = calculate_max_run(A_hat, D, timepoints)
        bootstrap_relaxation_time = calculate_relaxation_time(A_hat, D, timepoints, True, relaxation_cutoff=relaxation_cutoff)
        
        bootstrap_autocorrelation_pvalue = calculate_pvalue_from_list(bootstraps_autocorrelation, bootstrap_autocorrelation)
        bootstrap_max_run_pvalue = calculate_pvalue_from_list(bootstraps_max_run, bootstrap_max_run)
        bootstrap_relaxation_time_pvalue = calculate_pvalue_from_list(bootstraps_relaxation_time, bootstrap_relaxation_time)
        
        bootstrap_combined_pvalue = calculate_composite_statistic([bootstrap_autocorrelation_pvalue, bootstrap_max_run_pvalue, bootstrap_relaxation_time_pvalue])

        #print (bootstrap_combined_pvalue, observed_combined_pvalue)

        if bootstrap_combined_pvalue <= observed_combined_pvalue:
            

            n_greater += 1

    #calculate final Pvalue 
    pvalue = (n_greater+1)/(current_draw+1)
    
    #decrease uncertainty for pvalue near boundary
    if 0.01 < pvalue < 0.05:
        if current_draw < max_extra_bootstraps:
            current_draw += 1
        A_hat = simulate_trajectory(A, D, timepoints)
       
        
        bootstrap_autocorrelation = calculate_autocorrelation(A_hat,D, timepoints)
        bootstrap_max_run = calculate_max_run(A_hat, D, timepoints)
        bootstrap_relaxation_time = calculate_relaxation_time(A_hat, D, timepoints, True, relaxation_cutoff=relaxation_cutoff)
        
        bootstrap_autocorrelation_pvalue = calculate_pvalue_from_list(bootstraps_autocorrelation, bootstrap_autocorrelation)
        bootstrap_max_run_pvalue = calculate_pvalue_from_list(bootstraps_max_run, bootstrap_max_run)
        bootstrap_relaxation_time_pvalue = calculate_pvalue_from_list(bootstraps_relaxation_time, bootstrap_relaxation_time)
        
        bootstrap_combined_pvalue = calculate_composite_statistic([bootstrap_autocorrelation_pvalue, bootstrap_max_run_pvalue, bootstrap_relaxation_time_pvalue])

        #print (bootstrap_combined_pvalue, observed_combined_pvalue)
        if bootstrap_combined_pvalue <= observed_combined_pvalue:
            n_greater += 1
            
    pvalue = (n_greater+1)/(current_draw+1) 
    
    
    if i is not None:
        return pvalue,i
    else:
        return pvalue

def ArgumentParser():
    """Argument parser"""

    parser = argparse.ArgumentParser(prog = "error model", #changing this later
        formatter_class = argparse.RawTextHelpFormatter,
        description = "Error model to pvalues to identify true SNVs per population from LTEE data SNV" ) 

    parser.add_argument('-d','--depth-file', required=True, action="store",type=str, help="csv file with the depths per position per timepoint")
    parser.add_argument('-a','--allele-count-file', required=True, action="store", type=str, help="csv file with the alternate allele count per position per timepoint")
    parser.add_argument('-o','--output', required=False, action="store",type=str, help="output file name (default: ./pvalues.csv) ")
    parser.add_argument('-b','--bootstrap', required=False, action="store", type=int, default=10000, help="number of bootstrap simulations to perform error model statistics on")
    parser.add_argument('-t', '--threads', required=False, action="store", type=int, default=1, help="number of threads to use (default: 1)")
    parser.add_argument('--verbose',required=False,action='store_true',help="print stuff to screen")

    
    #check if arguments were entered 
    if len(sys.argv[1:])<1: 
        parser.print_help()
        sys.exit(-1)
    else:
        return parser.parse_args()
    
def main():
    args = ArgumentParser()

    verbose = True if args.verbose else False

    #check if depth file and allele count files exist
    depth_file = os.path.join(cwd, args.depth_file)
    if not os.path.isfile(depth_file):
        sys.stderr.write(F"Error: cannot find {depth_file}.\n Check if file is properly specified")
        sys.exit(-1)

    allele_count_file = os.path.join(cwd, args.allele_count_file)
    if not os.path.isfile(allele_count_file):
        sys.stderr.write(F"Error: cannot find {allele_count_file}.\n Check if file is properly specified")
        sys.exit(-1)

    depths = pd.read_csv(depth_file).set_index(["refname", "nucpos"])
    allele_count = pd.read_csv(allele_count_file).set_index(["refname", "nucpos", "read_base"])
    timepoints = depths.columns.tolist()
    
    #get depth and allele count matrix for population
    Dmt = depths[timepoints].values
    Amt = allele_count[timepoints].values

    #check if Dmt and Amt are the same shape
    if Dmt.shape != Amt.shape:
        sys.stderr.write(F"Error: Dmt shape and Amt shape are not equal")
        sys.exit(-1)

    bootstraps = args.bootstrap
    threads = args.threads

    #get raw  p-value 
    if not args.output:
        p_value_file = os.path.join(cwd, f"pvalues.csv")
    else:
        p_value_file = os.path.join(cwd, args.output)

    if not os.path.isfile(p_value_file):
        if verbose:
            sys.stdout.write(f"Calculating p-values for trajectories")
        pvalues = np.zeros(Amt.shape[0])

        with mp.Pool(processes=threads) as p:
            for result in tqdm([p.apply_async(calculate_final_pvalue, args=(Amt[i], Dmt[i], timepoints, bootstraps, 5, i)) for i in range(Amt.shape[0])]):
                pvalue, i  = result.get()
                pvalues[i] = pvalue
        
        #save raw pvalues
        pd.DataFrame(pvalues, columns=["p-value"]).to_csv(p_value_file, index=False)
    


if __name__ == "__main__":
    main()
