#!/usr/bin/env python3
#Script by Annelies de Rooij

import os, sys, snakemake, argparse,re
import pandas as pd, numpy as np
import multiprocessing


filedir = os.path.abspath(os.path.dirname(__file__)) #directory of this file
datadir = os.path.join(filedir, "data")
refdir = os.path.join(datadir, "references")
coding_json = os.path.join(datadir, "CDS_H3N2.json")
cwd = os.getcwd() #working dir


def ArgumentParser():
    """Argument parser"""

    parser = argparse.ArgumentParser(prog = "flu-tube", #changing this later
        formatter_class = argparse.RawTextHelpFormatter,
        description = "Analysis pipeline for raw NGS A/H3N2 data from the LTEE experiment." ) 

    parser.add_argument('-d','--data-dir', required=True, action="store",type=str, help="Directory where the either the raw NGS data files or partial flu-tube output can be found")
    parser.add_argument('-i', '--metadata', required=True, action="store",type=str, help="Metadata file containing the samples to be analyzed")
    parser.add_argument('-o','--output', required=False, action="store",type=str, help="Directory where the output file will be stored (default: ./output)")
    parser.add_argument('-r','--redo',required=False, action="store_true", help="Redo all steps if data already exists")
    parser.add_argument('-rs','--redo-step', required=False,action='store',choices=["RenameSamples","FASTQC", "AnalyzeFASTQC", "TrimReads", "SummarizeTrim", "PostTrimFASTQC", 
                                                                                    "AnalyzePostTrimFASTQC", "MapReads", "Sam2Bam", "PostMapFASTQC", "AnalyzePostMapFASTQC", 
                                                                                    "Deduplicate", "PostDeduplicationFASTQC", "AnalyzePostDeduplicationFASTQC", "ParseBam",
                                                                                    "ComputeCoverage", "CoveragePerSample", "CoveragePlots", "CallSNVsVariants"],
                        help="Redo from a specific step onwards")
    parser.add_argument('-t','--threads',required=False,action='store',type=int,help="Number of threads (default 4 (used for trimmomatic and bowtie,vcftools,minimap))")
    parser.add_argument('--verbose',required=False,action='store_true',help="print stuff to screen")
    
   
    #tool specs > parameters for tools used in snakemake 
    tool_specs = parser.add_argument_group("tool parameters", "parameters flag for tools used within the pipeline.\nall tools have default settings.")
    tool_specs.add_argument('-a', '--adapter', required=False, action="store", choices=["NexteraPE-PE", "TruSeq2-PE", "TruSeq2-SE", "TruSeq3-PE-2", "TruSeq3-PE", "TruSeq3-SE"], default="TruSeq3-PE-2",
                            help="adapter used during sequencing (default: TruSeq3-PE-2)")
    tool_specs.add_argument('-bqt', '--base-qual-threshold', required=False, action="store", type=int, default=30, help="Minimal accepted phred score for base quality (default=30)")
    tool_specs.add_argument('-ml', '--min-length', required=False, action="store", type=int, default=150, help="minimum required length of reads (default:150)")
    tool_specs.add_argument('-mq', '--min-mapq', required=False, action="store", type=int, default=20, help="minimum mapping quality of reads (default:20)")
    tool_specs.add_argument('-bet', '--base-error-tolerance', required=False, action="store", type=float, default=0.001, help="Max tolerated error rate with which variant is called instead of error (default=0.001)")
    tool_specs.add_argument('-mc', '--min-coverage', required=False, action="store", type=int, default=100, help="Minimum coverage nucleotide/amino acid per position for variant calling (default=100)")
     
    #check if arguments were entered 
    if len(sys.argv[1:])<1: #MIGHT CAUSE A PROBLEM LATER ON if wrapped in a setuptools package
        parser.print_help()
        sys.exit(-1)
    else:
        return parser.parse_args()
    
def main():
    args = ArgumentParser()

    #determine number of threads
    if args.threads:
        threads = args.threads
        #check threads are available 
        if threads <= multiprocessing.cpu_count():
            print (f"Number of threads being used is {threads}")
        else:
            threads = multiprocessing.cpu_count()
            print (f"Number of threads being used is {threads}")
    elif multiprocessing.cpu_count() >= 4:
        threads = 4
        print (f"Number of threads being used is {threads}. This is the default number of threads use '-t' to specify the desired number of threads.")
    else:
        threads = multiprocessing.cpu_count()
        print (f"Number of threads being used is {threads}.")


    #get input data folder 
    input_dir = os.path.join(cwd, args.data_dir)
    if not os.path.isdir(input_dir):
        sys.stderr.write(F"Error: cannot find {input_dir}.\nCheck if directory is properly specified")
        sys.exit(-1)

    #get output directory
    if args.output:
        output_dir = os.path.join(cwd, args.output)
    else:
        output_dir = os.path.join(cwd, "output")
    if not os.path.exists(output_dir):
        try:
            os.mkdir(output_dir)
        except:
            sys.stderr.write(f"Error: cannot create directory {output_dir}.\nPlease try creating the output directory manually.")
            sys.exit(-1)

    #find snakefiles
    snakefile = os.path.join(filedir, 'scripts', 'flu-tube.smk')
    if not os.path.exists(snakefile):
        sys.stderr.write(f"Cannot find correct snakefile at {os.path.join(filedir, 'scripts')}.")
        sys.exit(-1)

    #get metadata file and check if it exists
    metadata_file = os.path.join(cwd, args.metadata)
    if not os.path.exists(metadata_file):
        sys.stderr.write(f"Error: cannot find {metadata_file}.\nCheck if file is properly specified")
        sys.exit(-1)

    #read metadata for sample IDs
    metadata = pd.read_excel(metadata_file, engine='openpyxl', keep_default_na=False)
    metadata["sample"] = [f"P{row['Passage']}-{row['Donor']}{row['Replicate']}"  if row["Donor"]!='stock' else row['Donor'] for i, row in metadata.iterrows()]
    metadata["sample_trep"] = [f"P{row['Passage']}-{row['Donor']}{row['Replicate']}_{row['Technical replicate']}"  if row["Donor"]!='stock' else f"{row['Donor']}_{row['Technical replicate']}" for i, row in metadata.iterrows()]
    metadata["long_sample"] = [f"{row['GS Flow cell ID']}_P{row['Passage']}-{row['Donor']}{row['Replicate']}_{row['Technical replicate']}" if row["Donor"]!='stock' else f"{row['GS Flow cell ID']}_{row['Donor']}_{row['Technical replicate']}"  for i, row in metadata.iterrows()]

    original_names = metadata["GS full name"].unique().tolist()
    samples = metadata['sample'].unique().tolist()
    sample_treps = metadata["sample_trep"].unique().tolist()
    long_samples = metadata['long_sample'].unique().tolist()
    treps = metadata["Technical replicate"].unique().tolist()

    reads = ["R1", "R2"]
    #segments = ["1-PB2", "2-PB1", "3-PA", "4-HA", "5-NP", "6-NA", "7-M", "8-NS"]
    segments = ["4-HA"]

    #tool specifications 
    if  args.base_error_tolerance >=0 and args.base_error_tolerance <=1:
        bet = args.base_error_tolerance
    else:
        sys.stderr.write(f"Base error tolerance should be between 0 and 1.")
        sys.exit(-1)
    
    bqt = args.base_qual_threshold
    mc = args.min_coverage
    min_length = args.min_length
    min_mapq= args.min_mapq

    adapter = args.adapter

    #initialize snakemake config
    config = {
        "metadata":metadata_file,
        "original_names":original_names,
        "long_samples":long_samples,
        "reads":reads,
        "sample_treps":sample_treps,
        "samples":samples,
        "treps": treps,
        "segments":segments,
        "input":input_dir,
        "refdir":refdir,
        "datadir":datadir,
        "output":output_dir,
        "bet":bet,
        "bqt":bqt,
        "min_coverage":mc,
        "min_length":min_length,
        "min_mapq":min_mapq,
        "adapter":adapter,
        "verbose":args.verbose,
    }

    if args.redo_step:
        snakemake.snakemake(snakefile, printshellcmds=args.verbose, forceall=args.redo, forcerun=[args.redo_step], config=config,cores=threads,lock=False,quiet=args.verbose)
    else:
        snakemake.snakemake(snakefile, printshellcmds=args.verbose, forceall=args.redo, config=config,cores=threads,lock=False,quiet=args.verbose)


        
if __name__ == "__main__":
    main()