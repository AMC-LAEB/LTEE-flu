#/usr/bin/env python3
import re,gzip, json, matplotlib, pysam
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from Bio import SeqIO
from scipy.stats import binom
from scipy.special import loggamma
from collections import OrderedDict
import statistics as st
import itertools
from scipy.optimize import minimize

#majority of code has been copied and edited from scripts by Alvin X. Han

#using agg as matplotlib backend to avoid crashing during multi-treading
matplotlib.use('Agg')


NL_08_coding_regions = {"1-PB2":{"start":28, "end":2308}, "2-PB1":{"start":25, "end":2298}, "3-PA":{"start":25, "end":2175}, "4-HA":{"start":30, "end":1730},
                        "5-NP":{"start":46, "end":1542}, "6-NA":{"start":20, "end":1429}, "7-M":{"start":26, "end":1010}, "8-NS":{"start":27, "end":868}}

codon_table = {'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 'AAC':'N',
               'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R', 'CTA':'L', 'CTC':'L',
               'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 'CAC':'H', 'CAT':'H', 'CAA':'Q',
               'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
               'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G',
               'GGC':'G', 'GGG':'G', 'GGT':'G', 'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 
               'TTA':'L', 'TTG':'L', 'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*', 'TGC':'C', 'TGT':'C', 'TGA':'*',
               'TGG':'W'}
all_bases = ['A', 'T', 'G', 'C', '-', 'N']
all_amino_acids = list(set(codon_table.values())) + ["X"]

cc =  ["refname", "protein", "nucpos", "aapos", "read_base", "snv", "read_codon", "read_aa", "aa_mut", "non_synonymous"]
snv_columns = ['refname', 'nucpos', 'coverage', 'read_base', 'nuc_count', 'nuc_freq', 'snv',]



class NpEncoder(json.JSONEncoder):
    "JSON encoder object > transforms np data types to regular python datatypes"
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super(NpEncoder, self).default(obj)

def analyze_fastqc(fastqcf,bqt=30):
    """Analyze FASTQC results"""
    with open(fastqcf) as f: #readlines
        lines = [line for line in f]

    mod = ''
    df_dict = {}
    for l in lines:
        if re.search('>>Per base sequence quality', l):
            mod= 'pbsq'
            continue 
        elif re.search('>>Adapter Content', l):
            mod = 'ac'
            continue
        elif re.search('Per sequence quality score', l):
            mod = 'psqs'
            continue 
        
        if mod != '':
            if re.search('>>END_MODULE', l):
                df_dict[mod] = pd.DataFrame(df)
                mod = ''
            elif re.search('^#', l):
                headers = re.sub('#', '', l.strip()).split('\t')
                df = {h:[] for h in headers}
            else:
                values = l.strip().split("\t")
                for _h,h in enumerate(headers):
                    df[h].append(values[_h])

        
    # % of reads with acceptable phred score > bqt
    df = df_dict['psqs']
    df['Quality'] = pd.to_numeric(df['Quality'], errors='coerce')
    df['Count'] = pd.to_numeric(df['Count'], errors='coerce')
    percent_abv_qualthres = sum(df[df['Quality']>=bqt]['Count'])/sum(df['Count'])

    # max position where median phred score > bqt
    df = df_dict['pbsq']
    df['Median'] = pd.to_numeric(df['Median'], errors='coerce')

    try:
        start_pos = int(re.search('^(\d+)-*', df[df['Median']>=bqt]['Base'].iloc[0]).group(1))
    except:
        start_pos = None
    try:
        end_pos = int(re.search('-*(\d+)$', df[df['Median']>=bqt]['Base'].iloc[-1]).group(1))
    except:
        end_pos = None


    df = df_dict['ac'].drop('Position', axis=1)
    ac = {}
    for adapter in list(df):
        val = np.sum(pd.to_numeric(df[adapter], errors='coerce'))/len(df[adapter])
        if adapter in ac:
            if val > ac[adapter]:
                ac[adapter] = val
        else:
            if val > 0.:
                ac[adapter] = val

    return percent_abv_qualthres, ac, start_pos, end_pos

def get_ref_pos_mapping(reffile, coding_regions=NL_08_coding_regions):
    """get reference mapping positions"""
    ref_seqs = {rec.id:rec.seq for rec in SeqIO.parse(reffile, "fasta")}

    ref_pos_mapping = {}
    for segment, seq in ref_seqs.items():
        start_pos = coding_regions[segment]["start"]-1 

        mapping = []
        count = 0
        for i in range(len(seq)): 
            mapping.append(i-start_pos)
        
        ref_pos_mapping[segment] = np.array(mapping)
    return ref_pos_mapping

def summarize_minimap(minimap_file, minimap_ref, segments, targets, bqt=30):
    """summarize minimap results and report non target hits"""
    non_target_reads = {}
    ref_gene_count = np.zeros((len(minimap_ref), len(segments)), dtype=np.int32)
    targets = ["stock|A/H3N2"]
    with gzip.open(minimap_file, "rb") as f:
        lines = f.readlines()

    for line in lines:
        line = line.decode('utf-8').strip()
        if line.startswith('@'):
            continue

		#read fields
        fields = line.split('\t')
        read_id, flag, refname, cigar = fields[0], int(fields[1]), fields[2], fields[5]
        read_seq = fields[9]
        read_qual = fields[10]

        read_qual= [ord(qual)-33 for qual in read_qual]

        #filter shit reads
        flag_set = {147, 163, 99, 83, 131, 179, 67, 115}
        read_pair_id = 2 if flag in {147, 163, 131, 179} else 1

        if flag not in flag_set or cigar == '*':
            continue

		#get strain name and segment
        sname  = '|'.join(refname.split('|')[2:4]).replace("_", "")
        segment = refname.split("|")[-1]
        ref_gene_count[minimap_ref.index(sname), segments.index(segment)] += 1

        if sname not in targets:
            #remove low quality bases
            read_seq = "".join([b if read_qual[i] >= bqt else "N" for i, b in enumerate(read_seq)])
            non_target_reads[f"{read_id}|{sname}|{segment}"] = read_seq
			

	#summarize per strain > df for overview     	
    results = []
    for i, sname in enumerate(minimap_ref):
        results.append([sname]+list(ref_gene_count[i]))

    results = pd.DataFrame.from_records(results, columns=["reference strain"]+segments)

    return results, non_target_reads

def parse_bam(bamfile, ref_pos_mapping, refname_dict,  bqt, mapq=20, non_target_reads=None):
    """parse proper concordant reads"""
    map_json = {}
    
    #read sam file
    with pysam.AlignmentFile(bamfile, "rb") as samfile:
    #print (samfile.header) 
        for read in samfile:
            fields = str(read).split("\t")
    
            read_id, flag, refname, start_pos, mq, cigar = fields[0], int(fields[1]), refname_dict[fields[2]], int(fields[3])-1, int(fields[4]), fields[5]
            if mq < mapq:
                continue
            #filtering out unwanted reads if specified
            if non_target_reads is not None:
                if read_id in non_target_reads:
                    continue

            read_seq, read_qual = fields[9], fields[10]
            read_qual = [int(i) for i in read_qual.split("[")[-1].split("]")[0].split(",")]
            
            flag_set = {147, 163, 99, 83, 131, 179, 67, 115}
            read_pair_id = 2 if flag in {147, 163, 131, 179} else 1

            if flag not in flag_set or cigar == '*':
                continue

            expanded_cigar =[i for l in [[m.group(2)]*int(m.group(1)) for m in re.finditer('(\d+)([MIDNSHP\=X])', cigar) if m.group(2) not in {"H", "P"}] for i in l]

            curr_aln_num, rbase_idx = -1, -1
            read_nuc_seq, base_quality, ref_positions, insertion_nuc_seq = [], [], [], {}


            #Iterate over cigar mapping pos and get mapped nucleotide sequence, base quality per pos, and the reference pos 
            for expr in expanded_cigar:
                # H/P consumes neither query or reference
                if expr in ["H", "P"]:
                    continue

                #I/S consumes only query:
                elif expr in ["I", "S"]: #insertion
                    rbase_idx +=1
                    
                    read_base = read_seq[rbase_idx]
                    if re.search("[^ATGC]", read_base):
                        read_base = "N"
                    else:
                        bq = read_qual[rbase_idx]#-33
                        if bq < bqt:
                            read_base = "N"
                    
                    if expr == "I":
                        ref_pos = int(ref_pos_mapping[refname][start_pos + curr_aln_num])
                        try:
                            insertion_nuc_seq[ref_pos].append(read_base)
                        except:
                            insertion_nuc_seq[ref_pos] = [read_base]
                    
                    
                else:  # add to reference_position for M=X and D/N
                    curr_aln_num += 1
                    ref_pos = ref_pos_mapping[refname][start_pos+ curr_aln_num]
                    ref_positions.append(ref_pos)

                    # D/N consumes only reference
                    if re.search('[DN]', expr):
                        if expr == "D":
                            read_nuc_seq.append('-') # deletion
                        else:
                            read_nuc_seq.append('N') # skipped region from reference
                        base_quality.append("J") # assume highest PHRED33 quality score
                    else:
                        rbase_idx += 1

                        # only tally ATGC
                        read_base = read_seq[rbase_idx]
                        # only tally bases > required quality threshold and isn't ambiguous
                        bq = read_qual[rbase_idx]#-33
                       

                        if (bq < bqt) or (re.search("[^ATGC]", read_base)):
                            read_base = "N"

                        # add read_base to read_nuc_seq
                        read_nuc_seq.append(read_base)
                        # save base quality
                        base_quality.append(chr(read_qual[rbase_idx]+33))
                    

            # save to json
            ref_positions = np.array(ref_positions)
            read_nuc_seq = np.array(read_nuc_seq)
            base_quality = np.array(base_quality)

            #get only part of the arrays where the actual reference sequence is covered
            read_nuc_seq = read_nuc_seq[ref_positions>-1]
            
            if len(read_nuc_seq) == 0:
                continue
            
            base_quality = base_quality[ref_positions>-1]
            ref_positions = ref_positions[ref_positions>-1]

            read_nuc_seq = "".join(read_nuc_seq)
            base_quality = "".join(base_quality)

            if len(insertion_nuc_seq) > 0:
                for pos, temp_list in insertion_nuc_seq.items():
                    insertion_nuc_seq[pos] = "".join(temp_list)

            temp_dict = {"refname":refname, "start_pos":ref_positions[0], "read_seq":read_nuc_seq, "base_qual":base_quality, "insertion_seq":insertion_nuc_seq, "map_qual":mq}

            try:
                map_json[read_id][f"R{int(read_pair_id)}"] = temp_dict
            except:
                map_json[read_id] = {f"R{int(read_pair_id)}":temp_dict}

    return map_json

def map_single_read(refname, seq_positions, protein_cds_aa, mapped_nucleotides, mapped_codons):

    positions = list(seq_positions.keys())
    cds = list(set([i for p in protein_cds_aa.keys() for i in protein_cds_aa[p]["nuc_cds_pos"]]))
    mask = np.isin(cds, positions)
    nucpos = np.array(cds)[mask]

    sequence = list(seq_positions.values())

    for nuc_pos, base in zip(nucpos, sequence):
        try: 
            mapped_nucleotides[(refname, nuc_pos)][base] += 1
        except:
            try:
                mapped_nucleotides[(refname, nuc_pos)][base] = 1
            except:
                mapped_nucleotides[(refname, nuc_pos)] = {base:1}
            

    #codon portion
    for protein, d in protein_cds_aa.items():
        nuc_cds_pos = d["nuc_cds_pos"]
        pro_aa_pos = d["pro_aa_pos"]
        
        mask = np.isin(nuc_cds_pos, positions)
        codon_nucpos = np.array(nuc_cds_pos)[mask]
        codon_aapos = np.array(pro_aa_pos)[mask]

        for i, aapos in enumerate(np.sort(np.unique(codon_aapos))):
            read_codon_pos = codon_nucpos[codon_aapos==aapos]
            codon = np.array(sequence)[np.isin(positions, read_codon_pos)]
            if len(codon) < 3 or "N" in codon or "-" in codon: 
                continue
            #read_aa = codon_table["".join(codon)]
            codon = "".join(codon)
            
            try:
                mapped_codons[(refname, protein,aapos)][codon] += 1
            except:
                try:
                    mapped_codons[(refname, protein,aapos)][codon] = 1
                except:
                    mapped_codons[(refname, protein,aapos)] = {codon:1}
                

    return mapped_nucleotides, mapped_codons

def map_reads(sample_json, coding_json, bqt=30, mapped_nucleotides=None, mapped_codons=None, non_target_reads=None):
    if mapped_nucleotides is None:
        mapped_nucleotides = {}
        mapped_codons = {}
    for read_id, reads in sample_json.items():

        #filtering out unwanted reads if specified
        if non_target_reads is not None:
            if read_id in non_target_reads:
                continue

        #both reads must be present
        if set(sample_json[read_id].keys()) != set(['R1', 'R2']):
            continue

        refname = reads["R1"]["refname"]

        #get info per read pair id 
        r1_positions = np.arange(reads["R1"]["start_pos"], reads["R1"]["start_pos"]+len(reads["R1"]["read_seq"])) #indices
        r1_sequence = np.array(list(reads["R1"]["read_seq"]))
        r1_bq = np.array([ord(t)-33 for t in list(reads["R1"]["base_qual"])])
        
        r2_positions = np.arange(reads["R2"]["start_pos"], reads["R2"]["start_pos"]+len(reads["R2"]["read_seq"])) #indices
        r2_sequence = np.array(list(reads["R2"]["read_seq"]))
        r2_bq = np.array([ord(t)-33 for t in list(reads["R2"]["base_qual"])])

        #determine overlap
        overlap = np.intersect1d(r1_positions, r2_positions)
        if len(overlap) > 0:

            # get overlapping regions for read pair
            r1_overlap_mask = np.isin(r1_positions, overlap)
            r2_overlap_mask = np.isin(r2_positions, overlap)

            r1_overlap_read = r1_sequence[r1_overlap_mask]
            r2_overlap_read = r2_sequence[r2_overlap_mask]

            if (r1_overlap_read==r2_overlap_read).all() != True:
                r1_overlap_positions = r1_positions[r1_overlap_mask]
                r2_overlap_positions = r2_positions[r2_overlap_mask]
                
                r1_overlap_bq = r1_bq[r1_overlap_mask]
                r2_overlap_bq = r2_bq[r2_overlap_mask]

                mask = r1_overlap_read != r2_overlap_read #get position where read bases are different
                for i, (r1_base, r2_base) in enumerate(zip(r1_overlap_read[mask], r2_overlap_read[mask])):
                    
                    if r1_base == "N" and re.search("[ATGC]", r2_base): 
                        r1_sequence[np.where(r1_positions==r1_overlap_positions[mask][i])] = r2_sequence[np.where(r2_positions==r2_overlap_positions[mask][i])]
                        r1_bq[np.where(r1_positions==r1_overlap_positions[mask][i])] = r2_bq[np.where(r2_positions==r2_overlap_positions[mask][i])]

                    elif re.search("[ATGC]", r1_base) and re.search("[ATGC]", r2_base):
                        if r1_overlap_bq[mask][i] < r2_overlap_bq[mask][i]: #better bq read_2
                            r1_sequence[np.where(r1_positions==r1_overlap_positions[mask][i])] = r2_sequence[np.where(r2_positions==r2_overlap_positions[mask][i])]
                        elif r1_overlap_bq[mask][i] == r2_overlap_bq[mask][i] and r2_base!=r1_base: #ambiguous
                            r1_sequence[np.where(r1_positions==r1_overlap_positions[mask][i])] = "N"
                    elif r1_base == "-" and r1_base != r2_base:
                        r1_sequence[np.where(r1_positions==r1_overlap_positions[mask][i])] = r2_sequence[np.where(r2_positions==r2_overlap_positions[mask][i])]
                        r1_bq[np.where(r1_positions==r1_overlap_positions[mask][i])] = r2_bq[np.where(r2_positions==r2_overlap_positions[mask][i])]


            # remove overlapped section in r2
            r2_positions = r2_positions[~r2_overlap_mask]
            r2_sequence = r2_sequence[~r2_overlap_mask]

        positions = dict(map(lambda i,j : (i+1,j) , r1_positions,r1_sequence))
        positions.update(dict(map(lambda i,j : (i+1,j) , r2_positions,r2_sequence)))
        positions = dict(OrderedDict(sorted(positions.items())))

        #get base quality
        bq = dict(map(lambda i, j: (i+1, j), r1_positions, r1_bq))
        bq.update(dict(map(lambda i, j: (i+1, j), r2_positions, r2_bq)))
        bq = dict(OrderedDict(sorted(bq.items())))

        #remove bases below the quality threshold 
        for p,  q in bq.items():
            if q < bqt:
                positions[p] = "N"

        mapped_nucleotides, mapped_codons = map_single_read(refname, positions, coding_json[refname], mapped_nucleotides, mapped_codons)

    return mapped_nucleotides, mapped_codons

# def map_reads(sample_json, coding_json, bqt=30, mapped_nucleotides=None, mapped_codons=None, non_target_reads=None):
#     if mapped_nucleotides is None:
#         mapped_nucleotides = {}
#         mapped_codons = {}
#     for read_id, reads in sample_json.items():
       
#         #filtering out unwanted reads if specified
#         if non_target_reads is not None:
#             if read_id in non_target_reads:
#                 continue

#         if set(reads.keys()) != set(['R1', 'R2']):
#             continue
    
#         #refname 
#         refname = reads["R1"]["refname"]

#         #get info per read pair id 
#         r1_positions = np.arange(reads["R1"]["start_pos"], reads["R1"]["start_pos"] + len(reads["R1"]["read_seq"]))
#         r1_sequence = np.frombuffer(reads["R1"]["read_seq"].encode(), dtype="S1").copy()
#         r1_bq = np.fromiter((ord(c) - 33 for c in reads["R1"]["base_qual"]), dtype=np.uint8)

#         r2_positions = np.arange(reads["R2"]["start_pos"], reads["R2"]["start_pos"] + len(reads["R2"]["read_seq"]))
#         r2_sequence = np.frombuffer(reads["R2"]["read_seq"].encode(), dtype="S1").copy()
#         r2_bq = np.fromiter((ord(c) - 33 for c in reads["R2"]["base_qual"]), dtype=np.uint8)

#         #determine overlap
#         overlap, r1_idx, r2_idx = np.intersect1d(r1_positions, r2_positions, return_indices=True)

#         if overlap.size > 0:

#             r1_overlap_read = r1_sequence[r1_idx]
#             r2_overlap_read = r2_sequence[r2_idx]
#             r1_overlap_bq = r1_bq[r1_idx]
#             r2_overlap_bq = r2_bq[r2_idx]

#             #get positions where bases are different
#             mismatch_mask = r1_overlap_read != r2_overlap_read
#             mismatch_indices = np.where(mismatch_mask)[0]
        
#             for i in mismatch_indices:
#                 r1_base = r1_overlap_read[i].decode()
#                 r2_base = r2_overlap_read[i].decode()

#                 if r1_base == "N" and re.search("[ATGC]", r2_base):
#                     r1_sequence[r1_idx[i]] = r2_sequence[r2_idx[i]]

#                 elif re.search("[ATGC]", r1_base) and re.search("[ATGC]", r2_base):
#                     if r1_overlap_bq[i] < r2_overlap_bq[i]: #better base qualitiy of read 2
#                         r1_sequence[r1_idx[i]] = r2_sequence[r2_idx[i]]
#                         r1_bq[r1_idx[i]] = r2_bq[r2_idx[i]]
#                     elif r1_overlap_bq[i] == r2_overlap_bq[i] and r1_base != r2_base: #ambiguous
#                         r1_sequence[r1_idx[i]] = b"N"
#                 elif r1_base == "-" and r2_base != "-":
#                     r1_sequence[r1_idx[i]] = r2_sequence[r2_idx[i]]
#                     r1_bq[r1_idx[i]] = r2_bq[r2_idx[i]]

#             # remove overlapped section in r2
#             r2_keep = np.ones(len(r2_sequence), dtype=bool)
#             r2_keep[r2_idx] = False
#             r2_positions = r2_positions[r2_keep]
#             r2_sequence = r2_sequence[r2_keep]
#             r2_bq = r2_bq[r2_keep]

#         all_positions = np.concatenate([r1_positions, r2_positions])
#         all_sequence = np.concatenate([r1_sequence, r2_sequence])
#         all_bq = np.concatenate([r1_bq, r2_bq])
        
#         all_sequence[all_bq < bqt] = b"N"
        
#         positions = {pos + 1: base.decode() for pos, base in zip(all_positions, all_sequence)}


#         mapped_nucleotides, mapped_codons = map_single_read(refname, positions, coding_json[refname], mapped_nucleotides, mapped_codons)
#     return mapped_nucleotides, mapped_codons

def coverage_plot(coverage, reffasta, plotname, log_yscale=False):
    """make coverage plot """
    #get length of the segments from reference 
    reflength = {seg.description:len(str(seg.seq))for seg in SeqIO.parse(reffasta,"fasta")}
    #sort by length 
    sorted_ref = sorted(reflength.keys(), key=reflength.get, reverse=True)

    
    with plt.style.context("default"):
        fig = plt.figure(figsize=(8.3, 0.2*11.8))
        gs = gridspec.GridSpec(nrows=1, ncols=len(sorted_ref), width_ratios=[reflength[ref] for ref in sorted_ref], wspace=0.2)
        
        if coverage["coverage"].max() > 1000:
            max_y = 1000*np.ceil(coverage["coverage"].max()/1000) 
        else:
            max_y = 100*np.ceil(coverage["coverage"].max()/100) 

        for i, ref in enumerate(sorted_ref):
            ax = fig.add_subplot(gs[i])
            
            reflen = reflength[ref]
            df = coverage[coverage["refname"]==ref]

            # print (df)
            # for pos in df["nucpos"]:
            #     df["nucpos"] = df["nucpos"].replace([pos], pos)
            # print (df)
             
            #df = df[(df["nucpos"]>0)&(df["nucpos"]<=coord_end)]

            #what is actually plotted
            ax.plot(df["nucpos"], df["coverage"], lw=1, c="k", zorder=100)
            
            #plot setting and cosmetics
            ax.axhline(y=50, lw=0.5, ls="-.", c="indianred", zorder=0)
            ax.axhline(y=100, lw=1, ls="--", c="indianred", zorder=0)
            if log_yscale:
                ax.set_yscale("log")
                ax.set_ylim(10, max_y)
            else:
                ax.set_yticks(np.linspace(0, max_y, 6))
                ax.set_ylim(0, max_y)

            ax.grid(axis='y', ls="--", lw=0.5, c="#929292", zorder=-100)
            for axis in ["left", "right"]:
                ax.spines[axis].set_visible(False)
            if i > 0:
                ax.set_yticklabels([])
            if reflen < 1000:
                ax.set_xticks([1, reflen])
            elif reflen < 2000:
                ax.set_xticks(np.linspace(1, reflen, 3))
            else:
                ax.set_xticks(np.linspace(1, reflen, 4))
            ax.tick_params(axis='x', labelsize=7, rotation=30)
            ax.tick_params(axis='y', labelsize=9)
            ax.set_title(ref, fontsize=9)
        plt.savefig(plotname,bbox_inches="tight")

def assess_variant(var_freq, depth, bqt, bet):
    """assess if called variant is a result of a base calling error """
    pe = 10**(-bqt/10)
    pErr = 1 - binom.cdf(var_freq - 1, depth, pe)

    if (np.isnan(pErr)):
        raise Exception(f'pErr = nan; check freq/coverage/q = {var_freq}, {depth}, {bqt}')
    
    if pErr < bet:
        return True
    else:
        return False

def call_snvs_and_variants(coverage, codon_coverage, coding_json, segment_nuc_cds, min_coverage, bqt, bet):
    coverage = coverage[coverage["coverage"]>= min_coverage]

    variants = []
    for i, row in coverage.iterrows():
        for b in all_bases[:-2]: #only A, T, G, C 
            s, p, c,  = row["refname"], row["nucpos"], row["coverage"]
            rb = segment_nuc_cds[s][p-1] #base at position in refenence
            if b!= rb: #if there is an actual mutation going on
                if row[b] > 0 and assess_variant(row[b], row["coverage"], bqt, bet):
                    #determine substitution
                    snv = f"{rb}{p}{b}"

                    #determine whether substitution is synonymous or non synonymous
                    for protein,d in coding_json[s].items():
                        cds, aas = d["nuc_cds_pos"], d["pro_aa_pos"]
                
                        
                        if p in cds:
                            cp = [cds[i] for i, aa in enumerate(aas) if aa==aas[cds.index(p)]].index(p) +1
                            #cp =3 if p%3==0 else 2 if p-1%3==0 else 1#position within codon
                            aap = aas[cds.index(p)]
                            pcodons = [codon for codon in codon_table.keys() if codon[cp-1]==b]
                            codons = codon_coverage.set_index(["refname", "aapos", "protein"]).loc[(s,aap,protein), pcodons].to_dict()
                            ccov = codon_coverage.set_index(["refname", "aapos", "protein"]).loc[(s,aap,protein), "codon_coverage"]

                            for codon, codon_count in codons.items():
                                if pd.isna(codon_count):
                                    continue
                                if assess_variant(codon_count, ccov, bqt, bet) or codon_count> min_coverage:
                                    codon_idx = [j for j, a in enumerate(aas) if a == aas[cds.index(p)]]
                                    rcodon = "".join(segment_nuc_cds[s][cds[j]-1] for j in codon_idx)

                                    try:
                                        raa = codon_table[str(rcodon).upper()]
                                    except:
                                        continue #no clear ref so can't call variants
                                    try:
                                        aa = codon_table[str(codon).upper()]
                                    except:
                                        aa = "X"
                                
                                    aamut = f"{raa}{aap}{aa}"

                                    ns = True if raa!=aa else False

                                    variants.append([s,protein, p, c, b, row[b], row[b]/c, snv, aap, ccov, codon,aa, codon_count, codon_count/ccov, aamut, ns])

    variants = pd.DataFrame.from_records(variants, columns = ["refname", "protein", "nucpos", "coverage", "read_base", "nuc_count", "nuc_freq", "snv", "aapos", 
                                                                "codon_coverage", "read_codon", "read_aa", "codon_count", "codon_freq", "aa_mut", "non_synonymous"])

    #get aa count and freq
    #remove duplicate codon entries > two snv can lead to the same codon
    df = variants[["protein", "aapos", "read_codon"]]
    df = df.drop_duplicates()
    df = variants.loc[df.index.tolist(),:].groupby(["protein", "aapos", "read_aa"])[ "codon_count"].sum().to_frame()

   
    #remove duplicate codon entries > two snvs leading to the same codon
    df2 = variants.set_index(["refname", "aapos", "read_codon"])[["protein", "codon_coverage", "read_aa", "codon_count", "codon_freq", "aa_mut"]].drop_duplicates().reset_index()
    df = df.join(df2.set_index(["protein", "aapos", "read_aa"])["codon_coverage"].to_frame())

    
    df = df.rename(columns={"codon_count":"aa_count", "codon_coverage":"aa_coverage"})
    df["aa_freq"] = df["aa_count"]/df["aa_coverage"]
    variants = variants.set_index(["protein", "aapos", "read_aa"]).join(df).reset_index()
    # variants = variants[["refname", "protein", "nucpos", "coverage", "read_base", "nuc_count", "nuc_freq", "snv", "aapos", 
    #                      "codon_coverage", "read_codon", "read_aa", "codon_count", "codon_freq", "aa_mut", "non_synonymous"]] .sort_values(["refname", "nucpos"])

    return variants.sort_values(["refname", "nucpos"])

def neg_log_likelihood(p, observations, trials):
        total_log_likelihood = 0
        for i, x in enumerate(observations):
            n = trials[i]
            total_log_likelihood += binom.logpmf(x, n, p)
        return -total_log_likelihood

def estimate_binomial_mean(count, depth):
    # calculate maximum likelihood estimate using Nelder-Mead for noisy measurements
    x0 = st.mean([count[i]/x for i, x in enumerate(depth)]) #initial guess 
    result = minimize(neg_log_likelihood, x0=x0, args=(count, depth), bounds=[(0, 1)], method="Nelder-Mead")
    if not result.success:
        raise Exception(f'optimizer was not successfull')
    p_mle = result.x[0]
    return p_mle

def map_reads_per_segment(sample_reads, segment, coding_json, bqt=30):
    mapped_reads = {}
    for read_id, reads in sample_reads.items():
        #both reads must be present
        if set(sample_reads[read_id].keys()) != set(['R1', 'R2']):
            continue

        refname = reads["R1"]["refname"]
        if refname != segment:
            continue

        #get info per read pair id 
        r1_positions = np.arange(reads["R1"]["start_pos"], reads["R1"]["start_pos"]+len(reads["R1"]["read_seq"])) #indices
        r1_sequence = np.array(list(reads["R1"]["read_seq"]))
        r1_bq = np.array([ord(t)-33 for t in list(reads["R1"]["base_qual"])])
        
        r2_positions = np.arange(reads["R2"]["start_pos"], reads["R2"]["start_pos"]+len(reads["R2"]["read_seq"])) #indices
        r2_sequence = np.array(list(reads["R2"]["read_seq"]))
        r2_bq = np.array([ord(t)-33 for t in list(reads["R2"]["base_qual"])])

        #determine overlap
        overlap = np.intersect1d(r1_positions, r2_positions)
        if len(overlap) > 0:

            # get overlapping regions for read pair
            r1_overlap_mask = np.isin(r1_positions, overlap)
            r2_overlap_mask = np.isin(r2_positions, overlap)

            r1_overlap_read = r1_sequence[r1_overlap_mask]
            r2_overlap_read = r2_sequence[r2_overlap_mask]

            if (r1_overlap_read==r2_overlap_read).all() != True:
                r1_overlap_positions = r1_positions[r1_overlap_mask]
                r2_overlap_positions = r2_positions[r2_overlap_mask]
                
                r1_overlap_bq = r1_bq[r1_overlap_mask]
                r2_overlap_bq = r2_bq[r2_overlap_mask]

                mask = r1_overlap_read != r2_overlap_read #get position where read bases are different
                for i, (r1_base, r2_base) in enumerate(zip(r1_overlap_read[mask], r2_overlap_read[mask])):
                    
                    if r1_base == "N" and re.search("[ATGC]", r2_base): 
                        r1_sequence[np.where(r1_positions==r1_overlap_positions[mask][i])] = r2_sequence[np.where(r2_positions==r2_overlap_positions[mask][i])]
                    elif re.search("[ATGC]", r1_base) and re.search("[ATGC]", r2_base):
                        if r1_overlap_bq[mask][i] < r2_overlap_bq[mask][i]: #better bq read_2
                            r1_sequence[np.where(r1_positions==r1_overlap_positions[mask][i])] = r2_sequence[np.where(r2_positions==r2_overlap_positions[mask][i])]
                        elif r1_overlap_bq[mask][i] == r2_overlap_bq[mask][i] and r2_base!=r1_base: #ambiguous
                            r1_sequence[np.where(r1_positions==r1_overlap_positions[mask][i])] = "N"
                    elif r1_base == "-" and r1_base != r2_base:
                        r1_sequence[np.where(r1_positions==r1_overlap_positions[mask][i])] = r2_sequence[np.where(r2_positions==r2_overlap_positions[mask][i])]

            # remove overlapped section in r2
            r2_positions = r2_positions[~r2_overlap_mask]
            r2_sequence = r2_sequence[~r2_overlap_mask]

        positions = dict(map(lambda i,j : (i+1,j) , r1_positions,r1_sequence))
        positions.update(dict(map(lambda i,j : (i+1,j) , r2_positions,r2_sequence)))
        positions = dict(OrderedDict(sorted(positions.items())))
        
        #get base quality
        bq = dict(map(lambda i, j: (i+1, j), r1_positions, r1_bq))
        bq.update(dict(map(lambda i, j: (i+1, j), r2_positions, r2_bq)))
        bq = dict(OrderedDict(sorted(bq.items())))

        #remove bases below the quality threshold 
        for p,  q in bq.items():
            if q < bqt:
                positions[p] = "N"

        mapped_reads[read_id] = positions
        
    mapped_reads= pd.DataFrame.from_dict(mapped_reads, orient="index")    
    cds = list(set([i for p in coding_json[segment].keys() for i in coding_json[segment][p]["nuc_cds_pos"]]))
    mapped_reads = mapped_reads[cds].drop_duplicates()
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
        positions = [nucpos for nucpos, allele in row.to_dict().items() if nucpos in possible_alleles.keys() and allele in possible_alleles[nucpos]]
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
    all_possible_haplotypes = list(itertools.product(*[list(nucpos_allele_count[nucpos].keys()) for nucpos in sorted_nucpos]))

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
        sample_freq, sample_haploset, _, _, _ = reconstruct_haplotype(partial_haplotypes, nuc_coverage,  max_iter, save_state, n_searches, overdispersion, sigma, verbose)


        # append sample_freq to distribution
        for j, haplotype in enumerate(ML_haploset):
            if haplotype in sample_haploset:
                ML_haplotype_freqdist[(i, j)] = sample_freq[sample_haploset.index(haplotype)]
        
    return ML_haplotype_freqdist
        