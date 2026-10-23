import os, re, gzip, json
import pandas as pd, numpy as np, multiprocessing as mp
import snakemake, ast

from utils import *

# Prepare inputs based on the condition
additional_variants = expand(f"{config['output']}/snvs/variant_frequencies_{{sample}}.csv", sample=config["samples"]) if len(config['treps']) > 1 else []
rule all:
    input:
        summary = f"{config['output']}/trimmomatic_stats.csv",
        posttrim_report = f"{config['output']}/posttrim_fastqc_results.csv",
        postmap_report = f"{config['output']}/postmap_fastqc_results.csv",
        postdeduplication_report = f"{config['output']}/postdeduplication_fastqc_results.csv",
        coverage_plot = expand(f"{config['output']}/plots/covplot_{{sample_trep}}.pdf", sample_trep=config["sample_treps"]),
        variants = expand(f"{config['output']}/snvs/variants_{{sample_trep}}.csv", sample_trep=config["sample_treps"]),
        estimated_variants = additional_variants,

##### RULES ########

rule RenameSamples:
    input:
        metadata = config["metadata"],
        reads = expand(f"{config['input']}/{{original_name}}_{{read}}.fastq.gz", original_name=config['original_names'], read=config['reads'])
    output:
        reads = expand(f"{config['input']}/{{long_sample}}_{{read}}.fastq.gz", long_sample=config['long_samples'], read=config['reads'])
    message: "renaming FASTQ files"
    threads: 1
    run:
        metadata = pd.read_excel(input.metadata, engine="openpyxl")
        metadata["sample"] = [f"P{int(row['Passage'])}-{row['Donor']}{row['Replicate']}"  if row["Donor"]!='stock' else row['Donor'] for i, row in metadata.iterrows()]
        metadata["sample_trep"] = [f"P{int(row['Passage'])}-{row['Donor']}{row['Replicate']}_{row['Technical replicate']}"  if row["Donor"]!='stock' else f"{row['Donor']}_{row['Technical replicate']}" for i, row in metadata.iterrows()]
        metadata["long_sample"] = [f"{row['GS Flow cell ID']}_P{int(row['Passage'])}-{row['Donor']}{row['Replicate']}_{row['Technical replicate']}" if row["Donor"]!='stock' else f"{row['GS Flow cell ID']}_{row['Donor']}_{row['Technical replicate']}"  for i, row in metadata.iterrows()]

        new_names = {row["GS full name"]:row["long_sample"] for i, row in metadata.iterrows()}

        for f in input.reads:
            fn = "_".join(f.split("/")[-1].split('.')[0].split("_")[:-1])
            #already renamed
            if fn in new_names.values():
                continue
            
            shell(f"mv {f} {f.replace(fn, new_names[fn])}")       

rule FASTQC:
    input:
        read= f"{config['input']}/{{long_sample}}_{{read}}.fastq.gz",
    params:
        fastqcdir = f"{config['output']}/fastqc",
    output:
        html = f"{config['output']}/fastqc/{{long_sample}}_{{read}}_fastqc.html",
        fastqc_zip = f"{config['output']}/fastqc/{{long_sample}}_{{read}}_fastqc.zip",

    message: "performing FASTQC analysis for {input.read}"
    threads: 1
    run:
        cmd = ['gunzip', '-c', input.read, '|', 'fastqc',f'stdin:{{wildcards.long_sample}}_{{wildcards.read}}',
                    f'--outdir={params.fastqcdir}']
        shell(" ".join(cmd))  

rule AnalyzeFASTQC:
    input:
        fastqc_zips = expand(f"{config['output']}/fastqc/{{long_sample}}_{{read}}_fastqc.zip", long_sample=config['long_samples'], read=config['reads'])
    params:
        output = config["output"],
        fastqcdir = f"{config['output']}/fastqc",
        bqt = config['bqt'] 
    output:
        report = f"{config['output']}/pretrim_fastqc_results.csv"
    message: "summarizing pretrim FASTQC results"
    threads:1
    run:
        adapters_to_rm = {}
        fastqc_results = []

        for zipfile in input.fastqc_zips:
            long_sample = "_".join(zipfile.split("/")[-1].split("_")[:-2])
            read = zipfile.split("_")[-2]
            
            report = zipfile.replace(".zip", "/fastqc_data.txt")
            #unzip only if not unzipped before
            if not os.path.isfile(report):
                shell(f"unzip {zipfile} -d {params.fastqcdir}")
            
            #get results
            percent_abv_qualthres, ac, start_pos, end_pos = analyze_fastqc(report, params.bqt)

            # add results 
            for k, v in ac.items():
                if k in adapters_to_rm and v > adapters_to_rm[k]:
                    adapters_to_rm[k] = v
                elif v > 0.:
                    adapters_to_rm[k] = v

            fastqc_results.append({"sampid":long_sample, "read":read, "percent_abv_qualthres":percent_abv_qualthres,
                                           "start_pos":start_pos, "end_pos":end_pos})

            #remove zip 
            shell(f'rm -r -d {zipfile.replace(".zip", "")}')

        if len(adapters_to_rm) > 0:
            print ('\n#-- Presence of adapter sequence (max. proportion of reads) --#')
            for k, v in adapters_to_rm.items():
                print (f"{k}: {round(v,2)}")

        fastqc_results = pd.DataFrame.from_dict(fastqc_results)
        fastqc_results.to_csv(f"{params.output}/pretrim_fastqc_results.csv", index=False)
    
rule TrimReads:
    input:
        fastq_report = rules.AnalyzeFASTQC.output.report,
        r1 = f"{config['input']}/{{long_sample}}_R1.fastq.gz",
        r2 = f"{config['input']}/{{long_sample}}_R2.fastq.gz"
    params:
        trimdir = f"{config['output']}/trimmed",
        adapter = f"{config['datadir']}/adapters/{config['adapter']}.fa"
    output:
        r1_paired = f"{config['output']}/trimmed/{{long_sample}}_R1.paired.fastq.gz",
        r1_unpaired = f"{config['output']}/trimmed/{{long_sample}}_R1.unpaired.fastq.gz",
        r2_paired = f"{config['output']}/trimmed/{{long_sample}}_R2.paired.fastq.gz",
        r2_unpaired = f"{config['output']}/trimmed/{{long_sample}}_R2.unpaired.fastq.gz",
        trimlog = temp(f"{config['output']}/trimmed/{{long_sample}}_log.txt") #deleting these files later on    
    threads: workflow.cores if workflow.cores < 4 else 4
    message: "Trimming PE reads for {wildcards.long_sample}"
    run:
        fastqcdf = pd.read_csv(input.fastq_report)
        sampiddf = fastqcdf[fastqcdf["sampid"]==wildcards.long_sample]
        
        cmd = ["trimmomatic", "PE", "-phred33", "-threads", str(threads), "-summary", output.trimlog, input.r1, input.r2,
                output.r1_paired, output.r1_unpaired, output.r2_paired, output.r2_unpaired]
        if params.adapter != "":
            cmd += [f"ILLUMINACLIP:{params.adapter}:2:3:10"]
        cmd += ["LEADING:3", "TRAILING:3", "MAXINFO:30:0.4",f"CROP:{int(sampiddf['end_pos'].max())}", "MINLEN:30"]

        shell(" ".join(cmd))

rule SummarizeTrim:
    input:
        logs = expand(f"{config['output']}/trimmed/{{long_sample}}_log.txt",long_sample=config["long_samples"])
    output:
        summary = f"{config['output']}/trimmomatic_stats.csv"
    message: "Summarizing trimmomatic results"
    threads:1
    run:
        trimmed_df = []
        for log in input.logs:
            sampid = log.split("/")[-1].rstrip("_log.txt") #cant assess wildcards so this will have to do
            with open(log, "r") as f:
                results = [line.rstrip("\n").split(" ")[-1] for line in f]
            trimmed_df.append({"sampid":sampid, "total_pairs":int(results[0]),"both":int(results[1]),
                                "forward":int(results[3]),"reverse":int(results[5]),"dropped":int(results[7])})


        trimmed_df = pd.DataFrame.from_dict(trimmed_df)
        trimmed_df["both_prop"] = trimmed_df["both"]/trimmed_df["total_pairs"]
        
        trimmed_df = trimmed_df[['sampid', 'total_pairs', 'both', 'both_prop', 'forward', 'reverse', 'dropped']]
        print (f"Overall mean prop of pairs retained: {round(np.mean(trimmed_df['both_prop']),3)} (SD: {round(np.std(trimmed_df['both_prop'], ddof=1),3)})")

        trimmed_df.to_csv(output.summary, index=False)

rule PostTrimFASTQC:
    input:
        read= f"{config['output']}/trimmed/{{long_sample}}_{{read}}.paired.fastq.gz",
    params:
        fastqcdir = f"{config['output']}/fastqc_posttrim",
    output:
        html = f"{config['output']}/fastqc_posttrim/{{long_sample}}_{{read}}_fastqc.html",
        fastqc_zip = f"{config['output']}/fastqc_posttrim/{{long_sample}}_{{read}}_fastqc.zip",
    message: "performing FASTQC analysis for {input.read}"
    threads: 1
    run: 
        cmd = cmd = ['gunzip', '-c', input.read, '|', 'fastqc',f'stdin:{{wildcards.long_sample}}_{{wildcards.read}}', f'--outdir={params.fastqcdir}']
        shell(" ".join(cmd))  

rule AnalyzePostTrimFASTQC:
    input:
        fastqc_zips = expand(f"{config['output']}/fastqc_posttrim/{{long_sample}}_{{read}}_fastqc.zip", long_sample=config['long_samples'], read=config['reads'])
    params:
        output = config["output"],
        fastqcdir = f"{config['output']}/fastqc_posttrim",
        bqt = config['bqt'] 
    output:
        report = f"{config['output']}/posttrim_fastqc_results.csv"
    message: "summarizing post trim FASTQC results"
    threads:1
    run:
        adapters_to_rm = {}
        fastqc_results = []

        for zipfile in input.fastqc_zips:
            long_sample = "_".join(zipfile.split("/")[-1].split("_")[:-2])
            read = zipfile.split("_")[-2]
            
            report = zipfile.replace(".zip", "/fastqc_data.txt")
            #unzip only if not unzipped before
            if not os.path.isfile(report):
                shell(f"unzip {zipfile} -d {params.fastqcdir}")
            
            #get results
            percent_abv_qualthres, ac, start_pos, end_pos = analyze_fastqc(report, params.bqt)

            #add results 
            for k, v in ac.items():
                if k in adapters_to_rm and v > adapters_to_rm[k]:
                    adapters_to_rm[k] = v
                elif v > 0.:
                    adapters_to_rm[k] = v

            fastqc_results.append({"sampid":long_sample, "read":read, "percent_abv_qualthres":percent_abv_qualthres,
                                           "start_pos":start_pos, "end_pos":end_pos})

            #remove zip 
            shell(f'rm -r -d {zipfile.replace(".zip", "")}')

        if len(adapters_to_rm) > 0:
            print ('\n#-- Presence of adapter sequence (max. proportion of reads) --#')
            for k, v in adapters_to_rm.items():
                print (f"{k}: {round(v,2)}")

        fastqc_results = pd.DataFrame.from_dict(fastqc_results)
        fastqc_results.to_csv(f"{params.output}/posttrim_fastqc_results.csv", index=False)
        
rule MapReads:
    input:
        r1 = rules.TrimReads.output.r1_paired,
        r2 = rules.TrimReads.output.r2_paired 
    params:
        mapdir = f"{config['output']}/align",
        ref = f"{config['refdir']}/ref_NL08.fasta",
        #mgl = 2300 #max gene length
    output:
        samzip = f"{config['output']}/align/{{long_sample}}.sam.gz"
    threads: workflow.cores if workflow.cores < 10 else 10
    resources:
        memory = (workflow.cores*4)*1024
    message: "Mapping reads for {wildcards.long_sample}"
    run:

        # cmd = ['bowtie2', '-x', params.ref, '-X', str(params.mgl), '-k', '2', '--very-sensitive-local', '-p', str(threads),
        #         '-1', input.r1, '-2', input.r2, '-S', f"{params.mapdir}/{wildcards.long_sample}.sam" ]
        cmd = ['bowtie2', '-x', params.ref, '-X', str(params.mgl), '--very-sensitive-local', '-p', str(threads),
                '-1', input.r1, '-2', input.r2, '-S', f"{params.mapdir}/{wildcards.long_sample}.sam" ]
        shell(" ".join(cmd))

        cmd = ['gzip', f"{params.mapdir}/{wildcards.long_sample}.sam"]
        shell(" ".join(cmd))

rule Sam2Bam:
    input:
        sam = rules.MapReads.output.samzip
    params:
        min_length = config['min_length'], #min length of reads
        min_map_q = config['min_mapq'], #minimum map quality should be at least 20
    output:
        bam =  temp(f"{config['output']}/align/{{long_sample}}.bam")
    message: "converting SAM file to BAM file for {wildcards.long_sample} and remove too short reads"
    threads:1 
    shell:
        """
        samtools view -e 'length(seq)>{params.min_length}' -q {params.min_map_q} -f 0x2 -S -b {input.sam} > {output.bam} 
        """

rule PostMapFASTQC:
    input:
        bam = rules.Sam2Bam.output.bam,
    params:
        fastqcdi = f"{config['output']}/fastqc_postmap",
    output:
        html = f"{config['output']}/fastqc_postmap/{{long_sample}}_fastqc.html",
        fastqc_zip = f"{config['output']}/fastqc_postmap/{{long_sample}}_fastqc.zip",
    message: "performing FASTQC analysis for {input.bam}"
    threads: 1
    run: 
        cmd = cmd = ['fastqc', input.bam , f'--outdir={params.fastqcdir}']
        shell(" ".join(cmd))  

rule AnalyzePostMapFASTQC:
    input:
        fastqc_zips = expand(f"{config['output']}/fastqc_postmap/{{long_sample}}_fastqc.zip", long_sample=config['long_samples'])
    params:
        output = config["output"],
        fastqcdir = f"{config['output']}/fastqc_postmap",
        bqt = config['bqt'] 
    output:
        report = f"{config['output']}/postmap_fastqc_results.csv"
    message: "summarizing post mapping FASTQC results"
    threads:1
    run:
        adapters_to_rm = {}
        fastqc_results = []

        for zipfile in input.fastqc_zips:
            long_sample = "_".join(zipfile.split("/")[-1].split("_")[:-2])
            read = zipfile.split("_")[-2]
            
            report = zipfile.replace(".zip", "/fastqc_data.txt")
            #unzip only if not unzipped before
            if not os.path.isfile(report):
                shell(f"unzip {zipfile} -d {params.fastqcdir}")
            
            #get results
            percent_abv_qualthres, ac, start_pos, end_pos = analyze_fastqc(report, params.bqt)

            #add results 
            for k, v in ac.items():
                if k in adapters_to_rm and v > adapters_to_rm[k]:
                    adapters_to_rm[k] = v
                elif v > 0.:
                    adapters_to_rm[k] = v

            fastqc_results.append({"sampid":long_sample, "read":read, "percent_abv_qualthres":percent_abv_qualthres,
                                           "start_pos":start_pos, "end_pos":end_pos})

            #remove zip 
            shell(f'rm -r -d {zipfile.replace(".zip", "")}')

        if len(adapters_to_rm) > 0:
            print ('\n#-- Presence of adapter sequence (max. proportion of reads) --#')
            for k, v in adapters_to_rm.items():
                print (f"{k}: {round(v,2)}")

        fastqc_results = pd.DataFrame.from_dict(fastqc_results)
        fastqc_results.to_csv(f"{params.output}/postmap_fastqc_results.csv", index=False)

rule Deduplicate:
    input:
        bam = rules.Sam2Bam.output.bam
    output:
        sortednames = temp(f"{config['output']}/sorted/{{long_sample}}.names.bam"),
        fixmate = temp(f"{config['output']}/mates/{{long_sample}}.bam"),
        sorted_bam = temp(f"{config['output']}/sorted-aligned/{{long_sample}}.bam"),
        bam = f"{config['output']}/deduplicated/{{long_sample}}.bam",
        #bam = temp(f"{config['output']}/deduplicated/{{long_sample}}.bam"),
    message: "removing PCR duplicates for {wildcards.long_sample}"
    threads: 1
    shell:
        """
        samtools sort -n -o {output.sortednames} -O BAM {input.bam}
        samtools fixmate -m {output.sortednames} {output.fixmate}
        samtools sort -o {output.sorted_bam} {output.fixmate}
        samtools markdup -r -s {output.sorted_bam} {output.bam}
        """

rule PostDeduplicationFASTQC:
    input:
        bam = rules.Deduplicate.output.bam
    params:
        fastqcdir = f"{config['output']}/fastqc_postdeduplication",
    output:
        html = f"{config['output']}/fastqc_postdeduplication/{{long_sample}}_fastqc.html",
        fastqc_zip = f"{config['output']}/fastqc_postdeduplication/{{long_sample}}_fastqc.zip",

    message: "performing FASTQC analysis for {input.bam}"
    threads: 1
    run: 
        cmd = cmd = ['fastqc', input.bam , f'--outdir={params.fastqcdir}']
        shell(" ".join(cmd))  

rule AnalyzePostDeduplicationFASTQC:
    input:
        fastqc_zips = expand(f"{config['output']}/fastqc_postdeduplication/{{long_sample}}_fastqc.zip", long_sample=config['long_samples'])
    params:
        output = config["output"],
        fastqcdir = f"{config['output']}/fastqc_postdeduplication",
        bqt = config['bqt'] 
    output:
        report = f"{config['output']}/postdeduplication_fastqc_results.csv"
    message: "summarizing post deduplication FASTQC results"
    threads:1
    run:
        adapters_to_rm = {}
        fastqc_results = []

        for zipfile in input.fastqc_zips:
            long_sample = "_".join(zipfile.split("/")[-1].split("_")[:-2])
            read = zipfile.split("_")[-2]
            
            report = zipfile.replace(".zip", "/fastqc_data.txt")
            #unzip only if not unzipped before
            if not os.path.isfile(report):
                shell(f"unzip {zipfile} -d {params.fastqcdir}")
            
            #get results
            percent_abv_qualthres, ac, start_pos, end_pos = analyze_fastqc(report, params.bqt)

            #add results 
            for k, v in ac.items():
                if k in adapters_to_rm and v > adapters_to_rm[k]:
                    adapters_to_rm[k] = v
                elif v > 0.:
                    adapters_to_rm[k] = v

            fastqc_results.append({"sampid":long_sample, "read":read, "percent_abv_qualthres":percent_abv_qualthres,
                                           "start_pos":start_pos, "end_pos":end_pos})

            #remove zip 
            shell(f'rm -r -d {zipfile.replace(".zip", "")}')

        if len(adapters_to_rm) > 0:
            print ('\n#-- Presence of adapter sequence (max. proportion of reads) --#')
            for k, v in adapters_to_rm.items():
                print (f"{k}: {round(v,2)}")

        fastqc_results = pd.DataFrame.from_dict(fastqc_results)
        fastqc_results.to_csv(f"{params.output}/postdeduplication_fastqc_results.csv", index=False)


rule MiniMap:
    input:
        r1 = rules.TrimReads.output.r1_paired,
        r2 = rules.TrimReads.output.r2_paired 
    params:
        minimap_dir =  f"{config['output']}/minimap",
        ref = f"{config['refdir']}/minimap_references.fa",
    output:
        minimap = temp(f"{config['output']}/minimap/raw_{{long_sample}}.sam.gz")
    threads: workflow.cores if workflow.cores < 4 else 4
    run:
        cmd = ["minimap2", "-ax", "sr", "-t", str(threads), "-c", params.ref, input.r1, input.r2, ">", f"{params.minimap_dir}/raw_{wildcards.long_sample}.sam"]
        shell(" ".join(cmd))

        cmd = ["gzip",f"{params.minimap_dir}/raw_{wildcards.long_sample}.sam" ]
        shell(" ".join(cmd))

rule MiniMapDeduplication:
    input:
        raw_minimap = rules.MiniMap.output.minimap,
    params:
        minimap_dir =  f"{config['output']}/minimap",
        min_length = config['min_length'], #min length of reads
        min_map_q = config['min_mapq'], 
    output:
        bam = temp(f"{config['output']}/minimap/{{long_sample}}.bam"),
        sortednames = temp(f"{config['output']}/minimap/{{long_sample}}.names.bam"),
        fixmate = temp(f"{config['output']}/minimap/fixmate_{{long_sample}}.bam"),
        sorted_bam = temp(f"{config['output']}/minimap/sorted_{{long_sample}}.bam"),
        dp_bam = temp(f"{config['output']}/minimap/deduplicated_{{long_sample}}.bam"),
        minimap = temp(f"{config['output']}/minimap/{{long_sample}}.sam.gz"),
    shell:
        """
        samtools view -e 'length(seq)>{params.min_length}' -q {params.min_map_q} -f 0x2 -S -b {input.raw_minimap} > {output.bam}
        samtools sort -n -o {output.sortednames} -O BAM {output.bam}
        samtools fixmate -m {output.sortednames} {output.fixmate}
        samtools sort -o {output.sorted_bam} {output.fixmate}
        samtools markdup -r -s {output.sorted_bam} {output.dp_bam}
        samtools view -h -o  {params.minimap_dir}/{wildcards.long_sample}.sam {output.dp_bam}
        gzip {params.minimap_dir}/{wildcards.long_sample}.sam
        """

rule MiniMapSummary:
    input:
        minimap = rules.MiniMapDeduplication.output.minimap,
    params:
        ref = f"{config['refdir']}/minimap_references.fa",
        h3_strains = ["stock|A/H3N2"],
        bqt = config["bqt"],
    output:
        summary = f"{config['output']}/minimap/{{long_sample}}.csv",
        non_target_reads = f"{config['output']}/minimap/nonH3reads_{{long_sample}}.fasta",
    threads:1
    run:
        #get ref
        minimap_ref = []
        segments = ['PB2','PB1','PA','HA','NP','NA','MP','NS']
        for r in SeqIO.parse(params.ref, "fasta"):
            strain = "|".join(r.id.split("|")[2:4]).replace("_", "")
            if strain not in minimap_ref:
                minimap_ref.append(strain)
        
        summary, non_target_reads = summarize_minimap(input.minimap, minimap_ref, segments, params.h3_strains, params.bqt)

        #write output 
        summary.to_csv(output.summary, index=False)
        
        with open(output.non_target_reads, "w") as fw:
            for rid, rseq in non_target_reads.items():
                fw.write(f">{rid}\n{rseq}\n")

rule MiniMapSummaryPerSample:
    input:
        summary_files = expand(f"{config['output']}/minimap/{{long_sample}}.csv", long_sample=config["long_samples"]),
    params:
        summary_dir = f"{config['output']}/minimap",
    output:
        summary = f"{config['output']}/minimap/{{sample_trep}}.csv"
    threads: 1
    run:
        minimap_summary = []
        for f in input.summary_files:
            if wildcards.sample_trep in f:

                res = pd.read_csv(f).set_index(["reference strain"])
                
                if isinstance(minimap_summary,list):
                    minimap_summary = res

                else:
                    for strain, row in res.iterrows():
                        for c in res.columns:
                            minimap_summary.loc[strain,c]  += row[c]

        minimap_summary.reset_index().to_csv(output.summary, index=False)

rule ParseBam:
    input:
        bam = rules.Deduplicate.output.bam,
        non_target_reads = rules.MiniMapSummary.output.non_target_reads
    params:
        ref =  f"{config['refdir']}/ref_NL08.fasta",
        sortdir = f"{config['output']}/sorted-aligned",
        bqt = config["bqt"], # base qual threshold
        min_map_q = config['min_mapq'], #extra check
    output:
        json = f"{config['output']}/sorted-aligned/{{long_sample}}.json.gz"
    threads: 1
    resources:
        memory = 4*1024
    message: "Parsing BAm file of {wildcards.long_sample} for proper concordant reads....."
    run:
        #get ref_pos_mapping
        ref_pos_mapping = get_ref_pos_mapping(params.ref)
        refname_dict = {f"#{i}":k for i,k in enumerate(list(ref_pos_mapping.keys()))}

        #get non target read ids 
        non_target_reads = [r.id.split("|")[0] for r in SeqIO.parse(input.non_target_reads, "fasta")]


        map_json = parse_bam(input.bam, ref_pos_mapping, refname_dict, params.bqt, params.min_map_q, non_target_reads)

        with gzip.open(output.json, 'w') as output:
            json_str = json.dumps(map_json, cls=NpEncoder) + "\n"
            json_bytes = json_str.encode('utf-8')
            output.write(json_bytes)    

rule ComputeCoverage: 
    input:
        json = rules.ParseBam.output.json,
        non_target_reads = rules.MiniMapSummary.output.non_target_reads

    params:
        codingjson = f"{config['datadir']}/CDS_H3N2.json",
        bqt = config['bqt'],
    output:
        nuc_coverage = f"{config['output']}/coverage/nuc_coverage_{{long_sample}}.csv",
        codon_coverage = f"{config['output']}/coverage/codon_coverage_{{long_sample}}.csv"
    message: "Computing nucleotide and codon coverage per segment per position for {wildcards.long_sample}"
    threads:1
    run:
        mapped_codons, mapped_nucleotides = {}, {}

        with open(params.codingjson) as json_file:
            coding_json = json.load(json_file)    

        with gzip.open(input.json) as json_file:
            sample_json = json.load(json_file)

        #get non target read ids 
        non_target_reads = [r.id.split("|")[0] for r in SeqIO.parse(input.non_target_reads, "fasta")]

        mapped_nucleotides, mapped_codons = map_reads(sample_json, coding_json, params.bqt, mapped_nucleotides, mapped_codons, non_target_reads)

        mapped_nucleotides = pd.DataFrame.from_dict(mapped_nucleotides, orient="index", columns=all_bases).sort_index().reset_index()
        mapped_nucleotides = mapped_nucleotides.rename(columns={"level_0":"refname","level_1":"nucpos"}).drop(columns = ["N"])
        nuc_coverage = pd.melt(mapped_nucleotides, id_vars=["refname", "nucpos"], value_name="coverage", var_name="read_base")
        nuc_coverage = nuc_coverage.groupby(['refname', "nucpos"]).sum()#.reset_index()
        mapped_nucleotides = pd.concat([mapped_nucleotides.set_index(["refname", "nucpos"]), nuc_coverage], axis=1).reset_index()

        mapped_codons = pd.DataFrame.from_dict(mapped_codons, orient="index", columns=codon_table.keys()).sort_index().reset_index()
        mapped_codons = mapped_codons.rename(columns={"level_0":"refname", "level_1":"protein", "level_2":"aapos"})
        codon_coverage = pd.melt(mapped_codons, id_vars=["refname", "protein", "aapos"], value_name="codon_coverage", var_name="codon")
        codon_coverage = codon_coverage.groupby(['refname', "protein", "aapos"]).sum()
        mapped_codons = pd.concat([mapped_codons.set_index(["refname", "protein", "aapos"]), codon_coverage], axis=1).reset_index()

        mapped_nucleotides.to_csv(output.nuc_coverage, index=False)
        mapped_codons.to_csv(output.codon_coverage, index=False)

rule CoveragePerSample: 
    input:
        nuc_coverage_files = expand(f"{config['output']}/coverage/nuc_coverage_{{long_sample}}.csv", long_sample=config["long_samples"]),
        codon_coverage_files = expand(f"{config['output']}/coverage/codon_coverage_{{long_sample}}.csv", long_sample=config["long_samples"]),
    params:
        coverage_dir = f"{config['output']}/coverage/",
        samples = config["samples"],
    output:
        nuc_coverage = f"{config['output']}/coverage/nuc_coverage_{{sample_trep}}.csv",
        codon_coverage = f"{config['output']}/coverage/codon_coverage_{{sample_trep}}.csv"
    message:  "Complining flowcell coverage results into a single file for {wildcards.sample_trep}"
    run:
        
        nuc_coverage = []
        for f in input.nuc_coverage_files:
            
            if wildcards.sample_trep in f:
                df = pd.read_csv(f, keep_default_na=False, na_values=['']).set_index(['refname', 'nucpos']).sort_index()
                if isinstance(nuc_coverage, list):
                    nuc_coverage = df
                else:
                    for i, row in df.iterrows():
                        for c in df.columns:
                            try:
                                nuc_coverage.loc[i, c] += row[c]  
                            except:
                                nuc_coverage.loc[i, c] = row[c]
                                if isinstance(nuc_coverage, pd.DataFrame):
                                    nuc_coverage = nuc_coverage.sort_index()

        codon_coverage = []
        for f in input.codon_coverage_files:
            if wildcards.sample_trep in f:
                df = pd.read_csv(f, keep_default_na=False, na_values=['']).set_index(["refname", "protein", "aapos"]).sort_index()
                if isinstance(codon_coverage, list):
                    codon_coverage = df
                else:
                    for i, row in df.iterrows():
                        if i not in codon_coverage.index.tolist():
                            codon_coverage.loc[i] = row.to_list()
                            codon_coverage = codon_coverage.sort_index()
                        
                        else:
                            for c in df.columns:
                                codon_coverage.loc[i, c] += row[c]
    
                
        nuc_coverage.sort_index().reset_index().to_csv(output.nuc_coverage, index=False)
        codon_coverage.sort_index().reset_index().to_csv(output.codon_coverage, index=False)

rule CoveragePlots:
    input:
        coverage = rules.CoveragePerSample.output.nuc_coverage
    params:
        reffasta = f"{config['refdir']}/stock_consensus_cds.fasta",
        logYscale=True,
    output:
        coverage_plot = f"{config['output']}/plots/covplot_{{sample_trep}}.pdf"
    message: "Computing coverage plot for {wildcards.sample_trep}"
    threads:1
    run: 
        coverage = pd.read_csv(input.coverage)
        coverage_plot(coverage, params.reffasta, output.coverage_plot, params.logYscale) 

rule CallSNVsVariants:
    input:
        nuc_coverage = rules.CoveragePerSample.output.nuc_coverage,
        codon_coverage = rules.CoveragePerSample.output.codon_coverage, 
    params:
        codingjson = f"{config['datadir']}/CDS_H3N2.json",
        ref =  f"{config['refdir']}/stock_consensus_cds.fasta",
        min_coverage = config['min_coverage'],
        bqt = config['bqt'], #default 30 #base quality thresholpd
        bet = config['bet'],#default = 0.001 #base error tolerance 
    output:
        variants = f"{config['output']}/snvs/variants_{{sample_trep}}.csv"
    message: "Call SNVs and variants for {wildcards.sample_trep}"
    threads: 1
    run:
        #load data
        mapped_codons = pd.read_csv(input.codon_coverage,keep_default_na=False, na_values=[''])#.set_index(["refname", "protein", "aapos"])
        mapped_nucleotides = pd.read_csv(input.nuc_coverage)#.set_index(["refname", "nucpos"])

        #prep ref info
        with open(params.codingjson) as json_file:
            coding_json = json.load(json_file)   
        #get CDS of reference 
        segment_nuc_cds = {r.id:r.seq for r in SeqIO.parse(params.ref, 'fasta')}
    
        variants = call_snvs_and_variants(mapped_nucleotides, mapped_codons, coding_json, segment_nuc_cds, params.min_coverage, params.bqt, params.bet)
        
        #write output
        variants.to_csv(output.variants, index=False)



rule EstimateVariantFrequency:
    input:
        variant_1 = f"{config['output']}/snvs/variants_{{sample}}_1.csv",
        nuc_coverage_1 = f"{config['output']}/coverage/nuc_coverage_{{sample}}_1.csv",
        codon_coverage_1= f"{config['output']}/coverage/codon_coverage_{{sample}}_1.csv",
        variant_2 = f"{config['output']}/snvs/variants_{{sample}}_2.csv",
        nuc_coverage_2 = f"{config['output']}/coverage/nuc_coverage_{{sample}}_2.csv",
        codon_coverage_2= f"{config['output']}/coverage/codon_coverage_{{sample}}_2.csv",
    output:
        variants = f"{config['output']}/snvs/variant_frequencies_{{sample}}.csv"
    message: "estimate true variant frequencies"
    threads: 1
    run:
        #get snv results
        variants_1 = pd.read_csv(input.variant_1, keep_default_na=False, na_values=['']).set_index(cc)
        variants_2 = pd.read_csv(input.variant_2,keep_default_na=False, na_values=['']).set_index(cc)
        #get nuc coverage and codon coverage 
        coverage_1 = pd.read_csv(input.nuc_coverage_1,keep_default_na=False, na_values=['']).set_index(["refname", "nucpos"])
        coverage_2 = pd.read_csv(input.nuc_coverage_2,keep_default_na=False, na_values=['']).set_index(["refname", "nucpos"])

        codon_coverage_1 = pd.read_csv(input.codon_coverage_1,keep_default_na=False, na_values=['']).set_index(["refname", "protein", "aapos"])
        codon_coverage_2 = pd.read_csv(input.codon_coverage_2,keep_default_na=False, na_values=['']).set_index(["refname", "protein", "aapos"])

        #merge dfs 
        coverage = coverage_1.join(coverage_2, how="outer", lsuffix="|tr1", rsuffix="|tr2")
        codon_coverage = codon_coverage_1.join(codon_coverage_2, how="outer", lsuffix="|tr1", rsuffix="|tr2")
        variants = variants_1.join(variants_2, how="outer", lsuffix="|tr1", rsuffix="|tr2")

        #fill missing values for variants
        for (s, prot, npos, aapos, rb, snv, read_codon, read_aa, aa_mut, ns), row in variants.iterrows():
            for c in variants.columns:
                if not pd.isna(row[c]):
                    continue
                tr = c.split("|")[-1]
                if c.startswith("coverage"):
                    variants.loc[(s, prot, npos, aapos, rb, snv, read_codon, read_aa, aa_mut, ns), c] = coverage.loc[(s, npos), c]
                elif c.startswith("codon_coverage") or c.startswith("aa_coverage"):
                    #codon and aa coverage are the same
                    variants.loc[(s, prot, npos, aapos, rb, snv, read_codon, read_aa, aa_mut, ns), c] = codon_coverage.loc[(s, prot, aapos), f"codon_coverage|{tr}"]
                else: 
                    variants.loc[(s, prot, npos, aapos, rb, snv, read_codon, read_aa, aa_mut, ns), c] = 0
        
        #reset index
        coverage = coverage.reset_index()
        codon_coverage = codon_coverage.reset_index()           
        variants = variants.reset_index().drop_duplicates()

        #get technical replicate columns
        trc = [c for c in variants.columns if "|tr" in c]

        #using maximum likelihood estimation of the binomial proportions estimate the nucleotide, codon, and amino acid frequency
        for i, row in variants.iterrows():
            #nuc frequency 
            nuc_count = row[[c for c in trc if "nuc_count" in c]].to_list()
            nuc_depth = row[[c for c in trc if c.startswith("coverage")]].to_list()
            variants.loc[i, "nuc_freq"] = estimate_binomial_mean([int(nc) for nc in nuc_count], [int(nd) for nd in nuc_depth])

            #codon frequency
            codon_count = row[[c for c in trc if "codon_count" in c]].to_list()
            codon_depth = row[[c for c in trc if c.startswith("codon_coverage")]].to_list()
            variants.loc[i, "codon_freq"] = estimate_binomial_mean([int(cc) for cc in codon_count], [int(cd) for cd in codon_depth])

            #aa frequency
            aa_count = row[[c for c in trc if "aa_count" in c]].to_list()
            #same of as codon 
            aa_depth = row[[c for c in trc if c.startswith("aa_coverage")]].to_list()
            variants.loc[i, "aa_freq"] = estimate_binomial_mean([float(ac) for ac in aa_count],[float(ad) for ad in aa_depth])

        variants.to_csv(output.variants, index=False)



            


    



