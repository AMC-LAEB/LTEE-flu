# LTEE-flu

## Long-term experimental evolution of seasonal influenza A/H3N2 virus reveals diverse evolutionary potential

Codes to reproduce De Rooij, et al.

## Requirements

All codes were tested/ran on macOS Tahoe 26.5 (25F71) and Red Hat Enterprise Linux 8.10 Oopta.

Python and dependencies to run all code with the expection of [simulations.py](./scripts//simulations.py) are listed in [ltee-general.yml](./environments/ltee-general.yml) and can be installed with `minicoda`.  See [https://www.anaconda.com/docs/getting-started/miniconda/main](https://www.anaconda.com/docs/getting-started/miniconda/main) for detailed instructions.  

Python and dependencies to run [simulations.py](./scripts//simulations.py) are listed in [ltee-sim.yml](./environments/ltee-sim.yml) which can be installed with `minicoda` as well.

### Data

Raw sequence reads can be downloaded from the [Sequence Read Archive (SRA)](https://www.ncbi.nlm.nih.gov/sra) under [BioProject Accession number PRJNA1549747](https://www.ncbi.nlm.nih.gov/bioproject/1549747). For more detailed on sequencing processed see [raw read processing](#raw-read-processing-via-flutubepy) below and the manuscript. 

(Intermediate) results directory can be obtained from [zenodo](https://zenodo.org/) under [link](). However most of intermediate results to run majority of the analyses can be directly obtained within this repository; see [results](./results/).

To compared our experimental evolution with global surveillance data we donwloaded all available human A/H3N2 sequence data from GISAID that was collected between 1968 and 2025-12-12 and submitted before May 13th 2026 from the [GISAID EpiFlu database](https://gisaid.org/). We can not share this data but if can be downloaded from GISAID directly.

Once again, we would like to thank all authors, submitting laboratories and contributors for sharing their influenza sequence data.

## Raw read processing via flutube.py

Library preparation and NGS of our viral RNA samples were performed by GenomeScan B.V.  
Raw sequencing reads as obtained directly from GenomeScan B.V. are deposited in the NCBI Sequence Read Archive (SRA) at [link].

Raw sequencing reads can be processed using the custom flu-tube pipeline.  
We processed the raw reads using the following command:

```python flutube.py -d <data_dir> -i <metadata_file>```

- where `<data_dir>` is the directory contained the raw sequence reads (as uploaded to SRA)
- and `<metadata_file>` is the [sample key](./data/complete_sample_key.csv)
All other input parameters were kept at default settings

flutube.py processes individual flowcells results independently and pools computed coverage and variant calls in the end.  

### final flutube output required for analysis

To replicate most of the analysis in you will need flu-tube.py's `{output}/coverage` as well as `{output}/sorted-aligned` (only for analysis involving haplotype reconstrucution), where `{output}` is the directory name supplied to `--output` (./output by default)

if you don't want to generate the flutube.py output yourself:

- coverage output files required to reproduce the analysis can be found in the [./results](./results) folder.
- sorted-aligned output can be downloaded a zip from [LINK] and should be stored in [./results](./results/) as folder [./results/sorted-aligned/_read_jsons](./results/sorted_aligned_read_jsons/)

## Calling SNVs 

Although the [flu-tube](./flu-tube.py) pipeline calls SNVs, we called our SNV once more in the [Coverage_to_snv_frequencies](notebooks/Coverage_to_snv_frequencies.ipynb) juypter notebook to assure reproducibility and to evaluate correlation between technical replicates.

We also accessed whether to implemented an error correction model in the [Error_model_test](./notebooks/Error_model_test.ipynb) juypter notebook, and decided **not to** implemented for calling SNVs.

## SNV frequency trajectories and SNV dynamics

SNV frequency trajectories as well a HMM states assignments are visualized and assessed within the [SNV_frequency_HMM_state_assignment](./notebooks/SNV_frequencies_HMM_state_assigment.ipynb) juypter notebook.
Statistical analysis of the timing to first fixation event (as identified by the HMM) can be found in [SNV_frequencies_HMM_state_assigment](./notebooks/SNV_frequencies_HMM_state_assigment.Rmd) R notebook.  

SNV frequency dynamics investigation as well as phylogenetic reconstruction and assessment of the nucleotide diverisity and total derived allele frequency over time can be found in the respective[Phylogenetic_reconstruction](./notebooks/Phylogenetic_reconstruction.ipynb) [Dynamics_of_observed_molecular_evolution](./notebooks/Dynamics_of_observed_molecular_evolution.ipynb) juypter notebooks.  
Statistical analysis investigating the rate of SNV accumulation (via total derived allele frequency) and nucleotide diversity can be found in the [Dynamics_of_observed_molecular_evolution](./notebooks/Dynamics_of_observed_molecular_evolution.Rmd) R notebook.  

## Identification of non-drift-like trajectories

To identify non-drift-like trajectories, we need the SNV trajectories and the max frequency and transition time associated with each SNV observation and calculated/determine in the [SNV_frequency_HMM_state_assignment](./notebooks/SNV_frequencies_HMM_state_assigment.ipynb) juypter notebook.  
To generate the null-model drift-like trajectories we ran [simulations.py](./scripts/simulations.py) as specified in [simulation_commands.sh](./scripts/bash/simulation_commands.sh). We then determined whether our observed trajectories deviated significantly from our simulated trajectories in the [Inference_of_non-drift-like_trajectories](./notebooks/Inference_of_non-drift-like_trajectories.ipynb) juypter notebook.

## Haplotype reconstruction

We reconstucted segment haplotype for sequence reads for all non-drift-like SNVs at timepoints at which their frequency was above 5%. Complete analysis can be found in the [Haplotype_reconstruction](./notebooks/Haplotype_reconstruction.ipynb) juypter notebook. This notebook also contains annotated versions of the muller plot per population.

## Investigation of confounding factors in non-drift-like trajectories

We assessed evidence for mutation bias and codon usage bias in the [Mutation_bias_and_CUB](./notebooks/Mutation_bias_and_CUB.ipynb) juypter notebook. Further statistical analysese into experiment-wide mutation bias can be found in the [Mutation_bias](./notebooks/Mutation_bias.Rmd) R notebook.  
We assessed genetic linkage between non-drift-like SNVS in the [Hitchhiker_and_cooccurence_identfication](./notebooks/Hitchhiker_and_cooccurence_identfication.ipynb) juypter notebook.
We generated the final summary figure for this analysis in the [Potential_confounding_factors_summary](./notebooks/Potential_confounding_factors_summary.ipynb) juypter notebook.

## Parallel evolution investigation

We investigated signatures of parallel evolution in the [Parallel_evolution_investigation](./notebooks/Parallel_evolution_investigation.ipynb).

## Comparison with global trajectories

global data was downloaded from [GISAID EpiFlu Database](https://app1.epicov.org/epi3/frontend#4c92d3) and has be download as we are not allowed to reshare this data. We downloaded all available human lineage influenza A/H3N2 virus sequences collected from human hosts between 1968 and December 31st 2025 and submitted before May 13th 2026. Preperation of global as well as subsequent comparison with experimental data can be found in the [Global_observation_frequencies](./notebooks/Global_observation_frequencies.ipynb) juypter notebook.

## Additional analyses

- Bottleneck input calculations can be found in the [Bottleneck_size_calculation](./notebooks/Bottleneck_size_calculation.ipynb) juypter notebook.
- Sensitivity analysis of the apical release fraction and wash recovery fraction parameters for the drift-like simulations can be found in the [Simulation_wash_apical_sensitivity](./notebooks/Simulation_wash_apical_sensitivity.ipynb) juypter notebook.

## Authors

When using this code or an adaptated version of this code, or any of the data generated in the experiment/analysis please cite: De Rooij, et al. .

All code published in this repository is written by Annelies J. H. de Rooij, a.j.h.derooij@amsterdamumc.nl.  
Code of simulations.py, haplotype_reconstruction.py and flu-tube.py are adapted from codes written by Alvin X. Han.

## Acknowledgements

We would like to thank all authors and contributors of the GISAID EpiFlu database for sharing their influenza A/H3N2 sequence data.
