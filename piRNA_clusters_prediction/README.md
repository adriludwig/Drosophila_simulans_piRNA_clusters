#  Prediction of piRNA clusters (piCs)

Two approaches were used to predict piRNA clusters: a coverage-based approach and the proTRAC approach.

## 1 - Short reads pre-processing

Short reads were trimmed with Cutadapt v1.2.1 [(Martin, 2011)](https://journal.embnet.org/index.php/embnetjournal/article/view/200) to remove adapters and select reads of 23-29 nucleotides.
Reads mapping to annotated rRNA, tRNAs, miRNA and miscRNAs were eliminated using Bowtie v1.3.1 (-v 2 -k1 -y -un) [(Langmead et al., 2009)](https://pmc.ncbi.nlm.nih.gov/articles/PMC2690996/).
The different sets of RNAs were downloaded from Flybase for *D. melanogaster*.
The remaining reads were then mapped to *D. simulans* mRNA sequences (GCF_016746395.2_Prin_Dsim_3.1_rna.fna) with [Bowtie](https://github.com/BenLangmead/bowtie) v1.3.1 (-n 1 -l 12 -a --best --strata -y).
Reads mapping the sense strand of genes were removed to eliminate potential small RNAs derived from degraded RNAs.
Because many TEs are located within introns, small RNAs mapping to TEs (using [BWA aln](https://pmc.ncbi.nlm.nih.gov/articles/PMC2705234/)) were added back to the dataset.
TBr2_collapse and TBr2_duster from the NGS toolbox [(Rosenkranz & Zischler, 2012)](https://pmc.ncbi.nlm.nih.gov/articles/PMC3293768/) were used to collapse the remaining reads to include only unique piRNAs and to remove low-complexity reads.

## 2 - **Coverage-based protocol**

Filtered reads from each sample were mapped to their respective genome assemblies, allowing only unique mapping.
To maximise the number of unique mappings, we performed two rounds of mapping using [Bowtie](https://github.com/BenLangmead/bowtie) v1.3.1.
In the first round, no mismatches were allowed (bowtie -n 0 -l 12 -a -m 1).
In the second round, the unmapped reads from the first round (recovered with the fastq function of samtools v1.15.1) were remapped, allowing for one mismatch (bowtie -n 1 -l 12 -a -m 1).
Finally, we used [samtools](https://github.com/samtools) merge to combine the SAM files.
Windows of 500 bp were created with bamCoverage v3.5.5 of deepTools ([Ramírez et al., 2014](https://academic.oup.com/nar/article/42/W1/W187/2435511?guestAccessKey=)) to calculate normalised counts (CPM).
Windows with an expression of at least 8 CPM were filtered and merged if they were located within 10 kb of each other using samtools.
RPKM values were then calculated for the resulting regions, and those with an RPKM value of at least 25 were selected.
To analyse the 1U bias, genomic regions were extracted using bedtools getfasta, and the reads were remapped to these regions with the same parameters.
The BAM files were subsequently divided into separate files, and we used [FastQC](https://github.com/s-andrews/FastQC) v0.12.1 to produce a quality control report for each of them from which the 1U proportion of each region was extracted. Clusters with 1U frequencies above 60% were maintained. Clusters that did not map with TEs were removed.

### **2.1 - Step 1**

**Script:** [01_coverage_based_approach.sh](scripts/01_coverage_based_approach.sh)

Mapping reads and getting the coverage, extracting potential piCs (CPM \>8) and estimating the percentage of 1U for each potential piC.

##### Inputs:

-   Deduplicated small RNA FASTQ files for all replicates of each strain: `*_dedup.fastq`
  
-   Strain-specific genome assemblies: `*_final_assembly.fasta`

##### **Main Outputs:**

- Final BAM file for each replicate containing uniquely mapping small RNAs: `*_combined_unique_sorted2.bam`

- Potential piC regions for each replicate: `*_500bp_8CPM_merge10kb.bed`

- Potential piC regions for each replicate with their percentage of 1U: [all_samples_1U_coverage_based.txt](results/all_samples_1U_coverage_based.txt)

### 2.2 - Step 2

**Script:** [02_coverage_based_approach.sh](scripts/02_coverage_based_approach.sh)

For each sample, it computes total mapped reads, raw coverage, RPM, and RPKM, writes a per-sample table, and then concatenates all per-sample results into a single combined file.

##### Input:

- Final BAM generated in Step 1: `*_combined_unique_sorted2.bam`

- BED file generated in Step 1 : `*_500bp_8CPM_merge10kb.bed`

##### **Output:** 

- Potential piC regions for each replicate with coverage measures: [all_samples_CPM8_coverage_based.bed](results/all_samples_CPM8_coverage_based.bed)

### 2.3 - Step 3

**Script:** [03_coverage_based_approach_filtering_pics.r](scripts/03_coverage_based_approach_filtering_pics.r)

This R script performs the filtering and final selection of candidate piRNA clusters based on the coverage metrics and sequence composition calculated in the previous steps.
\>60%1U, \>25 RPKM, \>2kb

##### Inputs:

- Table produced in Step 2 containing coverage and normalised expression values (RPM and RPKM) for candidate regions: `all_samples_CPM8_coverage_based.bed`

- Table produced in Step 1 containing the proportion of reads with a 1U bias for each candidate region: `all_samples_1U_coverage_based.txt`
- TE annotation:  `all_samples_TE_annot_filtered.bed`

##### Output:

- Summary of all predicted piCs predicted by strain using coverage-based approach, after filtering: [piCs_coverage_based.txt](results/piCs_coverage_based.txt)
