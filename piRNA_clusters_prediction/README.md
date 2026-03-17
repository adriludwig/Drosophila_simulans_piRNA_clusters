#  Prediction of piRNA clusters (piCs)

Two approaches were used to predict piRNA clusters: a coverage-based approach and the proTRAC approach.

## Table of Contents

1. [Short reads pre-processing](#1---short-reads-pre-processing)  
2. [Coverage-based protocol](#2---coverage-based-protocol)  
3. [proTRAC protocol](#3---protrac-protocol)  
4. [Merge coverage-based and proTRAC piC annotations within strains](#4---merge-coverage-based-and-proTRAC-piC-annotations-within-strains)  


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

**Script:** [2.1_coverage_based_approach.sh](scripts/2.1_coverage_based_approach.sh)

Mapping reads and getting the coverage, extracting potential piCs (CPM \>8) and estimating the percentage of 1U for each potential piC.

##### Inputs:

-   Deduplicated small RNA FASTQ files for all replicates of each strain: `*_dedup.fastq`
  
-   Strain-specific genome assemblies: `*_final_assembly.fasta`

##### **Main Outputs:**

- Final BAM file for each replicate containing uniquely mapping small RNAs: `*_combined_unique_sorted2.bam`

- Potential piC regions for each replicate: `*_500bp_8CPM_merge10kb.bed`

- Potential piC regions for each replicate with their percentage of 1U: [all_samples_1U_coverage_based.txt](results/all_samples_1U_coverage_based.txt)

### 2.2 - Step 2

**Script:** [2.2_coverage_based_approach.sh](scripts/2.2_coverage_based_approach.sh)

For each sample, it computes total mapped reads, raw coverage, RPM, and RPKM, writes a per-sample table, and then concatenates all per-sample results into a single combined file.

##### Input:

- Final BAM generated in Step 1: `*_combined_unique_sorted2.bam`

- BED file generated in Step 1: `*_500bp_8CPM_merge10kb.bed`

##### **Output:** 

- Potential piC regions for each replicate with coverage measures: [all_samples_CPM8_coverage_based.bed](results/all_samples_CPM8_coverage_based.bed)

### 2.3 - Step 3

**Script:** [2.3_coverage_based_approach_filtering_pics.r](scripts/2.3_coverage_based_approach_filtering_pics.r)

This R script performs the filtering and final selection of candidate piRNA clusters based on the coverage metrics and sequence composition calculated in the previous steps.
\>60%1U, \>25 RPKM, \>2kb

##### Inputs:

- Table produced in Step 2 containing coverage and normalised expression values (RPM and RPKM) for candidate regions: `all_samples_CPM8_coverage_based.bed`

- Table produced in Step 1 containing the proportion of reads with a 1U bias for each candidate region: `all_samples_1U_coverage_based.txt`
- TE annotation:  `all_samples_TE_annot_filtered.bed`

##### Output:

- Summary of all predicted piCs predicted by strain using coverage-based approach, after filtering: [piCs_coverage_based.txt](results/piCs_coverage_based.txt)

## 3 - **proTRAC protocol**

Filtered reads from each sample were mapped to their respective genome assemblies, allowing multi-mapping with [Bowtie](https://github.com/BenLangmead/bowtie) v1.3.1 (bowtie -n 1 -l 12 -a --best --strata -y --chunkmbs 1024 -S).
Next, proTRAC_2.4.3 [(Rosenkranz & Zischler, 2012)](https://pmc.ncbi.nlm.nih.gov/articles/PMC3293768/) was run on sam files (-swsize 500 -swinger 100 -size 500 -1Tor10A 0.6 -min 23 -chitin 10 -pens 0.05 -pti).
Clusters that did not map with TEs were removed.

### 3.1 - Step 1

Mapping filtered reads against the respective genome assemblies and running proTRAC.

**Script:** [3.1_protrac_approach.sh](script/3.1_protrac_approach.sh)

##### Inputs:

-   Deduplicated small RNA FASTQ files for all replicates of each strain: `*_dedup.fastq`
-   
-   Strain-specific genome assemblies: `*_final_assembly.fasta`

##### **Main Outputs:**

-   SAM files for each replicate containing multiple mapping small RNAs: `*_Protrac.sam`

-   proTRAC output directories containing predicted piRNA clusters and associated statistics: `proTRAC_*`

### 3.2 Step2

Parsing results and getting RPKM

**Script:** [3.2_protrac_approach.sh](scrips/3.2_protrac_approach.sh)

##### Inputs:

-   proTRAC cluster predictions generated in Step 1: `proTRAC_*/clusters.gtf`

-   SAM files generated in Step 1: `*_Protrac.sam`

##### **Outputs:**

-   parsed proTRAC result for each strain: `*_parsed_clusters_simple.txt`

-   proTRAC final result from each strain obtained with the command below: [all_samples_protrac_clusters.txt](results/all_samples_protrac_clusters.txt)

    ```{bash eval = FALSE}
    first=1
    for file in *_parsed_clusters_simple.txt; do
        sample=$(basename "$file" _parsed_clusters_simple.txt)
        if [[ $first -eq 1 ]]; then
            awk -v sample="$sample" 'BEGIN {OFS="\t"} NR==1 {print $0, "Sample"} NR>1 {print $0, sample}' "$file"
            first=0
        else
            awk -v sample="$sample" 'BEGIN {OFS="\t"} NR>1 {print $0, sample}' "$file"
        fi
    done > all_samples_protrac_clusters.txt
    ```

    ### 3.3 Step 3

Merging piCs and eliminating those that do not intersect with TEs

**Script:** [3.3_proTRAC_approach.r](scripts/3.3_proTRAC_approach.r)
##### Input: 

-   proTRAC final result: `all_samples_protrac_clusters.txt`

-   TE annotation: `all_samples_TE_annot_filtered.bed`

##### **Output**:

Summary of all piCs predicted by strain using proTRAC approach: [piCs_protrac.txt](results/piCs_protrac.txt)

# 4 - Merge coverage-based and proTRAC piC annotations within strains

Summary: Predicted piRNA clusters from both approaches were combined and merged using bedtools.
Clusters intersecting with genes were split to eliminate gene regions.
Clusters with a size less than 2kb were removed.
More details are described below.

### 4.1 Concatenating results and splitting by strain

The files `piCs_coverage_based.txt` and `piCs_protrac.txt` were obtained previously.

```{bash eval = FALSE}
cat piCs_coverage_based.txt piCs_protrac.txt |sort -k1,1V -k2,2n >piCs_coverage_protrac.bed
awk -F'\t' '{ out="./piCs_per_strain/" $4 "_coverage_protrac.bed"; print > out }' piCs_coverage_protrac.bed
```

### 4.2 Merging coverage-based and proTRAC results per strain

```{bash eval = FALSE}
cd ./tables/piCs_per_strain
for file in *_coverage_protrac.bed; do
    strain="${file%%_coverage_protrac.bed}"
    out="${strain}_coverage_protrac_merged.bed"
	  bedtools merge -i $file -c 4,5,6,7 -o distinct,collapse,collapse,collapse >$out
done
```

### 4.3 Eliminating gene regions from piC annotation

Genic regions are frequently eliminated from piC annotations.
For this purpose, we first prepared a filtered reference gene annotation and lifted the coordinates to each strain genome assembly, and later, we transferred the coordinates.

#### 4.3.1 Eliminating lncRNAs and TEs from REF GFF

The *D. simulans* reference annotation ([GCF_016746395.2_Prin_Dsim_3.1_genomic.gff](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/746/395/GCF_016746395.2_Prin_Dsim_3.1/)) was downloaded from NCBI.
Because some lncRNAs may correspond to piCs, and some TEs are annotated as genes, we filter them out of the reference annotation.
lncRNAs were removed by string searching, and TEs were removed after blasting against exons.
AGAT was used to filter the GFF.

**Script:** [4.3.1_filtering_REF_GFF.sh](scripts/4.3.1_filtering_REF_GFF.sh)

##### Inputs:

-   Reference genome: `GCF_016746395.2_Prin_Dsim_3.1_genomic.fna`

-   Reference gene annotation: `GCF_016746395.2_Prin_Dsim_3.1_genomic.gff`

-   TE library: `Manual_Curated_Dsim_TE_lib.fasta`

##### **Output:**

-   Filtered GFF for the reference: `GCF_filtered_no_lncRNA_noTE.gff`

### 4.3.2 Liftoff gene annotation

This script projects the filtered reference gene annotation to each strain genome assembly using [Liftoff](https://github.com/agshumate/Liftoff).\
**Script:** [4.3.2_liftoff_REF_GFF.sh](scripts/4.3.2_liftoff_REF_GFF.sh)

##### Inputs:

-   Reference genome: `GCF_016746395.2_Prin_Dsim_3.1_genomic.fna`

-   Filtered reference annotation: `GCF_filtered_no_lncRNA_noTE.gff`

-   Genome assemblies: `*_final_assembly.fasta`

##### **Output:**

-   Lifted gene annotations: `*`\_final_assembly_liftoff_noTE_nolncRNA.gff.bed


### 4.3.3 Splitting piCs based on the coordinates of the genes

**Script**: [4.3.3_spliting_pics_intersecting_genes.sh](scripts/4.3.3_spliting_pics_intersecting_genes.sh)

##### Inputs:

-   Predicted piRNA clusters:`*_coverage_protrac_merged.bed`

-   Gene annotations generated in the previous step (transformed to bed using [gff_to_bed.py](scripts/gff_to_bed.py): `*_final_assembly_liftoff_noTE_nolncRNA.gff.bed`

-   TE annotations:`*_final_TE_annotation_morestringent.bed`

##### Output:

-   piRNA cluster annotation: `*_pic_outside_gene_2kb_final.bed`

### 4.3.4 Getting final standardised piC results - before liftover

##### Input: 

- piRNA cluster annotation: `pic_outside_gene_2kb_final_all.bed`

##### Output:

- Final piC annotation per strain before liftover: [summary_piCs_before_liftover.tsv](results/summary_piCs_before_liftover.tsv)

