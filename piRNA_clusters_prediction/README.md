#  Prediction of piRNA clusters (piCs)

Two approaches were used to predict piRNA clusters: a coverage-based approach and the proTRAC approach.

## Table of Contents

1. [Short reads pre-processing](#1---short-reads-pre-processing)  
2. [Coverage-based protocol](#2---coverage-based-protocol)  
3. [proTRAC protocol](#3---protrac-protocol)  
4. [Merge coverage-based and proTRAC piC annotations within strains](#4---merge-coverage-based-and-protrac-pic-annotations-within-strains)
5. [Cross-strain comparison of piRNA clusters - coordinates transfer](#5---cross-strain-comparison-of-pirna-clusters---coordinates-transfer)
6. [Intersecting piCs across Strains](#6---intersecting-pics-across-strains)
7. [Getting final non-lifted piCs](#7---getting-final-non-lifted-piCs)
8. [P-element, Shellder and Spoink insertions into piRNA clusters](#8---p-element-,shellder-and-poink-insertions-into-piRNA-clusters)

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

Mapping reads and getting the coverage, extracting potential piCs (CPM \>8) and estimating the percentage of 1U for each potential piC.


**Script:** [2.1_coverage_based_approach.sh](scripts/2.1_coverage_based_approach.sh)


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

**Script:** [2.3_coverage_based_approach_filtering_pics.R](scripts/2.3_coverage_based_approach_filtering_pics.R)

This R script performs the filtering and final selection of candidate piRNA clusters based on the coverage metrics and sequence composition calculated in the previous steps.
\>60%1U, \>25 RPKM, \>2kb

##### Inputs:

- Table produced in Step 2 containing coverage and normalised expression values (RPM and RPKM) for candidate regions: `all_samples_CPM8_coverage_based.bed`

- Table produced in Step 1 containing the proportion of reads with a 1U bias for each candidate region: `all_samples_1U_coverage_based.txt`
- TE annotation:  `all_samples_TE_annot_filtered.bed`

##### Output:

- Summary of all predicted piCs predicted by strain using coverage-based approach, after filtering: [piCs_coverage_based.txt](results/piCs_coverage_based.txt)

[Back to Table of Contents](#table-of-contents)

## 3 - **proTRAC protocol**

Filtered reads from each sample were mapped to their respective genome assemblies, allowing multi-mapping with [Bowtie](https://github.com/BenLangmead/bowtie) v1.3.1 (bowtie -n 1 -l 12 -a --best --strata -y --chunkmbs 1024 -S).
Next, proTRAC_2.4.3 [(Rosenkranz & Zischler, 2012)](https://pmc.ncbi.nlm.nih.gov/articles/PMC3293768/) was run on sam files (-swsize 500 -swinger 100 -size 500 -1Tor10A 0.6 -min 23 -chitin 10 -pens 0.05 -pti).
Clusters that did not map with TEs were removed.

### 3.1 - Step 1

Mapping filtered reads against the respective genome assemblies and running proTRAC.

**Script:** [3.1_protrac_approach.sh](script/3.1_protrac_approach.sh)

##### Inputs:

-   Deduplicated small RNA FASTQ files for all replicates of each strain: `*_dedup.fastq`
  
-   Strain-specific genome assemblies: `*_final_assembly.fasta`

##### **Main Outputs:**

-   SAM files for each replicate containing multiple mapping small RNAs: `*_Protrac.sam`

-   proTRAC output directories containing predicted piRNA clusters and associated statistics: `proTRAC_*`

### 3.2 - Step 2

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

    ### 3.3 - Step 3

Merging piCs and eliminating those that do not intersect with TEs

**Script:** [3.3_proTRAC_approach.R](scripts/3.3_proTRAC_approach.R)
##### Input: 

-   proTRAC final result: `all_samples_protrac_clusters.txt`

-   TE annotation: `all_samples_TE_annot_filtered.bed`

##### **Output**:

Summary of all piCs predicted by strain using proTRAC approach: [piCs_protrac.txt](results/piCs_protrac.txt)

[Back to Table of Contents](#table-of-contents)

# 4 - Merge coverage-based and proTRAC piC annotations within strains

Summary: Predicted piRNA clusters from both approaches were combined and merged using bedtools.
Clusters intersecting with genes were split to eliminate gene regions.
Clusters with a size less than 2kb were removed.
More details are described below.

### 4.1 - Concatenating results and splitting by strain

The files `piCs_coverage_based.txt` and `piCs_protrac.txt` were obtained previously.

```{bash eval = FALSE}
cat piCs_coverage_based.txt piCs_protrac.txt |sort -k1,1V -k2,2n >piCs_coverage_protrac.bed
awk -F'\t' '{ out="./piCs_per_strain/" $4 "_coverage_protrac.bed"; print > out }' piCs_coverage_protrac.bed
```

### 4.2 - Merging coverage-based and proTRAC results per strain

```{bash eval = FALSE}
cd ./tables/piCs_per_strain
for file in *_coverage_protrac.bed; do
    strain="${file%%_coverage_protrac.bed}"
    out="${strain}_coverage_protrac_merged.bed"
	  bedtools merge -i $file -c 4,5,6,7 -o distinct,collapse,collapse,collapse >$out
done
```

### 4.3 - Eliminating gene regions from piC annotation

Genic regions are frequently eliminated from piC annotations.
For this purpose, we first prepared a filtered reference gene annotation and lifted the coordinates to each strain genome assembly, and later, we transferred the coordinates.

#### 4.3.1 - Eliminating lncRNAs and TEs from REF GFF

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

### 4.3.2 - Liftoff gene annotation

This script projects the filtered reference gene annotation to each strain genome assembly using [Liftoff](https://github.com/agshumate/Liftoff).\
**Script:** [4.3.2_liftoff_REF_GFF.sh](scripts/4.3.2_liftoff_REF_GFF.sh)

##### Inputs:

-   Reference genome: `GCF_016746395.2_Prin_Dsim_3.1_genomic.fna`

-   Filtered reference annotation: `GCF_filtered_no_lncRNA_noTE.gff`

-   Genome assemblies: `*_final_assembly.fasta`

##### **Output:**

-   Lifted gene annotations: `*`\_final_assembly_liftoff_noTE_nolncRNA.gff.bed


### 4.3.3 - Splitting piCs based on the coordinates of the genes

**Script**: [4.3.3_spliting_pics_intersecting_genes.sh](scripts/4.3.3_spliting_pics_intersecting_genes.sh)

##### Inputs:

-   Predicted piRNA clusters:`*_coverage_protrac_merged.bed`

-   Gene annotations generated in the previous step (transformed to bed using [gff_to_bed.py](scripts/gff_to_bed.py)): `*_final_assembly_liftoff_noTE_nolncRNA.gff.bed`

-   TE annotations:`*_final_TE_annotation_morestringent.bed`

##### Output:

-   piRNA cluster annotation: `*_pic_outside_gene_2kb_final.bed`

### 4.3.4 - Getting final standardised piC results - before liftover

**Script:** [4.3.4_getting_final_standardised_results_before_liftover.R](scripts/4.3.4_getting_final_standardised_results_before_liftover.R)
##### Input: 

- piRNA cluster annotation: `pic_outside_gene_2kb_final_all.bed`

##### Output:

- Final piC annotation per strain before liftover: [summary_piCs_before_liftover.tsv](results/summary_piCs_before_liftover.tsv)

[Back to Table of Contents](#table-of-contents)

## 5 - Cross-strain comparison of piRNA clusters - coordinates transfer

Once potential piRNA clusters (piCs) were identified in each strain, the next step was to determine whether these clusters are strain-specific or shared across strains.
To enable this comparison, cluster coordinates were first transferred to the reference genome coordinates and then intersected across strains

### 5.1 - Transferring coordinates

For each strain, piC coordinates were transferred to Prin_Dsim_3.1 genome using the halLiftover tool [(Hickey et al., 2013)](https://pmc.ncbi.nlm.nih.gov/articles/PMC3654707/) using the hierarchical alignment (HAL) obtained by the pangenome analysis.
Around 70% of coordinates were transferred with this approach.
Additional efforts to obtain more transferred coordinates involved: 1) Using halLiftover for 10 kb upstream and downstream piC regions, which recovered an additional 381 coordinates. 2) Using Liftoff, we recovered an additional 109 coordinates, keeping only entries for larger piCs (\>10 kb) where the transferred coordinates were in the same chromosome with close proximity (\<100 kb).
Each step is described below:

#### 5.1.1 - Transferring coordinates from Pangenome alignment

Several file preparations are needed before running halLiftover:

#### 5.1.1.1 - Softmasking the genomes

Genomes need to be softmasked to run Cactus pangenome.
The genomes were masked using bedtools and their respective TE annotation.

**Script**: [5.1.1.1_getting_soft_masked_genomes.sh](scripts/5.1.1.1_getting_soft_masked_genomes.sh)

##### Inputs:

-   Genome assemblies for each strain: `*_final_assembly.fasta`

-   TE annotations for each strain: `*_final_TE_annotation_morestringent_new.bed`

##### **Output:**

-   Soft-masked genomes: `*final_assembly.softmasked.fasta`

#### **5.1.1.2 - Cleaning genome fasta headers**

Genome FASTA headers were simplified because the downstream program did not accept complex identifiers.

```{bash eval = FALSE}
for file in *_final_assembly.softmasked.fasta; do
  awk '/^>/ {
    header = substr($0, 2);
    sub(/_[^A-Za-z0-9]+$/, "", header);  # Remove ONLY non-alphanumeric suffixes after last _
    print ">" header;
    next
  } { print }' "$file" > "${file%.fasta}.cleaned.fasta"
done
```

##### **Output:**

-   Soft-masked genomes with edited headers: **`*`**`_final_assembly.softmasked.cleaned.fasta`

#### **5.1.1.3 - Cleaning piC BED file headers**

Because contig names were simplified in the genome FASTA files, the same modification had to be applied to the piC BED files to maintain coordinate compatibility.

##### **Input:**

-   piC annotation obtained for each strain: `*_pic_outside_gene_2kb_final.bed`

```{bash eval = FALSE}
for bed in *_pic_outside_gene_2kb_final.bed; do
  sample=$(basename "$bed" _pic_outside_gene_2kb_final.bed)
  awk '{
    chrom = $1;
    sub(/_[^A-Za-z0-9]+$/, "", chrom);  # Remove ONLY non-alphanumeric suffixes 
    $1 = chrom;
    print
  }' OFS='\t' "$bed" > "${sample}.cleaned.bed"
done
```

##### **Output:**

-   piC annotation with edited coordinates for each strain: `*_cleaned.bed`

#### **5.1.1.4 - Edit piC bed files to simplify the information**

To facilitate coordinate transfer and tracking of piRNA clusters, the BED files were simplified to contain only the genomic coordinates and a unique identifier in the fourth column.
Each piRNA cluster (piC) was assigned a sequential identifier (piC1, piC2, etc.), allowing clusters to be tracked after coordinate conversion.

**Script**: [5.1.1.4_simplify_piC_bed.py](scripts/5.1.1.4_simplify_piC_bed.py)

##### Input:

-   piC annotation with edited coordinates for each strain: `*_cleaned.bed`

##### **Output**:

-   Simplified piC annotation: `*_simplified.bed`

-   Reference table linking the original BED entries to the newly assigned piC identifiers: `*_reference_with_piC.txt`

### 5.1.2 - Running pangenome

The [cactus](https://github.com/ComparativeGenomicsToolkit/cactus) pangenome pipeline was run to construct a graph-based genome alignment.
The analysis was executed within the Cactus Singularity container.

```{bash eval = FALSE}
source /pub63/adriana/bin/cactus-bin-v2.9.8/venv-cactus-v2.9.8/bin/activate

singularity exec ~/bin/cactus-singularity/cactus_v2.9.7.sif cactus-pangenome ./jobStore  Dsim_genomes_new_pan.txt --outDir pangenome_new  --outName pangenome_new --reference REF --vcf --giraffe --gfa --gbz 
```

##### **Input:**

-   Configuration file listing the softmasked genomes to be included in the pangenome analysis and specifying the reference genome (REF): `Dsim_genomes_new_pan.txt`

-   Softmasked genomes assemblies: `*_final_assembly.softmasked.cleaned.fasta`

##### **Output:**

-   Directory containing all pangenome files output: `pangenome_new/`

### 5.1.3 - Transferring coordinates using halLiftover

#### 5.1.3.1 - Run halLiftover

This step transfers piC coordinates from each strain to the reference coordinates using the multiple genome alignment (.hal file) generated by Cactus pangenome.

```{bash eval = FALSE}
for bed in *_simplified.bed; do
  sample=$(basename "$bed" _simplified.bed)
  echo "Running halLiftover for $sample"
  halLiftover ../pangenome_new/pangenome_new.full.hal "$sample" "$bed" REF "../halLiftover/${sample}_liftover.bed"
done
```

##### Inputs:

-   Simplified piRNA cluster coordinates for each strain: `*_simplified.bed`

-   Multiple genome alignment generated by the Cactus pangenome: `pangenome_new.full.hal`

##### Output:

-   BED file for each strain containing piRNA cluster coordinates converted to the REF coordinates (a single piC may be represented by multiple adjacent intervals, reflecting alignment blocks in the HAL genome alignment.): `*_liftover.bed`
-   

#### 5.1.3.2 - Getting start and end for each piC

The output of halLiftover may contain multiple entries for the same piRNA cluster because the coordinate transfer follows the genome alignment blocks.
As a result, a single piC may be split into several adjacent segments in the reference genome.
To reconstruct the full cluster coordinates, these fragments were merged by retaining the first start position and the last end position for each piC.

```{bash eval = FALSE}
e.g. piC2 coordinates will be NC_052520.2:13508030-13508216
NC_052520.2     13508030        13508166        piC2 #first start:13508030
NC_052520.2     13508166        13508175        piC2
NC_052520.2     13508175        13508201        piC2
NC_052520.2     13508201        13508216        piC2 #last end:13508216
```

**Script**: [5.1.3.2_merge_transfered_pic_ranges.py](scripts/5.1.3.2_merge_transfered_pic_ranges.py)

```{bash eval = FALSE}
for bed in *_liftover.bed; do
  sample=${bed%_liftover.bed}
  python ./merge_transfered_pic_ranges.py "$bed" "${sample}_merged.bed"
done
```

##### Output: 

BED files containing a single interval per piRNA cluster in reference genome coordinates: `*_merged.bed`

## 5.2 - Transferring coordinates using flanking regions and pangenome

Some piRNA clusters could not be transferred to the reference genome using halLiftover.
To recover these clusters, an alternative strategy based on flanking regions and the pangenome alignment was applied.

### 5.2.1 - Getting non-transferred piCs

This step identifies piRNA clusters that were not successfully transferred by halLiftover by comparing the original piC identifiers with those present in the lifted coordinates.
Clusters missing from the liftover output are extracted for downstream coordinate recovery.

**Script:** [5.2.1_extract_nonlifted_piCs_1.sh](scripts/5.2.1_extract_nonlifted_piCs_1.sh)

##### **Output:**

`/non_lifted/* _non_lifted.bed`

### 5.2.2 - Getting flanking coordinates - 10kb

For piRNA clusters that were not transferred by halLiftover, this step extracts the 10 kb upstream and downstream coordinates.
These flanking intervals are used to try to recover the reference coordinates using halLiftover.

**Script**: [5.2.2_getting_flanking_non_lifted.py](scripts/5.2.2_getting_flanking_non_lifted.py)

##### **Output:** 

-   10 kb upstream coordinates for each non-transferred piC: `lift_remaining_flanking/*_non_lifted_upstream.bed`

-   10 kb downstream coordinates for each non-transferred piC: `lift_remaining_flanking/*_non_lifted_downstream.bed`

### 5.2.3 - Using halLiftover - flanking regions

#### Run halLiftover for flanking regions

```{bash eval = FALSE}
for bed in *_non_lifted_downstream.bed; do
  sample=$(basename "$bed" _non_lifted_downstream.bed)
  echo "Running halLiftover for $sample..."
  halLiftover ../../../pangenome_new/pangenome_new.full.hal "$sample" "$bed" REF "${sample}_downstream_liftover.bed"
done

for bed in *_non_lifted_upstream.bed; do
  sample=$(basename "$bed" _non_lifted_upstream.bed)
  echo "Running halLiftover for $sample..."
  halLiftover ../../../pangenome_new/pangenome_new.full.hal "$sample" "$bed" REF "${sample}_upstream_liftover.bed"
done
```

#### Getting start and end for each piC

**Script**: same used before [5.1.3.2_merge_transfered_pic_ranges.py](scripts/5.1.3.2_merge_transfered_pic_ranges.py)

```{bash eval = FALSE}
for bed in *_upstream_liftover.bed; do
  sample=${bed%_upstream_liftover.bed}
  python ./merge_transfered_pic_ranges.py "$bed" "${sample}_upstream_liftover_merged.bed"
done
for bed in *_downstream_liftover.bed; do
  sample=${bed%_downstream_liftover.bed}
  python ./merge_transfered_pic_ranges.py "$bed" "${sample}_downstream_liftover_merged.bed"
done
```

#### Merging up and down alignments of flanking regions to get the transferred coordinates

**Script**: [5.2.3_merging_up_and_down_flanking_lifted.py](scripts/5.2.3_merging_up_and_down_flanking_lifted.py)

##### Main output: 

Additional piCs (non-lifted with the first approach) with transferred coordinates for each strain: `*_non_lifted_lifted.bed`

#### 5.2.4 - Concatenating with previous lifted piCs

```{bash eval = FALSE}
for file in *_non_lifted_lifted.bed; do
  sample=$(basename "$file" _non_lifted_lifted.bed)
  
  # Concatenate and sort by chromosome (1st col) and start (2nd col)
  cat "../../${sample}_merged.bed" "${sample}_non_lifted_lifted.bed" \
    | sort -k1,1 -k2,2n \
    > "concat_lifted/${sample}_lifted.bed"
done
```
##### Output:

- `concat_lifted/*_lifted.bed`

## 5.3 - Transferring coordinates Using Liftoff

An additional attempt to transfer coordinates for the remaining piCs was done with Liftoff.

### 5.3.1 - Getting non-transferred piCs

**Script:** [5.3.1_extract_nonlifted_piCs_2.sh](scripts/5.3.1_extract_nonlifted_piCs_2.sh)

##### **Outputs:**

`lift_remaining_flanking/concat_lifted/non_lifted/*_non_lifted.bed`

### 5.3.2 - Running Liftoff

**Script:** [5.3.2_liftoff_remaining_piCs.sh](scripts/5.3.2_liftoff_remaining_piCs.sh)

##### **Main output:**

Additional piCs with transferred coordinates for each strain (non-lifted with the two previous approaches): `/lift_remaining_flanking/concat_lifted/non_lifted/liftoff_remaining/\*\_non_lifted_lifted.bed`

### 5.3.3 - Concatenating all transferred coordinates

All lifted piCs for each strain were concatenated in a final file: `*_final_lifted.bed`


## 5.4 - Recovering original information of piCs
This script appends the original information of the piC annotation to the lifted coordinates

Script: [5.4_recovering_original_piC_info.sh](scripts/5.4_recovering_original_piC_info.sh)

##### Output: `*_final_lifted_ed.bed`

Bed-like file to intersect
```{bash eval = FALSE}
for file in *_final_lifted_ed.bed; do
     sample=$(basename "$file" _final_lifted_ed.bed)
      awk 'OFS="\t"{print $1,$2,$3}' $file >${sample}_final_lifted_to_intersect.bed
done
```
[Back to Table of Contents](#table-of-contents)

## 6 - Intersecting piCs across Strains

### 6.1 - Concatenating piC annotation files

```{bash eval = FALSE}
#Select columns chrom, start, end, strain
>cat_final_lifted_to_intersect.bed
for file in *_final_lifted_ed.bed; do awk 'OFS="\t"{print $1,$2,$3,$8}' $file >>cat_final_lifted_to_intersect.bed; done
#sort
sort -k1,1 -k2,2n cat_final_lifted_to_intersect.bed >cat_final_lifted_to_intersect_sorted.bed
```

### 6.2 - Subtracting gene regions again

During coordinate transfer, some regions overlapping genes were reintroduced into the dataset. Although gene filtering had already been applied per strain earlier in the pipeline (which is important for accurate statistics within strains), this step is required to remove any gene-overlapping regions that reappeared after coordinates transfer. Therefore, gene filtering was applied again.

More info on how the reference GFF was filtered is here: [4.3.1- Eliminating gene regions from piC annotation](#431---eliminating-gene-regions%20from-pic)

```{bash eval = FALSE}
#subtract genes
bedtools subtract -a cat_final_lifted_to_intersect_sorted.bed -b GCF_filtered_no_lncRNA_noTE_ID_gene.bed |awk 'OFS="\t" {print $0,"cluster_"NR}'>cat_final_lifted_nogenes.tsv

#getting a bed file with only genes that intersect with the piCs
bedtools intersect -a GCF_filtered_no_lncRNA_noTE_ID_gene.bed -b cat_final_lifted_to_intersect_sorted.bed  -wa |awk 'OFS="\t"{print $1,$2,$3,$4, "gene"}'|sort|uniq > genes_intersecting
```

### 6.3 - Merging piCs across strains

```{bash eval = FALSE}
sort -k1,1 -k2,2n cat_final_lifted_nogenes.tsv|bedtools merge -i - -c 4 -o collapse >merged_lifted_final.bed
```

### 6.4 - Minimal cluster curation

We manually curated the clusters by inspecting their TE composition. Clusters without TEs were removed, and clusters corresponding to the same locus across strains were given the same ID (when eventually they initially had different ones).

Output: `merged_lifted_final_piCs_inspected.bed`

### 6.5 - Getting full information of piCs

To get the original information of piRNA clusters in the strains, we intersected the merged annotations with the original files.

```{bash eval = FALSE}
pwd:/data/adriana/TE_Invasions/Nanopore/cactus/cactus_singularity/halLiftover/final_piCs_lifted
for file in *_final_lifted_ed.bed; do
    sample=$(basename "$file" _final_lifted_ed.bed)
bedtools intersect -a merged_lifted_final_piCs_inspected.bed -b $file -wa -wb >${sample}_reverted_piCs
done

cat *_reverted_piCs >cat_reverted_piCs.txt
```

### 6.6 - Getting final standardised piC results - after liftover**

Script: [6.6_getting_final_standardised_results_after_liftover.R](scripts/6.6_getting_final_standardised_results_after_liftover.R)

##### Input:

-   piRNA cluster annotation: `cat_reverted_piCs.txt`

##### Output:

-   Final piC annotation per strain after liftover: [summary_piCs_after_liftover.tsv](results/summary_piCs_after_liftover.tsv)

[Back to Table of Contents](#table-of-contents)

## 7 - Getting final non-lifted piCs

**Script:** [7_extract_nonlifted_piCs_final.sh](scripts/7_extract_nonlifted_piCs_final.sh)

**Output:** \*\_non_lifted.bed;

Output for all strains: [all_final_non_lifted_pics.bed](results/all_final_non_lifted_pics.bed)

[Back to Table of Contents](#table-of-contents)

# 8 - P-element, Shellder and Spoink insertions into piRNA clusters

We analysed the long-read sequencing data using [Sniffles2](https://github.com/fritzsedlazeck/Sniffles) to identify structural variants, followed by a pipeline to extract all insertion sequences and identify those that are P-element, Spoink and Shellder. The TSD sequences were also extracted.

See [\*\*TE insertion detection from long reads](**%3Chttps://github.com/adriludwig/Drosophila_simulans_piRNA_clusters/tree/main/TE_insertion_detection_from_long_reads%3E**).\*\* 
Later, we intersected the Sniffles insertions with the TE genome annotation to identify the assembled insertions. We also intersect the insertion coordinates with the piRNA cluster annotation. The results were reviewed and edited as needed. TE copies located within or near low-complexity regions introduced during assembly polishing (e.g., Racon and Medaka) were occasionally missed due to local sequence distortions and minor coordinate discrepancies. Only insertions confirmed by sniffles or by visual inspection in IGV were considered.

## 8.1 - Filtering and checking insertions

Some filtering was applied to eliminate insertions with low support (support \<2).

Insertions in piRNA clusters with less than 5 reads of support (total) or less than 3 reads spanning the insertions were flagged to be inspected in IGV. Only the more confident ones were maintained in the dataset.

**Script:** [Filtering&Checking_Insertions.Rmd](scripts/Filtering&Checking_Insertions.Rmd)

##### Inputs:

-   Edited dataset of all insertions of P, Shellder and Spoink found using Sniflles2 crossmatched with the TE annotation from the assembled genomes and the piRNA cluster annotation: [P_Shellder_Spoink_all_strains_insertions_ed.xlsx](results/P_Shellder_Spoink_all_strains_insertions_ed.xlsx)

-   piRNA cluster annotation (to recover complete information of clusters):[summary_piCs_after_liftover.tsv](piRNA_clusters_prediction/results/summary_piCs_after_liftover.tsv)

##### Outputs:

-   [piCs_containing_P_insertions.tsv](results/piCs_containing_P_insertions.tsv)

-   [piCs_containing_Spoink_insertions.tsv](results/piCs_containing_P_insertions.tsv)

-   [piCs_containing_Shellder_insertions.tsv](results/piCs_containing_Shellder_insertions.tsv)


[Back to Table of Contents](#table-of-contents)

