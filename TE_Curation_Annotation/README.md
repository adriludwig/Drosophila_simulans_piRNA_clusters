# TE curation and annotation

## Data availability

***Drosophila simulans* curated TE library**: The *Drosophila simulans* curated TE library obtained in this work is available at [Manual_Curated_Dsim_TE_lib.fasta](TE_libraries/Manual_Curated_Dsim_TE_lib.fasta)

***Drosophila simulans* strains TE annotation**: the final TE annotation for each strain is available at [TE_annotation](TE_annotation)

# 1 - *D. simulans* TE library curation

## **1.2 - Pipeline overview**

![**Workflow of TE library curation**](https://www.dropbox.com/scl/fi/w0n79kfyqgzuwgq1x2kd4/Workflow_TE_lib_curation.png?rlkey=nx8sgv0ivhmhpnit7sq9k6vho&raw=1)

RepeatModeler v2.0.3 (Flynn et al., 2020) was used to predict *de novo* TEs in the reference *D. simulans* genome Prin_Dsim_3.1 (GCF_016746395.2). A number of TE consensus sequences were correctly curated either with MCHelper v1.7.0 ([https://github.com/GonzalezLab/MCHelper)](https://github.com/GonzalezLab/MCHelper)) with default parameters and fully automated analysis or TE-trimmer v1.4.0 (<https://github.com/qjiangzhao/TEtrimmer>) default parameters.

Sequences that were clearly not correct after checking TE-aid plots ([https://github.com/clemgoub/TE-Aid)](https://github.com/clemgoub/TE-Aid) were manually curated. Manual curation was performed according to the TE manual curation guide ([Goubert et al., 2022](https://link.springer.com/article/10.1186/s13100-021-00259-7)). In some cases, an online BLAST self-alignment (using the same sequence as both query and subject, with the parameters “align two or more sequences” and “somewhat similar sequences”) was used to help determine the boundaries of TE families, either by identifying their terminal repeats or target-site duplications.

Additional TE sequences from another 22 *D. simulans* samples were added to the curated TE library. RepeatModeler was run for each additional genome, followed by MCHelper. The MCHelper consensus sequences for each sample were concatenated into a single multifasta file, and redundancy was removed using cd-hit (-d 0 -aS 0.8 -c 0.8 -G 0 -g 1 -b 500). The remaining sequences were then blasted against the first version of the curated library (-outfmt "6 qseqid sseqid pident length mismatch qstart qend qlen slen" -evalue 1e-20 -reward 2 -penalty -3 -gapopen 5 -gapextend 2 -max_target_seqs 1). Sequences that did not meet the 80% identity criterion over 80% query sequence coverage were inspected and curated as described above.

Moreover, four new TE sequences were added to the library after running [pantera](https://link.springer.com/article/10.1186/s13100-024-00323-y), which identifies polymorphic TEs from pangenome data (Sierra & Durbin, 2024). The final library consists of 224 sequences.

Several library entries correspond to incomplete or highly degenerated TE fragments for which a full consensus sequence could not be reconstructed. These sequences, named fragX, were retained in the library because they remain informative for detecting TE-derived sequences, particularly within piRNA clusters.

## 1.3 - TE families nomenclature

TEs were named following a stepwise decision process after checking their similarity to known TEs in several resources: 1) the online Repbase database; 2) the *D. melanogaster* MCTE library ([Rech et al., 2022](https://www.nature.com/articles/s41467-022-29518-8)) a multispecies *Drosophila* library ([Chakraborty et al., 2021](https://genome.cshlp.org/content/31/3/380.long)) a set of new LTR sequences identified by ([Bargues & Lerat, 2017](https://pubmed.ncbi.nlm.nih.gov/28465726).

When a TE showed high similarity (\>95%) to a previously described *D. simulans* element, the existing name was retained. If no such match was found, the closest homolog in other *Drosophila* species was considered. If this name was not already in use in *D. simulans*, we adopted that name (e.g., G5_Ds, which is closer to G5 from *D. melanogaster).*

If the closest name was already used in *D. simulans*, or if no closely related element was identified in any database, a new name was assigned using a sequential numbering system within the corresponding TE superfamily (e.g., Copia3_Ds, following Copia2_Ds).

In a few cases, the same name has been assigned to distinct TEs in different databases (e.g., Gypsy-12_DSim from Repbase and GYPSY12_Dsi from Bargues & Lerat, 2017). In such cases, we assigned distinct names (e.g., Gypsy12_Ds and Gypsy12A_Ds) to avoid ambiguity.

For consistency, sequential numbering was written without spaces between the TE family name and the number. Subfamilies were indicated using a dash, followed by the suffix \_Ds, and the TE classification was provided after \# according to the standard format (e.g., Gypsy2-2_Ds#LTR/Gypsy).

Additional details on the naming rationale are provided in the file [Manual_curated_Dsim_lib_Info.xlsx](TE_libraries/Manual_curated_Dsim_lib_Info.xlsx)

# 2 - *D. simulans* strains TE annotation

## **2.1 - Pipeline overview**

![**Workflow of the TE annotation**](https://www.dropbox.com/scl/fi/sh8umoi05tnx3oyc1ima2/Workflow_TE_annot.png?rlkey=ef7e7cfaka8utx1zfyyx6uhhp&raw=1)

Genome assemblies were initially annotated using RepeatMasker with the curated *D. simulans* TE library. However, this approach revealed two sources of inconsistency in TE annotation. First, some sequences containing low-complexity regions were not effectively masked by RepeatMasker, leading to spurious TE annotations. Second, small genomic regions were sometimes incorrectly annotated as fragments of TE families known to be absent from a given strain (e.g., Spoink in the reference genome). Manual inspection indicated that these regions corresponded to remnants of other TE families that were more accurately classified when using a broader TE library.

To address these issues, we developed a dedicated pipeline that combines [RepeatMasker](https://www.repeatmasker.org/) and [BLAST](https://blast.ncbi.nlm.nih.gov/doc/blast-help/references.html) to enhance the accuracy and reliability of TE annotations. The annotation process is divided into two main modules: first, we use the curated *D. simulans* TE library to annotate the more conserved TE families. In the second module, we employ the combined *Drosophila* TE library to annotate the more divergent TE fragments.

**Module 1:** Genome assemblies were initially annotated using RepeatMasker v4.1.1 with the curated *D. simulans* TE library, employing -div 20 to recover more conserved sequences and -nolow to avoid fragmentation of complete copies. Potential TE fragment sequences were retrieved from the respective genome using getfasta from [BEDtools](https://bedtools.readthedocs.io/en/latest/)v2.30.0 based on the RepeatMasker annotation. To retain only reliable TE fragments, they were further masked for simple repeats and low-complexity regions using RepeatMasker (-noint) and a blastn search was performed against the curated *D. simulans* TE library (-e-value 1e^-10^; -word-size 7; -max_target_seqs 1; -max_hsps 1) eliminating sequences with no significant matches. Finally, potential split fragments were merged using merge from BEDtools v2.30.0 with a window size of 50 (-d 50) and fragments less than 80 bp were excluded. An additional blastn step is performed against the curated library to obtain complete TE classification, generating a high-confidence annotation of *D. simulans* more conserved TEs.

**Module 2:** The annotated TEs from Module 1 were masked in the genome assemblies using bedtools maskfasta. The masked genome was submitted to RepeatMasker (-nolow) to annotate the remaining more divergent TE fragments using the combined TE library. The same steps from Module 1 were performed to obtain a high-confidence annotation of the remaining TE fragments. Finally, the annotations from the two modules were combined.

## 2.2 - Code availability

Main and helper scripts required by this workflow:

-   [te_annotation_general.sh](scripts/te_annotation_general.sh)

-   [rename_headers_numbers.py](scrips/rename_headers_numbers.py)

-   [recode_headers_from_numbers.py](recode_headers_from_numbers.py)

### 2.2.1 Dependencies

The following programs must be installed and available in the system PATH:

-   RepeatMasker
-   BLAST+
-   BEDtools
-   Python3

### 2.2.2 Required inputs

-f: text file with one sample ID per line

-p: directory containing genome assemblies

-g: genome filename suffix

-t: curated D. simulans TE library [Manual_Curated_Dsim_TE_lib.fasta](Manual_Curated_Dsim_TE_lib.fasta)

-c: combined Drosophila TE library [General_Drosophila_TE_final_lib.fasta](General_Drosophila_TE_final_lib.fasta)

### 2.2.3 Usage

```{bash eval =FALSE}
bash te_annotation_general.sh \
  -f <sample_list.txt> \
  -p /path/to/genomes \
  -g _final_assembly.fasta \
  -t Manual_Curated_Dsim_TE_lib.fasta \
  -c General_Drosophila_TE_final_lib.fasta
```

### 2.2.4 Main outputs

For each sample, the script generates:

**Primary output files**

`${sample}_final_TE_annotation.bed`: combined TE annotation from both modules

`${sample}_final_TE_annotation_morestringent.bed`: filtered annotation retaining fragments longer than 79 bp

**Intermediate files**

The pipeline also produces multiple intermediate FASTA, BED, BLAST, and RepeatMasker files, as well as the directories `RM\_${sample}_nolow/` and `RM_${sample}_lib_others/`, which are used during TE filtering, validation, and re-annotation.

# 3 - Filtering and standardisation of TE annotations

Because the annotation included the use of a general *Drosophila* TE library from multiple sources, the classification labels are not fully standardised across elements. As a result, different naming conventions appear for the same TE categories(e.g. LINE elements may be annotated either as LINE or as Non-LTR_retrotransposon).

To ensure consistency in the downstream analysis, TE classification names were standardised after annotation. This step was necessary in order to correctly summarise the TE composition by superfamily, since inconsistent labels would otherwise artificially split elements belonging to the same TE group into separate categories.

Moreover, the annotation was manually inspected to remove matches that were unlikely to represent genuine TE sequences. Although some weak matches contribute little to the overall genome coverage, retaining them would artificially introduce additional TE superfamilies into the analysis. Therefore, candidate elements were inspected individually and filtered based on sequence similarity, alignment quality, and biological plausibility.

## 3.1 - TE entries removed from the dataset

The following TE annotations were excluded because they correspond to very short genomic fragments and the sequence similarity is insufficient to confidently classify them as TE remnants:

`Other`,`MINIME_DN`,`RTE-X`,`DNA/Zator`, `L2`, `DNA3-1_DF`,`DNA8-1_DEl`

## 3.2 - Standardisation of TE superfamily names in the final annotation

To ensure consistency across TE classifications, several TE names were standardised:

`DM_Micropia2#LTR#Gypsy` → `DM_Micropia2#LTR/Gypsy`

`DM_HMS-Beagle/LTR/Gypsy` → `DM_HMS-Beagle#LTR/Gypsy`

`Non-LTR_retrotransposon` → `LINE`

`LTR/Gyspy` → `LTR/Gypsy`

`LTR/gypsy` → `LTR/Gypsy`

`RC/Helitron` → `DNA/Helitron`

`DNA/hAT-Ac, DNA/hAT-hobo, DNA/hAT-Pegasus` → `DNA/hAT`

`DNA/TcMar-Tc1, DNA/Tc1-Mariner, DNA/TcMar-Mariner` → `DNA/Tc1`

`LINE_484/Jockey, LINE/LINE` → `LINE/Jockey`

`DNA/MULE-NOF` → `DNA/MuDR`

`DNA/PIF-Harbinger` → `DNA/Harbinger`

`DNA/CMC-Transib` → `DNA/Transib`
