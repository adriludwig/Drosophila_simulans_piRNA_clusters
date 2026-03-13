# TE curation and annotation

# *1. D. simulans* TE library curation

## **1.1 - Data availability**

The *Drosophila simulans* curated TE library obtained in this work is available at [Manual_Curated_Dsim_TE_lib.fasta](TE_libraries/Manual_Curated_Dsim_TE_lib.fasta)

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
