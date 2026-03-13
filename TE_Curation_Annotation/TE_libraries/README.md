# TE libraries

TE libraries used in this work:

## 1 - Manual_Curated_Dsim_TE_lib.fasta
[Manual_Curated_Dsim_TE_lib.fasta](./Manual_Curated_Dsim_TE_lib.fasta)
*Drosophila simulans* curated library obtained in this work. Please
check **xxxxx** for information on how it was obtained and
[Manual_curated_Dsim_lib_Info.xlsx](./Manual_curated_Dsim_lib_Info.xlsx) for more information on how the sequences were named.

## 2 - Dmel_MCTE_edited.fasta
[Dmel_MCTE_edited.fasta](Dmel_MCTE_edited.fasta)
A modified version of the *Drosophila melanogaster* TE library from ([Rech
et al., 2022](https://www.nature.com/articles/s41467-022-29518-8)).
Several sequences were edited or deleted. The library was also
supplemented with 5 additional sequences from Repbase23.11 and the BDGP
TE library. All changes made are described in [Dmel_MCTE_info.xlsx](./Dmel_MCTE_info.xlsx)

## 3 - DrosREPlibrary_Chakraborty_filtered.fasta
[DrosREPlibrary_Chakraborty_filtered.fasta](DrosREPlibrary_Chakraborty_filtered.fasta)
A library of 2,316 TE sequences from diverse *Drosophila* species was
obtained from ([Chakraborty et al.,
2021](https://genome.cshlp.org/content/31/3/380.long)), inspected to
eliminate satellites, rDNA, Low-Complexity sequences, and clear
chimaeras or non-TEs. P-element from *D. simulans* was added as the
library was used as a standalone library in another study.

Filterings:

-All sequences assigned as Satellite rDNA and Low_complexity were
eliminated

Extra filtering to eliminate chimaeras and not real TE families:
noCat_dsim_final-B-P113.23-Map3_novel_repeat#Other/de_novo;
noCat_dsim_final-B-P81.43-Map3_novel_repeat#Other/de_novo;
noCat_dsim_final-B-P139.100-Map3_novel_repeat#Other/de_novo;
noCat_dsech_final-B-P78.38-Map3_novel_repeat#Other/de_novo;
noCat_dsech_final-B-P110.14-Map3_novel_repeat#Other/de_novo;
NTS_DM#Other/IGS; TransibN2_DP#DNA/CMC-Transib

## 4- General_Drosophila_TE_final_lib.fasta
[General_Drosophila_TE_final_lib.fasta](General_Drosophila_TE_final_lib.fasta)
This is a redundant TE library resulting from the concatenation of the
three libraries above.
