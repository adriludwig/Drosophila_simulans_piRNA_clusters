# Genome assembly

This repository documents the pipeline used to assemble the 22 Drosophila simulans genomes.

## 1 - Data availability

The sequencing data and assemblies generated in this study are available in public repositories:

Nanopore long reads: NCBI [**xxx**]{.underline}

Illumina reads: NCBI [**xxx**]{.underline}

Final genome assemblies: [**link**]

## 2 - Pipeline overview

Base-calling of Nanopore reads was performed using Guppy v3.2.10 (Oxford Nanopore Technologies PLC) or Dorado v0.7.2 (Oxford Nanopore Technologies PLC). Adapters were trimmed with Porechop v0.2.4 (<https://github.com/rrwick/Porechop>) and Filtlong v0.2.1 (<https://github.com/rrwick/Filtlong>) was used for quality trimming (–min_length 500 –min_mean_q 10). NanoPlot v1.42.0 ([De Coster et al., 2018](https://academic.oup.com/bioinformatics/article/34/15/2666/4934939)) was used to access the statistical summary of reads. De novo reads assembly was performed using Flye v2.9.2 ([Kolmogorov et al., 2019](%3Chttps://www.nature.com/articles/s41587-019-0072-8)) with two polishing iterations. Additional polishing was performed using the nanopore reads with three rounds of Racon v1.6.0 ([Vaser et al., 2017](https://genome.cshlp.org/content/27/5/737.long))%20and%20Medaka%20v1.11.3%20(%3Chttps://github.com/nanoporetech/medaka). For Lps5, Hin17, Cro18 and SGA27, Pilon v1.24 ([Walker et al., 2014](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0112963)) was used to polish the genomes with Illumina reads (1 round). Craq v1.0.9-0 ([Li et al., 2023](https://www.nature.com/articles/s41467-023-42336-w)) was used to detect structural errors and break chimeric contigs. Kraken2 v2.1.3 ([Wood et al., 2019](https://link.springer.com/article/10.1186/s13059-019-1891-0)) was then used to identify and eliminate contaminants from contaminant sources. For scaffolding the new assemblies, we used RagTag v2.1.0 (<https://github.com/malonge/RagTag>) and Dsim_3.1 (GCF_016746395.2) genome assembly as a reference. Several tools were used to assess the quality of genome assemblies: 1) BUSCO v5.7.1 ([Seppey et al., 2019](https://link.springer.com/protocol/10.1007/978-1-4939-9173-0_14)) was used to assess the gene content using the diptera_odb10 database. 2) Quast v5.2.0 ([Gurevich et al., 2013](https://pmc.ncbi.nlm.nih.gov/articles/PMC3624806/)) computed the basic assembly statistics such as the total length, contig numbers and N50. 3) Craq v1.0.9-0 was used to quantify clip-based regional and structural errors.

## 2.1 - Base-calling

Base-calling of Nanopore reads was performed using Guppy v3.2.10 using default parameters (Oxford Nanopore Technologies PLC) or Dorado v0.7.2 (Oxford Nanopore Technologies PLC).

Dorado parameters:

```{bash eval = FALSE}
dorado basecaller --verbose -x cuda:0 -b 64 --emit-fastq -c 5000 fast <pod5_directory>
```

## **2.2 - Trim adapters with Porechop**

Adapters were trimmed with Porechop v0.2.4 (<https://github.com/rrwick/Porechop>).

##### Input:

-Nanopore reads: `<reads.fastq>`

```{bash eval = FALSE}
porechop -i <reads.fastq> -t 72 -v 2 --format fastq -o reads_trimmed.fastq --discard_middle
```

## **2.3 - Quality filtering with Filtlong**

Filtlong v0.2.1 (<https://github.com/rrwick/Filtlong>) was used for quality filtering using a minimum read length of 500 bp and a minimum mean quality score of Q10.

```{bash eval = FALSE}
filtlong --min_length 500 --min_mean_q 10 <trimmed.fastq> > <trimmed_filtered.fastq>
```

## **2.4 - Genome assembly with Flye**

De novo genome assembly of Nanopore reads was performed using Flye v2.9.2 ([Kolmogorov et al., 2019](https://www.nature.com/articles/s41587-019-0072-8)) with two polishing iterations.

##### Input:

-Nanopore reads after trimming and filtering steps: `<trimmed_filt.fastq>`

```{bash eval = FALSE}
flye --nano-raw <trimmed_filt.fastq> --genome-size 150M --iterations 2 --threads 72 -o ./flye_out 2>&1 | tee flye.output.log
```

## **2.5 - Genome polishing with Racon**

Genome polishing was performed using Racon v1.6.0 ([Vaser et al., 2017](https://pubmed.ncbi.nlm.nih.gov/28100585/)) with Nanopore reads. Reads were aligned to the assembly using minimap2 (<https://github.com/lh3/minimap2>), and polishing was performed in three iterative rounds.

##### Input:

Nanopore reads: `<reads.fastq>`

Flye assembly: `<assembly.fasta>`

```{bash eval = FALSE}
minimap2 -t 72 -x map-ont <assembly.fasta> <reads.fastq> > alignment.paf

racon -t 72 -m 8 -x -6 -g -8 -w 500 <reads.fastq> alignment.paf <assembly.fasta> > polished.fasta
```

## **2.6 - Genome polishing with Medaka**

Medaka v1.11.3 (<https://github.com/nanoporetech/medaka> ) was used for final polishing. The genome assembly after medaka is consensus.fasta

##### Input:

-Nanopore reads: `<reads_trimmed.fastq>`

-Assembly after Racon polishing: `<polished_round3.fasta>`

```{bash eval = FALSE}
medaka_consensus -i reads_trimmed.fastq -d polished_round3.fasta -o medaka_out_folder -t 72 -m r941_min_high_g344
```

## 2.7 - Genome polishing with Pilon

For Lps5, Hin17, Cro18 and SGA27, Pilon v1.24 (Walker et al., 2014) was used to polish the genomes with Illumina paired-end reads (1 round). Short reads were aligned to the assembly using BWA-mem, and the resulting alignments were sorted and indexed using samtools. Genome polishing was then performed with Pilon.

##### Input:

-Assembly after Medaka polishing: `<consensus.fasta>`

-Illumina paired-end reads: `<reads_R1.fastq.gz>` and `<reads_R2.fastq.gz>`

```{bash eval = FALSE}
bwa index consensus.fasta
bwa mem -t 72 consensus.fasta <reads_R1.fastq.gz> <reads_R2.fastq.gz> > alignment.sam
samtools view -@ 72 -b alignment.sam > alignment.bam
samtools sort -@ 72 alignment.bam -o alignment_sorted.bam
samtools index alignment_sorted.bam

java -Xmx600G -jar pilon-1.24.jar --genome <consensus.fasta> \
--bam alignment_sorted.bam --output <pilon_round1 --changes --vcf
```

## 2.8 - CRAQ contig breaking

Craq v1.0.9-0 ([Li et al., 2023](https://www.nature.com/articles/s41467-023-42336-w)) was used to detect structural errors and split potentially misassembled contigs.

##### Input:

\- Genome assembly after polishing:`<assembly.fasta>`\
- Nanopore reads: `<reads.fastq>`

```{bash eval = FALSE}
perl craq -g <assembly.fasta> -sms <reads.fastq> -b T -t <threads> -pl T --output_dir <craq_output>
```

## 2.9 Contamination screening with Kraken2

Kraken2 was used to taxonomically classify assembled contigs in order to detect potential contamination. Contigs assigned to taxa outside the genus *Drosophila* were considered contaminants and were removed from the assembly.

##### Input:

\- Genome assembly after contig splitting by Craq: `<assembly.fasta>`\
- Kraken database: `<kraken2_database>`

```{bash eval = FALSE}
kraken2 --db <kraken2_database> --threads 48 --use-names --report kraken2_report.txt --output kraken2_output.txt --classified-out kraken2.classified --unclassified-out kraken2.unclassified <assembly.fasta>
```

## 2.10 - Reference-guided scaffolding with RagTag

Contigs were scaffolded against the *Drosophila simulans* reference genome Prin_Dsim_3.1 (GCF_016746395.2) using RagTag v2.1.0

##### Input:

\- Genome assembly: `<assembly.fasta>`\
- Reference genome: `<reference.fasta`\>

```{bash eval = FALSE}
ragtag.py correct reference.fasta consensus.fasta -o ragtag_output
ragtag.py scaffold reference.fasta ragtag_output/ragtag.correct.fasta -o ragtag_output
```
