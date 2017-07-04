# IntronGetting

IntronGetting: Python scripts for getting ortholog EPIC sequences for primer design or phylogenomic analysis.

The input requierements for IntronGetting include:

(1): The absolute path of a folder which contains exon sequences of all species. 
For each file in the folder, the name of file consisted of genus and species, linked by "_", and end with "_exon". (e.g. Homo_sapiens_exon).
For each exon sequences, the head information format should be: Gene stable ID|Gene name|Chromosome/scaffold name|Exon stable ID|Genomic coding start|Genomic coding end|strand
																											
(2): The absolute path of a folder which contains genome sequences of all species. 
For each file in the folder, the name of file should be with a .fa suffix, and start with genus and species linked by "_", and followed by a dot. (e.g. Homo_sapiens.GRCh38.dna_sm.toplevel.fa)
																									
(3): The name of reference for two-way BLAST analysis. The name of genus and species, linked by "_". (e.g. Homo_sapiens)
																									
(4): Cpu_threads. Describing how many threads you would like to use in the processing.
																														

Usage: python get_EPIC.py exon_dir genome_dir reference_name cpu_threads
Note:  all directories should be provided with its absolute path.


The processing includes seven steps:

(1): Sort exons of each species according to their genomic position.

(2): Find orthologous exons using mutual best hit (MBH) in BLAST, selecting one species as the reference species. e.g. we selected Homo_sapiens as reference. Note: Add BLAST to your PATH environment variable.

(3): Covert BLAST results to tabular and filter low quality results. You can set coverage and identity cutoff to remove low quality result by changing EPIC.py.

(4): Analyze result files and find two-way BLAST ortholog exons. 

(5): Find ortholog exons across species.

(6): Filter candidates by taxon coverage, intron_length, intron/exon length ratio.

(7): Retrieve EPIC sequences. The step requires large memory as it recruits genome sequences.

##Disclaimer:

This software is experimental, in active development and comes without warranty. More detailed documentation is in preparation.
IntronGetting scripts were developed and tested using python 2.7 on Linux (Fedora) and Windows.

##Citation

Chen MY, Liang D, Zhang P. 2017. Phylogenomic Resolution of the Phylogeny of Laurasiatherian Mammals: Exploring Phylogenetic Signals within Coding and Noncoding Sequences.
