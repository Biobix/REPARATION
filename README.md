# REPARATION
REPARATION (Ribosome Profiling Assisted (Re-)Annotation of Bacterial genomes) is a pipeline that uses ribosome profiling data for a de novo open reading frame delineation in prokaryotic (bacterial) genomes.

The software was created at VIB-UGent Center for Medical Biotechnology and Lab of Bioinformatics and Computational Genomics (Biobix), University of Gent, belgium.

The REPARATION workflow consist of:

1) Genome wide ribosome profiling profile based on choosen occupancy position.

2) Traverse the genome and generate all possible open reading frames.
	- Generate traning set for random forest model.
	
3) Train a random forest to learn and use it for ORF prediction.

4) Start Site selection.


# Requirements
The REPARATION software is primarily written in perl, python2.7 and R. For each dataset the complete REPARATION workflow (from the sam files to final results) took under 45 minutes.
It is also dependant on a set of tool binaries which should all be installed on your system before the pipeline can execute all of its commands.

	Perl packages:
	--------------
		- Perl + v.5.10
		- Getopt::Long
		- Cwd	
		- BioPerl
		- POSIX
		
	R packages:
	--------------
		- ggplot2
		- randomForest
		- ROCR
		- SiZer

	Tool binaries:
	--------------
		- prodigal[*] https://github.com/hyattpd/Prodigal
		- USEARCH[*] (http://www.drive5.com/usearch/download.html)
		- R (http://www.r-project.org/)
		
		[*] The usearchX.Y.Z and prodigalX.Y.Z executatble should be renamed to "usearch" and "prodigal" respectively.
	
	The tool binary paths should be included in the $PATH variable or within the script/bin folder of the tool.

	
# Install
The tool does not rerquire any special installation outside the prerequisites. 

# Usage

Usage: ./reparation.pl -b Ribo_samfile -g genome_DNA_fasta_file -sdir scripts_directory -db Curated_bacterial_protein_sequences

The tool requires 4 mandatory input varaibles. A complete list of optional varibles are within the reparation.pl script.

NB: the Ribo-seq must be aligned to the same genome to be analyzed!


