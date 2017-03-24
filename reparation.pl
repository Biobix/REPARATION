#!/usr/bin/perl -w

#####################################
##	REPARARTION: Ribosome Profiling Assisted (Re-) Annotation of Bacterial genomes
##
##	Copyright (C) 2017 Elvis Ndah
##
##	This program is free software: you can redistribute it and/or modify
##	it under the terms of the GNU General Public License as published by
##	the Free Software Foundation, either version 3 of the License, or
##	(at your option) any later version.
##	
##	This program is distributed in the hope that it will be useful,
##	but WITHOUT ANY WARRANTY; without even the implied warranty of
##	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##	GNU General Public License for more details.
##
##	You should have received a copy of the GNU General Public License
##	along with this program.  If not, see <http://www.gnu.org/licenses/>.
##
##	contact: elvis.ndah@gmail.com
#####################################

use strict;
use warnings;
use diagnostics;

use Getopt::Long;
use Cwd;
use File::stat;


########################
##	Usage: perl ./reparation.pl -a alignment_file -sdir script_directory -g genome_fasta_file -db reference_fasta_protein_db
########################


#----------------------------------------------------
#			VARIBLES
#----------------------------------------------------

# Mandatory variables
my $genome;								# Prokaryotic genome file in fasta format
my $Ribo_alignment;						# Ribosome profiling alignment file [sam formate]
my $blastdb;							# protein blast database in fasta format
my $script_dir = "scripts";				# Directory where the script files are stored (defaults to the current directory)

# Input variables
my $threads = 4;						# # number of threads used by the USEARCH tool
my $workdir;							# working directory to store files (defaults to current directoy)
my $experiment;							# Experiment name
my $gtf;								# Genome annotation file [if avialable]
my $positive_set;						# Tab delimited file of positive examples (format [orf_id	strand], orf_id=region:start-stop) 
my $occupancy = 3;						# p site for the reads (3 = 3 prime end of read  and 5 = 5 prime end of read)
my $min_read_len = 22;					# Minimum RPF read length
my $max_read_len = 40;					# Maximum RPF read length
my $MINORF = 30;						# Minimum ORF length
my $MINREAD = 3;						# Only ORFs with at least this number of RPF reads within the start regions
my $OFFSET_START= 45;					# Offset at the start of the ORF
my $OFFSET_STOP	= 21;					# Offset at stop of the ORF
my $OFFSET_SD = 5;						# distance uptream of start codon to start search for SD sequence and position
my $SEED = "GGAGG";						# The seed shine dalgano sequence
my $USESD = 1;							# Flag to determine if RBS energy is included in the predictions [1= use RBS, 0=do not use RBS]
my $identity = 0.75;					# identity threshold for comparative psotive set selection
my $evalue = 1e-5;						# e value threshold for comparative psotive set selection
my $pgm = 1;							# Program to generate positive set prodigal=1, glimmer=2, provided_by_user=3
my $start_codons = "ATG,GTG,TTG";		# Comma seperated list of possible start codons
my $start_codon_nset = "CTG";			# Start codon for the negative set
my $start_codon_pset = "ATG,GTG,TTG";	# Start codon for the positive set. Defualts to start codons set


# track processing time
my $startRun = time();

# Get command line arguments
GetOptions(
	'g=s'=>\$genome,
	'gtf=s'=>\$gtf,
	'a=s'=>\$Ribo_alignment,
	'sdir=s'=>\$script_dir,
	'wdir=s'=>\$workdir,
	'en=s'=>\$experiment,
	'p=i'=>\$occupancy,
	'mn=i'=>\$min_read_len,
	'mx=i'=>\$max_read_len,
	'db=s'=>\$blastdb,
	'mo=i'=>\$MINORF,
	'mr=i'=>\$MINREAD,
	'ost=i'=>\$OFFSET_START,
	'osp=i'=>\$OFFSET_STOP,
	'osd=i'=>\$OFFSET_SD,
	'seed=s'=>\$SEED,
	'sd=i'=>\$USESD,
	'id=f'=>\$identity,
	'ev=i'=>\$evalue,
	'pg=i'=>\$pgm,
	'pset=s'=>\$positive_set,
	'cdn=s'=>\$start_codons,
	'ncdn=s'=>\$start_codon_nset,
	'pcdn=s'=>\$start_codon_pset
);


# Check for mandatory variables;
my %params = (
	g => $genome,
	b => $Ribo_alignment,
	sdir => $script_dir,
	db => $blastdb
);

#check if mandatory parameters are properly initialized
my @invalid = grep uninitialized_param($params{$_}), keys %params;	
die "Not properly initialized: @invalid\n" if @invalid;


# ensure start codons are uppercase
$start_codons = uc($start_codons);
$start_codon_nset = uc($start_codon_nset);
unless($start_codon_pset) {$start_codon_pset = $start_codons}

$experiment = ($experiment) ? $experiment."_": "";

# check if script directoryis properly initialized and contains all scripts
if ($script_dir) {
	$script_dir = (substr($script_dir, -1) eq '/') ? $script_dir: $script_dir."/"; # check and add forward slash

	unless (-e $script_dir."Ribo_seq_occupancy.py") {
		print "Script 'Ribo_seq_occupancy.py' not found in script directory.\nEnsure the directory  '$script_dir' contains all required file (see readme).\n";
		exit(1);
	}

	unless (-e $script_dir."generate_all_ORFs.pl") {
		print "Script 'generate_all_ORFs.pl' not found in script directory.\nEnsure the directory '$script_dir' contains all required file (see readme).\n";
		exit(1);
	}

	unless (-e $script_dir."profiles.pl") {
		print "Script 'profiles.pl' not found in script directory.\nEnsure the directory '$script_dir' contains all required file (see readme).\n";
		exit(1);
	}

	unless (-e $script_dir."plot_profile.R") {
		print "Script 'plot_profile.R' not found in script directory.\nEnsure the directory '$script_dir' contains all required file (see readme).\n";
		exit(1);
	}

	unless (-e $script_dir."Random_Forest.R") {
		print "Script 'Random_Forest.R' not found in script directory.\nEnsure the directory '$script_dir' contains all required file (see readme).\n";
		exit(1);
	}

	unless (-e $script_dir."post_processing.pl") {
		print "Script 'post_processing.pl' not found in script directory.\nEnsure the directory '$script_dir' contains all required file (see readme).\n";
		exit(1);
	}

	# Check prerequisits
	if ($pgm == 1) {
		check_if_pgm_exist('prodigal');
	} elsif ($pgm == 2) {
		check_if_pgm_exist('glimmer3');
	}
} else {
	print "Script directory not properly initialized. Please ensure the script directory is properly defined\n";
	exit(1);
}

# Check if working directory exist
my ($work_dir, $tmp_dir) = check_working_dir($workdir);

# Check prerequisits
if ($pgm == 1) {
	$positive_set = $work_dir."positive_set.txt";
} elsif ($pgm == 2) {
	$positive_set = $work_dir."positive_set.txt";
} elsif ($pgm == 3) {
	if ($positive_set) {
		if (-e $positive_set) {
			system("mv $positive_set ".$work_dir."positive_set.txt");
		} else {
			print "Positive '".$positive_set."'set is not a file.\nPlease provide a file of positive examples.\n";
			exit(1);
		}
	} else {
		print "Positive set not defined. If flag pg = 3 then provide an appropriate file for flag pset\n";
		exit(1);
	}
}

# print input variables
if ($identity <= 0 or $identity > 1) {
	print "-id can only take a value greater than zero and less than or equal to 1\n";
	exit;
}

# Generate occupancy file
print "Generating ribosome occupancy file..\n";
my $bedgraphS = $work_dir.$experiment."Ribo-seq_Sense_".$occupancy."prime";
my $bedgraphAS = $work_dir.$experiment."Ribo-seq_AntiSense_".$occupancy."prime";
my $occupancyFile = $work_dir.$experiment."Ribo-seq_".$occupancy."_occupancy.txt";
my $cmd_occupancy = "python ".$script_dir."Ribo_seq_occupancy.py $Ribo_alignment $occupancy $min_read_len $max_read_len $bedgraphS $bedgraphAS $occupancyFile";
print "$cmd_occupancy\n";
system($cmd_occupancy);
print "Done.\n\n";


# Generate all possible ORFs
print "Generate all possible ORFs...\n";
my $codons = $start_codons.",".$start_codon_nset;	# combine the start codons
my $ORF_file = $work_dir."tmp/all_valid_ORFs.txt";
my $cmd_orf_gen = "perl ".$script_dir."generate_all_ORFs.pl $genome $blastdb $occupancyFile $positive_set $ORF_file $MINORF $OFFSET_START $OFFSET_STOP $OFFSET_SD $SEED $identity $evalue $codons $start_codon_pset $pgm $work_dir $script_dir $threads";
print "$cmd_orf_gen\n";
system($cmd_orf_gen);
print "Done.\n\n";

# Generate metagenic profile
print "Meta-genic plots\n";
my $cmd_meta = "perl ".$script_dir."profiles.pl $positive_set $occupancyFile $MINREAD $work_dir $script_dir";
print "$cmd_meta\n";
system($cmd_meta);
print "Done.\n\n";

# ORF prediction
print "Performing ORF prediction analysis..\n";
my $RF_command = "Rscript ".$script_dir."Random_Forest.R $ORF_file $positive_set $work_dir $start_codons $start_codon_nset $USESD";
print "$RF_command\n";
system($RF_command);
print "Done.\n\n";

# Post processing
print "Cleaning up RF predictions..\n";
my $threshold = $work_dir."tmp/threshold.txt";
my $RF_prediction = $work_dir."tmp/RF_predicted_ORFs.txt";
my $output_prefix = $work_dir.$experiment."Predicted_ORFs";
my $processing_cmd = "perl ".$script_dir."post_processing.pl $RF_prediction $genome $occupancyFile $output_prefix $threshold $OFFSET_START $MINREAD $gtf";
print "$processing_cmd\n";
system($processing_cmd);
print "Done.\n\n";


timer($startRun);	# Get total Run time



##################
##	SUBS
##################

sub check_if_pgm_exist {

	my $pgm = $_[0];

	my $search = `which $pgm 2>&1`;
	chomp($search);
	if ($search =~ /^which: no/) {
		if ($pgm eq "prodigal") {
			unless (-e $script_dir."bin/prodigal") {	# if prodigal not install
				print "Could not locate ' $pgm '. Please ensure it is installed and executable\n";
				exit(1);
			}
		} elsif ($pgm eq "glimmer3") {
			unless (-e $script_dir."bin/glimmer/glimmer3" and -e $script_dir."bin/glimmer/bin/build-icm") {
				print "Could not locate '$pgm'. Please ensure it is installed and executable\n";
				exit(1);
			}
		}
	} 
}

sub uninitialized_param {
	my ($v) = @_;
	not ( defined($v) and length $v );
}


sub check_working_dir {

	my $work_dir = $_[0];
	my $tmp_dir;

	my @months = qw( Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec );
	my @days = qw(Sun Mon Tue Wed Thu Fri Sat Sun);
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();

	if ($work_dir) {
		$tmp_dir = $work_dir."/tmp";
		if (!-d $work_dir) {
			system("mkdir -p $work_dir" or die "Couldn't create '$work_dir.\n");
			system("mkdir -p $tmp_dir" or die "Couldn't create '$tmp_dir'.\n");
		} else {
			system("rm -rf $work_dir" or die "Can delete '$work_dir': $!\n");
			system("mkdir -p $work_dir" or die "Couldn't create '$tmp_dir'.\n");
			system("mkdir -p $tmp_dir" or die "Couldn't create '$tmp_dir'.\n");
		}
	} else {
		$work_dir = getcwd();
		$work_dir = $work_dir."/".$experiment."reparation_".$months[$mon].$mday;
		$tmp_dir = $work_dir."/tmp";

		if (!-d $work_dir) {
			system("mkdir -p $work_dir");
			system("mkdir -p $tmp_dir");
		} 
	}

	$work_dir = $work_dir."/";
	$tmp_dir = $tmp_dir."/";
	return $work_dir, $tmp_dir;
}


sub timer {

	my $startRun = shift;

	my $endRun 	= time();
	my $runTime = $endRun - $startRun;

	printf("\nTotal running time: %02d:%02d:%02d\n\n", int($runTime / 3600), int(($runTime  % 3600) / 60), int($runTime % 60));

}

