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

use Bio::SeqIO;
use Getopt::Long;
use POSIX;

# Get command line arguments
my ($genome,$blastdb,$positive_set,$min_read_len,$max_read_len,$MINORF,$identity,$evalue,$pcodons,$pgm,$work_dir,$script_dir,$threads, $seedBYpass, $genetic_code) = @ARGV;

my $upstream = 50;
#----------------------------------------------------
#			EXECUTION
#----------------------------------------------------

my %translationHash = 
	(GCA => "A", GCG => "A", GCT => "A", GCC => "A", GCN => "A",
     TGC => "C", TGT => "C",
     GAT => "D", GAC => "D",
     GAA => "E", GAG => "E",
     TTT => "F", TTC => "F",
     GGA => "G", GGG => "G", GGC => "G", GGT => "G", GGN => "G",
     CAT => "H", CAC => "H",
     ATA => "I", ATT => "I", ATC => "I",
     AAA => "K", AAG => "K",
     CTA => "L", CTG => "L", CTT => "L", CTC => "L", CTN => "L", TTA => "L", TTG => "L",
     ATG => "M",
     AAT => "N", AAC => "N",
     CCA => "P", CCT => "P", CCG => "P", CCC => "P", CCN => "P",
     CAA => "Q", CAG => "Q",
     CGA => "R", CGG => "R", CGC => "R", CGT => "R", CGN => "R",
     AGA => "R", AGG => "R",
     TCA => "S", TCG => "S", TCC => "S", TCT => "S", TCN => "S",
     AGC => "S", AGT => "S",
     ACA => "T", ACG => "T", ACC => "T", ACT => "T", ACN => "T",
     GTA => "V", GTG => "V", GTC => "V", GTT => "V", GTN => "V",
     TGG => "W",
     TAT => "Y", TAC => "Y",
     TAG => "*", TAA => "*", TGA => "*");

# Stop codons
my $stop_codons = {TAG => "*", TAA => "*", TGA => "*"};

# valid positive set codons
my $positive_codons = {};
my @pcodons = split /,/, $pcodons;
foreach my $codon(@pcodons) {
	$positive_codons->{$codon} = 1;
}


# read genome fasta file
my $genomes = read_fasta($genome);


# generate positive set 
my $positive_set_gtf = $work_dir."/tmp/positive.gtf";
if ($pgm == 1) { 
	# postive set generated with prodigal
	print "Generating positive samples training set using Prodigal...\n";
	my $positive_fasta = generate_prodigal_ORFs();
	create_positive_file($positive_set,$positive_fasta);
	print "Set of Positive examples created.\n";

} elsif ($pgm == 2) { 
	# postive set generated with Glimmer
	print "Generating positive samples training set using glimmer...\n";
	my $positive_fasta = generate_glimmer_ORFs();
	create_positive_file($positive_set,$positive_fasta);
	print "Set of Positive examples created.\n";

}



##################
##	SUBS
##################

sub create_positive_file {

	my $positive_set = $_[0];
	my $positive_fasta = $_[1];

	my $ORF_positive = {};

	my $blast_rpt = $work_dir."tmp/positive_set.bls";
	my $search = `which usearch 2>&1`;
	chomp($search);
	if ($search =~ /^which: no/) {
		system($script_dir."/bin/usearch -ublast $positive_fasta -db $blastdb -id $identity -evalue $evalue -maxhits 1 -blast6out $blast_rpt -threads $threads 2> /dev/null");
	} else {
		system("usearch -ublast $positive_fasta -db $blastdb -id $identity -evalue $evalue -maxhits 1 -blast6out $blast_rpt -threads $threads 2> /dev/null");
	}

	# GTF file
	open (G, ">".$positive_set_gtf) or die "Error creating file $positive_set_gtf\n";
	print G "#!positive set gtf file\n";
	print G "#!positive set gtf file\n";

	open (OUT, ">".$positive_set) or die "Error creating file $positive_set\n";
	print OUT "orf_id\tstrand\tstart_codon\n";
	open (F, $blast_rpt) or die  "Error reading file: $blast_rpt\n";
	while (<F>) {
		next if (/^Query_id/);
		chomp $_;

		my @line = split '\t', $_;
		my $query_id = $line[0];
		my $target_id = $line[1];
		my $perc_ident = $line[2];

		$query_id =~ s/\s+//g;
		my ($region)= $query_id=~ /^(.*):/;
		my ($start)= $query_id =~ /:(\d+)-/;
		my ($stop)= $query_id =~ /-(\d+)#/;
		my ($strand)= $query_id =~ /#(.)\|/;
		my ($start_codon)= $query_id =~ /\|(.*)$/;

		my $length = $stop - $start + 1;

		if (exists $positive_codons->{$start_codon} and $perc_ident >= $identity and $length >= $MINORF) {

		# ensure only valid ORFs are in positive set
			($query_id) = $query_id =~ /^(.*)#/;
			print OUT "$query_id\t$strand\t$start_codon\n";
			$ORF_positive->{$query_id} = $strand;

			my $gene = ($strand eq "+") ? $region.":FAM+".$stop: $region.":FAM-".$start;
			my $start_gene = ($strand eq '+') ? $start - $upstream: $start - 3;
			my $stop_gene = ($strand eq '+') ? $stop + 3: $stop + $upstream;
			my $start_cdn_strt = ($strand eq '+') ? $start: $stop - 2;
			my $start_cdn_stp = ($strand eq '+') ? $start + 2: $stop;
			my $stop_cdn_strt = ($strand eq '+') ? $stop_gene - 2 : $start_gene;
			my $stop_cdn_stp = ($strand eq '+') ? $stop_gene: $start_gene + 2;

			print G "$region\tprodigal\tgene\t$start_gene\t$stop_gene\t.\t$strand\t.\tgene_id \"$gene\"; gene_version \"1\"; gene_name \"$gene\"; gene_source \"prodigal\"; gene_biotype \"protein_coding\";\n";
			print G "$region\tprodigal\ttranscript\t$start_gene\t$stop_gene\t.\t$strand\t.\tgene_id \"$gene\"; gene_version \"1\"; transcript_id \"$query_id\"; transcript_version \"1\"; gene_name \"$gene\"; gene_source \"prodigal\"; gene_biotype \"protein_coding\"; transcript_name \"$gene"."-1"."\"; transcript_source \"prodigal\"; transcript_biotype \"protein_coding\";\n";
			print G "$region\tprodigal\texon\t$start_gene\t$stop_gene\t.\t$strand\t.\tgene_id \"$gene\"; gene_version \"1\"; transcript_id \"$query_id\"; transcript_version \"1\"; exon_number \"1\"; gene_name \"$gene\"; gene_source \"prodigal\"; gene_biotype \"protein_coding\"; transcript_name \"$gene"."-1"."\"; transcript_source \"prodigal\"; transcript_biotype \"protein_coding\"; exon_id \"$gene"."-e1"."\"; exon_version \"1\";\n";
			print G "$region\tprodigal\tCDS\t$start\t$stop\t.\t$strand\t.\tgene_id \"$gene\"; gene_version \"1\"; transcript_id \"$query_id\"; transcript_version \"1\"; exon_number \"1\"; gene_name \"$gene\"; gene_source \"prodigal\"; gene_biotype \"protein_coding\"; transcript_name \"$gene-1\"; transcript_source \"prodigal\"; transcript_biotype \"protein_coding\"; protein_id \"$query_id"."-p1"."\"; protein_version \"1\";\n";
			print G "$region\tprodigal\tstart_codon\t$start_cdn_strt\t$start_cdn_stp\t.\t$strand\t.\tgene_id \"$gene\"; gene_version \"1\"; transcript_id \"$query_id\"; transcript_version \"1\"; exon_number \"1\"; gene_name \"$gene\"; gene_source \"prodigal\"; gene_biotype \"protein_coding\"; transcript_name \"$gene-1\"; transcript_source \"prodigal\"; transcript_biotype \"protein_coding\";\n";
			print G "$region\tprodigal\tstart_codon\t$stop_cdn_strt\t$stop_cdn_stp\t.\t$strand\t.\tgene_id \"$gene\"; gene_version \"1\"; transcript_id \"$query_id\"; transcript_version \"1\"; exon_number \"1\"; gene_name \"$gene\"; gene_source \"prodigal\"; gene_biotype \"protein_coding\"; transcript_name \"$gene-1\"; transcript_source \"prodigal\"; transcript_biotype \"protein_coding\";\n";

		}
	}
	close OUT;
	close G;

	return $ORF_positive;
}


# prodigal positive
sub generate_prodigal_ORFs {

	my $positive_gff = $work_dir.'tmp/positive_set.gff';

	my $search = `which prodigal 2>&1`;
	chomp($search);

    my $prodigal_cmd = "";
	if ($search =~ /^which: no/) {
		$prodigal_cmd = $script_dir."/bin/prodigal -i $genome -f gff -o $positive_gff -g $genetic_code -q";
	} else {
		$prodigal_cmd = "prodigal -i $genome -f gff -o $positive_gff  -g $genetic_code -q";
	}

    if ($seedBYpass eq "Y") {$prodigal_cmd = $prodigal_cmd." -n"}
    system($prodigal_cmd);

	my $positive_fasta = $work_dir.'tmp/positive_set.fasta';
	open (OUT, ">".$positive_fasta) or die "error creating file $positive_fasta\n";
	open (IN, $positive_gff) or die "Error reading file: $positive_gff\n";
	while (<IN>) {

		next if (/^#/);
		chomp $_;
		my @line = split '\t', $_;
		my $region = $line[0];
		my $strand = $line[6];

		my $start = ($strand eq '+') ? $line[3]: $line[3] + 3;
		my $stop = ($strand eq '+') ? $line[4] - 3: $line[4];
		my $orf = $region.":".$start."-".$stop;
		
		my ($aa_seq, $start_codon) = translate($start,$stop,$strand,$region);
        next unless ($aa_seq);

		print OUT ">$orf"."#".$strand."|".$start_codon."\n$aa_seq\n";

	}
	close IN;
	close OUT;

	return $positive_fasta;
}


# generate using glimmer
sub generate_glimmer_ORFs {

	# Run glimmer
	my $search = `which glimmer3 2>&1`;
	chomp($search);

	my $glimmer_icm = $work_dir."/tmp/glimmer.icm";
	my $positive_predict = $work_dir."/tmp/positive_set";
	my $start_cdns = lc($pcodons);

    #print $script_dir."/bin/glimmer/glimmer3 $genome $glimmer_icm $positive_predict -A $start_cdns -g $MINORF -o 300\n"; exit;

	if ($search =~ /^which: no/) {
		system($script_dir."/bin/glimmer/build-icm $glimmer_icm< $genome");
		system($script_dir."/bin/glimmer/glimmer3 $genome $glimmer_icm $positive_predict -A $start_cdns -g $MINORF -o 300");
	} else {
		system("build-icm $glimmer_icm< $genome");
		system("glimmer3 $genome $glimmer_icm $positive_predict -A $start_cdns -g $MINORF -o 300");
	}

	my $positive_fasta = $work_dir."/tmp/positive_set.fasta";
	open (OUT, ">".$positive_fasta) or die "error creating file $positive_fasta\n";
	open (IN, $positive_predict.".predict") or die  "Error reading file: ".$positive_predict.".predict\n";
	my $region;
	while (<IN>) {
		chomp $_;
		if (/^>/) {
			$_ =~ s/>//;
			$region = (split ' ', $_)[0];
		} else {
			my @line = split '\s+', $_;
			my $strand = (split '', $line[3])[0];

			my $start = ($strand eq '+') ? $line[1]: $line[2] + 3;
			my $stop = ($strand eq '+') ? $line[2] - 3: $line[1];
			my $orf = $region.":".$start."-".$stop;

			my ($aa_seq, $start_codon) = translate($start,$stop,$strand,$region);

            next unless ($aa_seq);

			print OUT ">$orf"."#".$strand."|".$start_codon."\n$aa_seq\n";
		}
	}
	close OUT;
	close IN;

	return $positive_fasta;
}


sub translate {

	my $start = $_[0];
	my $stop = $_[1];
	my $strand = $_[2];
	my $region = $_[3];

	my $length = $stop - $start + 1;
	my $dna_seq = ($strand eq '+') ? substr($genomes->{$region}, $start - 1,$length) : revdnacomp(substr($genomes->{$region}, $start-1,$length));

    unless ($dna_seq) {return "", ""}

	my $start_codon = uc(substr($dna_seq, 0, 3));
	my $aa_seq = "";
	for (my $i = 0; $i <= (length($dna_seq) - 3); $i = $i + 3) {
		my $codon = uc(substr($dna_seq, $i, 3));
		last if ($codon eq '*');
        my $aa = $translationHash{$codon};

        # if amino acid doesn't exist then 
        unless ($aa) {$aa_seq = ""; last}
		$aa_seq = $aa_seq.$aa;
	}

	return $aa_seq, $start_codon;

}


sub revdnacomp {
    my $dna = shift;
    my $revcomp = reverse($dna);
    $revcomp =~ tr/ACGTacgt/TGCAtgca/;
    return $revcomp;
}


sub read_fasta {
	my $file = $_[0];
	my $cdna = {};

	my $in  = Bio::SeqIO->new(-file => $file, -format => "fasta");
	while(my $seqs = $in->next_seq) {
		my $id  = $seqs->display_id;	
		my $seq = $seqs->seq;
		my $desc = $seqs->desc;
		$cdna->{$id} = uc($seq);
	}
	return $cdna;
}
