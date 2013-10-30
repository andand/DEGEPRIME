#!/usr/bin/perl -w


=head1 NAME

TrimAlignment.pl - Trims and formats a multiple sequence alignment for use as input to DegePrime.pl 

=head1 USAGE

perl TrimAlignment.pl -i ALIGNMENT_FILE -o OUTPUT_FILE [-min MIN_OCCUPANCY] [-ref REF_SEQUENCE] [-trailgap] [-h]

Specify MIN_OCCUPANCY if you want to remove columns of low occupancy, or REF_SEQUENCE if you want to remove columns that are not occupied by the reference sequence.

The output alignment file will contain upper and lower case nucleotide letters, the lower case letters indicate that the succeding column(s) have been removed. (If the  input alignment file contains any lower case letters these will first be changed to upper case.) Any U:s in the input file will also be replaced by T:s in the output.

If neither MIN_OCCUPANCY nor REF_SEQUENCE is specified, all columns will be kept (and all nucleotide letters will be upper case). We recommend that you process your alignment file like this even if you don't wish to remove columns, in order to remove unwanted lower case letters.



=head1 POSITIONAL ARGUMENTS

-i ALIGNMENT_FILE			Specify alignment file (fasta format)

-o OUTPUT_FILE				Specify output file name (overwrites existing file with same name)


=head1 OPIONAL ARGUMENTS

-min MIN_OCCUPANCY			Specify minimum fraction of sequences (in the range 0 to 1) with a nucleotide in a column for it to be included, default 0

-ref REF_SEQUENCE			Specify id of a reference sequence in the alignment file. Only columns where the reference has a nucleotide will be included

-trailgap					Indicate that you want to count trailing sequences as gaps, default mode is to not do so 

-h							Print this help message


=cut

use Getopt::Long;

$alignment_file = undef;
$outfile = undef;
$cutoff = undef;
$reference = undef;
$count_trailing_as_gap = undef;

&GetOptions('i=s' => \$alignment_file, 'o=s' => \$outfile, 'min=f' => \$cutoff, 'ref=s' => \$reference, 'trailgap!' => \$count_trailing_as_gap, 'h!' => \$help);
if (!$alignment_file or !$outfile or $help) {
	system ('perldoc', $0);
	exit;
}
if ($cutoff) {
	if ($cutoff < 0 || $cutoff > 1) {	
		print "\n  Error: allowed minimum occupancy (-min) range: 0 to 1 \n";
		print "\n  perl degeprime.pl -h for more help \n\n";
		exit;
	}
}
if ($cutoff and $reference) {
	print "\n  Error: not possible to specify both minimum occupancy and reference sequence \n";
	print "\n  perl degeprime.pl -h for more help \n\n";
	exit;
}
if (!$cutoff and !$reference) {
	$cutoff = 0; # the default value of $cutoff
}

#$reference = "S001099426"; # E. coli K12  GenBank J01695
#$reference = "S000497753"; # Sulfolobus solfataricus P2; AE006720
#$aligned_file = "eukaryotic/arb-silva.de_2011-06-01_id2791.fasta";
#$aligned_file = "23S_bact/arb-silva.de_2011-06-20_id5497.fasta";
####################

print"\nRunning TrimAlignment\n";
if ($reference) {
	print"Getting positions occupied by reference sequence\n";
	&get_ref_positions;
} else {
	print"Calculating occupancies\n";
	&get_rich_positions;
}
print"Printing trimmed alignment\n";
&get_short_seqs;
print"Finnished TrimAlignment succesfully\n\n";

####################

sub get_rich_positions {
	open (INFILE, $alignment_file) || die ("Can't open $alignment_file");
	$id = "NA";
	$seq = "";
	$total = 0;
	@counts = ();
	while (<INFILE>) {
		chomp;
		$row = $_;
		if (substr($_, 0, 1) eq ">") {
			if ($seq ne "") {
				$total++;
				$seq =~ s/\s+//g;
				$l = length($seq);
				$lengths{$l} = 1;
				if ($count_trailing_as_gap) {
					$seq =~ s/\./\-/g;
				}
				for ($i = 0; $i < length($seq); $i++) {
					if (  substr($seq, $i, 1) ne "-" ) { # all
						$counts[$i]++;
					}
				}	
     			}
				substr($row, 0, 1) = "";
	           	@fields = split(/\s+/, $row);
    	       	$id = $fields[0];
			$seq = "";
		} else {					
			$seq = $seq.$_;
		}
	}
	$total++;
	$seq =~ s/\s+//g;
	$l = length($seq);
	$lengths{$l} = 1;
	if ($count_trailing_as_gap) {
		$seq =~ s/\./\-/g;
	}
	for ($i = 0; $i < length($seq); $i++) { 
		if (  substr($seq, $i, 1) ne "-" ) { # all
			$counts[$i]++;
		}
	}
	close(INFILE);

	if ((keys %lengths) > 1) {
		print "\nError: Aligned sequences have different lengths! \n\n"; die;
	}

	@comp_pos = ();
	for ($i = 0; $i < @counts; $i++) {
		if (($counts[$i]/$total) >= $cutoff) { 
			push(@comp_pos, $i);
		}
	}
}

sub get_ref_positions {
	open (INFILE, $alignment_file) || die ("Can't open $alignment_file");
	$id = "NA";
	$seq = "";
	$ref_found = 0;
	while (<INFILE>) {
		chomp;
		$row = $_;
		if (substr($_, 0, 1) eq ">") {
			if ($seq ne "") {
				$seq =~ s/\s+//g;
				if ($id eq $reference) {
					$ref_found = 1;
					#print"$reference\n"; die;
					@comp_pos = ();
					for ($i = 0; $i < length($seq); $i++) {
						#if (  substr($seq, $i, 1) =~ m/[A-Z]/ ) { # Only capital
						if (  substr($seq, $i, 1) =~ m/[A-Z]/i ) { # All
							push(@comp_pos, $i);
						}
					}
					last;
				}			
    	 	}
			substr($row, 0, 1) = "";
	        @fields = split(/\s+/, $row);
    	    $id = $fields[0];
			$seq = "";
		} else {					
			$seq = $seq.$_;
		}
	}
	$seq =~ s/\s+//g;
	if ($id eq $reference) {
		$ref_found = 1;
		@comp_pos = ();
		for ($i = 0; $i < length($seq); $i++) { 
			#if (  substr($seq, $i, 1) =~ m/[A-Z]/) { # Only capital
			if (  substr($seq, $i, 1) =~ m/[A-Z]/i ) { # All
				push(@comp_pos, $i);
			}
		}
	}			
	close(INFILE);
	if ($ref_found == 0) {
		print "\n Error: The reference sequence was not found. Note: if the sequence ids are split by spaces it's only the part before the first space that is considered. \n\n";
		die;
	}
}

sub get_short_seqs {
	open (INFILE, $alignment_file) || die ("Can't open $alignment_file");
	open (OUT, ">$outfile");
	$id = "NA";
	$seq = "";
	while (<INFILE>) {
    	chomp;
		$row = $_;
		if (substr($row, 0, 1) eq ">") {
			if ($seq ne "") {
				$seq =~ s/\s+//g;
				$seq =~ tr/[a-z]/[A-Z]/;
				$seq =~ tr/U/T/;
				$seqshort = "";
				for ($i = 0; $i < (@comp_pos - 1); $i++) {
					$nt = substr($seq, $comp_pos[$i], 1);
					$gap = substr($seq, ($comp_pos[$i] + 1),  ($comp_pos[$i+1] - $comp_pos[$i] - 1));
					if ($gap =~ m/\w/) {
						#print"$id pos:$i $gap\n";
						$nt =~ tr/[A-Z]/[a-z]/;
					}
					$seqshort = $seqshort.$nt;
				}
				$nt = $nt = substr($seq, $comp_pos[-1], 1);
				$seqshort = $seqshort.$nt;
				print OUT ">$id\n$seqshort\n";
	     	} 
			substr($row, 0, 1) = "";
        	@fields = split(/\s+/, $row);
	        $id = $fields[0];
			$seq = "";
		} else {
			$seq = $seq.$_;
		}
	}
	$seq =~ s/\s+//g;
	$seq =~ tr/[a-z]/[A-Z]/;
	$seq =~ tr/U/T/;
	$seqshort = "";
	for ($i = 0; $i < (@comp_pos - 1); $i++) {
		$nt = substr($seq, $comp_pos[$i], 1);
		$gap = substr($seq, ($comp_pos[$i] + 1),  ($comp_pos[$i+1] - $comp_pos[$i] - 1));
		if ($gap =~ m/\w/) {
			#print"$id pos:$i $gap\n";
			$nt =~ tr/[A-Z]/[a-z]/;
		}
		$seqshort = $seqshort.$nt;
	}
	$nt = $nt = substr($seq, $comp_pos[-1], 1);
	$seqshort = $seqshort.$nt;
	print OUT ">$id\n$seqshort\n";
	close(INFILE);
	close(OUT);
}
