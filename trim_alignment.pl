# Anders Andersson 2008 - 2012

use Getopt::Long;

$aligned_file = undef;
$cutoff = undef;
$count_trailing_as_gap = "F";
$cutoff = 0; # the default value of $cutoff

&GetOptions('i=s' => \$aligned_file, 'c=f' => \$cutoff);
if (!$aligned_file) {
  print "\n Error: No alignment file name given. Should be either:\n\n   perl trim_alignment.pl -i align_file > trimmed_align_file\n\n or, if you want to remove columns of low occupancy:\n\n   perl trim_alignment.pl -i align_file -c cutoff > trimmed_align_file \n\n (0 < cutoff <= 1)\n\n";
  die;
}
#if (!$cutoff) {
#  print "\n  Error: No minimum alignment column occupancy given (0 <= cutoff <= 1). Should be: \n  perl print_short_seqs.pl -i seq_file -c cutoff > outfile \n\n";
#  die;
#}

#$reference = "S001099426"; # E. coli K12  GenBank J01695
#$reference = "S000497753"; # Sulfolobus solfataricus P2; AE006720
#$aligned_file = "eukaryotic/arb-silva.de_2011-06-01_id2791.fasta";
#$aligned_file = "23S_bact/arb-silva.de_2011-06-20_id5497.fasta";
#$cutoff = 0.90;

#&get_ref_positions;
&get_rich_positions;
&get_short_seqs;

####################

sub get_rich_positions { ###
open (INFILE, $aligned_file) || die ("Can't open");
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
			#print"$total\n";
			$seq =~ s/\s+//g;
			#print"$seq\n\n";
			$l = length($seq);
			#print"$l\n";
			$lengths{$l} = 1;
			if ($count_trailing_as_gap eq "T") {
				$seq =~ s/\./\-/g;
			}
			for ($i = 0; $i < length($seq); $i++) {
				if (  substr($seq, $i, 1) ne "-" ) { # all
					@counts[$i]++;
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
#print"$seq\n\n";
$l = length($seq);
#print"$l\n";
$lengths{$l} = 1;
if ($count_trailing_as_gap eq "T") {
	$seq =~ s/\./\-/g;
}
for ($i = 0; $i < length($seq); $i++) { 
	if (  substr($seq, $i, 1) ne "-" ) { # all
		@counts[$i]++;
	}
}
#print"Total: $total\n\n";
close(INFILE);

if ((keys %lengths) > 1) {
	print "\nError: Aligned sequences have different lengths! $temp \n\n"; die;
}

@comp_pos = ();
for ($i = 0; $i < @counts; $i++) {
	if ((@counts[$i]/$total) >= $cutoff) { 
		push(@comp_pos, $i);
	}
}
} ###

sub get_ref_positions { ###
open (INFILE, $aligned_file) || die ("Can't open");
$id = "NA";
$seq = "";
while (<INFILE>) {
     	chomp;
	$row = $_;
	if (substr($_, 0, 1) eq ">") {
		if ($seq ne "") {
			$seq =~ s/\s+//g;
			if ($id eq $reference) {
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
	@comp_pos = ();
	for ($i = 0; $i < length($seq); $i++) { 
		#if (  substr($seq, $i, 1) =~ m/[A-Z]/) { # Only capital
		if (  substr($seq, $i, 1) =~ m/[A-Z]/i ) { # All
			push(@comp_pos, $i);
		}
	}
}			
close(INFILE);
} ###

sub get_short_seqs { ###
open (INFILE, $aligned_file) || die ("Can't open");
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
			print">$id\n$seqshort\n";
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
print">$id\n$seqshort\n";
close(INFILE);
} ###


