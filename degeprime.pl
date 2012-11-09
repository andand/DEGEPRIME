# Anders Andersson 2008 - 2012

#use strict;
use Storable qw<dclone>;
use Getopt::Long;
use Time::HiRes qw(time);

### default parameter setting ###

$skip_length = 20; # numb of bases in both ends of sequence that will not be considered in analysis.
$min_depth = 10; # minimum number of spanning sequences without in/dels to include in a window.
$number_iteratons = 100; # number of iterations in randomized summation

######################

$aligned_short_file = undef;
&GetOptions('i=s' => \$aligned_short_file, 'w=i' => \$window_length, 'd=i' => \$max_deg, 'skip=i' => \$skip_length, 'iter=i' => \$number_iteratons);
if (!$aligned_short_file) {
  print "\n Error: No alignment file name given. Should be:\n\n    perl degeprime.pl -i trimmed_align_file -w window_size -d max_degeneracy > output_file \n\n";
  die;
}
if (!$window_length) {
  print "\n  Error: No window length given. Should be: \n\n    perl degeprime.pl -i trimmed_align_file -w window_size -d max_degeneracy > output_file \n\n"; 
  die;
}
if (!$max_deg) {
  print "\n  Error: No maximum degeneracy given. Should be: \n\n    perl degeprime.pl -i trimmed_align_file -w window_size -d max_degeneracy > output_file \n\n";
  die;
}

$print_taxon_coverage = 0;
$taxon_level = 4;

#$hyden_outfile = "hyden_out.txt";
#$hyden_infile = "hyden_in.fasta";

#$taxonomy_file = "rdp_10_7_taxonomy.txt";
#$aligned_short_file = "arch_align_rdp_10_7_short.fasta";
#$aligned_short_file = "rdp_download_138807seqs_v9_1200nt_good_aligned_short_r2000.fna";
#$aligned_short_file = "polyG_mouse_all_slim.fa";
#$aligned_short_file = "rdp_download_138807seqs_v9_1200nt_good_aligned_short.fna";
#$aligned_short_file = "eukaryotic/arb-silva.de_2011-06-01_id2791_0.90.fasta";

#############

$totaltime = time;

&make_iupac;
if ($print_taxon_coverage == 1) {
	&get_id_taxonomy;
}
&read_short_seqs;
&checking_lengths;
&calc_coverage;

$totaltime = time - $totaltime;

#print"Total time:\t$totaltime\n";

#############

sub calc_coverage {
	#print"Pos\tTotal\tUnique\tExp_Deg\tExp_Matching\tCon_Deg\tCon_Matching\tSum_Deg\tSum_Matching\tRan_Deg\tRan_Matching\n";
	print"Pos\tTotalSeq\tUniqueMers\tEntropy\tPrimerDeg\tPrimerMatching\tPrimerSeq\n";
	for ($pos = 0; $pos <= ($seq_length - $window_length); $pos++) {
	#for ($pos = 385; $pos <= 385; $pos++) {
	#for ($pos = 470; $pos <= 519; $pos++) {
		$time0 = time;
		#($total_spanning, $zero_insertions) = &get_base_ranking($pos);
		($total_spanning, $zero_insertions, $entropy) = &get_mer_ranking($pos);
		if ($zero_insertions >= $min_depth) {
			$unique = @sorted_mer;
			print"$pos\t$total_spanning\t$unique\t$entropy";
			#&expansion;
			#&constriction;
			#&summation($pos);
			&random_summation($pos);
			#&hyden($pos);
			$time1 = time;
			$time = ($time1 - $time0);
			#print"\t$time\t$timeB\t$timeC";
			print"\n";	
		}
		
	}
	#print"\n\n  Elapsed Time: $time sec\n\n";
}

# for comparison with Hyden
sub hyden {
	local($startpos) = $_[0];
	local($total) = 0;
	local($deg) = undef;
	local($match) = 0;
	local($deg_primer) = "";
	local($new_primer) = "";
	local(%pos_n_primers) = ();
	local($pos) = 0;
	open (OUT, ">$hyden_infile") || die "Can't make $hyden_infile\n\n";
	foreach $id (keys %id_seq) {
		if (($start{$id} + $skip_length) <= $startpos) {
			if (($end{$id} - $skip_length) >= $startpos + $window_length) {
				$seq = $id_seq{$id};
				$mer = substr($seq, $startpos, $window_length);
				if ($mer !~ m/[^(A|T|C|G)]/) {
					$mer = $mer."A";
					print OUT ">$id\n$mer\n";
					$total++;
				}
			}
		}
	}
	$command = "hyden.exe -dna $hyden_infile -len5 $window_length -len3 1 -deg5 $max_deg -deg3 1 -from5 0 -to5 0 -from3 0 -to3 0 -mis5 0 -mis3 0 -nentropy $total -mprimers 1 > $hyden_outfile";
	system($command);
	($deg, $match, $deg_primer) = &get_results;
	print"\t$total\t$deg\t$match\t$deg_primer";
	#unlink($outfile);
	#die;
	#next;

	$pos_n_primers{-1}{""} = 1;
	$match = 0;
	for ($pos = 0; $pos < $window_length; $pos++) {
		$deg_nt = substr($deg_primer, $pos , 1);
		foreach $primer (keys %{$pos_n_primers{$pos - 1}}) {
			foreach $nt (keys %{$deg_n_nt{$deg_nt}}) {
				$new_primer = $primer.$nt;
				$pos_n_primers{$pos}{$new_primer} = 1;
			}
		}
	}
	foreach $primer (keys %{$pos_n_primers{$window_length - 1}}) {
		if (defined $mer_count{$primer}) {
			$match = $match + $mer_count{$primer};
		}				
	}
	print"\t$match";
}

sub get_results {
	local($stop) = 0;
	local($primer) = undef;
	local($match) = undef;
	local($total) = undef;
	local($deg) = undef;

	open (INFILE, $hyden_outfile) || die ("Can't open $outfile");
	while (<INFILE>) {
		chomp;
		if ($stop == 1) {
			@fields = split(/\s+/);
			$primer = $fields[1];
			last;
		}
		if (substr($_, 0, 6) eq "Best 5") {
			@fields = split(/\s+/);
			$match = $fields[4];
			$total = $fields[8];
			substr($total5, -1, 1) = "";
		}
		if (substr($_, 0, 7) eq "Primer:") {
			@fields = split(/\s+/);
			$deg = $fields[4];
			$stop = 1;
		}		
	}
	return($deg, $match, $primer);
	#print"$i\t$deg\t$match\t$primer\n";
}

sub get_id_taxonomy {
	$in = 0;
	open (INFILE, $taxonomy_file) || die ("Can't open");
	while (<INFILE>) {
		chomp;
		$row = $_;
		@fields = split(/\s+/, $row);
		$rdp_id = $fields[0];
		@tax = split(/;/, $fields[1]);
		$tax = $tax[0];
		for ($i = 1; $i < $taxon_level; $i++) {
			$tax = $tax.";".$tax[$i];
		}
		print"$rdp_id\t$tax\n";			
	}
}

sub read_short_seqs {
	open (INFILE, $aligned_short_file) || die ("Can't open");
	$seq = "";
	while (<INFILE>) {
		chomp;
		$row = $_;
		if (substr($row, 0, 1) eq ">") {
			if ($seq ne "") {
				$id_seq{$id} = $seq;
				$seq = "";
			}
			@fields = split(/\s+/, $row); # split id at space
			$id = $fields[0];
			substr($id, 0, 1) = "";
		} else {
			$seq = $seq.$row;
		}	
	}
	$id_seq{$id} = $seq;
	close (INFILE);
}

sub checking_lengths {
	foreach $id (keys %id_seq) {
		$seq = $id_seq{$id};
		$seq_length = length($seq);
		$length_id{$l} = $id;
		die if ($seq !~ m/^-*/); # ??
		$seq =~ m/(^\W*)/;
		#$seq =~ m/(^-*)/;
		$start{$id} = length($1); # position for first base
		die if ($seq !~ m/-*$/); # ??
		$seq =~ m/(\W*$)/;
		#$seq =~ m/(-*$)/;
		$end{$id} = length($seq) - length($1) - 1; # position for last base
	}
	if ((keys %length_id) > 1) {
		die ("Error: Not all aligned sequences have same length !!!\n");
	}
}

sub get_mer_ranking {
	local($startpos) = $_[0];
	local($mer);
	local($total_spanning) = 0;
	local($zero_gaps) = 0;
	local($entropy) = 0;
	local(%mer_count_entropy) = ();
	%mer_count = ();
	@sorted_mer = ();
	foreach $id (keys %id_seq) {
		if (($start{$id} + $skip_length) <= $startpos) {
			if (($end{$id} - $skip_length) >= $startpos + $window_length - 1) {
				$seq = $id_seq{$id};
				$mer = substr($seq, $startpos, $window_length);
				substr($mer, -1, 1) =~ tr/[a-z]/[A-Z]/;
				if ($mer !~ m/[^(A|T|C|G)]/) {
					$mer_count{$mer}++;
					$zero_gaps++;
				}
				$total_spanning++;
				
				#for entropy calc:
				$mer_count_entropy{$mer}++;	
			}
		}
	}
	#@sorted_mer = (sort {$mer_count{$b} <=> $mer_count{$a}} keys %mer_count); # original
	#@sorted_mer = (sort {$mer_count{$a} <=> $mer_count{$b}} keys %mer_count);
	@sorted_mer = keys %mer_count;
	# calculate entropy in window
	if ($total_spanning > 0) {
		foreach $mer (keys %mer_count_entropy) {
			$entropy = $entropy - ($mer_count_entropy{$mer}/$total_spanning)*log($mer_count_entropy{$mer}/$total_spanning)/log(2);
		}
	}
	return($total_spanning, $zero_gaps, $entropy);
}

sub get_base_ranking {
	local($startpos) = $_[0];
	local($pos);
	local($total_spanning) = 0;
	local($zero_gaps) = 0;
	@nt_matr = ();
	@freq_matr = ();
	%pos_n_base = ();
	foreach $id (keys %id_seq) {
		if (($start{$id} + $skip_length) <= $startpos) {
			if (($end{$id} - $skip_length) >= $startpos + $window_length - 1) {
				$seq = $id_seq{$id};
				$mer = substr($seq, $startpos, $window_length);
				substr($mer, -1, 1) =~ tr/[a-z]/[A-Z]/; 
				if ($mer !~ m/[^(A|T|C|G)]/) {
					for ($pos = 0; $pos < $window_length; $pos++) {
						$nt = substr($mer, $pos, 1);
						$pos_n_base{$pos}{$nt}++;
					}
					$zero_gaps++;
				}
				$total_spanning++;
			}
		}
	}
	if ($zero_gaps > 0) {
		for ($pos = 0; $pos < $window_length; $pos++) {
			@sorted_nt = (sort {$pos_n_base{$pos}{$b} <=> $pos_n_base{$pos}{$a}} keys %{$pos_n_base{$pos}});
			for ($i = 0; $i < 4; $i++) {
				if ($pos_n_base{$pos}{$sorted_nt[$i]} > 0) {
					$nt_matr[$pos][$i] = $sorted_nt[$i];
					$freq_matr[$pos][$i] = $pos_n_base{$pos}{$sorted_nt[$i]} / $total;
				} else {
					$freq_matr[$pos][$i] = 0;
					$nt_matr[$pos][$i] = "X";
				}		
			}
		}
	}
	return($total_spanning, $zero_gaps);
}

# just for comparison
sub expansion {
	local($pos);
	local($deg) = 1;
	local($matching);
	@rank = ();
	# set rank to most frequent base:
	for ($pos = 0; $pos < $window_length; $pos++) {
		$rank[$pos] = 0;
	}
	while ($deg < $max_deg) {
		$highest_freq = 0;
		for ($pos = 0; $pos < $window_length; $pos++) {
			if ($freq_matr[$pos][$rank[$pos] + 1] > $highest_freq) {
				$newdeg = $deg * ($rank[$pos] + 2) / ($rank[$pos] + 1);
				if ($newdeg <= $max_deg) {
					$highest_freq = $freq_matr[$pos][$rank[$pos] + 1];
					$highest_pos = $pos;
				}
			}
		}
		last if ($highest_freq == 0);
		$newdeg = $deg * ($rank[$highest_pos] + 2) / ($rank[$highest_pos] + 1);
		$deg = $newdeg;
		$rank[$highest_pos]++;
	}
	$matching = &check_matching;
	print"\t$deg\t$matching";
	#print"\t$matching";
	#&greedy($matching);
}

# just for comparison
sub constriction {
	local($pos);
	local($deg) = 1;
	local($matching);
	@rank = ();
	# set rank to base with lowest freq (above zero)
	for ($pos = 0; $pos < $window_length; $pos++) {
		for ($rank = 0; $rank < 4; $rank++) {
			if ($freq_matr[$pos][$rank] > 0) {
				$rank[$pos] = $rank;
			}
		}
	}
	for ($pos = 0; $pos < $window_length; $pos++) {
		$deg = $deg * ($rank[$pos] + 1);
	}
	while ($deg > $max_deg) {
		$lowest_freq = 10;
		for ($pos = 0; $pos < $window_length; $pos++) {
			if ($rank[$pos] > 0) {
				if ($freq_matr[$pos][$rank[$pos]] < $lowest_freq) {
					$lowest_freq = $freq_matr[$pos][$rank[$pos]];
					$lowest_pos = $pos;
				}
			}
		}
		last if ($lowest_freq == 10);
		$newdeg = $deg * ($rank[$lowest_pos]) / ($rank[$lowest_pos] + 1);
		$deg = $newdeg;
		$rank[$lowest_pos]--;
		last if ($newdeg < $max_deg);
	}
	$matching = &check_matching;
	print"\t$deg\t$matching";
	#print"\t$matching";
}

sub summation {
	local($pos);
	local($deg);
	local($matching) = 0;
	local(%pos_n_base) = (); #there is a global one as well!
	local(%temp_pos_n_base) = (); #there is a global one as well!
	local(%pos_n_primers) = ();

	#sum best primers while deg < max_deg:
	foreach $mer (@sorted_mer) {
		$deg = 1;
		%temp_pos_n_base = %{dclone(\%pos_n_base)};
		for ($pos = 0; $pos < $window_length; $pos++) {
			$nt = substr($mer, $pos, 1);
			$temp_pos_n_base{$pos}{$nt} = 1;
			$deg = $deg * (keys %{$temp_pos_n_base{$pos}});
		}
		if ($deg <= $max_deg) {
			%pos_n_base = %{dclone(\%temp_pos_n_base)};
		} else {
			last;
		}
	}

	#count matches:
	$pos_n_primers{-1}{""} = 1;
	for ($pos = 0; $pos < $window_length; $pos++) {
		foreach $primer (keys %{$pos_n_primers{$pos - 1}}) {
			foreach $nt (keys %{$pos_n_base{$pos}}) {		
				$new_primer = $primer.$nt;
				$pos_n_primers{$pos}{$new_primer} = 1;
			}
			
		}
	}
	$match = 0;
	foreach $primer (keys %{$pos_n_primers{$window_length - 1}}) {
		if (defined $mer_count{$primer}) { 
			$match = $match + $mer_count{$primer};
		}				
	}

	#make deg primer:
	$deg_primer = "";
	$deg = 1;
	for ($pos = 0; $pos < $window_length; $pos++) {
		$string = join("", keys %{$pos_n_base{$pos}});
		$deg_primer = $deg_primer.$deg{$string};
		$deg = $deg * (keys %{$pos_n_base{$pos}});
	}
	print"\t$deg\t$match\t$deg_primer";
	#print"\t$deg\t$match";
	#print"\t$match";
}

sub random_summation {
	local($pos);
	local($deg);
	local($matching);
	local($min_index);
	local($max_index);
	local(%pos_n_base);
	local(%temp_pos_n_base);
	local(%pos_n_primers);
	local(%already_chosen);
	local(@end_indices) = ();
	local(@temp_end_indices);
	local(@sorted_indices);
	local($bestmatching) = 0;
	local($bestmatch) = 0;
	local($bestdeg) = "";
	local($adjust_counts);
	local($adjust_index);

	# Make list from which to make random draws:
	$end_indices[0] = $mer_count{ $sorted_mer[0] };
	for ($i = 1; $i < @sorted_mer; $i++) {
		$end_indices[$i] = $end_indices[$i - 1] + $mer_count{ $sorted_mer[$i] };
	}
	
	# Iterate randomised summation procedure:
	$timeA = $timeB = $timeC = 0;
	for ($iter = 0; $iter < $number_iteratons; $iter++) {
		#print"$iter\n";
		@temp_end_indices = @end_indices;
		@temp_sorted_mer = @sorted_mer;
		$deg = 0;
		$trial = 0;
		$matching = 0;
		%pos_n_base = ();
		%already_chosen = ();
		while ($deg < $max_deg) {
			$trial++;
			#print" trials: $trial\n";
			last if ($trial > 100);
			#get random mer:
			$random = int(rand($temp_end_indices[-1]));
			$min_index = 0;
			$max_index = @temp_end_indices;
			$timeA = time;
			while ($max_index > $min_index) {
				$i = int( ($min_index + $max_index)/2 );
				#print"min: $min_index max: $max_index i: $i ival: $temp_end_indices[$i - 1] i-1val: $temp_end_indices[$i] r: $random\n";
				if ($temp_end_indices[$i] < $random) {
					$min_index = $i;
				} elsif ($temp_end_indices[$i - 1] >= $random) {
					$max_index = $i;
				} else {
					last;					
				}
			}
			$timeB = $timeB + time - $timeA;
			$selected_index = $i;	
			$mer = $temp_sorted_mer[$selected_index];

			if (defined $already_chosen{$selected_index}) {
				$trial--;
				$fails++;
				#print"\n  fails: $fails\n";
			} else {
				$already_chosen{$selected_index} = $selected_index;
				$fails = 0;

				# Calculate newdeg:
				$newdeg = 1;
				%temp_pos_n_base = %{dclone(\%pos_n_base)};
				for ($pos = 0; $pos < $window_length; $pos++) {
					$nt = substr($mer, $pos, 1);
					$temp_pos_n_base{$pos}{$nt} = 1;
					$newdeg = $newdeg * (keys %{$temp_pos_n_base{$pos}});
				}
				# Add mer if newdeg <= max_deg:
				if ($newdeg <= $max_deg) {
					$deg = $newdeg;
					$matching = $matching + $mer_count{$mer};
					%pos_n_base = %{dclone(\%temp_pos_n_base)};
				}
			}

			#update @temp_end_indices & @temp_sorted_mer
			if ($fails > 4) { # 4		
				$timeA = time;
				@sorted_indices = (sort {$already_chosen{$a} <=> $already_chosen{$b}} keys %already_chosen);
				$adjust_counts = 0;
				$adjust_index = 0;				
				for ($i = $sorted_indices[0]; $i < (@temp_end_indices - @sorted_indices - 1); $i++) {	
					while (defined $already_chosen{$i + $adjust_index}) {
						$adjust_counts = $adjust_counts + $mer_count{ $temp_sorted_mer[$i + $adjust_index] };
						$adjust_index++;
					}		
					$temp_end_indices[$i] = $temp_end_indices[$i + $adjust_index] - $adjust_counts;
				}
				#print"temp_end_indices: @temp_end_indices\n";
				for ($i = 0; $i < @sorted_indices; $i++) {
					splice(@temp_end_indices, -1, 1);			
				}
				#print"	\ntemp_sorted_mer BEFORE: @temp_sorted_mer\n";		
				#print"sorted_indices: @sorted_indices\n";
				for ($i = (@sorted_indices - 1); $i >= 0; $i--) {
					splice(@temp_sorted_mer, $sorted_indices[$i], 1);
					#print"temp_sorted_mer AFTER: @temp_sorted_mer\n";		
				}
				#print"\n";
				%already_chosen = ();
$timeC = $timeC + time - $timeA;
 			}
			# Stop if all mers have been tried:
			last if (@temp_end_indices == 0);
		}
		
		# Count matches:
		%pos_n_primers = ();
		$pos_n_primers{-1}{""} = 1;
		for ($pos = 0; $pos < $window_length; $pos++) {
			foreach $primer (keys %{$pos_n_primers{$pos - 1}}) {
				foreach $nt (keys %{$pos_n_base{$pos}}) {		
					$new_primer = $primer.$nt;
					$pos_n_primers{$pos}{$new_primer} = 1;
				}
			}	
		}
		$match = 0;
		foreach $primer (keys %{$pos_n_primers{$window_length - 1}}) {
			if (defined $mer_count{$primer}) { 
				$match = $match + $mer_count{$primer};
			}				
		}

		if ($match > $bestmatch) {
			$bestmatch = $match;
			$bestdeg = $deg;
			$deg_primer = "";
			for ($pos = 0; $pos < $window_length; $pos++) {
				$string = join("", keys %{$pos_n_base{$pos}});
				$deg_primer = $deg_primer.$deg{$string};
			}
			$bestprimer = $deg_primer;
		}
		#print"trial $trial iter $iter\n";
	}
	print"\t$bestdeg\t$bestmatch\t$bestprimer";
}

sub get_taxon_coverage  {
	local($overlapping) = 0;
	local($matching) = 0;
	local($perc) = 0;
	$matching = $overlapping = $perc = 0;
	if ($sorted_primers[$i] ne "") { 
		foreach $id (keys %{$seqfile_n_id{$seqfile}}) {
			next if ($start{$id} > $pos);
			next if ($end{$id} < $pos);
			$overlapping++;
			if ($id_rawseq{$id} =~ m/$sorted_primers[$i]/) {
				$matching++;
			}
		}
	}
	if ($overlapping > 0) {
		$perc = int(10000 * $matching/$overlapping)/100;
	}
	return($matching, $overlapping, $perc);
}

sub print_primer {
	for ($pos = 0; $pos < $window_length; $pos++) {
		print"$pos\t";
		for ($rank = 0; $rank <= $rank[$pos]; $rank++) {
			print" $nt_matr[$pos][$rank]";
		}
		print"\n";
	}
}

sub check_matching {
	local($match) = 0;
	local($primer) = "";
	local($new_primer) = "";
	local(%pos_n_primers) = ();
	local($pos) = 0;	
	$pos_n_primers{-1}{""} = 1;
	for ($pos = 0; $pos < $window_length; $pos++) {
		foreach $primer (keys %{$pos_n_primers{$pos - 1}}) {
			for ($rank = 0; $rank <= $rank[$pos]; $rank++) {
				$nt = $nt_matr[$pos][$rank];
				$new_primer = $primer.$nt;
				$pos_n_primers{$pos}{$new_primer} = 1;
			}
		}
	}
	foreach $primer (keys %{$pos_n_primers{$window_length - 1}}) {
		if (defined $mer_count{$primer}) {
			$match = $match + $mer_count{$primer};
		}				
	}
	return ($match);
}

sub make_iupac {
$deg{"A"} = "A";
$deg{"C"} = "C";
$deg{"G"} = "G";
$deg{"T"} = "T";

$deg{"AT"} = "W";
$deg{"TA"} = "W";

$deg{"CG"} = "S";
$deg{"GC"} = "S";

$deg{"AC"} = "M";
$deg{"CA"} = "M";

$deg{"GT"} = "K";
$deg{"TG"} = "K";

$deg{"AG"} = "R";
$deg{"GA"} = "R";

$deg{"CT"} = "Y";
$deg{"TC"} = "Y";

$deg{"CGT"} = "B";
$deg{"CTG"} = "B";
$deg{"TCG"} = "B";
$deg{"TGC"} = "B";
$deg{"GTC"} = "B";
$deg{"GCT"} = "B";

$deg{"AGT"} = "D";
$deg{"ATG"} = "D";
$deg{"TAG"} = "D";
$deg{"TGA"} = "D";
$deg{"GAT"} = "D";
$deg{"GTA"} = "D";

$deg{"ACT"} = "H";
$deg{"ATC"} = "H";
$deg{"TAC"} = "H";
$deg{"TCA"} = "H";
$deg{"CTA"} = "H";
$deg{"CAT"} = "H";

$deg{"ACG"} = "V";
$deg{"AGC"} = "V";
$deg{"GAC"} = "V";
$deg{"GCA"} = "V";
$deg{"CAG"} = "V";
$deg{"CGA"} = "V";

$deg{"ACGT"} = "N";
$deg{"ACTG"} = "N";
$deg{"ATCG"} = "N";
$deg{"ATGC"} = "N";
$deg{"AGTC"} = "N";
$deg{"AGCT"} = "N";

$deg{"CAGT"} = "N";
$deg{"CATG"} = "N";
$deg{"CTAG"} = "N";
$deg{"CTGA"} = "N";
$deg{"CGAT"} = "N";
$deg{"CGTA"} = "N";

$deg{"GACT"} = "N";
$deg{"GATC"} = "N";
$deg{"GTAC"} = "N";
$deg{"GTCA"} = "N";
$deg{"GCTA"} = "N";
$deg{"GCAT"} = "N";

$deg{"TACG"} = "N";
$deg{"TAGC"} = "N";
$deg{"TGAC"} = "N";
$deg{"TGCA"} = "N";
$deg{"TCAG"} = "N";
$deg{"TCGA"} = "N";

$deg_n_nt{"A"}{"A"} = 1;
$deg_n_nt{"C"}{"C"} = 1;
$deg_n_nt{"G"}{"G"} = 1;
$deg_n_nt{"T"}{"T"} = 1;

$deg_n_nt{"W"}{"A"} = 1;
$deg_n_nt{"W"}{"T"} = 1;

$deg_n_nt{"S"}{"C"} = 1;
$deg_n_nt{"S"}{"G"} = 1;

$deg_n_nt{"M"}{"A"} = 1;
$deg_n_nt{"M"}{"C"} = 1;

$deg_n_nt{"K"}{"G"} = 1;
$deg_n_nt{"K"}{"T"} = 1;

$deg_n_nt{"R"}{"A"} = 1;
$deg_n_nt{"R"}{"G"} = 1;

$deg_n_nt{"Y"}{"C"} = 1;
$deg_n_nt{"Y"}{"T"} = 1;

$deg_n_nt{"B"}{"C"} = 1;
$deg_n_nt{"B"}{"G"} = 1;
$deg_n_nt{"B"}{"T"} = 1;

$deg_n_nt{"D"}{"A"} = 1;
$deg_n_nt{"D"}{"G"} = 1;
$deg_n_nt{"D"}{"T"} = 1;

$deg_n_nt{"H"}{"A"} = 1;
$deg_n_nt{"H"}{"C"} = 1;
$deg_n_nt{"H"}{"T"} = 1;

$deg_n_nt{"V"}{"A"} = 1;
$deg_n_nt{"V"}{"C"} = 1;
$deg_n_nt{"V"}{"G"} = 1;

$deg_n_nt{"N"}{"A"} = 1;
$deg_n_nt{"N"}{"C"} = 1;
$deg_n_nt{"N"}{"G"} = 1;
$deg_n_nt{"N"}{"T"} = 1;

}
