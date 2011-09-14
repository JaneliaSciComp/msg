#!/usr/bin/perl -w
# updateRef.pl
# Example call:
# perl updateRef.pl dsim-all-chromosome-r1.3.fasta.msg update_reads-aligned-parent1.pileup 0 1 2 5

use strict;
use warnings;
use Getopt::Std;
use File::Basename;

my $src = dirname $0 ;

my ($refFile,$pileupFile,$add_indels,$minQV,$min_coverage,$max_coverage_stds) = @ARGV;
die "\n USAGE: updateRef.pl <reference_fasta_file> <pileup_file> <incorporate indels> <min QV (Sanger)> <min read coverage> <max read coverage std. dev.s>\n\n" if (@ARGV<4);

$min_coverage = 2 unless (defined $min_coverage);
#default for max_coverage_stds is undef, which means there is no max.

my %refReads = &readFasta($refFile);
my $outfile  = "$refFile.updated";
print "UPDATE REFERENCE FILE: $outfile\n";
print "\tmin read coverage: $min_coverage\n" .
		"\tmin qv threshold: $minQV\n";
if (defined $max_coverage_stds) {
    print "\tmax read coverage (in std. dev.s): $max_coverage_stds\n";
}
else {
    print "\tmax read coverage (in std. dev.s): NA\n";
}

unlink ("$outfile.fastq") if (-e "$outfile.fastq");

my $max_coverage;
if (defined $max_coverage_stds) {
    $max_coverage = get_max_coverage($pileupFile, $max_coverage_stds);
    print "\tmax read coverage value: $max_coverage\n";
}

my ($genome_total,$genome_uncovered,%coveredContigs) = &pileup2fq($pileupFile,"$outfile.fastq",\%refReads,
    $max_coverage); # replace regions in the reference with the solexa calls


### print out uncovered contigs
print "\nSummary by chromosome/contig:";
foreach my $contig (keys %refReads) {
	if (exists $coveredContigs{$contig}) {
		print "\n\t$contig (",length($refReads{$contig}),") $coveredContigs{$contig} length";
		print " *** lengths differ ***" if (length($refReads{$contig}) != $coveredContigs{$contig});

	} else {
		$genome_total += length($refReads{$contig});
		$genome_uncovered += length($refReads{$contig});
		
		open(OUT,">>$outfile.fastq") || die "ERROR: Can't write to $outfile.fastq: $!\n";
		my $printout;
		my $fake_qual ='!' x length($refReads{$contig});
		print OUT "\@$contig\n"; $printout = &p2q_print_str($contig,'seq',\$refReads{$contig}); print OUT "$printout\n";
		print OUT "+$contig\n"; $printout = &p2q_print_str($contig,'qual',\$fake_qual); print OUT "$printout\n";
		close OUT;
	}
}

print "\nTotal sites: $genome_total\n";
print "Total sites not covered: $genome_uncovered\n";
print "Fraction of reference not covered: ",$genome_uncovered/$genome_total,"\n";

system("$src/fastq_2_fasta.pl -q $minQV < $outfile.fastq > $outfile.fasta ");

exit;


####################################################################################################
####################################################################################################

sub sum_of_array {
    my ($values_ref) = @_;
    my $sum = 0;
    foreach (@$values_ref) { $sum += $_ }
    return $sum;
}

sub get_std_deviation {
    my ($values_ref) = @_;
    my $mean = sum_of_array($values_ref) / scalar(@$values_ref);
    my $total_for_var = 0;
    foreach my $val (@$values_ref) {
        $total_for_var += (($val - $mean)**2);
    }
    my $variance = $total_for_var / scalar(@$values_ref);
    my $std = sqrt($variance);
    return $std;
}

sub get_max_coverage {
    #The user parameter for max coverage specifies how many standard deviations
    #from the mean to allow.  This function returns what the actual number 
    #that would be.
    my ($pileup_file, $max_coverage_stds) = @_;

    my $max_coverage = 0;
    my @coverage_values;

    open(FILE,$pileup_file) || die "ERROR Can't open $pileup_file: $!\n";
    while (<FILE>) { my $line = $_; chomp $line;;
        #column 7 (0 based) is coverage
        my @t = split(/\s+/,$line);
        # Only look at reads with coverage
        if ($t[7] > 0) {
            push (@coverage_values, $t[7]);
        }
    } close FILE;
    my $mean = sum_of_array(\@coverage_values)/scalar(@coverage_values);
    my $std = get_std_deviation(\@coverage_values);
    $max_coverage = $mean + ($max_coverage_stds * $std);
    return sprintf "%.0f", $max_coverage; #round to int
}

sub readFasta {
	my ($file,$refReads) = @_;
	my %reads;

	my ($header,$seq,$keeper) = ('','',0);
	open (FILE,$file) || die "ERROR Can't read $file: $!\n";
	while (<FILE>) { chomp $_;
		if ($_ =~ /^>(\S+)/) { 
			$reads{$header} = $seq if ($seq);

			$header = $1;
			$seq = '';
			$keeper = 1;

		} elsif ($_ =~ />/) { $keeper = 0;
		} elsif ($keeper==1) { $seq .= $_; }
	} close FILE;
	$reads{$header} = $seq if ($seq);

	return %reads;
}


### ADAPTED FROM samtools.pl file
# Usage:   samtools.pl pileup2fq [options] <in.cns-pileup>
# Options: -d INT    minimum depth        [$opts{d}]
#          -D INT    maximum depth        [$opts{D}]
#          -Q INT    min RMS mapQ         [$opts{Q}]
#          -G INT    minimum indel score  [$opts{G}]
#          -l INT    indel filter winsize [$opts{l}]\n
sub pileup2fq {
  my ($file,$outfile,$refReads,$max_coverage) = @_;
  my %opts = (d=>1, D=>255, Q=>0, G=>0, l=>0);
#  my %opts = (d=>3, D=>255, Q=>25, G=>25, l=>10);
  getopts('d:D:Q:G:l:', \%opts);

  my (@gaps, $counter, $last_line, %contigs);
  my ($seq, $qual) = ('','');
  my $_Q = $opts{Q};
  my $_d = $opts{d};
  my $_D = $opts{D};

  my $last_pos = 0;
  my $skip_to_pos = -1;
  my $last_chr = '';
  my ($uncovered,$total) = (0,0);
  my ($genome_uncovered,$genome_total) = (0,0);

  open(FILE,$file) || die "ERROR Can't open $file: $!\n";
  while (<FILE>) { my $line = $_; chomp $line; $counter++;

   # 0-2   chromosome, coordinate, reference base, 
   # 3-6   consensus base, Phred-scaled consensus quality, SNP quality, RMS mapping quality
   # 7-9   coverage, read bases, read qualities and alignment mapping qualities
	my @t = split(/\s+/,$line);
	die "ERROR (line $counter): position -->$t[1]<-- not numeric in $pileupFile\n" unless ($t[1] =~ /\d+/);
	die "ERROR (line $counter): chromsome/scaffold -->$t[0]<-- not found in reference\n" unless (exists $$refReads{$t[0]});

	if ($last_chr ne $t[0]) { # new contig OR skip b/c part of indel

     if ($last_chr) {
	     ### TINA add to the end
        if ($last_pos < length($$refReads{$last_chr})) {
				$seq  .= uc(substr($$refReads{$last_chr},$last_pos,(length($$refReads{$last_chr}) - $last_pos))); ### TINA
				$qual .= '!' x (length($$refReads{$last_chr}) - $last_pos);

				$uncovered += (length($$refReads{$last_chr}) - $last_pos);
				$total += (length($$refReads{$last_chr}) - $last_pos);
        }

        &p2q_post_process($outfile, $last_chr, \$seq, \$qual, \@gaps, $opts{l});
		  $contigs{$last_chr} = length($seq);

		  $genome_uncovered += $uncovered;
		  $genome_total += $total;
		  print "...(line $counter) $last_chr $last_pos ::: total $total uncovered $uncovered ",$uncovered/$total," (genome $genome_uncovered $genome_total)\n";
     }
	  
	  $last_chr = $t[0]; die "ERROR (line $counter): read -->$t[0]<-- not found in $refFile\n" unless (exists $$refReads{$t[0]});
	  $last_pos = 0;
	  $skip_to_pos = -1;
	  $seq = ''; $qual = '';
	  @gaps = ();

	  ($uncovered,$total) = (0,0);
	}

	next unless ($t[1] > $skip_to_pos); ### deletion

	### PAD IN BETWEEN LINES
	if ($t[1] - $last_pos > 1) {

		### extends past reference length
		if ($t[1] > length($$refReads{$last_chr})) { $seq .= uc('n' x ($t[1] - $last_pos - 1)); }
		else { $seq .= uc(substr($$refReads{$last_chr},$last_pos,($t[1] - $last_pos - 1))); } ### TINA
		$qual .= '!' x ($t[1] - $last_pos - 1);
		
		$uncovered += ($t[1] - $last_pos - 1);
		$total += ($t[1] - $last_pos - 1);
	}


	if ($t[2] eq '*') {
	  push(@gaps, $t[1]);# if ($t[5] >= $opts{G});

	  if ($t[3] !~ /^\*\/\*/) {
	  
			my ($allele1,$allele2) = ($t[8],$t[9]);
			my ($n1,$n2) = ($t[10],$t[11]);
			
			my $indel_allele_num = 0; 
			$indel_allele_num += $n1 if ($allele1 ne '*');
			$indel_allele_num += $n2 if ($allele2 ne '*');

			if ($indel_allele_num/($n1+$n2) > .5) {
				my ($update_seq,$update_qual) = ('','');
				my @genotypes = split(/\//,$t[3]);

				my ($ins_or_del,$indel) = ('','');
				if ($genotypes[0] =~ /(\+|\-)(\S+)/) { ($ins_or_del,$indel) = ($1,$2); }
				elsif ($genotypes[1] =~ /(\+|\-)(\S+)/) { ($ins_or_del,$indel) = ($1,$2); }

				if ($add_indels!=0) {
					if ($ins_or_del eq '-')    { $update_seq = ''; $skip_to_pos = $t[1]+length($indel); }
					elsif ($ins_or_del eq '+') { $update_seq = $indel; }

					my $q = $t[4] + 33; $q = 126 if ($q > 126);
					$update_qual = chr($q) x length($update_seq);

					$seq  .= uc($update_seq);
	  				$qual .= $update_qual;
					$total += length($indel) if ($ins_or_del eq '-');

				### IGNORE INDELS, basically pad with N's if it's a deletion from the reference
				} else {
					if ($ins_or_del eq '-') {
						$seq  .= ('N' x length($indel)); $skip_to_pos = $t[1]+length($indel);
						$qual .= ('!' x length($indel));
						$total += length($indel);
					} 
				}
			}
		}

	# read bases are NOT all *
	} elsif ($t[8] !~ /^\*+$/) {

		### satisfies min_coverage AND are not all 0 (!) QVs
		if (($t[7]>=$min_coverage) && ($t[9] !~ /^!+$/) && 
            ((not (defined $max_coverage)) || $t[7]<=$max_coverage)) {
			$seq .= uc($t[3]);
			my $q = $t[4] + 33;
			$q = 126 if ($q > 126);
			$qual .= chr($q);
			$total++;
		### use default reference
		} else {
			$seq .= uc($t[2]);
			$qual .= '!';
			$total++;
		}

	} else { # read bases are all * --- represents a deletion
			$seq .= 'N';
			$qual .= '!';
			$total++;
	}

	die "ERROR updateRef ($last_chr last_pos $last_pos skip_to_pos $skip_to_pos line $counter): lengths not equal (seq ",length($seq),"!=",length($qual),")\nLINE: $line\n" if (length($seq)!=length($qual));

	$last_pos = $t[1];
	$last_pos = $skip_to_pos if ($skip_to_pos>$t[1]);
	$last_line = $line;
	die "ERROR updateRef ($last_chr last_pos $last_pos skip_to_pos $skip_to_pos line $counter): pos lengths not equal (seq ",length($seq),"!=$last_pos)\nLINE: $line\n" if (length($seq)!=$last_pos);

  } close FILE;


  if ($last_chr) {

     ### TINA add to the end
     if ($last_pos < length($$refReads{$last_chr})) {
			$seq .= uc(substr($$refReads{$last_chr},$last_pos,(length($$refReads{$last_chr}) - $last_pos))); ### TINA
			$qual .= '!' x (length($$refReads{$last_chr}) - $last_pos);

			$uncovered += (length($$refReads{$last_chr}) - $last_pos);
			$total += (length($$refReads{$last_chr}) - $last_pos);
     }
     &p2q_post_process($outfile, $last_chr, \$seq, \$qual, \@gaps, $opts{l});
     $contigs{$last_chr} = length($seq);
	  
	  $genome_uncovered += $uncovered;
	  $genome_total += $total;
	  print "...(line $counter) $last_chr $last_pos ::: total $total uncovered $uncovered ",$uncovered/$total," (genome $genome_uncovered $genome_total)\n";
	  ($uncovered,$total) = (0,0);
  } 

  return ($genome_total,$genome_uncovered,%contigs);
}



sub p2q_post_process {
  my ($outfile, $chr, $seq, $qual, $gaps, $l) = @_;
  die "p2q_post_process ERROR ($chr): lengths not equal (seq ",length($$seq),"!=",length($$qual),"\n"
    if (length($$seq)!=length($$qual));
  open (OUT,">>$outfile") || die "ERROR Can't create $outfile: $!\n";
  print OUT "\@$chr\n"; 
  my $printout = &p2q_print_str($chr,'seq',$seq); print OUT "$printout\n";

  print OUT "+$chr\n"; 
  $printout = &p2q_print_str($chr,'qual',$qual); print OUT "$printout\n";
  close OUT;
}



sub p2q_filter_gaps {
  my ($seq, $gaps, $l) = @_;

  for my $g (@$gaps) {
	my $x = $g > $l? $g - $l : 0;
	substr($$seq, $x, $l + $l) = lc(substr($$seq, $x, $l + $l));
  }
}



sub p2q_print_str {
  my ($chr,$str_qual,$s) = @_;
  my $printout = '';

  my $l = length($$s);
  for (my $i = 0; $i < $l; $i += 60) {
    $printout .= substr($$s, $i, 60);
  }

  return $printout;
}
