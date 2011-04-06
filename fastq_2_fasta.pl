#!/usr/bin/perl
# edited by Tina to print out 60 sequence chars per line
#
# This program converts Sanger formatted fastq files to fasta files
# The Sanger format encodes quality scores 0-93 as ASCII values 33-126
# It is not compatable with the Illumina fastq format.
#

use Getopt::Std;

getopt("q");

$usage = <<END;
	-q [0] quality score cutoff 

*** EDITED BY TINA 12.03.09
END

die $usage unless (defined $opt_q);

if ($opt_q > 0 && $opt_q <= 93) {
	$opt_q += 33;
}
else {
	$opt_q = 0;
}

while(<>) {
	if ( /^@(\S+)/ ) { ### TINA
		$id = $1;
		$seq = <>;
		chomp $seq;
		<>;
		$qual = <>;
		chomp $qual;
		#assert( length($seq) == length($qual) );
		if ( $opt_q ) {
			for ($i = 0; $i < length($qual); $i++) {
				if ( ord(substr($qual, $i, 1)) < $opt_q ) {
					substr($seq, $i, 1) = "N";
				}
			} 
		}
		if ( $opt_m ) {
			foreach ( @{ $masks{$id} } ) {
				substr($seq, $_->[0]+1, $_->[1]-$_->[0]+1) = "N" x ($_->[1]-$_->[0]+1);
			}
		}

		print ">$id\n" . &p2q_print_str(\$seq); ### TINA
	}
}
exit;


sub p2q_print_str {
  my ($s) = @_;
  my $printout = '';

  my $l = length($$s);
  for (my $i = 0; $i < $l; $i += 60) {
    $printout .= substr($$s, $i, 60);
    $printout .= "\n";
  }

  return $printout;
}
