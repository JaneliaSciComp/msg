#!/usr/bin/perl
### 04.29.10
### Tina
### outputs fasta file for samtools

use Getopt::Std;
getopt("io");

$usage = <<END;
	-i input fasta file
	-o output fasta file
END

die $usage unless (-e $opt_i && $opt_o);

my %reads = &readFasta($opt_i);

open OUT, ">$opt_o" || die "ERROR (reformatFasta4sam): Could not create $opt_o: $!\n";
foreach my $read (sort keys %reads) {
	print OUT ">$read\n" . &p2q_print_str(\$reads{$read});
} close OUT;


exit;


###########################################################################
###########################################################################
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


sub readFasta {
	my ($file) = @_;
	
	my %reads;
	my ($read,$seq);
	open(FILE,$file) || die "ERROR (reformatFasta2sam): Can't open $file: $!\n";
	while (<FILE>) { chomp $_;
		if ($_ =~ /^>(\S+)/) {
			$reads{$read} = $seq if ($seq);
			$read = $1;
			$seq = '';
		} else { $seq .= $_; }
	} close FILE;
	$reads{$read} = $seq if ($read);

	return %reads;  
}
