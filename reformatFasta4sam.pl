#!/usr/bin/perl
### 04.29.10
### Tina
### outputs fasta file for samtools

# For gzipped output instead; we don't need it right now
#use IO::Compress::Gzip qw(gzip $GzipError);

use lib qw(./msg .);
use Utils;

use Getopt::Std;
getopt("io");

$usage = <<END;
	-i input fasta file
	-o output fasta file
END

die $usage unless (-e $opt_i && $opt_o);

my %reads = &Utils::readFasta($opt_i, 0);

open OUT, ">$opt_o" || die "ERROR (reformatFasta4sam): Could not create $opt_o: $!\n";
# For gzipped output instead; we don't need it right now
#my $OUT = new IO::Compress::Gzip $opt_o 
#    or die "ERROR (reformatFasta4sam): Could not create $opt_o: IO::Compress::Gzip failed: $GzipError\n";

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
