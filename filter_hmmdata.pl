#!/usr/bin/env perl
use warnings;
use strict;
use Pod::Usage;
use Getopt::Long qw(GetOptions);
Getopt::Long::Configure qw(gnu_getopt);

#Original script by Molly Schumer
#Edited for automated integration into MSG by Patrick Reilly

=pod

=head1 NAME

filter_hmmdata.pl - Filter hmmdata file for MSG to exclude AIMs within a read length of each other
=head1 SYNOPSIS

filter_hmmdata.pl [options] <Input .hmmdata file>

 Options:
  --help,-h,-?          Display this help documentation
  --output_file,-o      Name for output filtered .hmmdata file
                        Default: [input file prefix].filtered.hmmdata
  --overwrite,-f        Force overwriting the input file with the output
  --read_length,-l      Read length, used as a static threshold distance between AIMs
                        Default: 100
                     
 Mandatory:
  Input .hmmdata file   Path to the input (un-filtered) .hmmdata file

=head1 DESCRIPTION

filter_hmmdata.pl omits sites from the .hmmdata file where two ancestry informative markers
(AIMs) are within a read length of each other. This is a rough heuristic for higher depth
read sets to avoid the bias induced by multiple AIM calls from the same read. The heuristic
built into MSG is overzealous on high depth read sets because it calls reads mapping within
5 bp of each other a single "contig" and only allows one AIM per "contig". Thus, a stretch
of 10 reads mapping within 5 bp of each other could span 1 kb, limiting resolution to 1 AIM
per kb in that region.

=cut

my $help = 0;
my $man = 0;
my $input_path = "";
my $output_path = "";
my $read_length = 100;
my $overwrite = 0;
GetOptions('output_file|o=s' => \$output_path, 'overwrite|f' => \$overwrite, 'read_length|l=i' => \$read_length, 'help|h|?' => \$help, man => \$man) or pod2usage(2);
pod2usage(-exitval => 1, -output => \*STDERR) if $help;
pod2usage(-exitval => 0, -verbose => 2, -output => \*STDERR) if $man;

#Read in the mandatory input .hmmdata file argument:
if (scalar @ARGV < 1 or ! -e $ARGV[0]) { #Not enough mandatory arguments
   print STDERR "Missing input .hmmdata file.\n";
   exit 2;
} else {
   $input_path = $ARGV[0];
}

if ($output_path eq "") {
   $output_path = $input_path;
   $output_path =~ s/\.hmmdata/\.filtered\.hmmdata/;
}

#Old calling information (no longer applies with GNU getopt ):
#example: perl filter_hmmdata.pl indiv0_ATGACA-group1.hmmdata 100 hmm_out 0/1
#enforces a distance of 1 read between AIMs
#argv 2  should be read length. For most recent TN5 data, this should be 150.
# 0 switch does not overwrite the original and 1 switch will overwrite the original with the new file

open IN, "<", $input_path or die "Can't open input .hmmdata file.\n";

my $position = 0;
my $site_current = 0;
my $last_site = 0; #Note: This rules out any AIMs within a read length of the start of the scaffold -- edge case

open OUT, ">", $output_path;

#Print a header for the .hmmdata file:
print OUT "pos\tref\tcons\treads\tquals\tA\tC\tG\tT\tN\tbad\tpar1ref\tpar2ref\n";

#Iterate over the sites in the input .hmmdata file:
while (my $line=<IN>) {
   chomp $line;

   my @elements = split(/\t/, $line);
   $position = $elements[0];
   my $par1 = $elements[11];
   my $par2 = $elements[12];

   #For non-N sites that differ between the parents (AIMs):
   if ($par1 ne $par2 and $par1 ne 'N' and $par2 ne 'N') {
      $site_current = $position;
      #print "$position\n";
      my $distance = $site_current - $last_site; #Off-by-one/fencepost?

      #If the AIMs are farther than a read's length apart, keep the site:
      if ($distance > $read_length){
         print OUT "$line\n";
      }
      $last_site = $site_current;
   }
}
close IN;
close OUT;

if ($overwrite == 1){
   system("mv -f $output_path $input_path");
}
