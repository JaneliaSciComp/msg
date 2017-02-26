#!/usr/bin/env perl
use warnings;
use strict;
use POSIX qw(ceil);
use File::Basename;
use File::Copy;
use lib qw(./msg .);
use Utils;

my $version = 1.0;

my $src = dirname $0;

my %default_params = (
   barcodes => 'barcode_file',
   priors => '1,0,0',
   chroms => 'all',
   sexchroms => 'X',
   chroms2plot => 'all',
   deltapar1 => '0.1',
   deltapar2 => '0.1',
   recRate => '3',
   rfac => '1',
   theta => '1',
   gff_thresh_conf => '0.95',
   one_site_per_contig => '1',
   pepthresh => '',
   logdir => qq(optimize_logs.$$),
   cluster => '0',
   barcodes_per_job => 1,
   msgrun2_submit_options => '',
   submit_cmd => q(sbatch -J ${jobname} -o $params{'logdir'}/${jobname}.stdout -e $params{'logdir'}/${jobname}.stderr),
   array_job_arg => '--array=',
   array_job_variable => '$SLURM_ARRAY_TASK_ID',
   array_job_depend_arg => '-d afterok:'
);

#Read in the configuration file:
my %params = %{Utils::parse_config('optimize.cfg', \%default_params)};
Utils::validate_config(\%params, qw(barcodes));

#Create the log file directory:
`mkdir -p $params{'logdir'}`;

#Output the input parameters as a record in STDOUT:
print "optimizeParameters version:  $version\n";
print `date`;
print "\n";
print "barcodes file:          ", $params{'barcodes'}, "\n";
print "priors:                 ", $params{'priors'}, "\n";
print "chroms:                 ", $params{'chroms'}, "\n";
print "sexchroms:              ", $params{'sexchroms'}, "\n";
print "chroms2plot:            ", $params{'chroms2plot'}, "\n";
print "deltapar1:              ", $params{'deltapar1'}, "\n";
print "deltapar2:              ", $params{'deltapar2'}, "\n";
print "recRate:                ", $params{'recRate'}, "\n";
print "rfac:                   ", $params{'rfac'}, "\n";
print "theta:                  ", $params{'theta'}, "\n";
print "gff_thresh_conf:        ", $params{'gff_thresh_conf'}, "\n";
print "one_site_per_contig:    ", $params{'one_site_per_contig'}, "\n";
print "pepthresh               ", $params{'pepthresh'}, "\n";
print "logfile_directory:      ", $params{'logdir'}, "\n";
print "cluster:                ", $params{'cluster'}, "\n";
print "barcodes_per_job:       ", $params{'barcodes_per_job'}, "\n";
print "msgrun2_submit_options: ", $params{'msgrun2_submit_options'}, "\n";
print "submit_cmd:             ", $params{'submit_cmd'}, "\n";
print "array_job_arg:          ", $params{'array_job_arg'}, "\n";
print "array_job_variable:     ", $params{'array_job_variable'}, "\n";
print "array_job_depend_arg:   ", $params{'array_job_depend_arg'}, "\n";

my $jobname = '';
my $logdir = $params{'logdir'}; #Quick bugfix for people that use $logdir in their msg.cfg

#Collect the individuals from the barcode file:
open(BARCODE, "<", $params{'barcodes'}) or die "Error: Could not open barcodes file.\n";
my @indivs = ();
while (my $line = <BARCODE>) {
   chomp $line;
   my ($index, $id) = split /\s+/, $line, 2;
   my $indiv = 'indiv' . $id . '_' . $index;
   push(@indivs, $indiv);
}
close(BARCODE);
#Generate a template script for performing the HMM fitting:
open(FIT, ">", "fitHMMpartitions.pl");
print FIT '#!/usr/bin/env perl', "\n",
'open(BARCODE, "<", "', $params{'barcodes'}, '") ',
'or die "Error: Could not open barcodes file.\n";', "\n",
'my @indivs = ();', "\n",
'while (my $line = <BARCODE>) {', "\n",
'   chomp $line;', "\n",
'   my ($index, $id) = split /\s+/, $line, 2;', "\n",
'   my $indiv = "indiv" . $id . "_" . $index;', "\n",
'   push(@indivs, $indiv);', "\n",
'}', "\n",
'close(BARCODE);', "\n";
if ($params{'cluster'}) {
   print FIT 'my $barcodes_per_job = ', $params{'barcodes_per_job'}, ';', "\n",
   'my $offset = (', $params{'array_job_variable'}, '-1)*${barcodes_per_job};', "\n";
} else {
   print FIT 'my $barcodes_per_job = ', scalar@indivs, ';', "\n",
   'my $offset = 0;', "\n";
}
print FIT 'for (my $i = 0; $i < $barcodes_per_job; $i++) {', "\n",
'   print "$indiv\n";', "\n",
'   `bash ', $src, '/fitParentalHMM.sh',
' -b ', $params{'barcodes'},
' -o hmm_data',
' -R hmm_fit',
' -i $indivs[$offset+$i]',
' -c ', $params{'chroms'},
' -x ', $params{'sexchroms'},
' -y ', $params{'chroms2plot'},
' -f ', $params{'deltapar1'},
' -g ', $params{'deltapar2'},
' -z ', $params{'priors'},
' -a ', $params{'recRate'},
' -r ', $params{'rfac'},
' -t ', $params{'theta'},
' -m ', $params{'gff_thresh_conf'},
' -u ', $params{'one_site_per_contig'},
' -j ', $params{'pepthresh'},
' > ', $params{'logdir'}, '/fitParentalHMM$indivs[$offset+$i].msg$$.stdout',
' 2> ', $params{'logdir'}, '/fitParentalHMM$indivs[$offset+$i].msg$$.stderr`;', "\n",
'}', "\n";
close(FIT);
#Run the template script, either using cluster commands, or locally:
$jobname = "fitHMMpartitions";
my ($jobidline, $jobid) = ('', '');
if ($params{'cluster'}) {
   my $array_size = ceil(@indivs/$params{'barcodes_per_job'});
   my $cmd = eval "qq($params{'submit_cmd'} $params{'msgrun2_submit_options'} $params{'array_job_arg'}1-${array_size} ./fitHMMpartitions.pl)";
   $jobidline = `$cmd 2>&1`;
   if ($jobidline =~ /\d+/) {
      print $jobidline;
      $jobid = $1;
      print "fitHMMpartitions.pl submitted as job ${jobid}.\n";
   } else {
      die "fitHMMpartitions.pl failed to submit\n";
   }
} else {
   `./fitHMMpartitions.pl`;
}

#Create the template script for the rest of the tasks:
open(SEGERR, ">", "plotSegCalcError.pl");
print SEGERR '#!/usr/bin/env perl', "\n",
'print "Plotting segregation for both parents:";', "\n",
'`Rscript ', $src, '/plotSegregation.R',
' -d hmm_fit',
' -c ', $params{'chroms'},
' -p ', $params{'chroms2plot'},
' 2> ', $params{'logdir'}, '/plotSegregation.msg$$.stderr',
' > ', $params{'logdir'}, '/plotSegregation.msg$$.stdout`;', "\n";
print SEGERR 'print "Calculating error_gamma:";', "\n",
'`perl ${src}/summary_mismatch.pl ', $params{'barcodes'}, '0',
' 2> ', $params{'logdir'}, '/summaryMismatch.msg$$.stderr',
' > ', $params{'logdir'}, '/summaryMismatch.msg$$.stdout`;';
close(SEGERR);

#Run the template script for the rest of the tasks:
$jobname = "plotSegCalcError";
if ($params{'cluster'}) {
   die "Error: Could not read job ID properly for fitHMMpartitions.pl, so did not submit plotSegCalcError.pl.\n" if $jobid eq '';
   my $cmd = eval "qq($params{'submit_cmd'} $params{'msgrun2_submit_options'} $params{'array_job_depend_arg'} ${jobid} ./plotSegCalcError.pl)";
   $jobidline = `$cmd 2>&1`;
   if ($jobidline =~ /\d+/) {
      $jobid = $1;
      print "plotSegCalcError.pl submitted as job ${jobid}.\n";
   } else {
      die "plotSegCalcError.pl failed to submit\n";
   }
} else {
   `./plotSegCalcError.pl`;
}
