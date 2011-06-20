#!/usr/bin/perl -w
### 04.27.10
### Tina
### run msg on cetus
use strict;

print "\nMSG\n";
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
printf "%4d-%02d-%02d %02d:%02d:%02d\n\n", $year+1900,$mon+1,$mday,$hour,$min,$sec;

### Default parameters
die "ERROR: Can't locate msg.cfg.\n" unless (-e 'msg.cfg');
my %params = (
         barcodes       => 'NULL',
         re_cutter      => 'MseI',
         linker_system  => 'Dros_SR_vII',
	      reads          => 'NULL',
	      parent1        => 'NULL',
	      parent2        => 'NULL',
	      chroms         => 'all',
	      sexchroms      => 'X',
	      chroms2plot    => 'all',
	      deltapar1      => '.01',
	      deltapar2      => '.01',
	      recRate        => '3',
	      rfac		      => '1',
	      thinfac	      => '1',
	      difffac	      => '.01',
	      priors         => '0,.5,.5',
	      bwaindex1      => 'bwtsw',
	      bwaindex2      => 'bwtsw',
	      pnathresh      => '0.03',
	      cluster        => '1',
	      threads        => '1',
	      queue          => '1day'
    );

open (IN,'msg.cfg') || die "ERROR: Can't open msg.cfg: $!\n";
while (<IN>) { chomp $_;
	next if ($_ =~ /^\#/);
	next unless ($_);
   $_ =~ s/^\s+//;
   $_ =~ s/\s+$//;
	my ($key,$val) = split(/=/,$_);
	$params{$key} = $val;
} close IN;
$params{'chroms2plot'} = $params{'chroms'} unless (defined $params{'chroms2plot'});
my $update_nthreads = $params{'threads'} if (defined $params{'threads'}); ## Number of BWA threads when updating genomes (must match msg.pl)


### check if all files exist
foreach my $param (qw( barcodes reads parent1 parent2 parent1_reads parent2_reads )) {
	if (exists $params{$param}) {
		die "Exiting from msgCluster: Missing file $params{$param}.\n" unless (-e $params{$param});
	}
}

### double check if the minimum exist
foreach my $key (sort keys %params) {
    die "ERROR (msgCluster): undefined parameter ($key) in msg.cfg.\n" unless ($params{$key} ne 'NULL');
    print "$key:\t$params{$key}\n" ;
}
print "\n" ;


### check if all the desired chroms are found in both parental files
### report their lengths also
my %par1_reads = &readFasta($params{'parent1'});
my %par2_reads = &readFasta($params{'parent2'});
my @chroms ;
if( $params{'chroms'} eq 'all') { 
	@chroms = keys %par1_reads ; 
} else { @chroms = split(/,/,$params{'chroms'}); }

my $numcontigs = length(@chroms) ;

open (OUT,'>msg.chrLengths') || die "ERROR (msgCluster): Can't create msg.chrLengths: $!\n";
print OUT "chr,length\n";
foreach my $chr (sort @chroms) { print OUT "$chr,$par1_reads{$chr}\n"; } 
close OUT;

if (exists $params{'parent1_reads'}) {
    $params{'update_minQV'} = 1 unless (exists $params{'update_minQV'});
    open (OUT,'>msgRun0-1.sh');
    print OUT "/bin/hostname\n/bin/date\n" .
	'perl msg/msg.pl ' .
	' --update ' .
	' --barcodes ' . $params{'barcodes'} .
	' --reads ' . $params{'reads'} . 
	' --parent1 ' . $params{'parent1'} . 
	' --update_minQV ' . $params{'update_minQV'} .
	' --parent1-reads ' . $params{'parent1_reads'} .
	' --threads ' . $params{'threads'} .
	' --bwaindex1 ' . $params{'bwaindex1'} .
	' --bwaindex2 ' . $params{'bwaindex2'} .
	" || exit 100\n";
    close OUT;
    system("chmod 755 msgRun0-1.sh");
    $params{'parent1'} .= '.msg.updated.fasta';
}


if (exists $params{'parent2_reads'}) {
    $params{'update_minQV'} = 1 unless (exists $params{'update_minQV'});
    open (OUT,'>msgRun0-2.sh');
    print OUT "/bin/hostname\n/bin/date\n" .
	'perl msg/msg.pl ' .
	' --update ' .
	' --barcodes ' . $params{'barcodes'} .
	' --reads ' . $params{'reads'} . 
	' --parent2 ' . $params{'parent2'} . 
	' --update_minQV ' . $params{'update_minQV'} .
	' --parent2-reads ' . $params{'parent2_reads'} .
	' --threads ' . $params{'threads'} .
	' --bwaindex1 ' . $params{'bwaindex1'} .
	' --bwaindex2 ' . $params{'bwaindex2'} .
	" || exit 100\n";
    close OUT;
    system("chmod 755 msgRun0-2.sh");
    $params{'parent2'} .= '.msg.updated.fasta';
}


####################################################################################################
### Parsing
open (OUT,'>msgRun1.sh');
print OUT "/bin/hostname\n/bin/date\n" .
    'perl msg/msg.pl ' .
    ' --barcodes ' . $params{'barcodes'} .
    ' --re_cutter ' . $params{'re_cutter'} .
    ' --linker_system ' . $params{'linker_system'} .
    ' --reads ' . $params{'reads'} . 
    ' --bwaindex1 ' . $params{'bwaindex1'} .
    ' --bwaindex2 ' . $params{'bwaindex2'} .
	 ' --parent1 ' . $params{'parent1'} .
	 ' --parent2 ' . $params{'parent2'} .
	 " --parse_or_map parse-only || exit 100\n";
close OUT;
system("chmod 755 msgRun1.sh");


### Mapping & Plotting
### qsub array: one for each line in the barcode file
chomp(my $tmp = `wc -l $params{'barcodes'}`); $tmp =~ s/^\s*//;
my @tmp2 = split(/\s+/,$tmp);
my $num_barcodes = shift @tmp2;

open (OUT,'>msgRun2.sh');

if ($params{'cluster'} != 0) {
   print OUT "#!/bin/bash\n/bin/hostname\n/bin/date\n" .
       "start=\$SGE_TASK_ID\n\n" .
       "let end=\"\$start + \$SGE_TASK_STEPSIZE - 1\"\n\n" .
       "for ((h=\$start; h<=\$end; h++)); do\n" .
#       "   sed -n '1,2p' $params{'barcodes'} > $params{'barcodes'}.\$h\n" .
#       "   sed -n \"\${h}p\" $params{'barcodes'} >> $params{'barcodes'}.\$h\n" .
       "   sed -n \"\${h}p\" $params{'barcodes'} > $params{'barcodes'}.\$h\n" .
       '   perl msg/msg.pl ' .
       ' --barcodes ' . $params{'barcodes'} . '.$h' .
       ' --reads ' . $params{'reads'} . 
       ' --parent1 ' . $params{'parent1'} . 
       ' --parent2 ' . $params{'parent2'} .
       ' --chroms ' . $params{'chroms'} .
       ' --sexchroms ' . $params{'sexchroms'} .
       ' --chroms2plot ' . $params{'chroms2plot'} .
       ' --parse_or_map map-only' .
       ' --deltapar1 ' . $params{'deltapar1'} .
       ' --deltapar2 ' . $params{'deltapar2'} .
       ' --recRate ' . $params{'recRate'} .
       ' --rfac ' . $params{'rfac'} .
       ' --priors ' . $params{'priors'} .
       " || exit 100\ndone\n" .
       "/bin/date\n";
} else {
   print OUT "#!/bin/bash\n/bin/hostname\n/bin/date\n" .
       '   perl msg/msg.pl ' .
       ' --barcodes ' . $params{'barcodes'} .
       ' --reads ' . $params{'reads'} . 
       ' --parent1 ' . $params{'parent1'} . 
       ' --parent2 ' . $params{'parent2'} .
       ' --chroms ' . $params{'chroms'} .
       ' --sexchroms ' . $params{'sexchroms'} .
       ' --chroms2plot ' . $params{'chroms2plot'} .
       ' --parse_or_map map-only' .
       ' --deltapar1 ' . $params{'deltapar1'} .
       ' --deltapar2 ' . $params{'deltapar2'} .
       ' --recRate ' . $params{'recRate'} .
       ' --rfac ' . $params{'rfac'} .
       ' --priors ' . $params{'priors'} .
       "\n";
    }
close OUT;
system("chmod 755 msgRun2.sh");


####################################################################################################
mkdir "msgOut.$$" unless (-d "msgOut.$$");
mkdir "msgError.$$" unless (-d "msgError.$$");

### Run jobs!
if (exists $params{'parent1_reads'} and exists $params{'parent2_reads'}) {
   if ($params{'cluster'} != 0) {
      system("qsub -N msgRun0-1.$$  -pe batch $update_nthreads -cwd -l excl=true -b y -V -sync n ./msgRun0-1.sh") ;
      system("qsub -N msgRun0-2.$$  -pe batch $update_nthreads -cwd -l excl=true -b y -V -sync n ./msgRun0-2.sh") ;
      system("qsub -N msgRun1.$$ -hold_jid msgRun0-1.$$,msgRun0-2.$$ -cwd -b y -V -sync n ./msgRun1.sh") ;
   } else {
      system("./msgRun0-1.sh > msgRun0-1.$$.out 2> msgRun0-1.$$.err") ;
      system("./msgRun0-2.sh > msgRun0-2.$$.out 2> msgRun0-2.$$.err") ;
      system("./msgRun1.sh > msgRun1.$$.out 2> msgRun1.$$.err") ;
   }

} elsif ( exists $params{'parent1_reads'} ) {
   if ($params{'cluster'} != 0) {
      system("qsub -N msgRun0-1.$$ -pe batch $update_nthreads -cwd -l excl=true -b y -V -sync n ./msgRun0-1.sh") ;
      system("qsub -N msgRun1.$$ -hold_jid msgRun0-1.$$ -cwd -b y -V -sync n ./msgRun1.sh") ;
   } else {
      system("./msgRun0-1.sh > msgRun0-1.$$.out 2> msgRun0-1.$$.err") ;
      system("./msgRun1.sh > msgRun1.$$.out 2> msgRun1.$$.err") ;
   }

} elsif ( exists $params{'parent2_reads'} ) {
   if ($params{'cluster'} != 0) {
      system("qsub -N msgRun0-2.$$ -pe batch $update_nthreads -cwd -l excl=true -b y -V -sync n ./msgRun0-2.sh") ;
      system("qsub -N msgRun1.$$ -hold_jid msgRun0-2.$$ -cwd -b y -V -sync n ./msgRun1.sh") ;
   } else {
      system("./msgRun0-2.sh > msgRun0-2.$$.out 2> msgRun0-2.$$.err") ;
      system("./msgRun1.sh > msgRun1.$$.out 2> msgRun1.$$.err") ;
   }
}
else { 
   if ($params{'cluster'} != 0) { system("qsub -N msgRun1.$$ -cwd -b y -V -sync n ./msgRun1.sh") ; }
   else { system("./msgRun1.sh > msgRun1.$$.out 2> msgRun1.$$.err") ; }
}

if ($params{'cluster'} != 0) {
   system("qsub -N msgRun2.$$ -hold_jid msgRun1.$$ -cwd -l excl=true -b y -V -sync n -t 1-${num_barcodes}:1 ./msgRun2.sh");
   #system("qsub -N msgRun2.$$ -hold_jid msgRun1.$$ -cwd -b y -V -sync n -t 3-${num_barcodes}:1 ./msgRun2.sh");
   system("qsub -N msgRun3.$$ -hold_jid msgRun2.$$ -cwd -l excl=true -b y -V -sync n Rscript msg/summaryPlots.R -c $params{'chroms'} -p $params{'chroms2plot'} -d hmm_fit -t $params{'thinfac'} -f $params{'difffac'} -b $params{'barcodes'} -n $params{'pnathresh'}");
   system("qsub -N msgRun4.$$ -hold_jid msgRun3.$$ -cwd -b y -V -sync n perl msg/summary_mismatch.pl $params{'barcodes'} 0");
   system("qsub -N msgRun5.$$ -hold_jid msgRun4.$$ -cwd -b y -V -sync n \"mv -f msgRun*.${$}.e** msgError.$$; mv -f msgRun*.${$}.o* msgOut.$$; rm -f $params{'barcodes'}.*\"");
} else { 
   system("./msgRun2.sh > msgRun2.$$.out 2> msgRun2.$$.err");
   system("Rscript msg/summaryPlots.R -c $params{'chroms'} -p $params{'chroms2plot'} -d hmm_fit -t $params{'thinfac'} -f $params{'difffac'} -b $params{'barcodes'} -n $params{'pnathresh'} > msgRun3.$$.out 2> msgRun3.$$.err");
   system("perl msg/summary_mismatch.pl $params{'barcodes'} 0");
   system("mv -f msgRun*.${$}.e** msgError.$$; mv -f msgRun*.${$}.o* msgOut.$$; rm -f $params{'barcodes'}.*");
}

print "\nNOTE: Output and error messages are located in: msgOut.$$ and msgError.$$ \n\n";
exit;


####################################################################################################
####################################################################################################
sub readFasta {
    my ($file) = @_;
    
    my %reads;
    my ($read,$seq);
    open(FILE,$file) || die "ERROR (msgCluster): Can't open $file: $!\n";
    while (<FILE>) { chomp $_;
		     if ($_ =~ /^>(\S+)/) {
			 $reads{$read} = length($seq) if ($seq);
			 $read = $1;
			 $seq = '';
		     } else { $seq .= $_; }
    } close FILE;
    $reads{$read} = length($seq) if ($read);

    return %reads;  
}
