#!/usr/bin/perl -w
### 04.27.10
### Tina
### run msg on cetus
use strict;
use lib qw(./msg .);
use Utils;

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
        rfac		      => '0.00001',
        thinfac	      => '1',
        difffac	      => '.01',
        priors         => '0,.5,.5',
        bwaindex1      => 'bwtsw',
        bwaindex2      => 'bwtsw',
        pnathresh      => '0.03',
        cluster        => '1',
        threads        => '8',
        theta        => '1',
        min_coverage => '2',
        max_coverage_exceeded_state => 'N',
        addl_qsub_option_for_exclusive_node => '',
        addl_qsub_option_for_pe => '',
        bwa_alg => 'aln',
        bwa_threads => '1',
        use_stampy => '0',
        stampy_premap_w_bwa => '1',
        stampy_pseudo_threads => '0',
        quality_trim_reads_thresh => '0',
        quality_trim_reads_consec => '30',
        indiv_stampy_substitution_rate => '0.001',
        parent1_stampy_substitution_rate => '0.001',
        parent2_stampy_substitution_rate => '0.001',
        indiv_mapq_filter => '0',
        parent1_mapq_filter => '0',
        parent2_mapq_filter => '0',
        index_file => '',
        index_barcodes => '',
    );

open (IN,'msg.cfg') || die "ERROR: Can't open msg.cfg: $!\n";
while (<IN>) { chomp $_;
	next if ($_ =~ /^\#/);
	next unless ($_);
	my ($key,$val) = split(/=/,$_,2);
	$params{Utils::strip($key)} = Utils::strip($val);
} close IN;

### Configure some parameters ###
$params{'chroms2plot'} = $params{'chroms'} unless (defined $params{'chroms2plot'});
my $update_nthreads = $params{'threads'} if (defined $params{'threads'}); ## Number of qsub slots when running pe option
#add space after qsub options so we can insert into commands, add thread/slot count to -pe option
if (defined $params{'addl_qsub_option_for_exclusive_node'} && $params{'addl_qsub_option_for_exclusive_node'}) {
    #example: go from user msg.cfg entered "-l excl=true" to "-l excl=true "
    $params{'addl_qsub_option_for_exclusive_node'} = $params{'addl_qsub_option_for_exclusive_node'}.' ';
}
else {
    $params{'addl_qsub_option_for_exclusive_node'} = '';
}
if (defined $params{'addl_qsub_option_for_pe'} && $params{'addl_qsub_option_for_pe'}) {
    #example: go from user msg.cfg entered "-pe batch" to "-pe batch 8 "
    $params{'addl_qsub_option_for_pe'} = $params{'addl_qsub_option_for_pe'}." $update_nthreads ";
}
else {
    $params{'addl_qsub_option_for_pe'} = '';
}

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
my %par1_reads = &Utils::readFasta($params{'parent1'}, 1);
my %par2_reads = &Utils::readFasta($params{'parent2'}, 1);
my @chroms ;
if( $params{'chroms'} eq 'all') { 
	@chroms = keys %par1_reads ; 
} else { @chroms = split(/,/,$params{'chroms'}); }

my $numcontigs = length(@chroms) ;

open (OUT,'>msg.chrLengths') || die "ERROR (msgCluster): Can't create msg.chrLengths: $!\n";
print OUT "chr,length\n";
foreach my $chr (sort @chroms) { print OUT "$chr,$par1_reads{$chr}\n"; } 
close OUT;

#Trim barcoded reads for quality if required, set reads.filtered paramater which all downstream calculations will look at.
# (skip this step for illumina indexed reads.  We will do it after splitting the indexes.)
my $do_run_reads_trim = 0;
if (($params{'quality_trim_reads_thresh'} > 0) && (!$params{'index_file'}) && (!$params{'index_barcodes'})) {
    $params{'reads.filtered'} = $params{'reads'} . '.trim.fastq.gz';
    if (-e $params{'reads.filtered'}) {
        print "Skipping quality trim step since $params{'reads.filtered'} already exists\n";
    }
    else {
        open (OUT,'>msgRun0-0.sh');
        print OUT "python msg/TQSfastq.py" . " -f " . $params{'reads'} . " -t " . $params{'quality_trim_reads_thresh'} .
                " -c " . $params{'quality_trim_reads_consec'} . " -q " . " -o " . $params{'reads'} . " -z";
        close OUT;
        &Utils::system_call("chmod 755 msgRun0-0.sh");
        $do_run_reads_trim = 1;
    }
}
else {
    #make new name same as old name
    $params{'reads.filtered'} = $params{'reads'};
}

if (exists $params{'parent1_reads'}) {
    $params{'update_minQV'} = 1 unless (exists $params{'update_minQV'});
    open (OUT,'>msgRun0-1.sh');
    print OUT "/bin/hostname\n/bin/date\n" .
	'perl msg/msg.pl ' .
	' --update ' .
	' --barcodes ' . $params{'barcodes'} .
	' --reads ' . $params{'reads.filtered'} . 
	' --parent1 ' . $params{'parent1'} . 
	' --update_minQV ' . $params{'update_minQV'} .
	' --min_coverage ' . $params{'min_coverage'} .
	' --max_coverage_exceeded_state ' . $params{'max_coverage_exceeded_state'} .    
	' --parent1-reads ' . $params{'parent1_reads'} .
	' --bwaindex1 ' . $params{'bwaindex1'} .
	' --bwaindex2 ' . $params{'bwaindex2'} .
	' --bwa_alg ' . $params{'bwa_alg'} .
	' --bwa_threads ' . $params{'bwa_threads'} .
    ' --use_stampy ' . $params{'use_stampy'} .
    ' --stampy_premap_w_bwa ' . $params{'stampy_premap_w_bwa'} .
    ' --stampy_pseudo_threads ' . $params{'stampy_pseudo_threads'} .
    ' --cluster ' . $params{'cluster'} .
    ' --quality_trim_reads_thresh ' . $params{'quality_trim_reads_thresh'} .
    ' --quality_trim_reads_consec ' . $params{'quality_trim_reads_consec'} .
    ' --parent_stampy_substitution_rate ' . $params{'parent1_stampy_substitution_rate'} .
    ' --parent_mapq_filter ' . $params{'parent1_mapq_filter'}
    ;
    #Add on optional arguments
    if (defined $params{'max_coverage_stds'}) {
        print OUT ' --max_coverage_stds ' . $params{'max_coverage_stds'};
    }
    if ($params{'addl_qsub_option_for_pe'}) {
        print OUT ' --addl_qsub_option_for_pe pe';
    }
	print OUT " || exit 100\n";
    close OUT;
    &Utils::system_call("chmod 755 msgRun0-1.sh");
    $params{'parent1'} .= '.msg.gz.updated.fasta';
}


if (exists $params{'parent2_reads'}) {
    $params{'update_minQV'} = 1 unless (exists $params{'update_minQV'});
    open (OUT,'>msgRun0-2.sh');
    print OUT "/bin/hostname\n/bin/date\n" .
	'perl msg/msg.pl ' .
	' --update ' .
	' --barcodes ' . $params{'barcodes'} .
	' --reads ' . $params{'reads.filtered'} . 
	' --parent2 ' . $params{'parent2'} . 
	' --update_minQV ' . $params{'update_minQV'} .
	' --min_coverage ' . $params{'min_coverage'} .
	' --max_coverage_exceeded_state ' . $params{'max_coverage_exceeded_state'} .  
	' --parent2-reads ' . $params{'parent2_reads'} .
	' --bwaindex1 ' . $params{'bwaindex1'} .
	' --bwaindex2 ' . $params{'bwaindex2'} .
	' --bwa_alg ' . $params{'bwa_alg'} .
	' --bwa_threads ' . $params{'bwa_threads'} .
    ' --use_stampy ' . $params{'use_stampy'} .
    ' --stampy_premap_w_bwa ' . $params{'stampy_premap_w_bwa'} .
    ' --stampy_pseudo_threads ' . $params{'stampy_pseudo_threads'} .
    ' --cluster ' . $params{'cluster'} .
    ' --quality_trim_reads_thresh ' . $params{'quality_trim_reads_thresh'} .
    ' --quality_trim_reads_consec ' . $params{'quality_trim_reads_consec'} .
    ' --parent_stampy_substitution_rate ' . $params{'parent2_stampy_substitution_rate'} .
    ' --parent_mapq_filter ' . $params{'parent2_mapq_filter'}
    ;
    #Add on optional arguments
    if (defined $params{'max_coverage_stds'}) {
        print OUT ' --max_coverage_stds ' . $params{'max_coverage_stds'};
    }
    if ($params{'addl_qsub_option_for_pe'}) {
        print OUT ' --addl_qsub_option_for_pe pe';
    }    
	print OUT " || exit 100\n";
    close OUT;
    &Utils::system_call("chmod 755 msgRun0-2.sh");
    $params{'parent2'} .= '.msg.gz.updated.fasta';
}


####################################################################################################
### Parsing
open (OUT,'>msgRun1.sh');
print OUT "/bin/hostname\n/bin/date\n" .
    'perl msg/msg.pl ' .
    ' --barcodes ' . $params{'barcodes'} .
    ' --re_cutter ' . $params{'re_cutter'} .
    ' --linker_system ' . $params{'linker_system'} .
    ' --reads ' . $params{'reads.filtered'} . 
    ' --bwaindex1 ' . $params{'bwaindex1'} .
    ' --bwaindex2 ' . $params{'bwaindex2'} .
    ' --bwa_alg ' . $params{'bwa_alg'} .
    ' --bwa_threads ' . $params{'bwa_threads'} .
    ' --use_stampy ' . $params{'use_stampy'} .
    ' --stampy_premap_w_bwa ' . $params{'stampy_premap_w_bwa'} .
    ' --parent1 ' . $params{'parent1'} .
    ' --parent2 ' . $params{'parent2'} .
    ' --indiv_stampy_substitution_rate ' . $params{'indiv_stampy_substitution_rate'} .
    ' --indiv_mapq_filter ' . $params{'indiv_mapq_filter'} .
    ' --quality_trim_reads_thresh ' . $params{'quality_trim_reads_thresh'} .
    ' --quality_trim_reads_consec ' . $params{'quality_trim_reads_consec'} .
    " --parse_or_map parse-only";
if ($params{'index_file'} && $params{'index_barcodes'}) {
    print OUT ' --index_file ' . $params{'index_file'} . ' --index_barcodes ' . $params{'index_barcodes'};
}
print OUT " || exit 100\n";
    
close OUT;
&Utils::system_call("chmod 755 msgRun1.sh");

### Replace barcodes file if using Illumina indexing since we will now have num indexes * num barcodes 
### barcoded individuals from parsing step
if ($params{'index_file'} && $params{'index_barcodes'}) {
    &Utils::system_call(
        "python msg/barcode_splitter.py --make_indexed_msg_barcodes_file --msg_barcodes " . $params{'barcodes'} .
        " --bcfile " . $params{'index_barcodes'});
    $params{'barcodes'} = $params{'barcodes'} . '.after.index.parsing';
}

### Mapping & Plotting
### qsub array: one for each line in the barcode file
my $num_barcodes = 0;
open(FILE,$params{'barcodes'}) || die "ERROR (msgCluster): Can't open $params{'barcodes'}: $!\n";
while (<FILE>) { chomp $_;
	     if ($_ =~ /^\S+\t.*$/) {
            $num_barcodes ++;
	     }
} close FILE;

print "num barcodes is $num_barcodes!\n";

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
        ' --reads ' . $params{'reads.filtered'} . 
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
        ' --theta ' . $params{'theta'} .
        ' --bwa_alg ' . $params{'bwa_alg'} .
        ' --bwa_threads ' . $params{'bwa_threads'} .
        ' --use_stampy ' . $params{'use_stampy'} .
        ' --stampy_premap_w_bwa ' . $params{'stampy_premap_w_bwa'} .
        ' --indiv_stampy_substitution_rate ' . $params{'indiv_stampy_substitution_rate'} .
        ' --indiv_mapq_filter ' . $params{'indiv_mapq_filter'} .
        " || exit 100\ndone\n" .
        "/bin/date\n";
} else {
   print OUT "#!/bin/bash\n/bin/hostname\n/bin/date\n" .
        '   perl msg/msg.pl ' .
        ' --barcodes ' . $params{'barcodes'} .
        ' --reads ' . $params{'reads.filtered'} . 
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
        ' --theta ' . $params{'theta'} .
        ' --bwa_alg ' . $params{'bwa_alg'} .
        ' --bwa_threads ' . $params{'bwa_threads'} .
        ' --use_stampy ' . $params{'use_stampy'} .
        ' --stampy_premap_w_bwa ' . $params{'stampy_premap_w_bwa'} .
        ' --indiv_stampy_substitution_rate ' . $params{'indiv_stampy_substitution_rate'} .
        ' --indiv_mapq_filter ' . $params{'indiv_mapq_filter'} .
       "\n";
    }
close OUT;
&Utils::system_call("chmod 755 msgRun2.sh");


####################################################################################################
mkdir "msgOut.$$" unless (-d "msgOut.$$");
mkdir "msgError.$$" unless (-d "msgError.$$");

### Run jobs!

#Trim barcoded reads for quality
my $hold_jid_for_trim_text = '';
if ($do_run_reads_trim) {
    if ($params{'cluster'} != 0) {
        &Utils::system_call("qsub -N msgRun0-0.$$ $params{'addl_qsub_option_for_pe'}-cwd $params{'addl_qsub_option_for_exclusive_node'}-b y -V -sync n ./msgRun0-0.sh") ;
        $hold_jid_for_trim_text = "-hold_jid msgRun0-0.$$ ";
    }
    else {
        &Utils::system_call("./msgRun0-0.sh > msgRun0-0.$$.out 2> msgRun0-0.$$.err");
    }
}

if (exists $params{'parent1_reads'} and exists $params{'parent2_reads'}) {
   if ($params{'cluster'} != 0) {
      &Utils::system_call("qsub -N msgRun0-1.$$ $hold_jid_for_trim_text $params{'addl_qsub_option_for_pe'}-cwd $params{'addl_qsub_option_for_exclusive_node'}-b y -V -sync n ./msgRun0-1.sh") ;
      &Utils::system_call("qsub -N msgRun0-2.$$ $hold_jid_for_trim_text $params{'addl_qsub_option_for_pe'}-cwd $params{'addl_qsub_option_for_exclusive_node'}-b y -V -sync n ./msgRun0-2.sh") ;
      &Utils::system_call("qsub -N msgRun1.$$ -hold_jid msgRun0-1.$$,msgRun0-2.$$ -cwd -b y -V -sync n ./msgRun1.sh") ;
   } else {
      &Utils::system_call("./msgRun0-1.sh > msgRun0-1.$$.out 2> msgRun0-1.$$.err") ;
      &Utils::system_call("./msgRun0-2.sh > msgRun0-2.$$.out 2> msgRun0-2.$$.err") ;
      &Utils::system_call("./msgRun1.sh > msgRun1.$$.out 2> msgRun1.$$.err") ;
   }

} elsif ( exists $params{'parent1_reads'} ) {
   if ($params{'cluster'} != 0) {
      &Utils::system_call("qsub -N msgRun0-1.$$ $hold_jid_for_trim_text $params{'addl_qsub_option_for_pe'}-cwd $params{'addl_qsub_option_for_exclusive_node'}-b y -V -sync n ./msgRun0-1.sh") ;
      &Utils::system_call("qsub -N msgRun1.$$ -hold_jid msgRun0-1.$$ -cwd -b y -V -sync n ./msgRun1.sh") ;
   } else {
      &Utils::system_call("./msgRun0-1.sh > msgRun0-1.$$.out 2> msgRun0-1.$$.err") ;
      &Utils::system_call("./msgRun1.sh > msgRun1.$$.out 2> msgRun1.$$.err") ;
   }

} elsif ( exists $params{'parent2_reads'} ) {
   if ($params{'cluster'} != 0) {
      &Utils::system_call("qsub -N msgRun0-2.$$ $hold_jid_for_trim_text $params{'addl_qsub_option_for_pe'}-cwd $params{'addl_qsub_option_for_exclusive_node'}-b y -V -sync n ./msgRun0-2.sh") ;
      &Utils::system_call("qsub -N msgRun1.$$ -hold_jid msgRun0-2.$$ -cwd -b y -V -sync n ./msgRun1.sh") ;
   } else {
      &Utils::system_call("./msgRun0-2.sh > msgRun0-2.$$.out 2> msgRun0-2.$$.err") ;
      &Utils::system_call("./msgRun1.sh > msgRun1.$$.out 2> msgRun1.$$.err") ;
   }
}
else { 
   if ($params{'cluster'} != 0) { &Utils::system_call("qsub -N msgRun1.$$ $hold_jid_for_trim_text -cwd -b y -V -sync n ./msgRun1.sh") ; }
   else { &Utils::system_call("./msgRun1.sh > msgRun1.$$.out 2> msgRun1.$$.err") ; }
}

if ($params{'cluster'} != 0) {
   &Utils::system_call("qsub -N msgRun2.$$ -hold_jid msgRun1.$$ -cwd $params{'addl_qsub_option_for_exclusive_node'}-b y -V -sync n -t 1-${num_barcodes}:1 ./msgRun2.sh");
   #&Utils::system_call("qsub -N msgRun2.$$ -hold_jid msgRun1.$$ -cwd -b y -V -sync n -t 3-${num_barcodes}:1 ./msgRun2.sh");
   &Utils::system_call("qsub -N msgRun3.$$ -hold_jid msgRun2.$$ -cwd $params{'addl_qsub_option_for_exclusive_node'}-b y -V -sync n Rscript msg/summaryPlots.R -c $params{'chroms'} -p $params{'chroms2plot'} -d hmm_fit -t $params{'thinfac'} -f $params{'difffac'} -b $params{'barcodes'} -n $params{'pnathresh'}");
   &Utils::system_call("qsub -N msgRun4.$$ -hold_jid msgRun3.$$ -cwd -b y -V -sync n perl msg/summary_mismatch.pl $params{'barcodes'} 0");
   #Run a simple validation
   &Utils::system_call("qsub -N msgRun5.$$ -hold_jid msgRun4.$$ -cwd -b y -V -sync n python msg/validate.py $params{'barcodes'}");
   #Cleanup - move output files to folders, remove barcode related files
   &Utils::system_call("qsub -N msgRun6.$$ -hold_jid msgRun5.$$ -cwd -b y -V -sync n \"mv -f msgRun*.${$}.e** msgError.$$; mv -f msgRun*.${$}.o* msgOut.$$; mv -f *.trim.log msgOut.$$; rm -f $params{'barcodes'}.*\"");
} else { 
   &Utils::system_call("./msgRun2.sh > msgRun2.$$.out 2> msgRun2.$$.err");
   &Utils::system_call("Rscript msg/summaryPlots.R -c $params{'chroms'} -p $params{'chroms2plot'} -d hmm_fit -t $params{'thinfac'} -f $params{'difffac'} -b $params{'barcodes'} -n $params{'pnathresh'} > msgRun3.$$.out 2> msgRun3.$$.err");
   &Utils::system_call("perl msg/summary_mismatch.pl $params{'barcodes'} 0");
   #Run a simple validation
   &Utils::system_call("python msg/validate.py $params{'barcodes'} > msgRun.validate.$$.out 2> msgRun.validate.$$.err");
   #Cleanup - move output files to folders, remove barcode related files
   &Utils::system_call("mv -f msgRun*.${$}.e** msgError.$$; mv -f msgRun*.${$}.o* msgOut.$$; mv -f *.trim.log msgOut.$$; rm -f $params{'barcodes'}.*");
}

print "\nNOTE: Output and error messages are located in: msgOut.$$ and msgError.$$ \n\n";
exit;
