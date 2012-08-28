#!/usr/bin/perl -w
use strict;
use lib qw(./msg .);
use Utils;

print "\nMSG - Update Parentals\n";
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
printf "%4d-%02d-%02d %02d:%02d:%02d\n\n", $year+1900,$mon+1,$mday,$hour,$min,$sec;

### Make sure all required dependencies are installed
&Utils::test_dependencies();

### Default parameters
### All of these parameters are required
my %default_params = (
        parent1        => 'NULL',
        parent2        => 'NULL',
        parent1_reads  => 'NULL',
        parent2_reads  => 'NULL',
        bwaindex1      => 'bwtsw',
        bwaindex2      => 'bwtsw',
        cluster        => '1',
        threads        => '8',
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
        parent1_stampy_substitution_rate => '0.001',
        parent2_stampy_substitution_rate => '0.001',
        parent1_mapq_filter => '0',
        parent2_mapq_filter => '0',
        debug => '0',
    );

my $params = Utils::parse_config('update.cfg', \%default_params);
Utils::validate_config($params, qw( parent1 parent2 parent1_reads parent2_reads));
my %params = %$params;

####################################################################################################

if (exists $params{'parent1_reads'}) {
    $params{'update_minQV'} = 1 unless (exists $params{'update_minQV'});
    open (OUT,'>msgRunU0-1.sh');
    print OUT "/bin/hostname\n/bin/date\n" .
    'perl msg/msg.pl ' .
    ' --update ' .
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
    ' --parent_mapq_filter ' . $params{'parent1_mapq_filter'} .
    ' --debug ' . $params{'debug'}
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
    &Utils::system_call("chmod 755 msgRunU0-1.sh");
    $params{'parent1'} .= '.msg.updated.fasta';
}


if (exists $params{'parent2_reads'}) {
    $params{'update_minQV'} = 1 unless (exists $params{'update_minQV'});
    open (OUT,'>msgRunU0-2.sh');
    print OUT "/bin/hostname\n/bin/date\n" .
    'perl msg/msg.pl ' .
    ' --update ' .
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
    ' --parent_mapq_filter ' . $params{'parent2_mapq_filter'} .
    ' --debug ' . $params{'debug'}
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
    &Utils::system_call("chmod 755 msgRunU0-2.sh");
    $params{'parent2'} .= '.msg.updated.fasta';
}

####################################################################################################
mkdir "msgUpdateParentalsOut.$$" unless (-d "msgUpdateParentalsOut.$$");
mkdir "msgUpdateParentalsError.$$" unless (-d "msgUpdateParentalsError.$$");

### Run jobs!

if (exists $params{'parent1_reads'} and exists $params{'parent2_reads'}) {
   if ($params{'cluster'} != 0) {
        &Utils::system_call("qsub -N msgRunU0-1.$$ $params{'addl_qsub_option_for_pe'}-cwd $params{'addl_qsub_option_for_exclusive_node'}-b y -V -sync n ./msgRunU0-1.sh") ;
        &Utils::system_call("qsub -N msgRunU0-2.$$ $params{'addl_qsub_option_for_pe'}-cwd $params{'addl_qsub_option_for_exclusive_node'}-b y -V -sync n ./msgRunU0-2.sh") ;
        #Cleanup - move output files to folders
        &Utils::system_call("qsub -N msgRunU6.$$ -hold_jid msgRunU0-1.$$,msgRunU0-2.$$ -cwd -b y -V -sync n \"mv -f msgRunU*.${$}.e** msgUpdateParentalsError.$$; mv -f msgRunU*.${$}.pe** msgUpdateParentalsError.$$; mv -f msgRunU*.${$}.o* msgUpdateParentalsOut.$$; mv -f msgRunU*.${$}.po* msgUpdateParentalsOut.$$; mv -f *.trim.log msgUpdateParentalsOut.$$\"");
   } else {
        &Utils::system_call("./msgRunU0-1.sh > msgRunU0-1.$$.out 2> msgRunU0-1.$$.err") ;
        &Utils::system_call("./msgRunU0-2.sh > msgRunU0-2.$$.out 2> msgRunU0-2.$$.err") ;
        #Cleanup - move output files to folders
        &Utils::system_call("mv -f msgRunU*.${$}.e** msgUpdateParentalsError.$$; mv -f msgRunU*.${$}.o* msgUpdateParentalsOut.$$; mv -f *.trim.log msgUpdateParentalsOut.$$");
   }
} elsif ( exists $params{'parent1_reads'} ) {
    if ($params{'cluster'} != 0) {
        &Utils::system_call("qsub -N msgRunU0-1.$$ $params{'addl_qsub_option_for_pe'}-cwd $params{'addl_qsub_option_for_exclusive_node'}-b y -V -sync n ./msgRunU0-1.sh") ;
        #Cleanup - move output files to folders
        &Utils::system_call("qsub -N msgRunU6.$$ -hold_jid msgRunU0-1.$$ -cwd -b y -V -sync n \"mv -f msgRunU*.${$}.e** msgUpdateParentalsError.$$; mv -f msgRunU*.${$}.pe** msgUpdateParentalsError.$$; mv -f msgRunU*.${$}.o* msgUpdateParentalsOut.$$; mv -f msgRunU*.${$}.po* msgUpdateParentalsOut.$$; mv -f *.trim.log msgUpdateParentalsOut.$$\"");
    } else {
        &Utils::system_call("./msgRunU0-1.sh > msgRunU0-1.$$.out 2> msgRunU0-1.$$.err") ;
        #Cleanup - move output files to folders
        &Utils::system_call("mv -f msgRunU*.${$}.e** msgUpdateParentalsError.$$; mv -f msgRunU*.${$}.o* msgUpdateParentalsOut.$$; mv -f *.trim.log msgUpdateParentalsOut.$$");
    }
} elsif ( exists $params{'parent2_reads'} ) {
   if ($params{'cluster'} != 0) {
        &Utils::system_call("qsub -N msgRunU0-2.$$ $params{'addl_qsub_option_for_pe'}-cwd $params{'addl_qsub_option_for_exclusive_node'}-b y -V -sync n ./msgRunU0-2.sh") ;
        #Cleanup - move output files to folders
        &Utils::system_call("qsub -N msgRunU6.$$ -hold_jid msgRunU0-2.$$ -cwd -b y -V -sync n \"mv -f msgRunU*.${$}.e** msgUpdateParentalsError.$$; mv -f msgRunU*.${$}.pe** msgUpdateParentalsError.$$; mv -f msgRunU*.${$}.o* msgUpdateParentalsOut.$$; mv -f msgRunU*.${$}.po* msgUpdateParentalsOut.$$; mv -f *.trim.log msgUpdateParentalsOut.$$\"");
   } else {
        &Utils::system_call("./msgRunU0-2.sh > msgRunU0-2.$$.out 2> msgRunU0-2.$$.err") ;
        #Cleanup - move output files to folders
        &Utils::system_call("mv -f msgRunU*.${$}.e** msgUpdateParentalsError.$$; mv -f msgRunU*.${$}.o* msgUpdateParentalsOut.$$; mv -f *.trim.log msgUpdateParentalsOut.$$");
   }
}

print "\nNOTE: Output and error messages are located in: msgUpdateParentalsOut.$$ and msgUpdateParentalsError.$$ \n\n";
if (exists $params{'parent1_reads'}) { print "Output file is $params{'parent1'}\n";}
if (exists $params{'parent2_reads'}) { print "Output file is $params{'parent2'}\n";}

exit;
