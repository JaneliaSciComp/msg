#!/usr/bin/env perl
use warnings;
use strict;
use File::Basename;
use FindBin;
use lib $FindBin::Bin;
use Utils;

sub trim { #Derived from Princeton CSES code
    my $s = $_[0];
    ($s) = $s =~ /^[a-zA-Z ]*<?(\d+[\w[\].-]*)>?/;
    return $s;
}

# get msg folder
my $src = dirname $0;

print "\nMSG\n";
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
printf "%4d-%02d-%02d %02d:%02d:%02d\n\n", $year+1900,$mon+1,$mday,$hour,$min,$sec;

### Determine the (relative or absolute, it doesn't really matter) path to the "msg" directory:
my $pathToMSG = dirname(__FILE__);
### Error out if some typical files can't be found in the above-determined directory
### since if they can't be found, the MSG install or placement of msgCluster.pl is atypical.
die "Something went wrong while trying to find MSG.  We tried this directory: ${pathToMSG}/, but it didn't work." unless (-e "${pathToMSG}/msg.pl" and -e "${pathToMSG}/parse_and_map.py" and -e "${pathToMSG}/parent1or2-hmm.sh" and -e "${pathToMSG}/summaryPlots.R");

### Make sure all required dependencies are installed
&Utils::test_dependencies($pathToMSG);

### Default parameters
### All of these parameters are required
my %default_params = (
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
        recRate        => '0',
        rfac          => '0.00001',
        thinfac	      => '1',
        difffac	      => '.01',
        priors         => '0,.5,.5',
        bwaindex1      => 'bwtsw',
        bwaindex2      => 'bwtsw',
        pnathresh      => '0.03',
        cluster        => '1',
        threads        => '8',
        theta        => '1',
        addl_qsub_option_for_exclusive_node => '',
        addl_qsub_option_for_pe => '',
        custom_qsub_options_for_all_cmds => '',
        bwa_alg => 'aln',
        bwa_threads => '1',
        use_stampy => '0',
        stampy_premap_w_bwa => '1',
        stampy_pseudo_threads => '0',
        quality_trim_reads_thresh => '0',
        quality_trim_reads_consec => '30',
        indiv_stampy_substitution_rate => '0.001',
        indiv_mapq_filter => '0',
        index_file => '',
        index_barcodes => '',
        email_host => '',
        notify_emails => '',
        debug => '0',
        gff_thresh_conf => '.95',
        new_parser => '0',
        new_parser_offset => '0',
        new_parser_filter_out_seq => '',
        pepthresh => '',
        one_site_per_contig => '1',
        full_summary_plots => '1',
        plot_lod_matrix => '1',
        max_mapped_reads => '0',
        filter_hmmdata => '0',
        read_length => '100',
        AS_XS_threshold => '6',
        samtools_path => 'samtools',
        barcodes_per_job => 1, #Allows for parallelization of msgRun2
        submit_cmd => 'qsub -N $jobname -cwd $params{"addl_qsub_option_for_exclusive_node"}$params{"custom_qsub_options_for_all_cmds"}-b y -V -sync n', #Default is for SGE, but can be customized, e.g. for SLURM
        #Example submit_cmd for SLURM: sbatch -J $jobname -o $logdir/$jobname.%j.stdout -e $logdir/$jobname.%j.stderr
        #Example submit_cmd for LSF: bsub -J ${jobname} -o ${logdir}/${jobname}.stdout -e ${logdir}/${jobname}.stder
        depend_arg => '-hold_jid $prev_jobname', #Default is for SGE
        #Example for SLURM: -d afterok:${prev_jobid}
        #Example for LSF: -w "done(${prev_jobid})"
        array_job_depend_arg => '-hold_jid $prev_jobname', #Not sure if this is redundant with depend_arg
        array_job_arg => '-t 1-${array_size}', #Default is for SGE
        #For SLURM: -a 1-${array_size}
        #For LSF, things are a bit wonky: -J "${jobname}[1-${array_size}]"
        #But this may conflict with -J specified in submit_cmd, so not tested
        array_job_variable => '$SGE_TASK_ID', #Default is for SGE
        #For SLURM: $SLURM_ARRAY_TASK_ID
        #For LSF: $LSB_JOBINDEX
        default_submit_options => '', #Leave blank for qsub, use with Slurm, applies to msgRun1, 2a, 2b, and 4-7
        #PFR might use on SLURM: -N 1 --ntasks-per-node=1 -t 1:00:00 --mem=4000
        #Stern tested on LSF with: -n 1 -W 60 -M 4000
        msgrun2_submit_options => '', #Leave blank for SGE, use with Slurm
        #Be sure to match up # cores (SLURM --cpus-per-task) with # threads above
        msgrun3_submit_options => '', #Leave blank for SGE, use with Slurm
        #msgRun3 is single-threaded (R generally is), so set # cores to 1 here
        #Example options for SLURM: -N 1 --ntasks-per-node=1 -t 1:00:00 --mem=4000
        #This example allocates one node, runs one process per node, with a timeout of 1 hour, and 4 GB RAM.
        verbose => 0,
        msgRun1 => 1, #Optionally run msgRun1, set to 0 to skip
        msgRun2 => 1, #Optionally run msgRun2, set to 0 to skip
        msgRunOther => 1 #Optionally run other steps of MSG, set to 0 to skip
    );

my $params = Utils::parse_config('msg.cfg', \%default_params);
Utils::validate_config($params, qw( barcodes reads parent1 parent2 ));
my %params = %$params;

#Fill in the rest of array_job_arg if only flag provided:
#e.g. if array_job_arg=-a for SLURM (by checking if it lacks whitespace),
# append " 1-${array_size}" to array_job_arg
#Note: Intended for backwards compatibility, and does not apply to LSF
$params{'array_job_arg'} .= ' 1-${array_size}' unless $params{'array_job_arg'} =~ /\s+/ or $params{'array_job_arg'} eq "";

#Fix to account for use of -J in LSF array jobs:
my $array_submit_cmd = $params{'submit_cmd'};
$array_submit_cmd =~ s/-J\s+\S+// if $params{'array_job_arg'} =~ /\s+/;

#Create a logfile directory for detailed logging:
&Utils::system_call("mkdir -p logs.$$");
my $logdir = "logs.$$"; #This variable can be used in msg.cfg within any of the cluster submission option strings
#e.g. For Slurm, submit_cmd might be "sbatch -J $jobname -o $logdir/$jobname.%j.stdout -e $logdir/$jobname.%j.stderr"
#This would output the msgRun* STDOUT and STDERR logs into the $logdir directory.
#Detailed logs are always placed in this directory though.

### check if all the desired chroms are found in both parental files
### report their lengths also
my %par1_reads = &Utils::readFasta($params{'parent1'}, 1);
my %par2_reads = &Utils::readFasta($params{'parent2'}, 1);
my @chroms;
#when people ask to analyze only part of the genome, then the size of the genome
#changes in the analysis! So, rather than trying to explain all this to end users,
#we decided to see if we can simply trick msg into always using the full genome size.
@chroms = keys %par1_reads; #<- This is that trick

my $numcontigs = scalar(@chroms);

open (OUT, '>', 'msg.chrLengths') || die "ERROR (msgCluster): Can't create msg.chrLengths: $!\n";
print OUT "chr,length\n";
foreach my $chr (sort @chroms) { print OUT "$chr,$par1_reads{$chr}\n"; } 
close OUT;

####################################################################################################
### Parsing
open (OUT, ">", "msgRun1.$$.sh");
print OUT "#!/bin/bash\n/bin/hostname\n/bin/date\n",
    "perl $src/msg.pl ",
    ' --barcodes ', $params{'barcodes'},
    ' --re_cutter ', $params{'re_cutter'},
    ' --linker_system ', $params{'linker_system'},
    ' --reads ', $params{'reads'},
    ' --bwaindex1 ', $params{'bwaindex1'},
    ' --bwaindex2 ', $params{'bwaindex2'},
    ' --bwa_alg ', $params{'bwa_alg'},
    ' --bwa_threads ', $params{'bwa_threads'},
    ' --use_stampy ', $params{'use_stampy'},
    ' --stampy_premap_w_bwa ', $params{'stampy_premap_w_bwa'},
    ' --parent1 ', $params{'parent1'},
    ' --parent2 ', $params{'parent2'},
    ' --indiv_stampy_substitution_rate ', $params{'indiv_stampy_substitution_rate'},
    ' --indiv_mapq_filter ', $params{'indiv_mapq_filter'},
    ' --quality_trim_reads_thresh ', $params{'quality_trim_reads_thresh'},
    ' --quality_trim_reads_consec ', $params{'quality_trim_reads_consec'},
    ' --new_parser ', $params{'new_parser'},
    ' --new_parser_offset ', $params{'new_parser_offset'},
    ' --parse_or_map parse-only',
    ' --samtools_path ', $params{'samtools_path'};

if ($params{'new_parser_filter_out_seq'}) {
    print OUT ' --new_parser_filter_out_seq ', $params{'new_parser_filter_out_seq'};
}
if ($params{'index_file'} && $params{'index_barcodes'}) {
    print OUT ' --index_file ', $params{'index_file'}, ' --index_barcodes ', $params{'index_barcodes'};
}
print OUT ' --logfile_directory ', $logdir;
print OUT " || exit 100\n";
    
close OUT;
&Utils::system_call("chmod 755 msgRun1.$$.sh");

### Replace barcodes file if using Illumina indexing since we will now have num indexes * num barcodes 
### barcoded individuals from parsing step
if ($params{'index_file'} && $params{'index_barcodes'}) {
    &Utils::system_call(
        "python $src/barcode_splitter.py --make_indexed_msg_barcodes_file --msg_barcodes " . $params{'barcodes'} .
        " --bcfile " . $params{'index_barcodes'} . " > logs.$$/barcode_splitter.$$.stdout 2> logs.$$/barcode_splitter.$$.stderr");
    $params{'barcodes'} = $params{'barcodes'} . '.after.index.parsing';
}

### Mapping & Plotting
### qsub array: one for each line in the barcode file
my $num_barcodes = 0;
open(FILE, "<", $params{'barcodes'}) || die "ERROR (msgCluster): Can't open $params{'barcodes'}: $!\n";
while (<FILE>) { chomp $_;
	     if ($_ =~ /^\S+\t.*$/) {
            $num_barcodes ++;
	     }
} close FILE;

print "num barcodes is $num_barcodes!\n";

# Note we include some parsing parameters here since the new style parser operates
# at the begining of msgRun2.
open (OUT,">", "msgRun2.$$.sh");

print OUT "#!/bin/bash\n/bin/hostname\n/bin/date\n";
if ($params{'cluster'}) {
   print OUT "let start=\"(($params{'array_job_variable'} - 1) * $params{'barcodes_per_job'}) + 1\"\n",
        "let max_arrayid=\"\$start + $params{'barcodes_per_job'} - 1\"\n",
        "end=\"\$((\$max_arrayid<$num_barcodes?\$max_arrayid:$num_barcodes))\"\n",
        "for ((h=\$start; h<=\$end; h++)); do\n",
        "   sed -n \"\${h} p\" $params{'barcodes'} > $params{'barcodes'}.$$.\$h\n", #Prevent multi-run collision of barcode files
        "   perl $src/msg.pl ",
        " --barcodes $params{'barcodes'}.$$.\$h"; #Prevent multi-run collision of barcode files
} else {
   print OUT "perl $src/msg.pl --barcodes $params{'barcodes'}";
}
   print OUT ' --reads ', $params{'reads'},
   ' --parent1 ', $params{'parent1'},
   ' --parent2 ', $params{'parent2'},
   ' --chroms ', $params{'chroms'},
   ' --sexchroms ', $params{'sexchroms'},
   ' --chroms2plot ', $params{'chroms2plot'},
   ' --parse_or_map map-only',
   ' --deltapar1 ', $params{'deltapar1'},
   ' --deltapar2 ', $params{'deltapar2'},
   ' --recRate ', $params{'recRate'},
   ' --rfac ', $params{'rfac'},
   ' --priors ', $params{'priors'},
   ' --theta ', $params{'theta'},
   ' --bwa_alg ', $params{'bwa_alg'},
   ' --bwa_threads ', $params{'bwa_threads'},
   ' --use_stampy ', $params{'use_stampy'},
   ' --stampy_premap_w_bwa ', $params{'stampy_premap_w_bwa'},
   ' --indiv_stampy_substitution_rate ', $params{'indiv_stampy_substitution_rate'},
   ' --indiv_mapq_filter ', $params{'indiv_mapq_filter'},
   ' --gff_thresh_conf ', $params{'gff_thresh_conf'},
   ' --new_parser ', $params{'new_parser'},
   ' --new_parser_offset ', $params{'new_parser_offset'},
   ' --re_cutter ', $params{'re_cutter'},
   ' --linker_system ', $params{'linker_system'},
   ' --quality_trim_reads_thresh ', $params{'quality_trim_reads_thresh'},
   ' --quality_trim_reads_consec ', $params{'quality_trim_reads_consec'},
   ' --one_site_per_contig ', $params{'one_site_per_contig'},
   ' --new_parser_filter_out_seq ', ($params{'new_parser_filter_out_seq'} || 'null'),
   ' --pepthresh ', ($params{'pepthresh'} || 'null'),
   ' --max_mapped_reads ', $params{'max_mapped_reads'},
   ' --filter_hmmdata ', $params{'filter_hmmdata'},
   ' --read_length ', $params{'read_length'},
   ' --repeat_threshold ', $params{'AS_XS_threshold'},
   ' --samtools_path ', $params{'samtools_path'},
   ' --logfile_directory ', $logdir;
if ($params{'cluster'}) {
   print OUT " || exit 100\ndone\n";
}
        
close OUT;
&Utils::system_call("chmod 755 msgRun2.$$.sh");


####################################################################################################
mkdir "msgOut.$$" unless (-d "msgOut.$$");
mkdir "msgError.$$" unless (-d "msgError.$$");

### Run jobs!
#Added parameters to generalize the cluster submission process
my ($jobname, $jobid, $prev_jobname, $prev_jobid);

#TODO: Clean up the following job submission code

if ($params{'msgRun1'}) { #Make msgRun1 optional
    if ($params{'cluster'} != 0) {
        $jobname = "msgRun1.$$";
        $jobid = &Utils::system_call(eval "qq{$params{'submit_cmd'} $params{'default_submit_options'} ./msgRun1.$$.sh}");
        $jobid = trim($jobid);
        if ($params{'verbose'}) { #Added minor verbosity
            print "Submitted $jobname (jobid $jobid)";
        }
        $prev_jobname = $jobname;
        $prev_jobid = $jobid;
    }
    else {
       &Utils::system_call("./msgRun1.$$.sh > logs.$$/msgRun1.$$.stdout 2> logs.$$/msgRun1.$$.stderr") ; 
    }
}

#The cluster submission generalization code here was written by members of Princeton CSES
if ($params{'msgRun2'}) { #Make msgRun2 optional
   if ($params{'cluster'} != 0) {
      my $array_size = int((${num_barcodes} + $params{'barcodes_per_job'} - 1) / $params{'barcodes_per_job'});
      $jobname = "msgRun2.$$";
      if ($params{'msgRun1'}) {
         $jobid = &Utils::system_call(eval "qq{$array_submit_cmd $params{'msgrun2_submit_options'} $params{'depend_arg'} $params{'array_job_arg'} ./msgRun2.$$.sh}");
      } else {
         $jobid = &Utils::system_call(eval "qq{$array_submit_cmd $params{'msgrun2_submit_options'} $params{'array_job_arg'} ./msgRun2.$$.sh}");
      }
      $jobid = trim($jobid);
      $prev_jobname = $jobname;
      $prev_jobid = $jobid;
   } else {
      &Utils::system_call("./msgRun2.$$.sh > logs.$$/msgRun2.$$.stdout 2> logs.$$/msgRun2.$$.stderr");
   }
}

if (!$params{'msgRunOther'}) {
   goto end;
}

if ($params{'cluster'} != 0) {
   $jobname = "msgRun2a.$$";
   &Utils::wrap_cmdline("$jobname.sh", "python $src/create_stats.py -i $params{'reads'} -b $params{'barcodes'}");
   $jobid = $params{'msgRun2'} ? &Utils::system_call(eval "qq{$params{'submit_cmd'} $params{'default_submit_options'} $params{'array_job_depend_arg'} ./${jobname}.sh}") : &Utils::system_call(eval "qq{$params{'submit_cmd'} $params{'default_submit_options'} ./${jobname}.sh}");
   $jobid = trim($jobid);
   if ($params{'verbose'}) {
      print "Submitted $jobname (jobid $jobid)\n";
   }
   if ($params{'pepthresh'} ne '') {
      $jobname = "msgRun2b.$$";
      &Utils::wrap_cmdline("$jobname.sh", "python $src/hmmprob_to_est.py -d hmm_fit -t $params{'pepthresh'} -o hmm_fits_est.csv");
      $jobid = $params{'msgRun2'} ? &Utils::system_call(eval "qq{$params{'submit_cmd'} $params{'default_submit_options'} $params{'array_job_depend_arg'} ./${jobname}.sh}") : &Utils::system_call(eval "qq{$params{'submit_cmd'} $params{'default_submit_options'} ./${jobname}.sh}");
      $jobid = trim($jobid);
      if ($params{'verbose'}) {
         print "Submitted $jobname (jobid $jobid)\n";
      }
   }
   $jobname = "msgRun3.$$";
   if ($params{'full_summary_plots'} == 1) {
      &Utils::wrap_cmdline("$jobname.sh", "Rscript $src/summaryPlots.R -c $params{'chroms'} -p $params{'chroms2plot'} -d hmm_fit -t $params{'thinfac'} -f $params{'difffac'} -b $params{'barcodes'} -n $params{'pnathresh'} -l $params{'plot_lod_matrix'}");
   }
   else {
      &Utils::wrap_cmdline("$jobname.sh", "python $src/combine.py -d hmm_fit");
   }
   $jobid = $params{'msgRun2'} ? &Utils::system_call(eval "qq{$params{'submit_cmd'} $params{'msgrun3_submit_options'} $params{'array_job_depend_arg'} ./${jobname}.sh}") : &Utils::system_call(eval "qq{$params{'submit_cmd'} $params{'msgrun3_submit_options'} ./${jobname}.sh}");
   $jobid = trim($jobid);
   if ($params{'verbose'}) {
      print "Submitted $jobname (jobid $jobid)\n";
   }
   $prev_jobname = $jobname;
   $prev_jobid = $jobid;
   
   $jobname = "msgRun4.$$";
   &Utils::wrap_cmdline("$jobname.sh", "perl $src/summary_mismatch.pl $params{'barcodes'} 0");
   $jobid = &Utils::system_call(eval "qq{$params{'submit_cmd'} $params{'default_submit_options'} $params{'depend_arg'} ./${jobname}.sh}");
   $jobid = trim($jobid);
   if ($params{'verbose'}) {
      print "Submitted $jobname (jobid $jobid)\n";
   }
   $prev_jobname = $jobname;
   $prev_jobid = $jobid;
   #Run a simple validation
   $jobname = "msgRun5.$$";
   &Utils::wrap_cmdline("$jobname.sh", "python $src/validate.py $params{'barcodes'}");
   $jobid = &Utils::system_call(eval "qq{$params{'submit_cmd'} $params{'default_submit_options'} $params{'depend_arg'} ./${jobname}.sh}");
   $jobid = trim($jobid);
   if ($params{'verbose'}) {
      print "Submitted $jobname (jobid $jobid)\n";
   }
   $prev_jobname = $jobname;
   $prev_jobid = $jobid;
   #Cleanup - move output files to folders, remove barcode related files
   $jobname = "msgRun6.$$";
   &Utils::wrap_cmdline("$jobname.sh", "mv -f msgRun*.${$}.e** msgError.$$; mv -f msgRun*.${$}.pe** msgError.$$; mv -f msgRun*.${$}.o* msgOut.$$; mv -f msgRun*.${$}.po* msgOut.$$; mv -f *.trim.log msgOut.$$; truncate -s0 temp.fq; rm -f $params{'barcodes'}.*");
   $jobid = &Utils::system_call(eval "qq{$params{'submit_cmd'} $params{'default_submit_options'} $params{'depend_arg'} ./${jobname}.sh}");
   $jobid = trim($jobid);
   if ($params{'verbose'}) {
      print "Submitted $jobname (jobid $jobid)\n";
   }
   $prev_jobname = $jobname;
   $prev_jobid = $jobid;
   #Notify users that MSG run has completed
   if ($params{'email_host'} && $params{'notify_emails'}) {
      $jobname = "msgRun7.$$";
      &Utils::wrap_cmdline("$jobname.sh", "python $src/send_email.py -e $params{'email_host'} -t $params{'notify_emails'} -s \"MSG Run $$ has completed\" -b \"NOTE: Output and error messages are located in msgOut.$$, msgError.$$, and logs.$$\"");
      $jobid = &Utils::system_call(eval "qq{$params{'submit_cmd'} $params{'default_submit_options'} $params{'depend_arg'} ./${jobname}.sh}"); #Not sure why depend_arg is used here instead of array_job_depend_arg
      $jobid = trim($jobid);
      if ($params{'verbose'}) {
         print "Submitted $jobname (jobid $jobid)\n";
      }
#        &Utils::system_call("qsub -N msgRun7.$$ -hold_jid msgRun6.$$ -cwd $params{'custom_qsub_options_for_all_cmds'}-b y -V -sync n python msg/send_email.py -e $params{'email_host'}" .
#            " -t $params{'notify_emails'} -s \\\"MSG Run has completed\\\"" .
#            " -b \\\"NOTE: Output and error messages are located in: msgOut.$$ and msgError.$$\\\""
#            );
   }
} else { 
   &Utils::system_call("./msgRun2.$$.sh > logs.$$/msgRun2.$$.stdout 2> logs.$$/msgRun2.$$.stderr");
   &Utils::system_call("python $src/create_stats.py -i $params{'reads'} -b $params{'barcodes'} > logs.$$/create_stats.$$.stdout 2> logs.$$/create_stats.$$.stderr");
   if ($params{'pepthresh'} ne '') {
       &Utils::system_call("python $src/hmmprob_to_est.py -d hmm_fit -t $params{'pepthresh'} -o hmm_fits_ests.csv > logs.$$/hmmprob_to_est.$$.stdout 2> logs.$$/hmmprob_to_est.$$.stderr");
   }
   if ($params{'full_summary_plots'} == 1) {
        &Utils::system_call("Rscript $src/summaryPlots.R -c $params{'chroms'} -p $params{'chroms2plot'} -d hmm_fit -t $params{'thinfac'} -f $params{'difffac'} -b $params{'barcodes'} -n $params{'pnathresh'} -l $params{'plot_lod_matrix'} > logs.$$/msgRun3.$$.stdout 2> logs.$$/msgRun3.$$.stderr");
   }
   else {
        &Utils::system_call("python $src/combine.py -d hmm_fit > logs.$$/combine.$$.stdout 2> logs.$$/combine.$$.stderr");
   }
   &Utils::system_call("perl $src/summary_mismatch.pl $params{'barcodes'} 0 > logs.$$/summary_mismatch.$$.stdout 2> logs.$$/summary_mismatch.$$.stderr");
   #Run a simple validation
   &Utils::system_call("python $src/validate.py $params{'barcodes'} $params{'full_summary_plots'} > logs.$$/validate.$$.stdout 2> logs.$$/validate.$$.stderr");
   #Cleanup - move output files to folders, remove barcode related files
   &Utils::system_call("mv -f msgRun*.${$}.e** msgError.$$; mv -f msgRun*.${$}.pe** msgError.$$; mv -f msgRun*.${$}.o* msgOut.$$; mv -f msgRun*.${$}.po* msgOut.$$; mv -f *.trim.log msgOut.$$; rm -f temp.fq; rm -f $params{'barcodes'}.*");
   #Notify users that MSG run has completed
   if ($params{'email_host'} && $params{'notify_emails'}) {
     &Utils::system_call("python $src/send_email.py -e $params{'email_host'}" .
        " -t $params{'notify_emails'} -s 'MSG Run has completed'" .
        " -b 'NOTE: Output and error messages are located in msgOut.$$, msgError.$$, and logs.$$ folders.'"
        );
   }
}

end:
print "\nNOTE: Output and error messages are located in msgOut.$$, msgError.$$, and logs.$$ folders.\n\n";
exit;
