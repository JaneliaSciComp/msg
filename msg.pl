#!/usr/bin/env perl
use File::Basename;
use File::Copy;
use Getopt::Long;
use lib qw(./msg .);
use Utils;

my $true = 1;
my $false = 0;
#Read in version from version file
my $version = do { local( @ARGV, $/ ) = './msg/version'; <> };

my $src = dirname $0;
my $update_genomes = $false;
my ($barcodes, $re_cutter, $linker_system, $raw_read_data, $parent1_genome, $parent2_genome, $parent1_reads, $parent2_reads, $update_minQV, $min_coverage, $max_coverage_stds, $max_coverage_exceeded_state,
    $parse_or_map, $priors, $chroms, $sexchroms, $chroms2plot, $deltapar1, $deltapar2, $recRate, $rfac, $bwaindex1, $bwaindex2, $theta, $bwa_alg, $bwa_threads, $use_stampy, $stampy_premap_w_bwa,
    $stampy_pseudo_threads, $cluster, $addl_qsub_option_for_pe, $quality_trim_reads_thresh, $quality_trim_reads_consec, $indiv_stampy_substitution_rate, $parent_stampy_substitution_rate,
    $indiv_mapq_filter, $parent_mapq_filter, $index_file, $index_barcodes, $debug, $gff_thresh_conf, $new_parser, $new_parser_offset, $new_parser_filter_out_seq, $custom_qsub_options_for_all_cmds,
    $one_site_per_contig, $pepthresh, $max_mapped_reads, $logdir, $filter_hmmdata, $read_length);

GetOptions(
    'barcodes|b=s' => \$barcodes,
    're_cutter=s' => \$re_cutter,
    'linker_system=s' => \$linker_system,
    'reads|i=s' => \$raw_read_data,
    'update|u' => \$update_genomes,
    'parent1=s' => \$parent1_genome,
    'parent2=s' => \$parent2_genome,
    'parent1-reads=s' => \$parent1_reads,
    'parent2-reads=s' => \$parent2_reads,
    'update_minQV=i' => \$update_minQV,
    'min_coverage=i' => \$min_coverage,
    'max_coverage_stds=i' => \$max_coverage_stds,
    'max_coverage_exceeded_state=s' => \$max_coverage_exceeded_state,
    'parse_or_map=s' => \$parse_or_map,
    'priors=s' => \$priors,
    'chroms=s' => \$chroms,
    'sexchroms=s' => \$sexchroms,
    'chroms2plot=s' => \$chroms2plot,
    'deltapar1=s' => \$deltapar1,
    'deltapar2=s' => \$deltapar2,
	'recRate=s' => \$recRate,    
    'rfac=s' => \$rfac,
    'bwaindex1=s' => \$bwaindex1,
    'bwaindex2=s' => \$bwaindex2,
    'theta=s' => \$theta,
    'bwa_alg=s' => \$bwa_alg,
    'bwa_threads=i' => \$bwa_threads,
    'use_stampy=i' => \$use_stampy,
    'stampy_premap_w_bwa=i' => \$stampy_premap_w_bwa,
    'stampy_pseudo_threads=i' => \$stampy_pseudo_threads,
    'cluster=i' => \$cluster,
    'addl_qsub_option_for_pe=s' => \$addl_qsub_option_for_pe,
    'quality_trim_reads_thresh=i' => \$quality_trim_reads_thresh,
    'quality_trim_reads_consec=i' => \$quality_trim_reads_consec,
    'indiv_stampy_substitution_rate=f' => \$indiv_stampy_substitution_rate,
    'parent_stampy_substitution_rate=f' => \$parent_stampy_substitution_rate,
    'indiv_mapq_filter=i' => \$indiv_mapq_filter,
    'parent_mapq_filter=i' => \$parent_mapq_filter,
    'index_file=s' => \$index_file,
    'index_barcodes=s' => \$index_barcodes,
    'debug=i' => \$debug,
    'gff_thresh_conf=s' => \$gff_thresh_conf,
    'new_parser=i' => \$new_parser,
    'new_parser_offset=i' => \$new_parser_offset,
    'new_parser_filter_out_seq=s' => \$new_parser_filter_out_seq,
    'custom_qsub_options_for_all_cmds=s' => \$custom_qsub_options_for_all_cmds,
    'one_site_per_contig=i' => \$one_site_per_contig,
    'pepthresh=s' => \$pepthresh,
    'max_mapped_reads=s' => \$max_mapped_reads,
    'logfile_directory=s' => \$logdir,
    'filter_hmmdata=i' => \$filter_hmmdata,
    'read_length=i' => \$read_length
    );

#### INTERNAL OPTIONS (for developers) #####

# Should MD tags be added to mapped SAM files when they are not included by default.
# It's slower to add them but they seem to affect the output.
my $GEN_MD_TAGS = $true;

############################################

print "msg version:  $version\n" ;
print `date` ;
print "\n" ;
print "parent1 genome: $parent1_genome\n" ;
print "parent2 genome: $parent2_genome\n" ;
print "update:         $update_genomes\n" ;
print "parse_or_map:   $parse_or_map\n" ;
print "chroms:         $chroms\n\n" ;
print "sexchroms:      $sexchroms\n\n" ;
print "max_coverage_stds $max_coverage_stds\n\n";
print "max_coverage_exceeded_state $max_coverage_exceeded_state\n\n";
print "bwaindex1 $bwaindex1\n\n";
print "bwaindex2 $bwaindex2\n\n";
print "bwa_alg $bwa_alg\n\n";
print "bwa_threads $bwa_threads\n\n";
print "use_stampy $use_stampy\n\n";
print "stampy_premap_w_bwa $stampy_premap_w_bwa\n\n";
print "stampy_pseudo_threads $stampy_pseudo_threads\n\n";
print "cluster $cluster\n\n";
print "addl_qsub_option_for_pe $addl_qsub_option_for_pe\n\n";
print "quality_trim_reads_thresh $quality_trim_reads_thresh\n\n";
print "quality_trim_reads_consec $quality_trim_reads_consec\n\n";
print "indiv_stampy_substitution_rate $indiv_stampy_substitution_rate\n\n";
print "parent_stampy_substitution_rate $parent_stampy_substitution_rate\n\n";
print "indiv_mapq_filter $indiv_mapq_filter\n\n";
print "parent_mapq_filter $parent_mapq_filter\n\n";
print "index_file $index_file\n\n";
print "index_barcodes $index_barcodes\n\n";
print "debug $debug\n\n";
print "gff_thresh_conf $gff_thresh_conf\n\n";
print "new_parser $new_parser\n\n";
print "new_parser_offset $new_parser_offset\n\n";
print "new_parser_filter_out_seq $new_parser_filter_out_seq\n\n";
print "one_site_per_contig $one_site_per_contig\n\n";
print "pepthresh $pepthresh\n\n";
print "max_mapped_reads $max_mapped_reads\n\n";
#only update parentals passes this in, so you won't see it on standard msg runs:
print "custom_qsub_options_for_all_cmds $custom_qsub_options_for_all_cmds\n\n";
print "logfile_directory $logdir\n\n";
print "filter_hmmdata $filter_hmmdata\n\n";
print "read_length $read_length\n\n";

if( $update_genomes ) {
	print "update genomes params:\n";
	print "   parent1 reads:	$parent1_reads\n" ;
	print "   parent2 reads:	$parent2_reads\n" ;
}
else{
	print "reads:		$raw_read_data\n" ;
	print "barcodes:	$barcodes\n" ;
	$parent1_genome or die "\nERROR: Must use --parent1 option to supply parent1 genome file" ;
	$parent2_genome or die "\nERROR: Must use --parent2 option to supply parent2 genome file" ;
}

my %genomes_fa;
my %genome_index;
$genomes_fa{'parent1'} = $parent1_genome if $parent1_genome;
$genomes_fa{'parent2'} = $parent2_genome if $parent2_genome;
$genome_index{'parent1'} = $bwaindex1;
$genome_index{'parent2'} = $bwaindex2;


sub post_update_cleanup {
    # When not in debug mode remove any files that aren't needed downstream and aren't needed for a re-run
    my ($sp, $parent_reads_for_updating_fq, $out_prefix) = @_;
    #examples from toy run:
    # $sp = parent2
    # $parent_reads_for_updating_fq = parent2_reads.fq.trim.fastq.gz
    # $out_prefix = update_reads-aligned-parent2
    if (!$debug) {
        if ($parent_reads_for_updating_fq =~ /trim\.fastq\.gz$/) {
            unlink($parent_reads_for_updating_fq);
        }
        unlink("$genomes_fa{$sp}.msg.updated.fastq");
        unlink("$genomes_fa{$sp}.msg");
        unlink glob("$out_prefix.*");
    }
}

sub run_stampy_on_cluster {
    # Run stampy on SGE cluster.  Note, design-wise, normally all qsub/cluster calls should be in 
    # msgCluster.pl but it wasn't practical in this case.
    
    my ($sp, $bwa_options, $out, $stampy_pseudo_threads, $reads_for_updating_fq_current_sp, $parent_stampy_substitution_rate) = @_;
    
    #Use 4 slots of the PE options is available since we need ~7GB of memory for each job.
    #We assume PE option is available if the user has specified the option in config file.
    if ($addl_qsub_option_for_pe) {$addl_qsub_option_for_pe = "-pe batch 4";} else {$addl_qsub_option_for_pe = "";}
    
    #write out a shell script since there are too many parameters to call qsub directly.
    open (OUT,">msgRun0-$sp-stampy.sh");
    print OUT "#!/bin/bash\n/bin/hostname\n/bin/date\n" .
        "start=\$SGE_TASK_ID\n\n" .
        "let end=\"\$start + \$SGE_TASK_STEPSIZE - 1\"\n\n" .
        "for ((h=\$start; h<=\$end; h++)); do\n" .
        "   stampy.py " . $bwa_options . 
        " --processpart=\${h}/$stampy_pseudo_threads" .
        " -g $sp.stampy.msg -h $sp.stampy.msg -M $reads_for_updating_fq_current_sp" .
        " --substitutionrate $parent_stampy_substitution_rate" .
        " -o $out.tmp.\${h}.sam" .
        "   || exit 100\n" .
        "   samtools view -bhS -o $out.tmp.\${h}.bam $out.tmp.\${h}.sam\n" .
        "done\n" .
        "/bin/date\n";
    close OUT;
    system("chmod 755 msgRun0-$sp-stampy.sh");
    #Note: using the uncompressed BAM option -u in the samtools view command above may give a speed
    #boost however it causes an "[bam_header_read] EOF marker is absent." error when merging.
    #I think this was fixed in later versions of samtools and this may be harmless here, but it's better
    #to be safe so I took it out.  (reference: http://seqanswers.com/forums/showthread.php?t=15363)
    
    #run it
    &Utils::system_call("qsub", "-t 1-${stampy_pseudo_threads}:1", $addl_qsub_option_for_pe, 
        "-N msgRun0-$sp-stampy.$$", "-V", "-sync y", "-cwd", &Utils::strip($custom_qsub_options_for_all_cmds),
        "-b y", "./msgRun0-$sp-stampy.sh");
    
    #we won't get here until qsub is done (-sync y above)
    #merge all of the temp files back together
    my @tmp_bam_file_names;
    for (1..$stampy_pseudo_threads) {push(@tmp_bam_file_names, "$out.tmp.$_.bam");}
    &Utils::system_call("samtools merge $out.stampy.tmp.bam " . join(" ", @tmp_bam_file_names));
    #convert back to SAM format
    &Utils::system_call("samtools view -h -o $out.sam $out.stampy.tmp.bam");
    
    #delete temp files and qsub rubbish
    if (!$debug) {
        unlink("$out.stampy.tmp.bam");
    }
    &Utils::system_call("rm -f $out.tmp.*");

    unless(-e 'stampy_logs' or mkdir 'stampy_logs') {
        die "Unable to create stampy_logs\n";
    }
    &Utils::system_call("mv -f msgRun0-$sp-stampy.* stampy_logs");
}

## -------------------------------------------------------------------------------
##
## UPDATE REFERENCES
my %reads_for_updating_fq;
my %reads_for_updating;
if( $update_genomes ) {
	$reads_for_updating_fq{'parent1'} = $parent1_reads if ($parent1_reads);
	$reads_for_updating_fq{'parent2'} = $parent2_reads if ($parent2_reads);
	#%reads_for_updating_fq = (parent1 => $parent1_reads, parent2 => $parent2_reads) ;

    #Make sure parents' reads files are specified and exist 
	for my $sp (keys %reads_for_updating_fq) {
		my $reads = $reads_for_updating_fq{$sp} ;
		$reads or die "Must supply --$sp-reads" ;
		-e $reads or die "No such file: $reads" ;
      my ($d, $s);
		($reads_for_updating{$sp}, $d, $s) = fileparse $reads, qr/\.[^.]*$/ ;
		$s eq ".fq" or $s eq ".fastq" or $s eq ".gz" or print "$sp reads extension is $s: expecting .fq or .fastq\n" ;
	}

}
else {
	$raw_read_data or die "Must use --reads option to supply raw data file" ;
	-e $raw_read_data or die "No such file: $raw_read_data" ;

	$barcodes or die "Must use --barcodes option to supply barcodes file" ;
	-e $barcodes or die "No such file: $barcodes" ;

	$parent1_reads and die "Must use --update option with --parent1-reads" ;
	$parent2_reads and die "Must use --update option with --parent2-reads" ;
}

if( $update_genomes ) {
    print "Updating genomes\n";
    foreach my $sp ( keys %reads_for_updating ) {
        if( -e $genomes_fa{$sp}.".msg.updated.fasta" ) {
            print "$sp: updated genome files already present\n" ;
            next ;
        }
        
        #Trim reads if required (you could run this conncurrently while trimming barcoded reads for a speed up in the future.)
        if ($quality_trim_reads_thresh > 0) {
            &Utils::system_call("python","$src/TQSfastq.py","-f",$reads_for_updating_fq{$sp},"-t",$quality_trim_reads_thresh,
                "-c",$quality_trim_reads_consec,"-q","-z","-o",$reads_for_updating_fq{$sp},"> $logdir/TQSfastq${sp}.msg$$.stdout 2> $logdir/TQSfastq${sp}.msg$$.stderr");
            $reads_for_updating_fq{$sp} = $reads_for_updating_fq{$sp} . '.trim.fastq.gz';
        }
    
        &Utils::system_call("perl", "$src/reformatFasta4sam.pl", "-i", $genomes_fa{$sp}, "-o", "$genomes_fa{$sp}.msg", "> $logdir/reformatFasta4sam${sp}.msg$$.stdout 2> $logdir/reformatFasta4sam${sp}.msg$$.stderr") ;
    
        my $out = "update_reads-aligned-$sp" ;
        unless( -e "$out.pileup" ) {
            #Always call index ( though you don't need it if $use_stampy = 1 and $stampy_premap_w_bwa = 0)
            &Utils::system_call("bwa", "index", "-a", $genome_index{$sp}, "$genomes_fa{$sp}.msg", "> $logdir/bwaIndex${sp}.msg$$.stdout 2> $logdir/bwaIndex${sp}.msg$$.stderr") 
                unless( -e "$genomes_fa{$sp}.msg.bwt" and -e "$genomes_fa{$sp}.msg.ann" );
            &Utils::system_call("samtools", "faidx", "$genomes_fa{$sp}.msg", "> $logdir/samtoolsFaidx${sp}.msg$$.stdout 2> $logdir/samtoolsFaidx${sp}.msg$$.stderr") ;
    
            unless (-e "$out.sam") {
                if ($use_stampy == 1) {
                    #Build stampy genome file
                    &Utils::system_call("stampy.py","-G", "$sp.stampy.msg", "$genomes_fa{$sp}.msg", "> $logdir/stampyG${sp}.msg$$.stdout 2> $logdir/stampyG${sp}.msg$$.stderr");
                    #Build stampy hash file
                    &Utils::system_call("stampy.py","-g", "$sp.stampy.msg", "-H", "$sp.stampy.msg", "> $logdir/stampyH${sp}.msg$$.stdout 2> $logdir/stampyH${sp}.msg$$.stderr");
                    #Call stampy mapping with or without bwa
                    my $bwa_options = "";
                    if ($stampy_premap_w_bwa == 1) {
                        $bwa_options = "--bwaoptions=\"-q10 $genomes_fa{$sp}.msg\"";
                    }
                    if (($stampy_pseudo_threads > 0) && ($cluster == 1)) {
                        run_stampy_on_cluster($sp, $bwa_options, $out, $stampy_pseudo_threads, $reads_for_updating_fq{$sp}, $parent_stampy_substitution_rate);
                    }
                    else {
                        #Standard Stampy run w/o qsub
                        &Utils::system_call("stampy.py", $bwa_options, 
                            "-g", "$sp.stampy.msg", "-h", "$sp.stampy.msg", 
                            "--substitutionrate", $parent_stampy_substitution_rate,
                            "-M", $reads_for_updating_fq{$sp}, "-o", "$out.sam", "> $logdir/stampyRun${sp}.msg$$.stdout 2> $logdir/stampyRun${sp}.msg$$.stderr");
                    }                    
                }
                elsif ($bwa_alg eq 'aln') {
                    &Utils::system_call("bwa", $bwa_alg, "-t", $bwa_threads, "$genomes_fa{$sp}.msg", $reads_for_updating_fq{$sp}, "> $out.sai", "2> $logdir/bwaAln${sp}.msg$$.stderr");
                    &Utils::system_call("bwa", "samse", "$genomes_fa{$sp}.msg", "$out.sai", $reads_for_updating_fq{$sp}, "> $out.sam", "2> $logdir/bwaSamse${sp}.msg$$.stderr");
                }
                elsif ($bwa_alg eq 'bwasw') {
                    &Utils::system_call("bwa", $bwa_alg, "-t", $bwa_threads, "$genomes_fa{$sp}.msg", $reads_for_updating_fq{$sp}, "> $out.sam", "2> $logdir/bwaSW${sp}.msg$$.stderr");
                }
                elsif ($bwa_alg eq 'mem') {
                    &Utils::system_call("bwa", $bwa_alg, "-t", $bwa_threads, "$genomes_fa{$sp}.msg", $reads_for_updating_fq{$sp}, "> $out.sam", "2> $logdir/bwaMem${sp}.msg$$.stderr");
                }
                else {die "Invalid bwa_alg parameter and not using stampy: $bwa_alg. Must be aln, mem, or bwasw (or set use_stampy=1)";}

            }
            
            # Filter out unmapped reads, etc
            &Utils::system_call("$src/filter-sam.py", "-i", "$out.sam", "-o", "$out.filtered.sam", "-a", $bwa_alg, "-s", $use_stampy, "> $logdir/filter-sam${sp}.msg$$.stdout 2> $logdir/filter-sam${sp}.msg$$.stderr") ;
            if ($parent_mapq_filter > 0) {
                &Utils::system_call("samtools", "view", "-bt", "$genomes_fa{$sp}.msg.fai", "-q $parent_mapq_filter", "-o $out.bam", "$out.filtered.sam", "> $logdir/samtoolsView${sp}.msg$$.stdout 2> $logdir/samtoolsView${sp}.msg$$.stderr");
            }
            else {
                &Utils::system_call("samtools", "view", "-bt", "$genomes_fa{$sp}.msg.fai", "-o $out.bam", "$out.filtered.sam", "> $logdir/samtoolsView${sp}.msg$$.stdout 2> $logdir/samtoolsView${sp}.msg$$.stderr");
            }
            &Utils::system_call("samtools", "sort", "$out.bam", "$out.bam.sorted", "> $logdir/samtoolsSort${sp}.msg$$.stdout 2> $logdir/samtoolsSort${sp}.msg$$.stderr");
            
            if (($GEN_MD_TAGS == $true) && ($bwa_alg eq 'bwasw' || $use_stampy == 1)) {
                if ($debug == $true) {
                    &Utils::system_call("cp","-f","$out.bam.sorted.bam","$out.bam.sorted.bam.beforecalmd.bam");
                }
                #Put back in MD tags, bwasw and stampy omit them:
                &Utils::system_call("samtools", "calmd", "-b", "$out.bam.sorted.bam", "$genomes_fa{$sp}.msg", "> $out.sorted.calmd.bam", "2> $logdir/samtoolsCalmd${sp}.msg$$.stderr");
                if ($debug == $true) {
                    #copy instead of move for debug mode so we can view each step
                    &Utils::system_call("cp","-f","$out.sorted.calmd.bam","$out.bam.sorted.bam");               
                }
                else {
                    #Can't output calmd/view to overwrite $out.sam since view starts overwriting as data is piped in, so do a move after the fact.
                    &Utils::system_call("mv","$out.sorted.calmd.bam","$out.bam.sorted.bam");
                }
            }
            &Utils::system_call("samtools", "index", "$out.bam.sorted.bam", "> $logdir/samtoolsIndex${sp}.msg$$.stdout 2> $logdir/samtoolsIndex${sp}.msg$$.stderr");
            #When possible, we should be able to switch between pileup and mpileup seamlessly
            &Utils::system_call("samtools", "pileup", "-f", "$genomes_fa{$sp}.msg", "$out.bam.sorted.bam", "-c", "> $out.pileup", "2> $logdir/samtoolsPileup${sp}.msg$$.stderr");

        }

        #Update the parent genome with the newly mapped reads
        #Note: updateRef.pl uses the consensus base from samtools pileup for dealing with indels
        #so unless we rewrite that section to use other data from a pileup generated by mpileup,
        #we won't be able to remove the version constraint of samtools
		&Utils::system_call("perl", "$src/updateRef.pl", "$genomes_fa{$sp}.msg", "$out.pileup", "0", $update_minQV, $min_coverage,
            $max_coverage_exceeded_state, $max_coverage_stds, "> $logdir/updateRef${sp}.msg$$.stdout 2> $logdir/updateRef${sp}.msg$$.stderr");
						
		## Re-run BWA against updated refs
		## Not doing that at the moment
		&post_update_cleanup($sp, $reads_for_updating_fq{$sp}, $out);
	}
	print "Finished updating genomes\n" ;
	exit;
}



## -------------------------------------------------------------------------------
##
## Use BWA to align raw reads against genomes, creating two sam files for each parent ref
my %genomes;
for my $sp (keys %genomes_fa) {
	-e $genomes_fa{$sp} or die "No such file: $genomes_fa{$sp}";

   my ($d, $s);
	($genomes{$sp}, $d, $s) = fileparse $genomes_fa{$sp}, qr/\.[^.]*$/;
	$s eq ".fa" or $s eq ".fasta" or $s eq ".gz" or print "$sp genome extension is $s: expecting .fa or .fasta\n";
	
	if ($parse_or_map eq 'parse-only') {
		### reformat for samtools (60 chars per line)
		&Utils::system_call("perl", "$src/reformatFasta4sam.pl", "-i", "$genomes_fa{$sp}", "-o", "${sp}_ref.fa", "> $logdir/reformatFasta4sam${sp}.msg$$.stdout 2> $logdir/reformatFasta4sam${sp}.msg$$.stderr") unless (-e "${sp}_ref.fa");
	}
}


%genomes_fa = (parent1 => 'parent1_ref.fa', parent2 => 'parent2_ref.fa') ;
$parse_or_map = "--$parse_or_map" if ($parse_or_map);

### BWA INDEXING
foreach my $sp ( keys %genomes ) {
	unless( -e "$genomes_fa{$sp}.bwt" and -e "$genomes_fa{$sp}.ann" ) {
		## 1. Align reads against existing genome -> sam file
		## Always call index ( though you don't need it if $use_stampy = 1 and $stampy_premap_w_bwa = 0)
		#@args = ("bwa", "index", "-a", $genome_index{$sp}, $genomes_fa{$sp}) ;
		print "Running index on reformatted fa:\n";
      #Replacing this system() call with &Utils::system_call() call
      &Utils::system_call("bwa", "index", "-a", $genome_index{$sp}, $genomes_fa{$sp}, "> $logdir/bwaIndex${sp}.msg$$.stdout 2> $logdir/bwaIndex${sp}.msg$$.stderr");
		#print "@args\n" ;
		#system("@args") == 0 or die "Error in @args: $?" ;
		#When mapping with stampy, perform the genome and hash building steps here
        if ($use_stampy == 1) {
            #Build stampy genome file
            &Utils::system_call("stampy.py","-G", "$genomes_fa{$sp}.stampy.msg", $genomes_fa{$sp}, "> $logdir/stampyG${sp}.msg$$.stdout 2> $logdir/stampyG${sp}.msg$$.stderr");
            #Build stampy hash file
            &Utils::system_call("stampy.py","-g", "$genomes_fa{$sp}.stampy.msg", "-H", 
                "$genomes_fa{$sp}.stampy.msg", "> $logdir/stampyH${sp}.msg$$.stdout 2> $logdir/stampyH${sp}.msg$$.stderr");
        }
	}
}


#Map each individual (when called by msgRun2 with map-only, or parse all individuals when
# called by msgRun1 with parse only)
my $samfiles_dir = basename($raw_read_data) . "_sam_files";
mkdir $samfiles_dir unless (-d $samfiles_dir);
&Utils::system_call('python', "$src/parse_and_map.py", '-i', $raw_read_data, '-b', $barcodes,
		 '--parent1', $genomes_fa{'parent1'}, '--parent2', $genomes_fa{'parent2'}, $parse_or_map,
       '--re_cutter', $re_cutter, '--linker_system', $linker_system, '--bwa_alg', 
        $bwa_alg, '--bwa_threads', $bwa_threads, 
        '--use_stampy', $use_stampy, '--stampy_premap_w_bwa', $stampy_premap_w_bwa,
        '--indiv_stampy_substitution_rate', $indiv_stampy_substitution_rate,
        '--indiv_mapq_filter', $indiv_mapq_filter, '--index_file', $index_file,
        '--index_barcodes', $index_barcodes, '--quality_trim_reads_thresh', $quality_trim_reads_thresh || '0',
        '--quality_trim_reads_consec', $quality_trim_reads_consec || '0','--dbg', $debug,
        '--new_parser', $new_parser || '0', '--new_parser_offset', $new_parser_offset || '0',
        '--new_parser_filter_out_seq', $new_parser_filter_out_seq || 'null',
        "> $logdir/parseAndMap.msg$$.stdout 2> $logdir/parseAndMap.msg$$.stderr");

## Strip species out of reference column

## Gzip

## Fit HMM to estimate ancestry

if ($parse_or_map eq '--map-only') {
   print "\nFit HMM to estimate ancestry...\n";

	open (BARCODE,$barcodes) || die "ERROR: Can't open $barcodes: $!\n";
   foreach my $bc_line (<BARCODE>) {
	   chomp $bc_line;

   	# fastq_file = 'indiv' + ind[1] + '_' + ind[0]
   	my @bc_bits = split(/\s+/,$bc_line);
   	my $indiv = 'indiv' . $bc_bits[1] . '_' . $bc_bits[0];
      print "\t$indiv\n";
   
   	&Utils::system_call('bash', "$src/parent1or2-hmm.sh",
   			 '-b', $barcodes,
   			 '-s', $samfiles_dir,
   			 '-o', 'hmm_data',
   			 '-R', 'hmm_fit',
   			 '-p', $genomes_fa{'parent1'},
   			 '-q', $genomes_fa{'parent2'},
   			 '-i', $indiv,
   			 '-c', $chroms,
   			 '-y', $chroms2plot,
   			 '-f', $deltapar1,
   			 '-g', $deltapar2,
             '-a', $recRate,             
   			 '-r', $rfac,
   			 '-x', $sexchroms,
   			 '-z', $priors,
   			 '-t', $theta,
             '-w', $bwa_alg,
             '-e', $use_stampy,
             '-m', $gff_thresh_conf,
             '-u', $one_site_per_contig,
             '-j', $pepthresh,
             '-n', $max_mapped_reads,
             '-v', $filter_hmmdata,
             '-l', $read_length,
   		  "> $logdir/parent1or2-hmm${indiv}.msg$$.stdout 2> $logdir/parent1or2-hmm${indiv}.msg$$.stderr");

   } close BARCODE;
}
