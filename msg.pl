#!/usr/bin/env perl
use File::Basename;
use File::Copy;
use Getopt::Long ;
$true = 1 ;
$false = 0 ;
$version = '0.0.1' ;
sub system_call {
	print "@_\n" ;
	system("@_") == 0 or die "Error in @_: $?" ;
}

$src = dirname $0 ;
$update_genomes = $false ;
$update_nthreads = 1 ; ## Number of BWA threads when updating genomes (must match msgCluster.pl)

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
	'threads|t=i' => \$update_nthreads,
	'parse_or_map=s' => \$parse_or_map,
	'priors=s' => \$priors,
	'chroms=s' => \$chroms,
	'sexchroms=s' => \$sexchroms,
	'chroms2plot=s' => \$chroms2plot,
	'deltapar1=s' => \$deltapar1,
	'deltapar2=s' => \$deltapar2,
	'rfac=s' => \$rfac,
	'bwaindex1=s' => \$bwaindex1,
	'bwaindex2=s' => \$bwaindex2
	) ;

print "msg version:  $version\n" ;
print `date` ;
print "\n" ;
print "parent1 genome: $parent1_genome\n" ;
print "parent2 genome: $parent2_genome\n" ;
print "update:         $update_genomes\n" ;
print "parse_or_map:   $parse_or_map\n" ;
print "chroms:         $chroms\n\n" ;
print "sexchroms:      $sexchroms\n\n" ;
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

$genomes_fa{'parent1'} = $parent1_genome if $parent1_genome;
$genomes_fa{'parent2'} = $parent2_genome if $parent2_genome;
$genome_index{'parent1'} = $bwaindex1;
$genome_index{'parent2'} = $bwaindex2;

## -------------------------------------------------------------------------------
##
## UPDATE REFERENCES
if( $update_genomes ) {
	%reads_for_updating_fq;
	$reads_for_updating_fq{'parent1'} = $parent1_reads if ($parent1_reads);
	$reads_for_updating_fq{'parent2'} = $parent2_reads if ($parent2_reads);
	#%reads_for_updating_fq = (parent1 => $parent1_reads, parent2 => $parent2_reads) ;

	for $sp (keys %reads_for_updating_fq) {
		$reads = $reads_for_updating_fq{$sp} ;
		$reads or die "Must supply --$sp-reads" ;
		-e $reads or die "No such file: $reads" ;
		($reads_for_updating{$sp}, $d, $s) = fileparse $reads, qr/\.[^.]*$/ ;
		$s eq ".fq" or $s eq ".fastq" or print "$sp reads extension is $s: expecting .fq or .fastq\n" ;
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
	foreach $sp ( keys %reads_for_updating ) {
		if( -e $genomes_fa{$sp}.".msg.updated.fasta" ) {
			print "$sp: updated genome files already present\n" ;
			next ;
		}
	
		&system_call("perl", "$src/reformatFasta4sam.pl", "-i", $genomes_fa{$sp}, "-o", "$genomes_fa{$sp}.msg") ;

		$out = "update_reads-aligned-$sp" ;
		unless( -e "$out.pileup" ) {

			&system_call("bwa", "index", "-a", $genome_index{$sp}, "$genomes_fa{$sp}.msg") 
				unless( -e "$genomes_fa{$sp}.msg.bwt" and -e "$genomes_fa{$sp}.msg.ann" );
			&system_call("samtools", "faidx", "$genomes_fa{$sp}.msg") ;

			unless (-e "$out.sam") {
				&system_call("bwa", "aln", "-t", $update_nthreads, "$genomes_fa{$sp}.msg", $reads_for_updating_fq{$sp}, "> $out.sai") ;
				&system_call("bwa", "samse", "$genomes_fa{$sp}.msg", "$out.sai", $reads_for_updating_fq{$sp}, "> $out.sam") ;
			}

			&system_call("$src/filter-sam.py", "-i", "$out.sam", "-o", "$out.filtered.sam") ;

			&system_call("samtools", "view", "-bt", "$genomes{$sp}.msg.fai", "-o $out.bam", "$out.filtered.sam") ;
			&system_call("samtools", "sort", "$out.bam", "$out.bam.sorted") ;
			&system_call("samtools", "index", "$out.bam.sorted.bam") ;
			&system_call("samtools", "pileup", "-f", "$genomes_fa{$sp}.msg", "$out.bam.sorted.bam", "-c", "> $out.pileup") ;

		}

		&system_call("perl", "$src/updateRef.pl", "$genomes_fa{$sp}.msg", "$out.pileup", "0", $update_minQV, $min_coverage);
		unlink("$genomes_fa{$sp}.msg");
						
		## Re-run BWA against updated refs
		## Not doing that at the moment
	}
	print "Finished updating genomes\n" ;
	exit;
}



## -------------------------------------------------------------------------------
##
## Use BWA to align raw reads against genomes, creating two sam files for each individual
for $sp (keys %genomes_fa) {
	-e $genomes_fa{$sp} or die "No such file: $genomes_fa{$sp}" ;

	($genomes{$sp}, $d, $s) = fileparse $genomes_fa{$sp}, qr/\.[^.]*$/ ;
	$s eq ".fa" or $s eq ".fasta" or print "$sp genome extension is $s: expecting .fa or .fasta\n" ;
	
	if ($parse_or_map eq 'parse-only') {
		### reformat for samtools (60 chars per line)
		&system_call("perl", "$src/reformatFasta4sam.pl", "-i", "$genomes_fa{$sp}", "-o", "${sp}_ref.fa") unless (-e "${sp}_ref.fa") ;
	}
}


%genomes_fa = (parent1 => 'parent1_ref.fa', parent2 => 'parent2_ref.fa') ;
$parse_or_map = "--$parse_or_map" if ($parse_or_map);

### BWA INDEXING
foreach $sp ( keys %genomes ) {
	unless( -e "$genomes_fa{$sp}.bwt" and -e "$genomes_fa{$sp}.ann" ) {
		## 1. Use BWA to align reads against existing genome -> sam file
		@args = ("bwa", "index", "-a", $genome_index{$sp}, $genomes_fa{$sp}) ;
		print "@args\n" ;
		system("@args") == 0 or die "Error in @args: $?" ;
	}
}

$samfiles_dir = basename($raw_read_data) . "_sam_files" ;
mkdir $samfiles_dir unless (-d $samfiles_dir);

print "\n\nRUNNING ", join(' ',('python', "$src/parse_and_map.py", '-i', $raw_read_data, '-b', $barcodes,
		 '--parent1', $genomes_fa{'parent1'}, '--parent2', $genomes_fa{'parent2'}, $parse_or_map,
       '--re_cutter', $re_cutter, '--linker_system', $linker_system )) ;
&system_call('python', "$src/parse_and_map.py", '-i', $raw_read_data, '-b', $barcodes,
		 '--parent1', $genomes_fa{'parent1'}, '--parent2', $genomes_fa{'parent2'}, $parse_or_map,
       '--re_cutter', $re_cutter, '--linker_system', $linker_system ) ;

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
   
   	&system_call('bash', "$src/parent1or2-hmm.sh",
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
   			 '-r', $rfac,
   			 '-x', $sexchroms,
   			 '-z', $priors
   		  ) ;

   } close BARCODE;
}
