#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use File::Basename;
use Switch;
$" = "\n";
my $msg_src = dirname $0 ;
my $verbose = 0;

# updated Greg Pinero April 2012
# - added Tn5-IonTorrent

# version 10 (02.25.10)
# - added more restriction enzymes
# - changed linker processing

# version 9.1 - reads are filtered according to a specified minimum quality score, runs of Ns at the end of a sequence are truncated, 
# also see USAGE for changes to input
# all reads with barcodes that don't exist are dumped to "bad_barcodes"
# all reads with unreadable barcodes are dumped to "unreadable_barcodes"
# all sequences with no insert are dumped to "linker"
# all sequences <20 bp are dumped to "junk"

# this works for an arbitrary number of barcodes of arbitrary length
# you need a separate file barcode file that lists each barcode on a separate line, a 2nd column denotes an additional tag to be added to each parsed sequence
# the program fq_all2std.pl is needed if you need to convert Illumina QVs to Sanger QVs - this should not be necessary much longer

# USAGE
my %opts = (e=>'null', p=>0, b=>'barcode_file', i=>'null', l=>'Dros_SR_vII', s=>0, m=>0, a=>0, c=>0);
getopts('e:p:b:i:l:s:m:a:c:', \%opts);
die(qq/
USAGE: parse_BCdata2BWA.pl [options] fastq_file

OPTIONS:
   -e restriction enzyme used for digestion [MseI|NdeI|Hpy188I|null]
   -p prepend bases to reads to complete the RE motif (useful for mapping reads) [0|1]
   -b barcode file [$opts{b}]
   -i file of barcodes to ignore [null]
   -l linker system used [Dros_SR_vI|Dros_SR_vII|Dros_SR_vIII_T_overhang|Dros_SR_vIII_blunt|ButtFinch_PE_vI|ButtFinch_PE_vII]
   -s strict (only retain reads with the correct RE motif) [0|1]
   -m minimum QV to mask from the 3' end [0]
   -a assign barcode if first n-1 bps are nonredundant [0]
   -c convert individual fastq files into fasta [0]
\n/) if (@ARGV == 0 && -t STDIN);

my ($infile, $prepend_enz, $prepend_flag, $barcode_file, $ignored_barcode_file, $linker_type, $strict, $minQV, $assignBC, $convert2fasta) =
   ($ARGV[$#ARGV], $opts{e} ,$opts{p} ,$opts{b} ,$opts{i} ,$opts{l} ,$opts{s} ,$opts{m} ,$opts{a} ,$opts{c});

print "Barcode file: $barcode_file\n";
print "RE: $prepend_enz\n";
print "Linker system: $linker_type\n\n";

my $opt_q = $minQV + 33 if ($minQV > 0 && $minQV <= 93);


### CLEANUP
####################################################################################################
my $dir = "${infile}_parsed"; mkdir $dir unless (-d $dir);
print "Cleaning up $dir...\n";
foreach my $file (qq(bad_barcodes bad_barcodes_addon linkers unreadable_barcodes junk RE_seq_mismatch RE_seq)) {
	unlink("$dir/$file") if (-e "$dir/$file");
}


### READ IN BARCODES
####################################################################################################
print "Reading in barcode file...\n";
my $BC_length=0;
my %BC_table = (); 
my %BC_table_chopLast = ();
open (BC,$barcode_file) || die "ERROR: Can't open $barcode_file: $!\n";
while (my $line = <BC>) { chomp $line;
print "$line\n";
	my ($key,$id)=split(/\s+/,$line);
	my $chopped_key = substr($key,0,length($key)-1);
	$BC_table{$key}=$id;
	$BC_table_chopLast{$chopped_key} = $key;
	$BC_length = length($key);

	unlink("$dir/indiv${id}_${key}") if (-e "$dir/indiv${id}_${key}");
	unlink("$dir/indiv${id}_${key}.$minQV.fa") if (-e "$dir/indiv${id}_${key}.$minQV.fa");
} close BC;
print "barcode length: $BC_length\nnumber of barcodes: " . scalar(keys(%BC_table)) . "\n";

if ($ignored_barcode_file ne 'null') {
	open (BC,$ignored_barcode_file) || die "ERROR: Can't open $ignored_barcode_file: $!\n";
	while (my $line = <BC>) { chomp $line;
		my ($key,$id)=split(/\s+/,$line);
		my $chopped_key = substr($key,0,length($key)-1);
		$key =~ s/N/[ATGCN]/g;
		$chopped_key =~ s/N/[ATGCN]/g;

		$BC_table{$key} = 'ignore';
		$BC_table_chopLast{$chopped_key} = $key;
	} close BC;
}

if ($assignBC == 1) {
	if (scalar(keys %BC_table)==scalar(keys %BC_table_chopLast)) {
		print "\nForce assigning barcodes based on chopping last position of barcode\n";
		foreach my $key (keys %BC_table_chopLast) {
			print "$key = $BC_table_chopLast{$key} = $BC_table{$BC_table_chopLast{$key}}\n"; 
		}
	} else { $assignBC = 0; }
}

### PRESETS
####################################################################################################
my $min_seq_length = 20;
my $seq_start = $BC_length;


### LINKER PROCESSING === barcode + overhang + reverseComplement_barcode + FC1
####################################################################################################
# Dros_Nico          NNT       Ann       TCGGAAGAGCTCGTATGCCGTCTTCTGCTTG
# Dros_SR_v1			NNNN   TA nnnn      TCGGAAGAGCTCGTATGCCGTCTTCTGCTTG	   TCGGAAGAG CTCGTATGCCGTCTTCTGCTTG 
# Dros_SR_vII			NNNNNG TA cnnnnn    TCGGAAGAGCTCGTATGCCGTCTTCTGCTTG	   TCGGAAGAG CTCGTATGCCGTCTTCTGCTTG 
# Dros_SR_vIII			NNNNNN    nnnnnn  GATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG	GA TCGGAAGAG CTCGTATGCCGTCTTCTGCTTG (blunt end ligation)
#							NNNNNN T  nnnnnn  GATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG	GA TCGGAAGAG CTCGTATGCCGTCTTCTGCTTG (if filled by Klenow) 
# ButtFinch_PE_vI		NNNG   TA cnnn    GATCGGAAGAGCGGTTCAGCAGGAATGCCGA		GA TCGGAAGAG CGGTTCAGCAGGAATGCCGA 
# ButtFinch_PE_vII	NNNi      nnn     GATCGGAAGAGCGGTTCAGCAGGAATGCCGAG		GA TCGGAAGAG CGGTTCAGCAGGAATGCCGAG
# Tn5-IonTorrent    NNNNNN AGATGTGTATAAGAGACAG some seq CTGTCTCTTATACACATCTATCACCGACTGCCCATAGAGAGGAAAGCGGAGGCGTAGTGG
#   barcode 19bp-transposon-end the-seqn-we-want  reverse-complement-of-linker-P1
#   i.e barcode barcode_addon  NO OVERHANG! FC1_linker_seq
#   NNNNNN AGATGTGTATAAGAGACAG some seq CTGTCTCTTATACACATCTATCACCGACTGCCCATAGAGAGGAAAGCGGAGGCGTAGTGG

### FC1_overhang === overhang selected for RE cut site
### overhang === overhang that gets sequenced (dependent on blunt ligation)
### FC1_linker_seq === FC1
my ($barcode_addon,$FC1_overhang,$overhang,$FC1_linker_seq);
switch ($linker_type) {
	case 'Dros_Nico'			{			$barcode_addon = '';  $FC1_overhang = 'TA'; $overhang = 'TA'; $FC1_linker_seq =   'TCGGAAGAGCTCGTATGCCGTCTTCTGCTTG'; }
	case 'Dros_SR_vI'			{			$barcode_addon = '';  $FC1_overhang = 'TA'; $overhang = 'TA'; $FC1_linker_seq =   'TCGGAAGAGCTCGTATGCCGTCTTCTGCTTG'; }
	case 'Dros_SR_vII'		{			$barcode_addon = ''; $FC1_overhang = 'TA'; $overhang = 'TA'; $FC1_linker_seq =   'TCGGAAGAGCTCGTATGCCGTCTTCTGCTTG'; }
	#case 'Dros_SR_vII'		{			$barcode_addon = 'G'; $FC1_overhang = 'TA'; $overhang = 'TA'; $FC1_linker_seq =   'TCGGAAGAGCTCGTATGCCGTCTTCTGCTTG'; }
	case 'Dros_SR_vIII_T_overhang'{	$barcode_addon = '';  $FC1_overhang = 'T';  $overhang = 'T';  $FC1_linker_seq = 'GATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG'; }
	case 'Dros_SR_vIII_blunt'		{	$barcode_addon = '';  $FC1_overhang = 'T';  $overhang = '';	  $FC1_linker_seq = 'GATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG'; }
	case 'ButtFinch_PE_vI'	{			$barcode_addon = 'G'; $FC1_overhang = 'TA'; $overhang = 'TA'; $FC1_linker_seq = 'GATCGGAAGAGCGGTTCAGCAGGAATGCCGA'; }
	case 'ButtFinch_PE_vII'	{			$barcode_addon = '';  $FC1_overhang = 'T';  $overhang = '';	  $FC1_linker_seq = 'GATCGGAAGAGCGGTTCAGCAGGAATGCCGAG'; } 
	case 'ButtFinch_PE_vII_ANA'	{	$barcode_addon = '';  $FC1_overhang = 'T';  $overhang = '';	  $FC1_linker_seq =   'TCGGAAGAGCGGTTCAGCAGGAATGCCGAG'; } # GA became part of barcode
	case 'Tn5-IonTorrent'	{		$barcode_addon = 'AGATGTGTATAAGAGACAG';  $FC1_overhang = '';  $overhang = ''; 
	    $FC1_linker_seq = 'CTGTCTCTTATACACATCTATCACCGACTGCCCATAGAGAGGAAAGCGGAGGCGTAGTGG'; }
	else							{			die "ERROR: Unknown linker_type: $linker_type\n"; }
}  

$seq_start = $BC_length + length($barcode_addon) if ($barcode_addon);
my $linker_search = 'TCGGAAGAG';
my $linker_length = length($FC1_linker_seq);
my $linker_strip_pos = 0; $linker_strip_pos = pos($FC1_linker_seq) if ($FC1_linker_seq =~ /$linker_search/g);
my $linker_search_pos = $BC_length + length($barcode_addon) + length($FC1_overhang) + $BC_length + $linker_strip_pos;

print "FC1 overhang for system $linker_type: ->$FC1_overhang<- (barcode addon ->$barcode_addon<-, sequenced overhang ->$overhang<-)\n";
print "self-ligated linker for $linker_type: $FC1_linker_seq ($linker_strip_pos, startpos $linker_search_pos)\n\n";


### RESTRICTION ENZYME STRICT SEARCH and PREPENDING
####################################################################################################
### padding between barcode and sequence read
my $RE_seq = ''; # RHS of cut
my $prepend_seq = ''; # LHS of cut
my $RE_overhang = 0;

# T added to 5'end to complete TTAA (MseI)
if ($prepend_enz eq 'MseI')	{			$prepend_seq = 'T';					$RE_seq = 'TAA';		$RE_overhang = -2;
} elsif ($prepend_enz eq 'NdeI')	{		$prepend_seq = 'CA';           	$RE_seq = 'TATG';		$RE_overhang = -2;
} elsif ($prepend_enz eq 'Hpy188I'){	$prepend_seq = 'TCn';          	$RE_seq = 'GA';		$RE_overhang = 1;
} elsif ($prepend_enz eq 'HphI'){		$prepend_seq = 'GGTGAnnnnnnnn';	$RE_seq = '';			$RE_overhang = 1;
} elsif ($prepend_enz eq 'MboII'){		$prepend_seq = 'GAAGAnnnnnnnn';	$RE_seq = '';			$RE_overhang = 1;
} elsif ($prepend_enz eq 'XcmI'){		$prepend_seq = 'CCAnnnnn';     	$RE_seq = 'nnnnTGG';	$RE_overhang = 1;
} elsif ($prepend_enz eq 'Hpy188III'){ $prepend_seq = 'TC';           	$RE_seq = 'nnGA';		$RE_overhang = -2;
} elsif ($prepend_enz eq 'AhdI') {		$prepend_seq = 'GACnnn';       	$RE_seq = 'nnGTC';	$RE_overhang = 1;
} elsif ($prepend_enz eq 'HpyAV') {		$prepend_seq = 'CCTTCnnnnnn';  	$RE_seq = '';			$RE_overhang = 1;
} elsif ($prepend_enz eq 'null') {		$prepend_seq = '';  	$RE_seq = '';			$RE_overhang = 0;
} else { die "ERROR: Unknown restriction enzyme ($prepend_enz)!\n"; }

### update RE_seq OR prepend_seq with the overhang
print "FC1_overhang $FC1_overhang\n";
print "BEFORE prepend_seq $prepend_seq.$RE_seq RE_seq\n";

### adjust RE_seq with FC1_overhang
if ($RE_overhang < 0) {
	$RE_seq = $FC1_overhang . substr($RE_seq,length($FC1_overhang)) if ($RE_seq);

### adjust prepend_seq with FC1_overhang
} elsif ($RE_overhang > 0) {
	$prepend_seq = substr($prepend_seq,0,length($prepend_seq)-length($FC1_overhang)) . $FC1_overhang if ($prepend_seq);
}
print "AFTER  prepend_seq $prepend_seq.$RE_seq RE_seq\n";

my ($append_bc_seq,$append_bc_seq_length,$check_overhang) = ('',length($RE_seq),0);
my $search_internal_RE = $prepend_seq . $RE_seq;
my $prepend_qual = 'I' x length($prepend_seq);
$RE_seq =~ s/n/[ATGCN]/g;
$search_internal_RE =~ s/n/[ATGCN]/g;
if ($RE_overhang < 0) {
	$append_bc_seq = $RE_seq;
	#$append_bc_seq = $barcode_addon . $RE_seq;
} elsif ($RE_overhang > 0) {
	$append_bc_seq = $FC1_overhang . $RE_seq;
	$append_bc_seq_length += length($FC1_overhang);
	$check_overhang = 1;
	#$append_bc_seq = $barcode_addon . $FC1_overhang . $RE_seq;
}


if ($prepend_flag == 0) { $prepend_seq = ''; $prepend_qual = ''; }
print "prepend seq for $prepend_enz (flag=$prepend_flag) = '$prepend_seq'\n";
print "internal RE search string = ->$prepend_seq<- . ->$RE_seq<- = $search_internal_RE\n";
print "RE_seq ($prepend_enz) ->$prepend_seq.$RE_seq<-\n";
print "append_bc_seq_length = $append_bc_seq_length ($append_bc_seq)\n";
print "FORMAT: barcode + barcode_addon($barcode_addon) + FC1_overhang($FC1_overhang) + RE_seq($RE_seq) --> barcode + $append_bc_seq + sequence\n";


### TRACK READS
####################################################################################################
my %counters_bad = qw(	total 0
								bad_codes 0
								unreadable_codes 0 
								forced_assigned 0 );
my %counters_bc;  
foreach my $key (keys %BC_table) {
	next if ($BC_table{$key} eq 'ignore');
	foreach my $type qw(
		linkers RE_seq_mismatch min_seq_length overhang_mismatch bad_codes_addon
		RE_seq_mismatch.min_seq_length small_insert internal_RE good) {
		$counters_bc{"$key.$BC_table{$key}"}{$type} = 0;
	}
}

my @linkers =();
my @bad_codes = ();
my @bad_codes_addon = ();
my @overhang_mismatch = ();
my @unreadable_codes = ();
my @junk = ();
my @RE_seq_mismatch = ();
my @RE_seqs = ();
my %BC_arrays = (); foreach my $key (keys %BC_table) { @{$BC_arrays{$BC_table{$key}}} = (); } # initialize

### Track quals of bad and good barcodes
my (@bad_barcode_quals,@good_barcode_quals,@mismatch_barcode_linkers,@unreadable_barcode_quals);

my $i=0;
my $k=0;
open (IN,$infile);
while (<IN>) { my $line = $_; chomp $line;
	
	####################################################################################################
	### WRITE OUT...
	if ($i==500000) {
		&printout($dir,\%BC_arrays,\%BC_table,\@bad_codes,\@bad_codes_addon,\@overhang_mismatch,\@unreadable_codes,\@linkers,\@junk,\@RE_seq_mismatch);

		foreach my $key (keys %BC_table) { @{$BC_arrays{$BC_table{$key}}} = (); }
		@bad_codes = ();
		@bad_codes_addon = ();
		@overhang_mismatch = ();
		@unreadable_codes = ();
		@linkers = ();
		@junk = ();
		@RE_seq_mismatch = ();

		$i=0;
		$k = $k+500000;
		print "completed: $k\n";
	}
	
	$i++;
	$counters_bad{'total'}++;

	# @GEN-SEQ-ANA:3:1:0:341#0/1
	# NGTGTATGGATAAAAAATCAAAGATTATCACATTTAAGGTTTGCTGATGACATAGCACTATTCAGTGATACGGCCC
	# +GEN-SEQ-ANA:3:1:0:341#0/1
	# DNY[[[[[[[[[[[[[[XVX[ZZZZZZXXXZY[ZZZ[YYXYYYVY[[XXYZXZYVSVYYZWVQXVSVYXUXSVSBB
	$line =~ s/-/_/g; # sequence header
	my $name = $line;
		
	##### BARCODE, ID and SEQUENCE
	$line = <IN>; chomp $line; my $orig_sequence = $line;
	my $bc       = substr($orig_sequence,0,$BC_length);
	my $identity = &matchBC($bc,%BC_table);
	my $barcode_addon_check = 0; $barcode_addon_check = 1 if (substr($orig_sequence,$BC_length,length($barcode_addon)) ne ($barcode_addon));
	my $overhang_check = 0; $overhang_check = 1 if (($check_overhang==1) &&
		(substr($orig_sequence,$BC_length+length($barcode_addon),length($FC1_overhang)) ne ($FC1_overhang)));

	### try to classify with the last position of the barcode chopped off
	if (($identity eq '0') && ($assignBC == 1)) {
		my $chopped_bc = substr($bc,0,length($bc)-1);
		if (exists $BC_table_chopLast{$chopped_bc}) {
			$bc = $BC_table_chopLast{$chopped_bc};
			$identity = &matchBC($bc,%BC_table);

			$counters_bad{'forced_assigned'}++;
		}
	}

	my $full_sequence = substr($line,$seq_start);
	my $unmasked_sequence = substr($line,$seq_start);
	$line = <IN>; $line = <IN>; chomp $line;
	my $bc_qual = substr($line,0,$seq_start);
	my $quality = substr($line,$seq_start);

	next if ($identity eq 'ignore');

	### MASK BASED ON QV
	print "\n*** $name ($bc $identity) ***\n" if ($verbose == 1);
	print "ORIGINAL\t$orig_sequence\nUNMASKED\t$unmasked_sequence\n" if ($verbose == 1);
	my $sequence = '';
	if ($opt_q) {
		for (my $q=0; $q<length($quality); $q++) {
			if (ord(substr($quality,$q,1)) < $opt_q ) { $sequence .= 'N'; } 
			else { $sequence .= substr($unmasked_sequence,$q,1); }
		}
	} else { $sequence = $unmasked_sequence; }
	print "MASKED\t\t$sequence\n" if ($verbose == 1);

	### STRIP OFF N's from the end
	if ($sequence =~ /(N+)$/g) {
		my $string_of_Ns = pos($sequence) - length($1);
		$sequence = substr($sequence,0,$string_of_Ns);
		$quality  = substr($quality,0,$string_of_Ns);
	}
	print "STRIPPED\t$sequence\n" if ($verbose == 1);

	### (1) CHECK BARCODE
	if ($identity eq "0") {

		# UNREADABLE barcode --- N in barcode section
		if ($bc =~/N/g){
			$counters_bad{'unreadable_codes'}++;

			push @unreadable_barcode_quals,$bc_qual;
			$name =~ s/GEN_SEQ_ANA/UNREADABLE_${bc}/g;
			push(@unreadable_codes, $name);
			push(@unreadable_codes, $orig_sequence);

		# BAD barcode --- not classified and doesn't have N's
		} else {
			$counters_bad{'bad_codes'}++;

			push @bad_barcode_quals,$bc_qual;
			$name =~ s/GEN_SEQ_ANA/BAD_${bc}/g;
			push(@bad_codes, $name);
			push(@bad_codes, $orig_sequence);
		}

	} elsif ($barcode_addon_check==1) {
		$counters_bc{"$bc.$identity"}{'bad_codes_addon'}++;

		$name =~ s/GEN_SEQ_ANA/BAD_ADDON_${bc}/g;
		push(@bad_codes_addon, $name);
		push(@bad_codes_addon, $orig_sequence);

	} elsif ($overhang_check==1) {
		$counters_bc{"$bc.$identity"}{'overhang_mismatch'}++;

		$name =~ s/GEN_SEQ_ANA/OVERHANG_MM_${bc}/g;
		push(@overhang_mismatch, $name);
		push(@overhang_mismatch, $orig_sequence);

	### NOW ALL OTHER CHECKS
	} else {
	
		### CASE 1a: search self-ligating linker products
		### CASE 1b: check for sequences too short (has linker) and take everything to the left
		my $self_linker_search_bc = '\w{0,' . length($FC1_overhang) . '}' . &revComp($bc . $barcode_addon);
		my $self_linker_search_pos = '\w{' .(length($bc)+length($barcode_addon)) . ',' . 
														(length($bc)+length($barcode_addon)+length($FC1_overhang)) . '}' . $FC1_linker_seq;
		print "SELF LIGATED SEARCH $self_linker_search_bc OR $self_linker_search_pos\n" if ($verbose == 1);

		if (($sequence =~ /^$self_linker_search_bc/) || ($sequence =~ /^$self_linker_search_pos/)) {
			$counters_bc{"$bc.$identity"}{'linkers'}++;

			push @linkers,($name,$orig_sequence); 
			print "...trashing (1a) self-ligated linker $sequence\n" if ($verbose == 1);

		} elsif ($sequence =~ /$linker_search/) {
			$counters_bc{"$bc.$identity"}{'small_insert'}++;

			$name =~ s/GEN_SEQ_ANA/JUNK_${bc}_Clone_too_short/g;
			push(@junk, $name);
			push(@junk, $orig_sequence);			
			print "...trashing (1b) short clone $sequence\n" if ($verbose == 1);

		} else {
			push @RE_seqs,substr($full_sequence,0,$append_bc_seq_length); ### KEEP TRACK OF RE_SEQ SITE

			# Establish if the reads make sense or not
			my ($junk_switch,$junk_name) = (0,'');

			if (($strict==1) && (substr($full_sequence,0,$append_bc_seq_length) !~ /^$append_bc_seq/)) { 
				$junk_switch = 1; 
				$junk_name = 'RE_seq_mismatch'; 
				print "$append_bc_seq != ",substr($full_sequence,0,$append_bc_seq_length),"\n" if ($verbose == 1);
			}
			if (length($sequence)<$min_seq_length) { 
				$junk_switch = 1; 
				if ($junk_name) { $junk_name .= '.min_seq_length'; }
				else { $junk_name = 'min_seq_length'; }
			}
	
			### CASE 3: keep the read
			if ($junk_switch==0){

				if (($strict==1) && ($sequence =~ /$search_internal_RE/) && ($search_internal_RE)) {
					print "...found internal RE $sequence\n" if ($verbose == 1);
					$counters_bc{"$bc.$identity"}{'internal_RE'}++;

					$junk_name = 'internal_RE_site'; 
					$name =~ s/GEN_SEQ_ANA/JUNK_${bc}_${junk_name}/g;
					push(@junk, $name);
					push(@junk, $orig_sequence);			

				} else {
					push @good_barcode_quals,$bc_qual;
					$counters_bc{"$bc.$identity"}{'good'}++;

					$name =~ s/GEN_SEQ_ANA/indiv_${identity}_${bc}/g;
					push @{$BC_arrays{$BC_table{$bc}}}, $name;
					push @{$BC_arrays{$BC_table{$bc}}}, $prepend_seq . $sequence;
						
               #$name =~ s/^@/+/;
					push @{$BC_arrays{$BC_table{$bc}}}, '+';
					push @{$BC_arrays{$BC_table{$bc}}}, $prepend_qual . $quality;
					print "...keeping (3) $prepend_seq . $sequence\n" if ($verbose == 1);
die "ERROR: how did this short read get past strict??? $name\n$sequence\n" if (length($sequence) < $min_seq_length);
				}

			# CASE 4: sequence starts out wrong is is therefore junk
			} else {
				$counters_bc{"$bc.$identity"}{$junk_name}++;
				if ($junk_name eq 'RE_seq_mismatch') {
					$name =~ s/GEN_SEQ_ANA/UNKNOWN_${bc}_${junk_name}/g;
					push(@RE_seq_mismatch, $name);
					push(@RE_seq_mismatch, $orig_sequence);			
					print "...trashing (4 $junk_name) $sequence\n" if ($verbose == 1);

				} else {
					$name =~ s/GEN_SEQ_ANA/JUNK_${bc}_${junk_name}/g;
					push(@junk, $name);
					push(@junk, $orig_sequence);
					print "...trashing (4 $junk_name) $sequence\n" if ($verbose == 1);
				}
			}
		}
	}
} close IN;

&printout($dir,\%BC_arrays,\%BC_table,\@bad_codes,\@bad_codes_addon,\@overhang_mismatch,\@unreadable_codes,\@linkers,\@junk,\@RE_seq_mismatch);
print "completed: ", $k+$i, "\n";



### convert to fasta and mask N's
####################################################################################################
if ($convert2fasta) {
	open (BC,$barcode_file) || die "ERROR: Can't open $barcode_file: $!\n";
	while (my $line = <BC>) { chomp $line;
		my ($key,$id)=split(/\s+/,$line);
		system("$msg_src/fastq_2_fasta.pl -i $dir/indiv${id}_${key} -o $dir/indiv${id}_${key}.$minQV.fa -q $minQV") if (-e "$dir/indiv${id}_${key}"); # server version
	} close BC;
}


### CREATE SEQLOGO OF THE RE_SEQ SITE
####################################################################################################
my @RE_seq_reps = &getFrequencies(\@RE_seqs);
open (RE,">$dir/RE_seq");
for (my $i=0; $i<=$#RE_seq_reps; $i++) { print RE ">$i\n$RE_seq_reps[$i]\n"; } close RE;
print "Creating RE_seq logo from ",scalar(@RE_seqs)," reads\n";

# [2:03pm-Mar.11.10] tina@yakuba--> ~/myBin/weblogo/seqlogo -F PDF -f RE_seq_mismatch_seq > Nde_vI_REseq_mismatch.pdf
system("$msg_src/seqlogo -F PDF -f $dir/RE_seq -x $prepend_enz -c -Y -n -e -S > $dir/RE_seq.$prepend_enz.pdf &") if (-e "$msg_src/seqlogo");



### SUMMARY
####################################################################################################
print "\n\nTOTAL BREAKDOWN:\n";
foreach my $type (sort keys %counters_bad) {
	print "$type,$counters_bad{$type}\n";
}

print "\nBY BARCODE\n";
foreach my $bc (sort keys %counters_bc) {
	foreach my $type (sort keys %{$counters_bc{$bc}}) {
		print "$bc,$type,$counters_bc{$bc}{$type}\n";
	}
}

exit;



####################################################################################################
####################################################################################################
####################################################################################################
sub printout {
	my ($dir,$BC_arrays,$BC_table,$bad_codes,$bad_codes_addon,$overhang_mismatch,$unreadable_codes,$linkers,$junk,$RE_seq_mismatch) = @_;

	foreach my $key (keys %{$BC_table}) {
		my $id = $$BC_table{$key};
		#Warning: parse_and_map.py expects these files to start with "indiv" so careful changing that
		open (OUT, ">>$dir/indiv${id}_${key}"); 
		print OUT "@{$BC_arrays{$BC_table{$key}}}\n" if (@{$BC_arrays{$BC_table{$key}}});
		@{$BC_arrays{$BC_table{$key}}} = ();
	}
			
	open (OUT, ">>$dir/bad_barcodes"); print OUT "@{$bad_codes}\n" if (@{$bad_codes});
	open (OUT, ">>$dir/bad_barcodes_addon"); print OUT "@{$bad_codes_addon}\n" if (@{$bad_codes_addon});
	open (OUT, ">>$dir/overhang_mismatch"); print OUT "@{$overhang_mismatch}\n" if (@{$overhang_mismatch});
	open (OUT, ">>$dir/unreadable_barcodes"); print OUT "@{$unreadable_codes}\n" if (@{$unreadable_codes});
	open (OUT, ">>$dir/linkers"); print OUT "@{$linkers}\n" if (@{$linkers});
	open (OUT, ">>$dir/junk"); print OUT "@{$junk}\n" if (@{$junk});
	open (OUT, ">>$dir/RE_seq_mismatch"); print OUT "@{$RE_seq_mismatch}\n" if (@{$RE_seq_mismatch});
	return;
}

sub revComp {
	my ($string) = @_;
	my $rev_string = reverse $string;
	$rev_string =~ tr/ATGC/TACG/;
	return $rev_string;
}


sub matchBC {
	my ($bc,%BC_table) = @_;

	my $identity;
	if (exists $BC_table{$bc}) { return $BC_table{$bc}; }
	else {
		foreach my $key (keys %BC_table) { return $BC_table{$key} if ($bc =~ /^$key$/); }
	}

	return '0';
}

sub getFrequencies {
	my ($seqs) = @_;

	### get counts of each base at each position
	my %pos_counts;
	foreach my $seq (@{$seqs}) {
		my @bases = split(//,$seq);
		for (my $i=0; $i<=$#bases; $i++) {
			$pos_counts{$i}{$bases[$i]}++;
		}
	}

	### create fake sequences representing the proportion of each base at each position
	my @fake_seqs;
	foreach my $pos (sort {$a<=>$b} keys %pos_counts) {
		my $bases = '';
		foreach my $allele (keys %{$pos_counts{$pos}}) {
			my $fraction = int(1000 * ($pos_counts{$pos}{$allele}/scalar(@{$seqs})));
			$bases .= $allele x $fraction;
		}
		$bases .= ('N' x (1000-length($bases))) if (length($bases) < 1000);

		for (my $i=0; $i<1000; $i++) { $fake_seqs[$i] .= substr($bases,$i,1); }
	}

	return @fake_seqs;
}
