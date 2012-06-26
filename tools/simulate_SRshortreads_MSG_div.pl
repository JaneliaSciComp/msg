#!perl -w
use strict;

# version 4 - models sequence reads with "divergence" in the form of snps and indels -- NOT a good model for sequence error because they occur uniformly 

# version 3 - this version restricts the reads to cloned library fragments in a specified size range (low, high)

# version 2 speeds things up 10-fold by searching for MseI sites more efficiently

# usage 
# perl simulate_SRshortreads_MSG2.pl reffile maxreadlength numreads low_size_lim high_size_lim


# note that this version starts reads from MseI sites to mimic short reads from an MSG libary

# simulates single-end reads from a recombinant genome specified in a simulated ancestry file
# to simulate PE reads use simulate_PEshortreads.pl

# no error model is currently implemented, all bases are given maximum QV (=I)
# RUI varition in read length is implemented: 75-$maxreadlength 

# NOTE - the genome reference file should stochiometrically represent a recombinant diploid or haploid individual
# example: F2 individual, there should be two recombinant genomes
# example: F1 backcross, the file should contain one recombinant genome and one parental genome
# example: haploid progeny, the file should contain only one recombinant genome

if (@ARGV < 4){print "\nusage: zcat ref.fa.gz | perl simulate_SRshortreads_MSG_div.pl maxreadlength numreads low_size_limit high_size_limit \n\n";exit;}


#################  SET MUTATION PARAMETERS HERE ################################################################################
my $mut_rate = 0.03;							# rate of divergence per base pair (0.03 for Dsim-Dsec)
my $indel_mu_scalar = 0.157;					# scalar to multiply indel divergence by (~0.157 for Drosophila)
my $indel_size = 1.2;							# mean of exponential that describes indel size distribution (~1.2 for Drosophila)
#################  SET MUTATION PARAMETERS HERE ################################################################################


use Math::Random qw(:all);
srand(time());

my $maxreadlength = shift(@ARGV); chomp $maxreadlength;			# read length from Illumina machine 
my $numreads = shift(@ARGV); chomp $numreads;					# number of reads

my $low = shift(@ARGV); chomp $low;  # low size limit
my $high = shift(@ARGV); chomp $high;  # high size limit

open READS, ">reads.fa";
open SAM, ">reads.sam";
print SAM "ReadID\tDirection\tStart\tCIGAR\tMDtag\n";


print "mutation parameters: mu: $mut_rate indel_mu_scalar: $indel_mu_scalar indel_size: $indel_size\n";


#**************** read in genome reference sequence and create "genome" string
my @data = ();
my @sequence_names = ();
my @sequence_temp =();  
my @data_temp =();	# needed to read in nexus files
my $genome_length =0; my $seq_count =0; my $sequence_temporary = ""; 
my $genome = "";
$genome = $genome.("N" x $maxreadlength); # initial padding with Ns ($maxreadlength) 
print "Reading in the reference genome ...\n";
while (<>) { my $line=$_;
		chomp $line;
		$line =~ s/\s+//g;			# remove all white spaces
		$line =~ s/type/\ttype/g;	# insert a tab
		if ($line =~ />/) 
			{
			$seq_count++;
			(my $name, my $other) = split(/\t/, $line); chomp $name;
			push(@sequence_names, $name);
			$sequence_temporary = join ("", @data_temp);
			push(@sequence_temp,$sequence_temporary); 
			@data_temp=();
			} 
		else {push(@data_temp, $line);}
		}
	$sequence_temporary = join ("", @data_temp);
	push(@sequence_temp,$sequence_temporary); 
	for my $x (0..(scalar(@sequence_names)-1))	# assemble "genome" string
		{
		$genome = $genome.$sequence_temp[$x+1];
		$genome = $genome.("N" x $maxreadlength); # padding with Ns ($maxreadlength) to separate chromosomes
		}
print "Done. Number of sequences read in: $seq_count\n";
#print OUT $genome;

my $length = length($genome); print "genome length: $length\n";


# catalog all MseI sites in the genome
my @MseIsites =();
while ($genome =~ /TTAA/g){push(@MseIsites, pos($genome)-4);}

print "number of MseI sites: ", scalar(@MseIsites), "\n";
# for my $k (0..(scalar(@MseIsites)-1)){print $MseIsites[$k], "\t";} print "\n";

# ****************** generate the short read data 
my $counter = 0;
while ($counter < $numreads){
	# for my $i (1..$numreads){
	my $pos =0; my $site = 0;
	my $diff_right = 0; my $diff_left = 0;
	my $pos_right =0; my $pos_left=0;
	
	my $read = ""; my $qv ="";  
	# $pos = random_uniform_integer(1, $maxreadlength, $length-$maxreadlength);	# start of forward read
		
	while (($pos==0) or ($pos<$maxreadlength) or ($pos>($length-$maxreadlength))){
		$site = random_uniform_integer(1, 2, (scalar(@MseIsites)-2));
		$pos = $MseIsites[$site];
		# print "search site: $site search pos: ", $pos-$maxreadlength+3, "\n";
		}		
	    # print "site: $site pos: ", $pos-$maxreadlength+3, "\n";
	
		my $polarity = random_uniform(1);												# direction of read and whether it is mated
		my $readlength = random_uniform_integer(1, 75, $maxreadlength);					# read length
		$pos_right = $MseIsites[$site+1]; 
		$pos_left  = $MseIsites[$site-1];
		$diff_right = $pos_right + (-1*$pos);
		$diff_left  = $pos + (-1*$pos_left);
	
	my $read_start=0;
	my $read_polarity=4;
	if (($polarity<=0.5) && ($diff_right>=$low) && ($diff_right<=$high)){			# forward read
			# print "read $i polarity (+) position: ", $pos-$maxreadlength+3, "\n";
			# print "forward read: ", $pos-$maxreadlength+3, " - ", $pos+$readlength-$maxreadlength+3, "\n";
			$read_polarity = 0;
			$read_start= $pos+6-100;
			$read = substr($genome, $pos+1, $readlength);
			$counter++;
			} # forward polarity
	elsif (($diff_left>=$low) && ($diff_left<=$high)) {								# reverse read
			# print "read $i polarity (-) position: ", $pos-$maxreadlength+3, "\n";
			# print "forward read: ", $pos-$readlength-$maxreadlength+3, " - ", $pos-$maxreadlength+3, "\n";
			$read_polarity = 16;
			$read_start= $pos-$readlength+9-101;
			$read = rev_comp(substr($genome, $pos-$readlength+3, $readlength));
			$counter++;
			} # reverse polarity
		
		$read =~ s/N//g;
		
		if ($read ne ""){
			my @mutated_read = mutate($read, $mut_rate, $indel_mu_scalar, $indel_size);		# mutate read
			# my @mutated_read = mutate($read, 0.03, 0.157, 1.2);							# mutate read
			# my @mutated_read = mutate($read, 0, 0, 1);									# DONT mutate read
			for my $k (0..(length($mutated_read[0])-1)){$qv = $qv."I";}						# read qualities
			print READS "\>read$counter\n$mutated_read[0]\n";
			print SAM "read$counter\t$read_polarity\t$read_start\t$mutated_read[1]\t$mutated_read[2]\n";
			# print "$read\t", length($read), "\n";
			}
		else{}	
		
	} # for numreads


#**************************************************************************************************************
# make reverse complement

sub rev_comp {
    
	my($sequence) = @_;
	
	$sequence =~ s/G/1/go; 
	$sequence =~ s/A/2/go; 
	$sequence =~ s/T/3/go; 
	$sequence =~ s/C/4/go; 
	$sequence =~ s/1/C/go; 
	$sequence =~ s/2/T/go; 
	$sequence =~ s/3/A/go; 
	$sequence =~ s/4/G/go; 
	
	$sequence = reverse($sequence);
	
	return $sequence;
}

# mutate read
# ****************************************************************************************************************

sub mutate {
	
	my ($sequence, $mutation_rate, $indel_mut_rate_scalar, $indel_size) = @_;
	# my $mutation_rate = 0.03;
	#my @array = ("G", "A", "T", "C");
	my @array = ("g", "a", "t", "c");
	my $mut_seq = $sequence;
	my $mu = $mutation_rate * length($sequence);
	my $num_snps=0; my $num_indels=0;
	if ($mu != 0){
		$num_snps = random_poisson(1, $mu);										# number of snp mutations
		$num_indels = random_poisson(1, $mu*$indel_mut_rate_scalar);			# number of indel mutations
		# if ($num_indels >1){$num_indels=1;}									# constrain to 1 indel
		}
	my @indel_positions = ();
	my @indel_length = ();
	
	if (($num_snps==0) && ($num_indels==0)){return ($sequence, length($sequence)."M", "MD:Z:".length($sequence));} # if no mutations at all
	else{
		# process snps here
		if ($num_snps !=0){
			my @snp_positions = random_uniform_integer($num_snps, 1, length($sequence)); # positions of snps
			for my $j (0..(scalar(@snp_positions)-1)){
				my $snp=random_permutation(@array);
				while (uc($snp) eq substr($sequence,$snp_positions[$j]-1,1)){$snp=random_permutation(@array);}
				#print "snp: $snp pos: $snp_positions[$j] mut: ", substr($sequence,$snp_positions[$j]-1,1), "\n";
				$mut_seq=substr($mut_seq,0,$snp_positions[$j]-1).$snp.substr($mut_seq,$snp_positions[$j]);
			} # for all snps
		} # if there is > 0 snp muts		

		# process indels here
		if ($num_indels !=0){
			my $check_indels=0;
			@indel_positions = ();
			@indel_length = ();
			while ($check_indels==0){
				my $flag = 0;
				@indel_positions = random_uniform_integer($num_indels, 1, length($sequence));  # positions of indels
				# indel lengths
				for my $k (0..(scalar(@indel_positions)-1)){
					my @polarity_array = (-1,1);
					my $indel_len=0; 
					while($indel_len==0){$indel_len= int(random_exponential(1, $indel_size)+0.5)*random_permutation(@polarity_array)};
					# print "indel: $indel_len bp pos: $indel_positions[$k]\n";
					$indel_length[$k]=$indel_len;
				}
				
				if ($num_indels<2){$check_indels=1; }# print "indel overlap is not an issue\n";}
				else{
					my %indel_sites =();
					for my $q (0..(scalar(@indel_positions)-1)){
						if ($indel_length[$q]>0){
							for my $r (0..$indel_length[$q]-1){
								my $position = $indel_positions[$q]+$r; 
								$indel_sites{$position}++;
								if ($indel_sites{$position}>1){$flag=1;}
							} # for length of indel
						} # if insertion
						else{   for my $r ($indel_length[$q]+1..0){
							my $position = $indel_positions[$q]+$r; 
							$indel_sites{$position}++;
							if ($indel_sites{$position}>1){$flag=1;}
						} # for length of indel
						} # if deletion
					} # for all indels
					if ($flag==0){$check_indels=1; }# print "multiple indels but non-overlapping\n"; }
					else{}# print "multiple indels but overlapping, redrawing ... \n";}
				} # if more than one indel
			} # while check_indels==0 
			
			# print "number of indels to process: ", scalar(@indel_positions), "\n";
			for my $t (0..(scalar(@indel_positions)-1)){
				# print OUT "$indel_length[$t]\n";
				# process deletions:
				if ($indel_length[$t]<0){
					while (((-1*$indel_length[$t])+$indel_positions[$t])>length($sequence)){$indel_positions[$t]+=-1; }# print "resetting indel position: $indel_positions[$t]\n";}
					while ($indel_positions[$t]<2){$indel_positions[$t]++; }# print "resetting indel position: $indel_positions[$t]\n";}
					if ($indel_length[$t]<0){
						$mut_seq=substr($mut_seq,0,$indel_positions[$t]-1).("_" x (-1*$indel_length[$t])).substr($mut_seq,$indel_positions[$t]-$indel_length[$t]-1);
						}
				} # if a deletion
				# process insertions:
				elsif ($indel_length[$t]>0){
					while ($indel_positions[$t]>length($sequence)){$indel_positions[$t]+=-1; } # print "resetting indel position: $indel_positions[$t]\n";}
					while ($indel_positions[$t]<2){$indel_positions[$t]++; }# print "resetting indel position: $indel_positions[$t]\n";}
					my $ins="";
					$sequence=substr($sequence,0,$indel_positions[$t]-1).("_" x $indel_length[$t]).substr($sequence,$indel_positions[$t]-1);
					$mut_seq=substr($mut_seq,0,$indel_positions[$t]-1).("x" x $indel_length[$t]).substr($mut_seq,$indel_positions[$t]-1);
				} # if insertion
			} # for all indels
		} # if there is > 0 indel muts	
	} # if there are any muts at all
	
	# print "$sequence\n";
	# print "$mut_seq\n";
	
	# test MD and CIGAR here
	# my $sequence = "GATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC";
	# my $mut_seq  = "tAxxTCGATCGATCGATCGATCGATCGATCG___GATCGATCctTCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATa";
	
	###### process cigar based on reference and mutant sequence created above
	my @cigar_elements=();
	my $M=0; my $D=0; my $I=0; my $offset = 0;
	my $converted_ref = $sequence; $converted_ref =~ s/[GATC]/N/g;
	my $converted_mutseq = $mut_seq; $converted_mutseq =~ s/[GATCgatc]/N/g;
	for my $m (0..length($converted_mutseq)-1){
		my $ref_base = substr($converted_ref, $m+$offset, 1);
		my $mut_seq_base = substr($converted_mutseq, $m, 1);
		# matches
		if ($mut_seq_base eq $ref_base){$M++;}#  print "M: $M\n"; }
		elsif ($M>0){push(@cigar_elements,$M."M"); $M=0; }#  print "x: reset M\n";}
		# deletions
		if ($mut_seq_base eq "_"){$D++;}#  print "D: $D\n"; }
		elsif ($D>0){push(@cigar_elements,$D."D"); $D=0; }#  print "D: reset D\n";}
		# insertions
		if ($mut_seq_base eq "x"){$I++; }# $offset+=-1; }#  print "I: $I\n"; }
		elsif ($I>0){push(@cigar_elements,$I."I"); $I=0; }#  print "I: reset I\n";
	} # for all sites
	if ($M>0){push(@cigar_elements,$M."M");}#  print "end\n";}  # finish up here
	# print cigar
	my $cigar=""; for my $n (0..scalar(@cigar_elements)-1){$cigar=$cigar.$cigar_elements[$n];}
	
	###### process MD based on reference and mutant sequence created above
	my @MD_elements=();
	$M=0; $D=""; $I=0; $offset = 0; 
	for my $m (0..length($mut_seq)-1){
		my $ref_base = substr($sequence, $m+$offset, 1);
		my $mut_seq_base = substr($mut_seq, $m, 1);	
		# insertions
		if ($mut_seq_base eq "x"){} # $offset+=-1;}
		# matches
		if (($mut_seq_base eq $ref_base) or ($mut_seq_base eq "x")){if ($mut_seq_base eq $ref_base){$M++;}}#  print "M: $M\n"; }}
		elsif ($M>0){push(@MD_elements,$M); $M=0; }#  print "x: reset M\n";}		
		# snps
		if (($ref_base ne $mut_seq_base) && ($mut_seq_base =~/[gatc]/)){push(@MD_elements,lc($ref_base));}		
		# deletions
		if ($mut_seq_base eq "_"){$D=$D.$ref_base;}#  print "D: $D\n"; }
		elsif ($D ne ""){push(@MD_elements,$D."D"); $D=""; }#  print "D: reset D\n";}		
		
	} # for all sites
	if ($M>0){push(@MD_elements,$M);  }#  print "end\n";}  # finish up here
	
	# construct initial MD
	my $MD=""; 
	for my $n (0..scalar(@MD_elements)-1){
		if ($MD_elements[$n] =~/D/){$MD=$MD."^".$MD_elements[$n];$MD =~s/D//;}
		else{$MD=$MD.$MD_elements[$n];}
	}
	# final processing of MD
	if (substr($MD,0,1)=~/[GATC]/){$MD="0".$MD;}
	if (substr($MD,length($MD)-1,1)=~/[GATC]/){$MD=$MD."0";}
	my $MD_new=$MD;
	for my $z (1..length($MD)-1){
		if ((substr($MD, $z,1)=~/[gatc]/) && (substr($MD, $z-1,1)=~/[gatc]/)){ 
			$MD_new=substr($MD, 0,$z)."0".substr($MD,$z);
		}
		else{$MD=$MD_new;}
	}			
	
	for my $w (0..length($mut_seq)-1){if (substr($mut_seq,$w,1) eq "x"){substr($mut_seq,$w,1)=random_permutation(@array);}}	# fill in insertions
	$mut_seq=~ s/\_//g;																										# remove deletions
	
	
	return (uc($mut_seq), $cigar, $MD);
	
	
} #end sub

