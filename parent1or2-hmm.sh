#!/usr/bin/env sh

usage () {
    echo usage: `basename $0` -b barcodes -s samdir -o outdir -R Routdir -p parent1 -q parent2 -i indiv -c chroms -w bwaalg
    exit 2
}

die () {
    echo "$1"
    exit ${2:-1}
}

src=$(dirname $0)

while getopts "b:s:o:R:p:q:i:c:x:y:f:g:z:r:t:h:w:e:" opt
do 
  case $opt in
      b) barcodes=$OPTARG ;;
      s) samdir=$OPTARG ;;
      o) outdir=$OPTARG ;;
      R) Routdir=$OPTARG ;;
      p) parent1=$OPTARG ;;
      q) parent2=$OPTARG ;;
      i) indiv=$OPTARG ;;
      c) chroms=$OPTARG ;;
      x) sexchroms=$OPTARG ;;
      y) chroms2plot=$OPTARG ;;
      f) deltapar1=$OPTARG ;;
      g) deltapar2=$OPTARG ;;
      z) priors=$OPTARG ;;
      r) rfac=$OPTARG ;;
      t) theta=$OPTARG ;;
      w) bwaalg=$OPTARG ;;
      e) usestampy=$OPTARG ;;
      *) usage ;;
  esac
done
shift $(($OPTIND - 1))

[ -n "$samdir" ] && [ -n "$outdir" ] && [ -n "$barcodes" ] && [ -n "$parent1" ] && [ -n "$parent2" ] && [ -n "$indiv" ]|| usage
[ -n "$Routdir" ] || usage
[ -d $outdir ] || mkdir -p $outdir

[ -n "$deltapar1" ] || deltapar1=.01
[ -n "$deltapar2" ] || deltapar2=$deltapar1
[ -n "$rfac" ] ||      rfac=.000001

date
echo "version 0.0"

plate=$(echo $indiv | perl -pe 's/^indiv([A-Z][0-9][0-9]?)_.+/$1/')
sex=$(perl -ne "print if /[ACGT]+\t$plate\t/" $barcodes | cut -f4)
[ -z "$sex" ] && sex=female

echo ; echo ; echo "---------------------------------------------------------------------" ; echo

echo "Processing INDIVIDUAL $indiv PLATE $plate SEX $sex DELTA $deltapar1,$deltapar2 RFAC $rfac"

indivdir=$outdir/$indiv
[ -d $indivdir ] || mkdir -p $indivdir

[ -e $indivdir/aln_${indiv}_par1-filtered.sam ] || [ -e $indivdir/aln_${indiv}_par1-filtered.sam.gz ] && \
    [ -e $indivdir/aln_${indiv}_par2-filtered.sam ] || [ -e $indivdir/aln_${indiv}_par2-filtered.sam.gz ] || {

    echo "Extracting reference allele information from SAM files for $indiv ($parent1 and $parent2)"
    echo "python $src/extract-ref-alleles.py -i $indiv -d $samdir -o $indivdir --parent1 $parent1 --parent2 $parent2 --chroms $chroms --bwa_alg $bwaalg --use_stampy $usestampy"
    python $src/extract-ref-alleles.py -i $indiv -d $samdir -o $indivdir --parent1 $parent1 --parent2 $parent2 --chroms $chroms --bwa_alg $bwaalg --use_stampy $usestampy || {
        echo "Error during extract-ref-alleles.py for $indiv"
    }
}

echo "Creating pileup for $indiv"
echo "bash $src/make-pileups.sh -i $indiv -d $indivdir -p $parent1 -q $parent2 2>&1 | grep -vF 'deleted'"
bash $src/make-pileups.sh -i $indiv -d $indivdir -p $parent1 -q $parent2 2>&1 | grep -vF 'deleted'


echo "Writing HMM input data for $indiv"
cmd="Rscript $src/write-hmm-data.R -i $indiv -d $indivdir -c $chroms"
echo $cmd
exec 3>&1; exec 1>&2; echo $cmd; exec 1>&3 3>&-
$cmd || {
    echo "Error during write-hmm-data.R for $indiv"
}


echo "Fitting HMM for $indiv"
Rindivdir=$Routdir/$indiv
[ -d $Rindivdir ] || mkdir -p $Rindivdir
cmd="Rscript $src/fit-hmm.R -d $outdir -i $indiv -s $sex -o $Routdir -p $deltapar1 -q $deltapar2 -r $rfac -c $chroms -x $sexchroms -y $chroms2plot -z $priors -t $theta"

exec 3>&1; exec 1>&2; echo $cmd; exec 1>&3 3>&-
echo $cmd
$cmd || {
    echo "Error during fit-hmm.R for $indiv"
}
echo "Done $indiv"
# block-22 ends here
