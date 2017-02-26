#!/usr/bin/env sh

usage () {
    echo usage: `basename $0` -b barcodes -o outdir -R Routdir -i indiv -c chroms
    exit 2
}

die () {
    echo "$1"
    exit ${2:-1}
}

src=$(dirname $0)

while getopts "b:o:R:i:c:x:y:f:g:z:a:r:t:m:u:j:" opt
do 
  case $opt in
      b) barcodes=$OPTARG ;;
      o) outdir=$OPTARG ;;
      R) Routdir=$OPTARG ;;
      i) indiv=$OPTARG ;;
      c) chroms=$OPTARG ;;
      x) sexchroms=$OPTARG ;;
      y) chroms2plot=$OPTARG ;;
      f) deltapar1=$OPTARG ;;
      g) deltapar2=$OPTARG ;;
      z) priors=$OPTARG ;;
      a) recRate=$OPTARG ;;      
      r) rfac=$OPTARG ;;
      t) theta=$OPTARG ;;
      m) gff_thresh_conf=$OPTARG ;;
      u) one_site_per_contig=$OPTARG ;;
      j) pepthresh=$OPTARG ;;
      *) usage ;;
  esac
done
shift $(($OPTIND - 1))

[ -n "$outdir" ] && [ -n "$barcodes" ] && [ -n "$indiv" ]|| usage
[ -n "$Routdir" ] || usage
if [[ ! -d $outdir ]]; then
   echo "Error: It appears you may not have run MSG on the parentals yet."
   echo "Please run MSG on the parentals one time before trying to optimize parameters,"
   echo "as optimizing the parameters requires some files made by MSG."
   exit 3
fi

[ -n "$deltapar1" ] || deltapar1=.01
[ -n "$deltapar2" ] || deltapar2=$deltapar1
[ -n "$recRate" ] ||   recRate=0
[ -n "$rfac" ] ||      rfac=.000001


#Barcode file structure: barcode indiv_id plate_id sex
#Column delimiter should still be tabs, but some text editors misinterpret the tab key as a fixed number of spaces.
#Thus we handle the general case of the delimiter being one or more whitespace characters (i.e. \s+).
#We also handle the general case for the individual's ID, as just being "indiv" followed by one or more non-whitespace characters.
indivnumber=$(echo $indiv | perl -pe 's/^indiv(\S+)_.+/$1/')
sex=$(perl -ne "print if /^[ACGT]+\s+$indivnumber\s+/" $barcodes | awk ' { print $4; } ')
plate=$(perl -ne "print if /^[ACGT]+\s+$indivnumber\s+/" $barcodes | awk ' { print $3; } ')
[ -z "$sex" ] && sex=female

echo "Processing INDIVIDUAL $indiv PLATE $plate SEX $sex DELTA $deltapar1,$deltapar2 RECRATE $recRate RFAC $rfac"

#Check for existence of the input file directory:
indivdir=$outdir/$indiv
if [[ ! -d $indivdir ]]; then
   echo "Error: It appears you may not have run MSG on the parentals yet."
   echo "Please run MSG on the parentals one time before trying to optimize parameters,"
   echo "as optimizing the parameters requires some files made by MSG."
   exit 3
fi
#Check for existence of the hmmdata files:
if [[ -z $(ls $indivdir | grep '\.hmmdata') ]]; then
   echo "Error: write-hmm-data.R failed to create any hmmdata files for $indiv."
   exit 4
fi


echo "Fitting HMM for $indiv"
Rindivdir=$Routdir/$indiv
[ -d $Rindivdir ] || mkdir -p $Rindivdir
cmd="Rscript $src/fit-hmm.R -d $outdir -i $indiv -s $sex -o $Routdir -p $deltapar1 -q $deltapar2 -a $recRate -r $rfac -c $chroms -x $sexchroms -y $chroms2plot -z $priors -t $theta -g $gff_thresh_conf -u $one_site_per_contig -j $pepthresh"

echo $cmd
$cmd || {
    echo "Error during fit-hmm.R for $indiv"
}
echo "Done $indiv"
