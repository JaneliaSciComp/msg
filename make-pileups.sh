#!/usr/bin/env sh

# [[file:~/Work/simsec/org/simsec.org::*Create%20sim%20pileup%20for%20each%20contig][block-26]]
#!/bin/bash -x

usage () {
    echo usage: `basename $0` -i indiv -d dir -a all -r range -p parent1 -q parent2
    exit 2
}

die () {
    echo "$1"
    exit ${2:-1}
}

all=1
#Defaults $all to true

while getopts "ai:d:p:q:" opt
do 
  case $opt in
      i) indiv=$OPTARG ;;
      d) dir=$OPTARG ;;
      a) all=$OPTARG ;;
      r) range=$OPTARG ;;
      p) parent1=$OPTARG ;;
      q) parent2=$OPTARG ;;
      *) usage ;;
  esac
done
shift $(($OPTIND - 1))

[ -n "$indiv" ] && [ -n "$dir" ] && [ -n "$parent1" ] && [ -n "$parent2" ] || usage

species=par1
[ -e $parent1 ] || die "$parent1 doesn't exist"

if [ $all -eq 1 ]; then
   refs=$(ls $dir/refs/par1 | grep '\-ref\.alleles' | perl -pe 's/-ref.alleles//')
   if [ -z $range ]; then
      refs=$(echo ${refs} | sed -n "${range} p")
   fi
fi

file=$dir/aln_${indiv}_${species}-filtered
echo "$file"

[ -e $parent1.fai ] || samtools faidx $parent1
[ -e $file-sorted.bam ] || {
    echo "samtools view -bt $parent1.fai $file.sam | samtools sort - $file-sorted"
    samtools view -bt $parent1.fai $file.sam | samtools sort - $file-sorted
}
[ -e $file-sorted.bam.bai ] || samtools index $file-sorted.bam

for ref in $refs ; do
    [ -e $file-$ref-sorted.pileup ] || {
        echo "Making pileup for $species contig $ref"
        samtools view -bu $file-sorted.bam $ref | samtools sort - $file-$ref-sorted
        samtools index $file-$ref-sorted.bam
        samtools pileup -Bcf $parent1 $file-$ref-sorted.bam > $file-$ref-sorted.pileup
        rm $file-$ref-sorted.bam $file-$ref-sorted.bam.bai
    }
done

## rm -f $file-sorted.bam $file-sorted.bam.bai
# block-26 ends here
