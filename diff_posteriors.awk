#!/bin/awk -f
#Expects 2 ancestry-prob-par*.tsv files from different MSG runs as arguments.
#The idea is to detect any differences between the files while ignoring
# permutations of rows and only emitting missing or additional AIMs.
#If you get no output, that means the inputs are effectively identical.
#This is particularly useful for validating MSG installations with example
# datasets like the toy, the Yang et al. (2019) Current Biology subset, etc.
BEGIN{
   FS="\t";
   OFS=FS;
}
#Establish a list and a hash of AIMs in the first file:
FNR==NR&&FNR==1{
   for (i=2; i<=NF; i++) {
      marker[i]=$i;
      markermap[$i]=i;
   };
}
#Store the individual IDs and posteriors from the first file:
FNR==NR&&FNR>1{
   indiv[$1]=FNR;
   for (i=2; i<=NF; i++) {
      posterior[$1,marker[i]]=$i;
   };
}
#Establish a list of AIMs in the second file, and output if any are
# missing or added relative to the first file:
#We don't detect or indicate if the mismatch is due to addition vs. omission.
#This can be inferred by the user, since markers are in sorted order.
FNR<NR&&FNR==1{
   for (i=2; i<=NF; i++) {
      markerb[i]=$i;
      if (length(markermap[$i]) == 0) {
         print "mismatch of markers "marker[i]" and "$i;
      };
   };
}
#Now check for matching posteriors for matched individuals.
FNR<NR&&FNR>1{
   indivb[$1]=FNR;
   for (i=2; i<=NF; i++) {
      if (posterior[$1,markerb[i]] != $i) {
         print "mismatch of posteriors for "$1" at "markerb[i]": "posterior[$1,markerb[i]]" vs. "$i;
      };
   };
}
#Finally, check for missing or added individuals:
#We do this rather inefficiently, with two loops acting as asymmetric
# set differences
END{
   for (id in indiv) {
      if (length(indivb[id]) == 0) {
         print "individual "id" missing from file 2";
      };
   };
   for (idb in indivb) {
      if (length(indiv[idb]) == 0) {
         print "individual "idb" missing from file 1";
      };
   };
}
