#!/usr/bin/env Rscript

args <- commandArgs()
dollar0 <- substring(args[grep("^--file=", args)], 8)

source(file.path(dirname(dollar0), "ded.R"))
options(error=quote(q("yes")))

decode.pileup.bases <- function(x, ref) {
    ## x is a vector of pileup base codes. Examples are
    ## c..c   .+2cc   ^fa^f,^f,
    ## ref is a vector containing the corresponding reference alleles
    ## http://samtools.sourceforge.net/pileup.shtml
    x <- gsub("\\^.", "", x) ## begin contiguous something...
    x <- gsub("\\$", "", x)  ## end contiguous something...

    ## get rid of indels for now
    ## x <- gsub("^[\\+-][ACGTNacgtn]+", "", x)
    x <- gsub("^[\\+-][ACGTNXMRWSYKVHDBacgtnxmrwsykvhdb]+", "", x)
    x <- gsub(".*\\*.*", "", x)

    len <- 1
    while(length(grep("[\\+-]", x)) > 0) {
        ## re <- paste("[\\+-]", len, paste(rep('[ACGTNacgtn]', len), collapse=""), sep="")
        re <- paste("[\\+-]", len, paste(rep('[ACGTNXMRWSYKVHDBacgtnxmrwsykvhdb]', len), collapse=""), sep="")
        x <- gsub(re, "", x)
        len <- len + 1
        
        ## We used to die if there were any indels > length 10 but we decided there's no harm in allowing them through.
        ##if(len > 10) {
        ##    print(x)
        ##    stop()
        ##}
    }
    ## x <- gsub("N", "", x)

    ## Replace [.,] with reference base
    cat("mapply...")
    x <- mapply(gsub, pattern=list("[\\.,]"), replacement=ref, x=x)
    x <- toupper(x)

    ## Count occurrences of alleles
    cat("countalleles...")
    tmpfile=tempfile()
    cat(x, file=tmpfile, sep="\n")
    pipa <- pipe(sprintf("%s/countalleles < %s", dirname(dollar0), tmpfile))
    tab <- matrix(scan(pipa, what=integer()), ncol=5, byrow=TRUE)
    close(pipa) ; file.remove(tmpfile)
    colnames(tab) <- allele.states
    cat("\n")
    tab
}



opts <- getopts()
indiv <- opts$i
dir <- opts$d

stopifnot(!is.null(indiv), !is.null(dir), !is.null(opts$c))

species <- c("par1", "par2")
pupsp <- "par1"
allele.states <- c("A","C","G","T","N")
if(opts$c == "all") {
    contigs <- system(sprintf("ls %s/refs/par1  | grep -- '-ref.alleles' | perl -pe 's/-ref.alleles//'", dir), intern=TRUE)
} else {
    contigs <- unlist(strsplit(opts$c,split=",")); ## "__CONTIG__"
}

cat("Writing HMM input data file for", length(contigs), "contigs.\n")

for(contig in contigs) {
    cat("Writing HMM input data file for", indiv, contig, "\n")

    outfile <- sprintf("%s/%s-%s.hmmdata", dir, indiv, contig, sep="-")
    if(file.exists(outfile)) {
        cat("File exists -- skipping:", outfile, "\n")
        next
    }
    pupfile <- sprintf("%s/aln_%s_%s-filtered-%s-sorted.pileup", dir, indiv, pupsp, contig)
    if(!file.exists(pupfile)) {
        cat("MISSING ",pupfile,"\n")
        next
    } else {
        ## indels in pileup result in lines with > 10 fields
        ## Even so, apparently read.table does not read these files correctly.
        ## It ends up with only 315140 lines out of over 2 million.
        ## pup <- read.delim(pipe(paste("cut -f1-10 <", pupfile)), header=FALSE, as.is=TRUE, na.strings=NULL)
        ## pup <- read.delim(pipe(paste("cut -f2,3,4,9 <", pupfile)), header=FALSE, as.is=TRUE, na.strings=NULL)
        ## colnames(pup) <- c("contig","pos","ref","cons","v","w","x","y","reads","z")
        pipa <- pipe(paste("cut -f2,3,4,9,10 <", pupfile))
        pup <- scan(pipa, what="", sep="\n")
        close(pipa)
        pup <- strsplit(pup, "\t")
        stopifnot(sapply(pup, length) == 5)
        pup <- matrix(unlist(pup), ncol=5, byrow=TRUE)
        pup <- data.frame(pos=as.integer(pup[,1]), ref=unlist(lapply(pup[,2],toupper)), cons=pup[,3], reads=pup[,4], quals=pup[,5])
        
        cat("\tStarting with a total of", nrow(pup), "positions.\n")
        ok <- pup$ref %in% allele.states
        cat("\tRemoving", sum(!ok), "positions at which ref is not [ACGT].\n")
        pup <- pup[ok,,drop=F]
                                        #pup <- pup[ok,]

        pup <- pup[grep("^\\*+$",as.vector(pup$reads),invert=T),,drop=F]
        cat("\tKeeping", nrow(pup), "non-indels in par1...")
	
        cat("\tDecoding read bases...")
        tab <- decode.pileup.bases(pup$reads, pup$ref)
	
        if(nrow(tab) != nrow(pup)) {
            print(str(tab))
            print(str(pup))
            stop("nrows differ after decoding bases")
        }
        pup <- cbind(pup, tab)
        pup$bad <- rep("", nrow(pup))
        cat("\n")

        for(sp in species) {
            reffile <- sprintf('%s/refs/%s/%s-ref.alleles', dir, sp, contig)
            ref <- read.delim(reffile, header=FALSE, as.is=TRUE, col.names=c("pos","allele"))

### remove positions that are indels in either par1 or par2 (based on ref)
            ## non_indels_pos <- pup[is.na(pup$pos)==F,]$pos
            non_indels_pos <- intersect(ref$pos,pup$pos)
            pup <- pup[pup$pos %in% non_indels_pos,,drop=F]
            ref <- ref[ref$pos %in% non_indels_pos,,drop=F]
            cat("\tKeeping", length(non_indels_pos), "non-indels...")

            if(nrow(ref) != nrow(pup)) {
                ## warning(paste("MINOR WARNING (write-hmm-data)",reffile,"nrow(ref) != nrow(pup)\nbodging\n"))
                map <- match(pup$pos, ref$pos)
                ref <- ref[map,]
            }
            stopifnot(nrow(ref) == nrow(pup), ref$pos == pup$pos)
            spref <- sprintf("%sref", sp)
            pup[[spref]] <- ref$allele
            
            ok <- !is.na(pup[[spref]])
            cat("\tRemoving", sum(!ok), "positions at which", sp, "ref was NA or N.\n")
            ## pup <- pup[ok,]
            pup$bad[!ok] <- paste(sp, "NA/N")
            
            ok <- pup[[spref]] %in% c("A","C","G","T")
            cat("\tRemoving", sum(!ok), "positions at which", sp, "ref is not [ACGT].\n")
            ## pup <- pup[ok,]
            pup$bad[!ok] <- paste(sp, "not ACGT")
        }

        ok <- !is.na(pup$ref)
    	cat("\tRemoving", sum(!ok), "positions at which par1 ref (pileup) is NA\n")
    	## pup <- pup[ok,]
    	pup$bad[!ok] <- "par1 ref NA"
    	
    	ok <- pup$ref %in% c("A","C","G","T")
    	cat("\tRemoving", sum(!ok), "positions at which par1 ref (pileup) is not [ACGT].\n")
    	## pup <- pup[ok,]
    	pup$bad[!ok] <- "par1 ref not ACGT"
    	
    	ok <- pup$ref == pup$par1ref
    	cat("\tRemoving", sum(!ok), "positions at which par1 ref (pileup) disagrees with par1 ref.\n")
    	## pup <- pup[ok,]
    	pup$bad[!ok] <- "par1 ref disagree"
    	
    	cat("\tAfter filtering, keeping", nrow(pup), "positions.\n")
    	write.table(pup, file=outfile, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
    }
}
