#!/usr/bin/env Rscript

args <- commandArgs()
dollar0 <- substring(args[grep("^--file=", args)], 8)

source(sprintf("%s/ded.R", dirname(dollar0)))
source(sprintf("%s/hmmlib.R", dirname(dollar0)))

## If you are running this code interactively (i.e. separately from a
## run of the MSG software) then don't evalue the following line.
opts <- getopts()

## Instead, if you have done a previous MSG run then code to create
## the 'opts' object will be in the output file for MSG run 3, as a
## result of this line:
cat("Input parameters for summaryPlots.R\nopts <-") ; dput(opts)

dir <- opts$d
contigs <- as.character(unlist(strsplit(opts$c,split=","))); #Require that scaffold names be strings
contigs2plot <- as.character(unlist(strsplit(opts$p,split=",")));

chrLengths <- read.csv("msg.chrLengths", row.names=1)
chrLengths <- structure(chrLengths$length, names=as.character(rownames(chrLengths)))

if (opts$c=="all") {
	contigs <- as.character(names(chrLengths));
} else {
	#contigs <- contigs[contigs %in% names(chrLengths)]
	contigs <- names(chrLengths)[match(contigs,names(chrLengths),nomatch="0")]
	contigs <- contigs[chrLengths[contigs] > 0]
}

if (opts$p=="all") contigs2plot <- as.character(names(chrLengths))


imagedir <- sprintf("%s_images", dir)
if(!file.exists(imagedir)) dir.create(imagedir)

get.ancestry.probs <- function(ancestry, thinfac, difffac, contigs2use, pna.thresh=1, type="all") {
    pp.all <- lapply(contigs2use, interpolate.probs, dir=dir, thinfac=thinfac, difffac=difffac, ancestry=ancestry)
    names(pp.all) <- contigs2use
    pp.pna <- lapply(lapply(pp.all, is.na), colMeans) ### what fraction of individuals is missing data

	 cat(type," ",ancestry,": Numbers of markers after removal for high missing data proportion (%NA<", pna.thresh, "):\n",sep="")
	 print(sapply(pp.pna, function(pna) paste(sum(pna < pna.thresh),"/",length(pna),sep="")))
    pp <- mapply(function(p, pp.pna) p[, pp.pna < pna.thresh, drop=FALSE], pp.all, pp.pna, SIMPLIFY=FALSE)
    names(pp) <- contigs2use
    pp
}

genomeplot <- function(y, ...) {
    plot(x=info$genomepos, y=y, type="l", bty="n", xaxt="n", xlab="", ...)
    abline(v=info$boundaries, lty=2, col="blue")
    axis(side=1, at=info$boundaries, labels=FALSE, tick=TRUE)
    axis(side=1, at=info$midpoints, labels=names(info$midpoints), tick=FALSE, las=2)
}

#Plot par1 segregation:
pp1 <- get.ancestry.probs("par1", thinfac=1, difffac=0, contigs2use=contigs2plot, type="plot") #pp1 never used
p1 <- do.call("cbind", pp1) #p1 never used

contig.lengths <- sapply(pp1, ncol)
contig.fac <- factor(rep(contigs, contig.lengths))
current_contigs <- levels(contig.fac)[match(contigs,levels(contig.fac),nomatch="0")]

pos <- as.integer(colnames(p1))
info <- contig.info(pos, contig.fac, chrLengths[current_contigs])
rm(pp1);
##
## Segregation proportions
##
pdf(file.path(imagedir, "segregation_par1.pdf"), width=10, height=5)
genomeplot(colMeans(p1, na.rm=TRUE), ylab="Mean probability", main="Average probability of par1 homozygosity (hemizygosity)")
abline(h=1/2, col="red", lty=2)
legend("bottomleft", legend="Expectation", col="red", lty=2, bty="n")
dev.off()
rm(p1);

#Plot par2 segregation:
pp2 <- get.ancestry.probs("par2", thinfac=1, difffac=0, contigs2use=contigs2plot, type="plot")
p2 <- do.call("cbind", pp2)

contig.lengths <- sapply(pp2, ncol)
contig.fac <- factor(rep(contigs, contig.lengths))
current_contigs <- levels(contig.fac)[match(contigs,levels(contig.fac),nomatch="0")]

pos <- as.integer(colnames(p2))
info <- contig.info(pos, contig.fac, chrLengths[current_contigs])
rm(pp2); #Reduce memory consumption
##
## Segregation proportions
##
pdf(file.path(imagedir, "segregation_par2.pdf"), width=10, height=5)
genomeplot(colMeans(p2, na.rm=TRUE), ylab="Mean probability", main="Average probability of par2 homozygosity (hemizygosity)")
abline(h=1/2, col="red", lty=2)
legend("bottomleft", legend="Expectation", col="red", lty=2, bty="n")
dev.off()
rm(p2);
cat("All done!");
