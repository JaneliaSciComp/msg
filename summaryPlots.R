#!/usr/bin/env Rscript

source(sprintf("%s/ded.R", "msg"))
source(sprintf("%s/hmmlib.R", "msg"))

plot.correlation.matrix <- TRUE
write.ancestry.probs <- TRUE
breakpoint.widths <- FALSE

## If you are running this code interactively (i.e. separately from a
## run of the MSG software) then don't evalue the following line.
opts <- getopts()

## Instead, if you have done a previous MSG run then code to create
## the 'opts' object will be in the output file for MSG run 3, as a
## result of this line:
cat("Input parameters for summaryPlots.R\nopts <-") ; dput(opts)

bc <- opts$b
dir <- opts$d
thinfac <- as.numeric(opts$t)
difffac <- as.numeric(opts$f)
pna.thresh <- as.numeric(opts$n)
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
    fname <- sprintf("ancestry-probs-%s-%f-%f-%f.%s.rda", ancestry, thinfac, difffac, pna.thresh, type)
    if(!file.exists(fname)) {
        pp.all <- lapply(contigs2use, interpolate.probs, dir=dir, thinfac=thinfac, difffac=difffac, ancestry=ancestry)
        names(pp.all) <- contigs2use
        save(pp.all, file=fname)
    } else {
		 cat(type," ",ancestry,": ",fname," file exists\n",sep="")
		 pp.all <- read.object(fname)

		 ### check if the contigs read in match contigs2use; if not, need to re-read the data in
		 if (sum(contigs2use %in% names(pp.all))!=length(contigs2use)) {
			pp.all <- lapply(contigs2use, interpolate.probs, dir=dir, thinfac=thinfac, difffac=difffac, ancestry=ancestry)
			names(pp.all) <- contigs2use
			save(pp.all, file=fname)
		 }
    }
    pp.pna <- lapply(lapply(pp.all, is.na), colMeans) ### what fraction of individuals is missing data

	 cat(type," ",ancestry,": Numbers of markers after removal for high missing data proportion (%NA<", pna.thresh, "):\n",sep="")
	 print(sapply(pp.pna, function(pna) paste(sum(pna < pna.thresh),"/",length(pna),sep="")))
    pp <- mapply(function(p, pp.pna) p[, pp.pna < pna.thresh, drop=FALSE], pp.all, pp.pna, SIMPLIFY=FALSE)
    names(pp) <- contigs2use
    pp
}

### full marker set, no thinning
cat("Extracting markers for the following contigs:\n",opts$c);
pp1 <- get.ancestry.probs("par1", thinfac=1, difffac=0, contigs2use=contigs, type="all")
pp <- pp2 <- get.ancestry.probs("par2", thinfac=1, difffac=0, contigs2use=contigs, type="all")
p1 <- do.call("cbind", pp1)
rm(pp1) #Reduce memory consumption
p <- p2 <- do.call("cbind", pp2)
rm(pp2); #Reduce memory consumption

contig.lengths <- sapply(pp, ncol)
rm(pp); #Reduce memory consumption
contig.fac <- factor(rep(contigs, contig.lengths))
#current_contigs <- contigs[contigs %in% levels(contig.fac)]
current_contigs <- levels(contig.fac)[match(contigs,levels(contig.fac),nomatch="0")]

pos <- as.integer(colnames(p))
rm(p); #Reduce memory consumption
info <- contig.info(pos, contig.fac, chrLengths[current_contigs])

if(write.ancestry.probs) {
    p1.table <- p1
    rm(p1); #Reduce memory consumption
    p2.table <- p2
    rm(p2); #Reduce memory consumption
    colnames(p1.table) <- colnames(p2.table) <- paste(contig.fac, pos, sep=":")
    p12.table <- 1 - (p1.table + p2.table)
    msg.write.table(round(p1.table, 6), file="ancestry-probs-par1.tsv")
    msg.write.table(round(p2.table, 6), file="ancestry-probs-par2.tsv")
    msg.write.table(round(p12.table, 6), file="ancestry-probs-par1par2.tsv")
    rm(p1.table, p2.table, p12.table); #Reduce memory consumption
} else {
    rm(p1, p2); #Reduce memory consumption
}

genomeplot <- function(y, ...) {
    plot(x=info$genomepos, y=y, type="l", bty="n", xaxt="n", xlab="", ...)
    abline(v=info$boundaries, lty=2, col="blue")
    axis(side=1, at=info$boundaries, labels=FALSE, tick=TRUE)
    axis(side=1, at=info$midpoints, labels=names(info$midpoints), tick=FALSE, las=2)
}

##
## Off-diagonal LOD profile
##
## TODO: should save these R objects into Robjects output dir
cat("\nOff-diagonal LOD profile...\n");
#pp1 <- get.ancestry.probs("par1", thinfac=1, difffac=0, contigs2use=contigs2plot, type="plot") #pp1 never used
pp <- pp2 <- get.ancestry.probs("par2", thinfac=1, difffac=0, contigs2use=contigs2plot, type="plot")
#p1 <- do.call("cbind", pp1) #p1 never used
p <- p2 <- do.call("cbind", pp2)
rm(p2, pp2); #Reduce memory consumption

contig.lengths <- sapply(pp, ncol)
rm(pp); #Reduce memory consumption
contig.fac <- factor(rep(contigs2plot, contig.lengths))
current_contigs <- levels(contig.fac)[match(contigs2plot,levels(contig.fac),nomatch="0")]

pos <- as.integer(colnames(p))
info <- contig.info(pos, contig.fac, chrLengths[current_contigs])

rhat.offdiag <- est.rf.p.profile(p=p, offdiag=TRUE, lod=FALSE)
save(rhat.offdiag, file="rhat-offdiag.rda")

rlod.offdiag <- est.rf.p.profile(p=p, offdiag=TRUE, lod=TRUE)
save(rlod.offdiag, file="rlod-offdiag.rda")

### save offdiag data
write.table(cbind(as.vector(contig.fac), pos, rhat.offdiag, rlod.offdiag),file="offdiagonal_data.tsv",sep="\t",quote=F,row.names=F,col.names=c("chrom", "pos", "rhat", "rlod"))

pdf(file.path(imagedir, "offdiagonal-lod.pdf"), width=10, height=5)
genomeplot(rlod.offdiag, ylab="LOD")
abline(h=25, col="red")
dev.off()

##
## Segregation proportions
##
pdf(file.path(imagedir, "segregation.pdf"), width=10, height=5)
genomeplot(colMeans(p, na.rm=TRUE), ylab="Mean probability", main="Average probability of par2 homozygosity (hemizygosity)")
abline(h=1/2, col="red", lty=2)
legend("bottomleft", legend="Expectation", col="red", lty=2, bty="n")
dev.off()

##
## Missing data proportions
##
pdf(file.path(imagedir, "missing.pdf"), width=10, height=5)
genomeplot(colMeans(is.na(p)), ylab="Missing proportion")
dev.off()

rm(p, rhat.offdiag, rlod.offdiag); #Reduce memory consumption


if(plot.correlation.matrix) {
    ##
    ## Thinned ancestry probabilities
    ##
    cat("\nThinning ancestry probabilities...\n");
    
    ## There is a problem currently in that the thinning procedure
    ## results in a different set of markers being selected for
    ## par1 and par2. While that is unresolved, we only output
    ## probabilities for par2.
    
    ## pp1.thin <- get.ancestry.probs("par1", thinfac=thinfac, difffac=difffac)
    pp.thin <- pp2.thin <- get.ancestry.probs("par2", thinfac=thinfac, difffac=difffac, contigs=contigs2plot, pna.thresh=pna.thresh, type="thinned_plot")
    ## p1.thin <- do.call("cbind", pp1.thin)
    p.thin <- p2.thin <- do.call("cbind", pp2.thin)
    rm(pp2.thin); #Reduce memory consumption
    contig.lengths.thin <- sapply(pp.thin, ncol)
    contig.fac.thin <- factor(rep(contigs2plot, contig.lengths.thin))
	 current_contigs <- levels(contig.fac.thin)[match(contigs2plot,levels(contig.fac.thin),nomatch="0")] ### Tina

	 if (length(current_contigs) > 0) {
	    pos.thin <- as.integer(colnames(p.thin))
	    info.thin <- contig.info(pos.thin, contig.fac.thin, chrLengths[current_contigs])
	
	    if(write.ancestry.probs) {
	        ## p1.thin.table <- p1.thin
	        p2.thin.table <- p2.thin
	        ## colnames(p1.thin.table) <- paste(contig.fac.thin, pos, sep=":")
	        colnames(p2.thin.table) <- paste(contig.fac.thin, pos.thin, sep=":")
	        ## p12.thin.table <- 1 - (p1.thin.table + p2.thin.table)
	        ## msg.write.table(round(p1.thin.table, 6), file="ancestry-probs-thinned-par1.tsv")
	        msg.write.table(round(p2.thin.table, 6), file="ancestry-probs-thinned-par2.tsv")
	        ## msg.write.table(round(p12.thin.table, 6), file="ancestry-probs-thinned-par1par2.tsv")
           rm(p2.thin.table); #Reduce memory consumption
	    }
	
	    ##
	    ## plot heatmap
	    ##
	    n <- nrow(pp.thin[[1]])
       rm(pp.thin); #Reduce memory consumption
	    lod.max <- n*log10(2)
	    lod.thin.fname <- sprintf("lod-thin-%f-%f-%f.rda",thinfac,difffac,pna.thresh)

	    if(!file.exists(lod.thin.fname)) {
	        lod.thin <- est.rf.p(p.thin, lod=TRUE, na.rm=TRUE)
	        save(lod.thin, file=lod.thin.fname)
	    } else {
	        lod.thin <- read.object(lod.thin.fname)
	    }
       rm(p.thin); #Reduce memory consumption
	
	    cat("Plotting heatmap...\n");
	    bitmap(file=file.path(imagedir, "lod-matrix.bmp"), width=100, height=100, bg="transparent");
	    plot.rf.ded(pos=info.thin$genomepos, lod.thin, zmax=lod.max, info=info.thin)
	    dev.off();
	} else {
		cat("Unable to create thinned lod-matrix.bmp plot: Try increasing pna.thresh (",pna.thresh,") or decreasing difffac (",difffac,")\n",sep="")
	}
}

if(breakpoint.widths) source(sprintf("%s/breakpoint-widths.R", "msg"))
