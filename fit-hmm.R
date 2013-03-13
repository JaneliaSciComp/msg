#!/usr/bin/env Rscript

## options(warn=2,error=recover);
options(error=quote(q("yes")))

args <- commandArgs()
dollar0 <- substring(args[grep("^--file=", args)], 8)

source(sprintf("%s/ded.R", dirname(dollar0)))
source(sprintf("%s/hmmlib.R", dirname(dollar0)))
library("R.methodsS3", lib.loc = dirname(dollar0))
library("R.oo", lib.loc = dirname(dollar0))

opts <- getopts()
indivs <- unlist(strsplit(opts$i,split=","))
sex <- opts$s
dir <- opts$d
outdir <- opts$o
deltapar1 <- as.numeric(opts$p)
deltapar2 <- as.numeric(opts$q)
rfac <- as.numeric(opts$r)
priors <- unlist(strsplit(opts$z,split=","))
theta <- as.numeric(opts$t)

stopifnot(!is.null(indivs), !is.null(dir), !is.null(outdir), length(indivs) == 1)

one.site.per.read <- TRUE
minCoverage <- 0;

contigLengths <- read.csv("msg.chrLengths",header=T,sep=",",as.is=T);
contigs <- sort(as.vector(contigLengths$chr))
rownames(contigLengths) <- contigLengths$chr

main.contigs <- unlist(strsplit(opts$c,split=","))
plot.contigs <- unlist(strsplit(opts$y,split=","))
if(opts$c == "all") main.contigs <- contigs;
if(opts$y == "all") plot.contigs <- contigs;
sex.chroms <- unlist(strsplit(opts$x,split=","))

aveSpace <- sum(as.numeric(contigLengths[contigLengths$chr %in% plot.contigs,]$length)) / length(plot.contigs)
plotPadding <- 10^(ceiling(log10(aveSpace))-2)

alleles <- c("A","C","G","T")


for(indiv in indivs) {
    cat(indiv, "\n")
    ## if(opts$c == "all")
    ##     all.contigs <- system(sprintf("ls %s/%s/refs/par1  | grep -F '-ref.alleles' | perl -pe 's/-ref.alleles//'", dir, indiv), intern=TRUE)
    ## else all.contigs <- contigs
    
    dataa <- list()
    hmmdata.file <- file.path(outdir, indiv, paste(indiv, "hmmprob.RData", sep="-"))
    if(file.exists(hmmdata.file)) {
        cat("HMM fit for indiv", indiv, "already exists\n")
        dataa <- read.object(hmmdata.file)
        ## next
    } 
    else {
        for(contig in main.contigs) {

            if(sex == "male" && contig %in% sex.chroms) {
                ploidy <- 1
                ancestries <- c("par1","par2")
					 phi <- rep(1/length(ancestries), length(ancestries))
            }
            else {
                ploidy <- 2
                ancestries <- c("par1/par1","par1/par2","par2/par2")
					 phi <- priors
            }
            cat("\t", contig, sex, ploidy, "\n")

            if (!file.exists(sprintf("%s/%s/%s-%s.hmmdata", dir, indiv, indiv, contig))) {
                cat("MISSING file for CONTIG ", contig, " INDIV ", indiv, "\n")
                cat(sprintf("%s/%s/%s-%s.hmmdata", dir, indiv, indiv, contig),"\n")
                next
            } 


            data <- read.data(dir, indiv, contig)
            data$read <- factor.contiguous(data$pos)
            total_sites <- length(unique(data$read));
            cat("\tRound 2: Total number of markers", total_sites, "\n")

            ok <- !is.na(data$bad) | !is.na(data$par1ref) & !is.na(data$par2ref) & !is.na(data$cons)
            cat("\tRound 2: Removing", sum(!ok), "sites at which par1/par2/cons allele unknown\n")
            data$bad[!ok] <- "par1/par2/cons unknown"

            ok <- data$A + data$C + data$G + data$T > 0
            cat("\tRound 2: Removing", sum(!ok), "sites at which cons allele is known but reads are unknown\n")
            data$bad[!ok] <- "reads unknown"

            ok <- !is.na(data$bad) | data$par1ref %in% alleles & data$par2ref %in% alleles
            data$bad[!ok] <- "par1/par2 not in ACGT"
                                        #data$bad[!ok] <- "par1/par2/cons not in ACGT"
            cat("\tRound 2: Removing", sum(!ok), "sites at which par1/par2 ref not %in% {", paste(alleles, collapse=", "), "}\n")
            ## data <- data[ok,,drop=F]
            
            badpos <- data$pos[!is.na(data$bad)]
            data <- data[is.na(data$bad),,drop=F]

            ok <- data$par1ref != data$par2ref
            cat("\tRemoving", sum(!ok), "sites at which par1 == par2\n")
            data <- data[ok,,drop=F]

			data$count <- data$A + data$C + data$G + data$T #+ data$N
			ok <- data$count >= minCoverage
            cat("\tRemoving", sum(!ok), "sites at where coverage is < ",minCoverage,"\n")
			data <- data[ok,,drop=F]

			if (nrow(data)==0) next
            
            if(one.site.per.read) {
                ## Sample one site per read
                data$read <- factor(data$read)
                ## ok <- 1:nrow(data) %in% sapply(levels(data$read),  function(x) sample(which(data$read == x), 1))
                ok <- !duplicated(data$read)
                cat("\tRemoving", sum(!ok), "sites from same reads\n")
                data <- data[ok,,drop=F]
                cat("\tNumber of informative markers:", nrow(data), "\n")
            }


            cat("\tFinal total of", nrow(data), "sites at which par1 != par2\n")
				#if (nrow(data)<20) next
			if (nrow(data)==0) next
			
            L <- nrow(data)
            K <- length(ancestries)

            ## Transition probabilities
            if(contig %in% main.contigs) {
                r <- 1 / contigLengths[contig,"length"]
				} else {
                cat("\tContig ", contig, " not found in main.contigs - defaulting to contig length of ", contigLengths[1,"chr"], "\n")
                r <- 1 / contigLengths[1,"length"] ## Arbitrarily use the first contig for unassembled contigs
				}

            d <- c(NA, diff(data$pos))
            p <- 1 - exp(-r*d*rfac)
            Pi <- array(dim=c(L,K,K), dimnames=list(NULL, ancestries, ancestries))
            if(ploidy == 2) {
                Pi[,"par1/par1","par1/par1"] <- Pi[,"par1/par2","par1/par2"] <- Pi[,"par2/par2","par2/par2"] <- 1-p
                Pi[,"par1/par1","par1/par2"] <- Pi[,"par1/par2","par1/par1"] <- Pi[,"par1/par2","par2/par2"] <- Pi[,"par2/par2","par1/par2"] <- p
                Pi[,"par1/par1","par2/par2"] <- Pi[,"par2/par2","par1/par1"] <- 0
            } else {
                Pi[,"par1","par1"] <- Pi[,"par2","par2"] <- 1-p
                Pi[,"par1","par2"] <- Pi[,"par2","par1"] <- p
            }
            Pi[1,,] <- NA
            
            ## Allele frequencies in parental backgrounds
            ppar1 <- ppar2 <- matrix(NA, nrow=4, ncol=4, dimnames=list(alleles, alleles))
            ppar1[] <- deltapar1/3
            diag(ppar1) <- 1-deltapar1
            ppar2[] <- deltapar2/3
            diag(ppar2) <- 1-deltapar2

            p1 <- ppar1[data$par1ref,,drop=F]
            p2 <- ppar2[data$par2ref,,drop=F]
            p12 <- array(c(p1,p2), dim=c(dim(p1),2))
            dimnames(p12) <- list(NULL, alleles, NULL)

		    ## Take all (<=50) reads for each site
		    N<-min(max(data$A+data$C+data$G+data$T+data$N),50) ## Total number of reads
			#theta <-1 ## quality value correction - defined by passed option
			eps<-paste('eps',seq(1,N,by=1),sep='')
			read<-paste('read',seq(1,N,by=1),sep='')
            y <- data[,c(alleles,"reads","quals","par1ref"),drop=F]
			y$selected_allele <- NA
			y[,eps]<-rep(0,N)
			y[,read]<-rep(5,N)
			fun<-function (x) {
				if (x=='A')	return (0)
				if (x=='C')	return (1)
				if (x=='G')	return (2)
				if (x=='T')	return (3)
				if (x=='N') return (5)
			}
            for(i in 1:nrow(y)) {
				total.reads <- unlist(strsplit(cleanupReadPileup(y[i,"reads"],y[i,"par1ref"]),''))
				y[i,"selected_allele"] <- total.reads[sample(length(total.reads),1)] ## Sample one read for plotting
				qual<-NULL
				qual_corrected<-NULL
				for (s in 1:min(length(total.reads),N)){
					y[i,read[s]] <-lapply(total.reads[s],fun)
					qual<-c(qual,(charToInt(unlist(strsplit(y[i,"quals"],''))[s])-33))
				}
						
				for (g in 1:length(qual)) {
					qual_corrected[g]<-qual[g]*(theta^(order(qual)[g]-1)) ## quality value correction
					y[i,eps[g]] <- 10^(-(qual_corrected[g])/10)
				}
			}

			data$read_allele <- as.vector(y[,"selected_allele"])
			
            ## Emission probabilities
            prob = Pr.y.given.z(y=y[,read,drop=F], p=p12, n=N, eps=y[,eps,drop=F], ploidy=ploidy, C=TRUE, dir=dirname(dollar0), chrom=contig, id=indiv)
   			colnames(prob) <- paste("Pr(y|", ancestries, ")")
            data <- cbind(data, prob)
            data$est <- apply(prob, 1, which.max)
            
            ## Posterior probability
            hmm <- forwardback.ded(Pi=Pi, delta=phi, prob=prob)
            #hmm <- forwardback.ded(Pi=Pi, delta=rep(1/K, K), prob=prob)
            Pr.z.given.y <- exp(hmm$logalpha + hmm$logbeta - hmm$LL)
            colnames(Pr.z.given.y) <- paste("Pr(", ancestries, "|y)")
            data <- cbind(data, Pr.z.given.y)
            attr(data, "badpos") <- badpos
            dataa[[contig]] <- data
        
        }
        cat("Saving data...")
        save(dataa, file=hmmdata.file)
        cat("OK\n")
    }
    
    
    contigLengths <- contigLengths[plot.contigs,]

    ## Track the width of breakpoints
    breakpoints <- {};
    matchMismatch <- {}

    cat("Plotting...")
    plotfile <- file.path(outdir, indiv, paste(indiv, "hmmprob.pdf", sep="-"))
    if(file.exists(plotfile)) { cat("plot already exists\n") ; next }
    pdf(file=plotfile, width=7, height=1.5)
    par(mar=c(2,2.5,0.5,0.5),bg="transparent",cex.main=.68,cex.lab=.8,font.lab=2,cex.axis=.38,mgp=c(1,.000001,0),xaxs="i")

	 plot(0, 0, xlab="", ylab="", col="transparent", xlim=c(1,sum(as.numeric(contigLengths$length)) + plotPadding*(length(plot.contigs)+1)), ylim=c(-1.01,1.01), axes=F) 

	 axis(side=2,at=c(-1,0,1),labels=c("","",""),col="gray38");
	 mtext(c("par2","par1"),side=2,line=.68,at=c(-1,1),font=2,cex=.8,col=c("blue","red"),las=2);
	 box(col="gray68"); 
    
	 current_start	<- plotPadding;
    for(contig in plot.contigs) {

        mtext(side=1,at=current_start,contig,font=2,cex=.8,line=1,xpd=T,adj=0)
        current_end <- current_start + contigLengths[contigLengths$chr == contig,"length"] - 1;

        if(sex == "male" && contig == "X") {
            ploidy <- 1
            ancestries <- c("par1","par2")
            par1homo_col <- 1;
            par2homo_col <- 2;
        }
        else {
            ploidy <- 2
            ancestries <- c("par1/par1","par1/par2","par2/par2")
            par1homo_col <- 1;
            par2homo_col <- 3;
        }

        if (sum(names(dataa) %in% contig)!=0) {
            contig_data <- dataa[[contig]];
            x <- contig_data$pos
            y <- contig_data[,paste("Pr(", ancestries, "|y)")]
	
				### divvy up homozygous and heterozygous blocks
            byBlocks <- breakpoint.width(x, y[,par1homo_col], y[,par2homo_col], indiv=indiv, contig=contig, conf1=.05 ,conf2=.95);
            if (is.null(byBlocks[["bps"]])==F) { breakpoints <- rbind(breakpoints,byBlocks[["bps"]]); }

            ### Output a gff file for geneious
            geneious_data <- {};
            #Calculate confidence thresholds
            gff_thresh <- .95; #!!TODO: make into parameter
            gff_thresh_inverse = 1-gff_thresh; #e.g., .05
            #Get breakpoints using both confidence thresholds
            gff_data <- breakpoint.width(x, y[,par1homo_col], y[,par2homo_col], indiv=indiv, 
                contig=contig, conf1=gff_thresh ,conf2=gff_thresh);
            gff_data_inverse <- breakpoint.width(x, y[,par1homo_col], y[,par2homo_col], 
                indiv=indiv, contig=contig, conf1=gff_thresh_inverse ,conf2=gff_thresh_inverse);
            #If we got any breakpoints, collect them in geneious_data
            #Data we get will have columns: x[end_index]-x[start_index]+1,indiv,contig,start_index,end_index,x[start_index],x[end_index],endpoints
            #looks like:
            #$blocks
            #   V1       V2              V3 V4  V5  V6      V7       V8
            #1 homozygous_par2 13707087 indivH12_TTGACG 2R 181 548 5855429 19562515
            #$bps
            #   V1              V2 V3  V4  V5      V6      V7               V8     V9
            #1 60682 indivH12_TTGACG 2R 180 181 5794748 5855429 5804668.36093989   5804668.36093989
            if (is.null(gff_data[["bps"]])==F) { geneious_data <- rbind(geneious_data,gff_data[["bps"]]); }
            if (is.null(gff_data_inverse[["bps"]])==F) { geneious_data <- rbind(geneious_data,gff_data_inverse[["bps"]]); }
            #Write out geneious_data to file in gff version 3 format (if there were any breakpoints found)
            #Testing:
            #Rscript msg/fit-hmm.R -d hmm_data -i indivH12_TTGACG -s male -o hmm_fit -p 0.01 -q 0.01 -r .000001 -c all -x X -y all -z 0.25,0.5,0.25 -t 1
            if (is.null(geneious_data)==F) {
                #Create a data frame with the following columns:
                #blank, "Geneious", run name, start, end, dot, dot, dot, attributes
                #(I'm not sure what version of gff this is but it matches what David is using)
                geneious_data <- geneious_data[6:7]; #Keep only bps columns (We want 6:7, 4:5 are markers?)
                geneious_data <- cbind(geneious_data, rep("", nrow(geneious_data))); #empty column to start
                geneious_data <- cbind(geneious_data, rep("Geneious", nrow(geneious_data))); #seqid
                geneious_data <- cbind(geneious_data, rep("msg_run", nrow(geneious_data))); #source
                geneious_data <- cbind(geneious_data, rep(".", nrow(geneious_data))); #
                geneious_data <- cbind(geneious_data, rep(".", nrow(geneious_data))); #
                geneious_data <- cbind(geneious_data, rep(".", nrow(geneious_data))); #
                geneious_data <- cbind(geneious_data, rep(paste("name",indiv, sep="="), nrow(geneious_data))); #attributes
                # Reorder columns to match geneious format
                reordered_df <- geneious_data[c(3,4, 5, 1, 2, 6, 7, 8, 9)];
                write.table(reordered_df,file=file.path(outdir, indiv, paste(indiv, contig, "breakpoints.gff", sep="-")),
                    append=F,quote=F,na="NA",row.names=F,col.names=F,sep="\t");
            }

				### plot
            like.par1 <- contig_data[contig_data$read_allele==contig_data$par1ref,]$pos;
            like.par2 <- contig_data[contig_data$read_allele==contig_data$par2ref,]$pos;
            plot.posterior(x+current_start, y, ancestries, like.par1+current_start, like.par2+current_start, bounds=c(1,contigLengths[contigLengths$chr==contig,]$length)+current_start-1, subtract=current_start, tickwidth=5*10^7)

            ### report mismatch fraction for homozygous regions
				if (nrow(byBlocks[["blocks"]])>0) {
					## plot fraction of par1/(par1+par2) among informative markers (between -1 and 1)
            	matchMismatch <- rbind(matchMismatch, reportCounts(contig_data, as.vector(byBlocks[["blocks"]][,"V1"]), as.vector(byBlocks[["blocks"]][,"V4"]), as.numeric(as.vector(byBlocks[["blocks"]][,"V7"])), as.numeric(as.vector(byBlocks[["blocks"]][,"V8"]))))

            }
        }

        current_start <- current_start + contigLengths[contigLengths$chr == contig,"length"] + plotPadding;
        if (contig != plot.contigs[length(plot.contigs)]) {
				abline(v=current_start-(plotPadding/2),col="gray68",lwd=1)
        }
    }

    etc <- ""
    main <- sprintf("%s (%s): delta=(%.0e, %.0e)", indiv, sex, deltapar1, deltapar2)
    dev.off()
    cat("OK\n")
    
    if (is.null(breakpoints)==F) {
        write.table(breakpoints,file=file.path(outdir, indiv, paste(indiv, "breakpoints.csv", sep="-")),append=F,quote=F,na="NA",row.names=F,col.names=F,sep=",");
    } 
    
    if (is.null(matchMismatch)==F) {
        write.table(as.data.frame(matchMismatch),file=file.path(outdir, indiv, paste(indiv, "matchMismatch.csv", sep="-")),
						append=F,quote=F,na="NA",row.names=F,col.names=T,sep=",");
   }
}
