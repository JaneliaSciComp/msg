#!/usr/bin/env Rscript

forwardback.ded <- function (x, Pi, delta, distn, pm, pn = NULL, fortran = FALSE, prob) {
    ## This is a slight modification of function forwardback in library HiddenMarkov
    ## - changed fortran default to FALSE
    ## - allow transition prob matrix to vary along the sequence
    ## - allow observation probs to be provided explicitly
    if(fortran) require(HiddenMarkov)

    m <- dim(Pi)[2]
    if(missing(prob)) { ## DED
        n <- length(x)
        dfunc <- makedensity(distn)
        prob <- matrix(as.double(0), nrow = n, ncol = m)
        for (k in 1:m) prob[, k] <- dfunc(x = x, getj(pm, k), pn, log = FALSE)
    }
    else { ## DED: inserted this option to provide observation probabilities
        n <- nrow(prob)
        stopifnot(dim(prob) == c(n,m))
    }
    phi <- as.double(delta)
    logalpha <- matrix(as.double(rep(0, m * n)), nrow = n)
    lscale <- as.double(0)
    cat("forward\n")
    if (fortran != TRUE) {
        for (i in 1:n) {
            if (i > 1) 
                phi <- phi %*% array(Pi[i,,],dim=dim(Pi)[2:3]) ## DED: allow transition matrix to vary along sequence
            phi <- phi * prob[i, ]
            sumphi <- sum(phi)
            phi <- phi/sumphi
            lscale <- lscale + log(sumphi)
            logalpha[i, ] <- log(phi) + lscale
            if(!i %% 1000) cat(i, "\r")
        }
        LL <- lscale
        cat("\n")
    }
    else {
        if (!is.double(Pi)) 
            stop("Pi is not double precision")
        if (!is.double(prob)) 
            stop("prob is not double precision")
        memory0 <- rep(as.double(0), m)
        loop1 <- .Fortran("loop1", m, n, phi, prob, Pi, logalpha, 
            lscale, memory0, PACKAGE = "HiddenMarkov")
        logalpha <- loop1[[6]]
        LL <- loop1[[7]]
    }
    logbeta <- matrix(as.double(rep(0, m * n)), nrow = n)
    phi <- as.double(rep(1/m, m))
    lscale <- as.double(log(m))
    cat("backward\n")
    if (fortran != TRUE) {
		if (n>1) {
        for (i in seq(n - 1, 1, -1)) {
            phi <- array(Pi[i+1,,],dim=dim(Pi)[2:3]) %*% (prob[i + 1, ] * phi) ## DED
            logbeta[i, ] <- log(phi) + lscale
            sumphi <- sum(phi)
            phi <- phi/sumphi
            lscale <- lscale + log(sumphi)
            if(!i %% 1000) cat(i, "\r")
        }
        cat("\n")
		}
    }
    else {
        memory0 <- rep(as.double(0), m)
        loop2 <- .Fortran("loop2", m, n, phi, prob, Pi, logbeta, 
            lscale, memory0, PACKAGE = "HiddenMarkov")
        logbeta <- loop2[[6]]
    }
    return(list(logalpha = logalpha, logbeta = logbeta, LL = LL))
}

Pr.y.given.z <- function(y, p, eps, n=50, ploidy=2, log=TRUE, C=FALSE, norm=FALSE, dir=".", chrom="none", id="id") {
    alleles <- c("A","C","G","T")
    J <- length(alleles)
    L <- nrow(y)
    if(C) {
        stopifnot(ploidy %in% c(1,2))
        #tmpfile <- "hmmprobs.input"
        tmpfile <- paste(tempfile(),".",id,sep="")

        stopifnot(dim(p) == c(L, J, 2))
        p1 <- matrix(p[,,1], nrow=L)
        p2 <- matrix(p[,,2], nrow=L)

        write.table(cbind(y, p1, p2, eps), file=tmpfile,
                    sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

        pipa <- pipe(sprintf("%s/hmmprobs -t %d -p %d -n %d < %s", dir, L, ploidy, n, tmpfile))
        #pipa <- pipe(sprintf("%s/hmmprobs -t %d -e %f -p %d < %s", dir, L, eps, ploidy, tmpfile))
		  on.exit({close(pipa) ; file.remove(tmpfile) })
		  #on.exit({close(pipa) ; })
        return(matrix(scan(pipa, quiet=TRUE), nrow=L, byrow=TRUE))
    }
    stopifnot(ploidy == 2)
    g <- make.genotypes(alleles)
    Pr.allele.given.genotype <- g/2 * (1 - 3*eps/4) + (1 - g/2) * eps/4
    ancestries <- c("0","1","2")
    
    Pr.genotype.given.ancestry <-
        sapply(ancestries,
               function(z) array(dim=c(J,J,L), dimnames=list(alleles, alleles, NULL)), simplify=FALSE)
    
    ## cat("Pr(g|z)")
    for(a1 in alleles)
        for(a2 in alleles) {
            ## cat("\n\t", a1, a2)
            Pr.genotype.given.ancestry[["0"]][a1,a2,] <- p[,a1,1] * p[,a2,1]
            Pr.genotype.given.ancestry[["1"]][a1,a2,] <- p[,a1,1] * p[,a2,2]
            Pr.genotype.given.ancestry[["2"]][a1,a2,] <- p[,a1,2] * p[,a2,2]
        }

    ## cat("\narray -> list")
    for(z in ancestries) {
        ## cat("\n\t", z)
        Pr.genotype.given.ancestry[[z]] <-
            lapply(1:L, function(l) Pr.genotype.given.ancestry[[z]][,,l])
    }
    
    ## cat("\nPr(y|g)")

    Pr.yl.given.genotype <-
        if(norm) function(yl)
            apply(Pr.allele.given.genotype, c(1,2), function(p) dmultinom(yl, prob=p))
        else function(yl)
            apply(Pr.allele.given.genotype, c(1,2), function(p) prod(p^yl))
    
    ## cat("\narray -> list")
    Pr.y.given.genotype <- lapply(1:L, function(l) { if(! l %% 1000) cat(l, "\r") ; Pr.yl.given.genotype(y[l,]) })

    ## cat("\nsum_g (Pr(y|g)Pr(g|z))")
    ## on.exit(cat("\n"))
    sapply(ancestries, function(z)
           sapply(mapply("*", Pr.y.given.genotype, Pr.genotype.given.ancestry[[z]], SIMPLIFY=FALSE), sum))
}

make.genotypes <- function(alleles) {
    J <- length(alleles)
    g <- array(0, dim=c(J,J,J), dimnames=rep(list(alleles), 3))
    for(a in alleles) {
        g[a,,a] <- 1
        g[,a,a] <- g[,a,a] + 1
    }
    g
}

read.data <- function(dir, indiv, contig) {
    f <- sprintf("%s/%s/%s-%s.hmmdata", dir, indiv, indiv, contig)
    if (file.exists(f)==T) {
		cat("Reading data from", f, "\n")
    	read.delim(f, header=TRUE, as.is=TRUE, na.strings="", quote="") ### need to eliminate quoting b/c of QVs
    } else { cat(paste("missing",f)); }
}

plot.posterior <- function(x, Pr.z.given.y, ancestries, like.par1, like.par2, bounds=NULL, subtract=0, tickwidth=10^6, ...) {
    K <- length(ancestries)
    colnames(Pr.z.given.y) <- ancestries
    col <- if(K == 3) c("red","purple","blue") else c("red","blue")
    #col <- if(K == 3) c("blue","purple","red") else c("blue","red")
    names(col) <- ancestries
    
    lines(bounds,c(0,0),lwd=5,col="gray68",lend=1);
    lines(x=x, y=Pr.z.given.y[,1], ylim=c(-1,1), lend=1, type="n", ...)

    plot_to  <- c(0:floor((bounds[2]-bounds[1]+1)/tickwidth));
    plot_to_labels <- rep("",length(plot_to));
    plot_indices <- seq(0,max(plot_to),by=1);
    #plot_indices <- seq(0,max(plot_to),by=5);
    plot_to_labels[plot_indices+1] <- plot_indices*5;
#	 if (length(plot_to_labels) > 1) { plot_to_labels[length(plot_to_labels)] <- ""; 
#	 } else { 
#		plot_to_labels <- c(plot_to_labels);
#		plot_to <- c(plot_to);
#	 }

    axis(side=1,at=(plot_to*tickwidth)+bounds[1],labels=plot_to_labels,tck=-.05,col="transparent",col.ticks="gray28",adj=1);

    rug(like.par1,side=3,col=col[1],lwd=.5);
    rug(like.par2,side=1,col=col[length(ancestries)],lwd=.5);
    #if(!missing(badpos)) rug(badpos, side=3, col="red")

    homozygous.ancestries <- if(K == 3) ancestries[c(1,3)] else ancestries
    Pr.z.given.y[,homozygous.ancestries[2]] <- -Pr.z.given.y[,homozygous.ancestries[2]]
    L <- length(x)
    for(anc in homozygous.ancestries)
        polygon(rbind(c(x[1],0),
                      cbind(x,Pr.z.given.y[,anc]),
                      c(x[L], 0)),
                col=col[anc])
}


findBreak <- function(indiv,contig,y,z1,z2,conf1,conf2,parentage) {
	z1.rle <- rle(z1)
	z2.rle <- rle(y>=conf2)

	blocks <- {}
	widths <- {}

   ### LOCATE BREAKPOINTS
   last_index <- 0;
	if (length(z1.rle$values)>1) {
		### breakpoint: transition
		for (i in 1:length(z1.rle$values)) {

			### find the first time the other condition is satisfied
			start_index <- last_index + z1.rle$lengths[i]
			if (i<length(z1.rle$values)) {
				end_index <- match(z1.rle$values[i],z2[(start_index+1):(start_index + z1.rle$lengths[i+1]-1)]) + start_index
				if (is.na(end_index)==F) {
					if ((y[start_index]<conf1 & y[end_index]>conf2) | (y[end_index]<conf1 & y[start_index]>conf2)) {
						endpoints <- approx(x=y[start_index:end_index], y=x[start_index:end_index], xout=c(conf1,conf2))$y;
					} else { endpoints <- c(NA,NA); }
					widths <- rbind(widths,c(x[end_index]-x[start_index]+1,indiv,contig,start_index,end_index,x[start_index],x[end_index],endpoints));
				}
			}
			last_index <- last_index + z1.rle$lengths[i]
		}
	} 

   ### LOCATE HOMOZYGOUS BLOCKS
   last_index <- 0;
	if (length(z2.rle$values)>0) {
		### breakpoint: transition
		for (i in 1:length(z2.rle$values)) {
			start_index <- last_index + 1;
			end_index <- start_index + z2.rle$lengths[i] -1;
         if (z2.rle$values[i]==TRUE) {
			   blocks <- rbind(blocks,c(parentage,x[end_index]-x[start_index]+1,indiv,contig,start_index,end_index,x[start_index],x[end_index]));
         }
         last_index <- last_index + z2.rle$lengths[i]
      }
   }

	list(blocks=as.data.frame(blocks),bps=as.data.frame(widths));
}

breakpoint.width <- function(x, y1, y2, indiv=NA, contig=NA, conf1=.05 ,conf2=.95) {
	width1 <- findBreak(indiv,contig,y1,y1>=conf1,y1<=conf2,conf1,conf2,'homozygous_par1')
	width2 <- findBreak(indiv,contig,y2,y2>=conf1,y2<=conf2,conf1,conf2,'homozygous_par2')

	list(blocks=rbind(width1[["blocks"]],width2[["blocks"]]),bps=rbind(width1[["bps"]],width2[["bps"]]))
}	


breakpoint.width.old <- function(x, y1, y2, indiv=NA, contig=NA, conf1=.05 ,conf2=.95) {
	z1 <- y1>=conf1 & y1<=conf2;
	z2 <- y2>=conf1 & y2<=conf2;
	z  <- z1 | z2;
	z.rle <- rle(z);

	blocks <- {};
	widths <- {};

	### find the endpoint indicies where z==TRUE
	last_index <- 1;
	for (i in 1:length(z.rle$values)) {

		### breakpoint: transition
		if (z.rle$values[i]==TRUE) {
			start_index <- last_index;
			end_index <- start_index + z.rle$lengths[i] -1;

			### must be flanked by FALSE on both sides
			if (start_index > 1 & end_index < (length(x)-1)) {
				if ((y1[start_index-1]<conf1 & y1[end_index+1]>conf2) | (y1[end_index+1]<conf1 & y1[start_index-1]>conf2)) {
					endpoints <- approx(x=y1[c(start_index-1):(end_index+1)], y=x[c(start_index-1):c(end_index+1)], xout=c(conf1,conf2))$y;
				} else if ((y2[start_index-1]<conf1 & y2[end_index+1]>conf2) | (y2[end_index+1]<conf1 & y2[start_index-1]>conf2)) {
					endpoints <- approx(x=y2[c(start_index-1):(end_index+1)], y=x[c(start_index-1):c(end_index+1)], xout=c(conf1,conf2))$y;
				} else { 
					endpoints <- c(NA,NA);
				}

				widths <- rbind(widths,c(x[end_index]-x[start_index]+1,indiv,contig,start_index,end_index,x[start_index],x[end_index],endpoints));
			}

		### block
		} else {
			start_index <- last_index + 1;
			end_index <- start_index + z.rle$lengths[i] -1 -1;

			if ((sum(y1[start_index:end_index]>conf2)>0) & (sum(y2[start_index:end_index]<conf1)>0)) { type <- "homozygous_par1" } 
			else if ((sum(y1[start_index:end_index]<conf1)>0) & (sum(y2[start_index:end_index]>conf2)>0)) { type <- "homozygous_par2" } 
			else { type <- "heterozygous"; }

			blocks <- rbind(blocks,c(type,x[end_index]-x[start_index]+1,indiv,contig,start_index,end_index,x[start_index],x[end_index]));
		}

		last_index <- last_index + z.rle$lengths[i];
	}

	list(blocks=as.data.frame(blocks),bps=as.data.frame(widths));
}


retrieve.data <- function(dir,globType) {
    ff <- list.files(dir, paste("\\-",globType,"$",sep=""), recursive=TRUE, full.names=TRUE)
    inames <- gsub(paste("-",globType,"$",sep=""), "", basename(ff))

    n <- length(ff)
    fi <- function(i) {
        cat(inames[i], "\r")
        res <- try(read.csv(ff[i], header=F,sep=","))
        if(!inherits(res, "try-error")) res
        #else recover() ## NULL
    }

    x <- lapply(1:n, fi)
    names(x) <- inames;
    x
}


reportCounts <- function(contig_data,block_types,contigs,window_starts,window_ends) {
	window_counts <- {};

	for (w in 1:length(window_starts)) {
		window_data <- contig_data[contig_data$pos>=window_starts[w] & contig_data$pos<=window_ends[w],];
		if (nrow(window_data)>0) {
			par1.par2 <- nrow(window_data[window_data$par1ref!=window_data$par2ref,]);
			like.par1 <- length(window_data[window_data$read_allele==window_data$par1ref,]$pos);
			like.par2 <- length(window_data[window_data$read_allele==window_data$par2ref,]$pos);
			#unlike.par1par2 <- length(window_data[window_data$read_allele!=window_data$par1ref & window_data$read_allele!=window_data$par2ref,]$pos);
			window_counts <- rbind(window_counts,c(block_types[w],contigs[w],window_starts[w],window_ends[w],par1.par2,like.par1,like.par2))#,unlike.par1par2));
		}
	}

	if (is.null(window_counts)==F) {
		window_counts <- as.data.frame(window_counts);
		names(window_counts) <- c("ancestry","contig","start","end","par1_par2_diff","like_par1","like_par2");
		window_counts;
	}
}


reportCounts.old <- function(contig_data,contig_start=1,contig_length,window_size=NA) {
	window_counts <- {};
	if (is.na(window_size)==T) { windows <- c(contig_start,contig_length); } 
	else { windows <- seq(from=contig_start,to=contig_length,by=window_size); }

	for (w in 2:length(windows)) {
		window_data <- contig_data[contig_data$pos>=windows[w-1] & contig_data$pos<windows[w],];
		if (nrow(window_data)>0) {
			par1.par2 <- nrow(window_data[window_data$par1ref!=window_data$par2ref,]);
			like.par1 <- length(window_data[window_data$read_allele==window_data$par1ref,]$pos);
			like.par2 <- length(window_data[window_data$read_allele==window_data$par2ref,]$pos);
			unlike.par1par2 <- length(window_data[window_data$read_allele!=window_data$par1ref & window_data$read_allele!=window_data$par2ref,]$pos);
			window_counts <- rbind(window_counts,c(windows[w-1],windows[w],par1.par2,like.par1,like.par2,unlike.par1par2));
			#window_counts <- rbind(window_counts,c((windows[w-1]+windows[w])/2,par1.par2,like.par1,like.par2,unlike.par1par2));
		}
	}

	if (is.null(window_counts)==F) {
		window_counts <- as.data.frame(window_counts);
		names(window_counts) <- c("start","end","par1.par2","like.par1","like.par2","unlike.par1par2");
		#names(window_counts) <- c("midpoint","par1.par2","like.par1","like.par2","unlike.par1par2");
		window_counts;
	}
}


map.indivs <- function(dir, f, indivs, pass.inames=FALSE, ...) {
    ff <- list.files(dir, "\\.RData$", recursive=TRUE, full.names=TRUE)
    inames <- gsub("-hmmprob.RData$", "", basename(ff))
    if(!missing(indivs)) {
        stopifnot(indivs %in% inames)
        ok <- inames %in% indivs
        ff <- ff[ok]
        inames <- inames[ok]
    }
    n <- length(ff)
    fi <- function(i) {
        ## cat(inames[i], "\r")
        if(pass.inames)
            res <- try(f(read.object(ff[i]), iname=inames[i], ...))
        else
            res <- try(f(read.object(ff[i]), ...))
        if(!inherits(res, "try-error")) res
        #else NULL
        else recover() ## NULL
    }
    x <- lapply(1:n, fi)
    names(x) <- inames
    x
}


map.indivs.hmmdata <- function(dir, contig, f=identity, indivs, pass.inames=FALSE, ...) {
    ff <- list.files(dir, sprintf("%s\\.hmmdata$", contig), recursive=TRUE, full.names=TRUE)
    inames <- gsub(sprintf("-%s.hmmdata$", contig), "", basename(ff))
    if(!missing(indivs)) {
        stopifnot(indivs %in% inames)
        ok <- inames %in% indivs
        ff <- ff[ok]
        inames <- inames[ok]
    }
    n <- length(ff)
    readf <- function(file) read.delim(file)
    fi <- function(i) {
        cat(inames[i], "\r")
        if(pass.inames)
            res <- try(f(read.object(ff[i]), iname=inames[i], ...))
        else
            res <- try(f(readf(ff[i]), ...))
        if(!inherits(res, "try-error")) res
        else recover() ## NULL
    }
    x <- lapply(1:n, fi)
    names(x) <- inames
    x
}
    

tabf <- function(data) {
    ## Proportion of >1 coverage sites with >1 allele
    f <- function(x) {
        y <- x[,c("A","C","G","T")]
        nreads <- rowSums(y)
        ispoly <- rowSums(y > 0) > 1
        mean(ispoly[nreads > 1])
    }
    sapply(data, f)
}

contig.info <- function(pos, fac, lengths) {
    boundaries <- c(0, cumsum(lengths))
    nsites <- as.data.frame(table(fac))
    contigs <- levels(fac)[match(names(lengths),levels(fac))]
    offsets <- rep(boundaries[-length(boundaries)], nsites[match(contigs,nsites$fac),]$Freq)
    #nsites <- table(fac)
    #contigs <- levels(fac)
    #offsets <- rep(boundaries[-length(boundaries)], nsites)
    names(offsets) <- fac
    midpoints <- boundaries[-length(boundaries)] + diff(boundaries)/2
    names(midpoints) <- names(boundaries)[-1]
    chrBreaks <- which(c(FALSE, diff(as.integer(fac)) != 0))
    #chrBreaks <- which(c(FALSE, diff(as.integer(fac)) == 1))
    names(chrBreaks) <- paste(contigs[-length(contigs)], contigs[-1], sep="/")
    list(pos=pos, fac=fac, boundaries=boundaries, midpoints=midpoints,
         offsets=offsets, genomepos=pos+offsets, chrBreaks=chrBreaks)
}

interpolate.probs <- function(contig, dir, thinfac=1, difffac=0, ancestry) {
    f <- function(d) {
        d <- d[[contig]]
        ancestries <- c(paste(ancestry, ancestry, sep="/"), ancestry)
        pnames <- paste("Pr(", ancestries, "|y)")
        names(pnames) <- ancestries
        if (!is.null(nrow(d))) {
            pname <- pnames[pnames %in% colnames(d)]
            #if(pname == pnames[ancestry]) cat("haploid\n")
            d[,c("pos", pname)]
        }
        else NULL
    }
    
    dd <- map.indivs(dir, f)

	 if (sum(unlist(lapply(dd,is.null)))==length(dd)) { #NULL
		 y <- t(as.matrix(t(rep(NA,length(dd)))))
		 rownames(y) <- names(dd)
		 colnames(y) <- "0"
		 y

	 } else {
	    ## retrieve the union of all positions (unique and sorted)
	    x <- sort(unique(unlist(lapply(dd, "[[", "pos"))))
	    n <- length(x)
		 
		 ## thin by taking somewhat uniformly spaced markers
		 if(thinfac < 1) {
		     thin <- seq(1, n, length.out=floor(n*thinfac))
		     x <- x[thin]
		 }
		 
		 ## determine linearly interpolated values for markers
		 #y <- sapply(dd, function(d) approx(d, xout=x)$y)
		 #y <- sapply(dd, function(d) { if (!is.null(d) && (nrow(d)>1)) { approx(d, xout=x)$y } else { rep(NA,n) } })
		 y <- sapply(dd, function(d) { 
		    if (!is.null(d) && (nrow(d)>1)) { 
		   	 interpolated <- approx(d, xout=x)$y 
		   	 if (sum(is.na(interpolated)==T)==n) { interpolated <- rep(mean(d[,2],na.rm=T),n) }
		    } else if (!is.null(d) && (nrow(d)==1)) { interpolated <- rep(d[[2]],n)
		    } else { interpolated <- rep(NA,n) } 
		    interpolated
		    })
	
		 if (n>1) { 
			if(difffac > 0) {
		 	    keep <- c(TRUE, rowSums(abs(diff(y)) > difffac, na.rm=TRUE) > 0)
		 	    y <- y[keep,,drop=FALSE]
		 	    x <- x[keep]
		 	}
		 } else { y <- as.matrix(t(y)) }
			 
		 rownames(y) <- x
		 t(y)
    }
}


est.rf.p <- function(p, lod=FALSE, na.rm=TRUE) {
    ## p[i,j] = PosteriorProb(individual i has homozygous ancestry at site j)
    ## The calculations below follow from
    ## b = E(#rec events in indiv i between sites j and k) =
    ##      1 x Pr(hom j)Pr(het k) + 1 x Pr(het j)Pr(hom k)
    ## Other outcomes either impossible in backcross or contribute
    ## zero recobinations to expectation.
    n <- nrow(p)
    sump <- colSums(p, na.rm=na.rm)
    if(!na.rm)
        cp <- crossprod(p)
    else {
        cov <- cov(p, use="pairwise.complete.obs")
        cp <- (n-1)*cov + sump %o% sump / n
    }
    Enr <- outer(sump, sump, "+") - 2*cp
    if(lod) {
        nrlogr <- Enr*log10(Enr)
        nrlogr[!is.finite(nrlogr)] <- 0

        n1rlog1r <- (n-Enr)*log10(n-Enr)
        n1rlog1r[!is.finite(n1rlog1r)] <- 0
        
        nrlogr - n*log10(n) + n1rlog1r + n*log10(2)
    }
    else
        Enr/n
}

est.rf.p.profile <- function(p, p1, offdiag=FALSE, lod=TRUE, na.rm=TRUE, logf=log10) {
    n <- nrow(p)
    L <- ncol(p)
    if(offdiag) {
        stopifnot(missing(p1))
        p1 <- p[,-1]
        p <- p[,-L]
    }
    rhat <- colMeans(p + p1 - 2*p*p1, na.rm=na.rm)
    if(lod) {
        nrlogr <- n*rhat*logf(rhat)
        nrlogr[!is.finite(nrlogr)] <- 0

        n1rlog1r <- n*(1-rhat)*logf(1-rhat)
        n1rlog1r[!is.finite(n1rlog1r)] <- 0
        
        ans <- nrlogr + n1rlog1r + n*logf(2)
    }
    else
        ans <- rhat
    if(offdiag) c(NA, ans) else ans
}

plot.rf.ded <- function (z, pos, type="imageplot", zmax=12, info, ...) {
    genomewide <- !missing(info)
    
    if(genomewide) {
        op <- par(xpd = TRUE, las = 1, lwd=6, cex.axis=.68, mar=c(8,8,.5,.5), adj=0)
        on.exit(par(op))
    }
    if(type != "n") {
        diag(z) <- zmax
        z[!is.na(z) & z > zmax] <- zmax
        z[is.na(z)] <- -1
    }
    if(type != "n") {
        spectrumf <- rev
        if ('gamma' %in% names(formals(rainbow))) { 
            col <- c("lightgray", spectrumf(rainbow(256, start = 0, end = 2/3, gamma = 0.6)))
        } else {
            col <- c("lightgray", spectrumf(rainbow(256, start = 0, end = 2/3)))
        }
        ## col <- c("lightgray", spectrumf(rainbow(256, start = 1/3, end = 1, gamma = 0.6)))
        br <- c(-1, seq(-1e-06, zmax, length = 257))
    }
    else { col <- NA ; br <- 1:2 }

    stopifnot(dim(z) == length(pos))
    d <- diff(pos)
    L <- length(pos)
    ## x <- c(pos[1]-d[1]/2, pos[-L] + d/2, pos[L] + d[L-1]/2)
    x <- c(0,pos) * 1e-6
    xaxt <- yaxt <- if(genomewide) "n" else "s"
    image(x, x, z, ylab = "", xlab = "", xaxt=xaxt, yaxt=yaxt, breaks = br, col = col)
    if(genomewide) {
        ## abline(v=info$boundaries[-1], h=info$boundaries[-1])
        ## abline(v=x, h=x, lty=2, lwd=0.05)
		  #bounds <- x[1+unlist(lapply(info$chrBreaks,function(t) mean(t+c(-1:0))))]
		  bounds <- x[info$chrBreaks] + 1e-6
		  abline(v=bounds, h=bounds, lty=1, xpd=FALSE)
        axis(1, at=c(0,bounds), labels=names(info$midpoints), tick=FALSE, line=-.8, las=3, hadj=1, padj=1)
        axis(2, at=c(0,bounds), labels=names(info$midpoints), tick=FALSE, line=-.5, las=1, hadj=1, padj=0)
    }
}
    
plot.hmm.fits <- function(dir, outdir=".", info, n=96) {
   plotf <- function(data, iname, contig) {
       data <- data[[contig]]
       etc <- ""
       if(contig == "X" && "Pr( par1 |y)" %in% colnames(data)) {
           etc <- "male"
           ancestries <- c("par1","par2")
       }
       else ancestries <- c("par1/par1","par1/par2","par2/par2")
       x <- data$pos
       y <- data[,paste("Pr(", ancestries, "|y)")]
       plot.posterior(x, y, ancestries, main=contig, ylab=paste(iname, etc), xlab="")
   }

   for(contig in c("2L","2R","3L","3R","4","X")) {
		 pdf(file=sprintf("%s/%s.pdf", outdir, contig), width=10, height=200)
       #png(sprintf("%s/%s.png", outdir, contig), width=1000, height=20000)
       par(mfrow=c(n, 1), mar=c(1,4,1,4), bg="transparent")
       map.indivs(dir, f=plotf, contig=contig, pass.inames=TRUE)
       dev.off()
   }
}

msg.write.table <- function(d, file)
    write.table(d, file=file, sep="\t", col.names=NA, quote=FALSE)


read.cross.msg <- function(file.AA, file.AB,
                           sex=rep("f", length(indivs)), phenos=NULL, na.rm=FALSE) {
    ## Read MSG ancestry probabilities from file and return an R/qtl
    ## "cross" object.

    ## MSG ancestry probabilities are contained in a simple tabular
    ## format with individual names in first column, and labels of the
    ## form "contig:position" in the first row

    require(qtl)
    genotypes <- c("A","H","B") ; alleles <- c("A","B")
    if(missing(file.AB)) backcross <- TRUE
    
    ## Read MSG probabilities
    prob.AA <- as.matrix(read.table(file.AA, header=TRUE, row.names=1,
                                    as.is=TRUE, check.names=FALSE))
    if(!backcross) {
        prob.AB <- as.matrix(read.table(file.AB, header=TRUE, row.names=1))
        stopifnot(dim(prob.AB) == dim(prob.AA))
        stopifnot(rownames(prob.AB) == indivs)
        stopifnot(colnames(prob.AB) == markers)
        if(na.rm) {
            bad <- colSums(is.na(prob.AA)) > 0 | colSums(is.na(prob.AB)) > 0
            cat("Removing", sum(bad), "markers with missing ancestry probabilities\n")
            prob.AA <- prob.AA[,!bad]
            prob.AB <- prob.AB[,!bad]
        }
    }
    else if(na.rm) {
        bad <- colSums(is.na(prob.AA)) > 0
        cat("Removing", sum(bad), "markers with missing ancestry probabilities\n")
        prob.AA <- prob.AA[,!bad]
    }
    indivs <- rownames(prob.AA)
    markers <- colnames(prob.AA)
    
    ## Construct data frame corresponding to R/qtl csv format
    col.labels <- strsplit(markers, ":")
    contigs <- sapply(col.labels, "[", 1)
    contigs.f <- factor(contigs)
    bp <- as.integer(sapply(col.labels, "[", 2))

    ## Get sex data (default is all female)
    if(is.null(sex)) sex <- rep(0, length(indivs))
    else {
        stopifnot(length(sex) == length(indivs))
        if(is.character(sex)) {
            sex <- factor(sex)
        }
        if(is.factor(sex)) {
            stopifnot(levels(sex) == c("female","male") || levels(sex) == c("f","m"))
            sex <- as.integer(sex) - 1
        }
    }
    stopifnot(sex %in% c(0,1))
    sex.column <- c("sex", "", sex)

    ## TODO: My inclination is to set the genotypes to all missing, so
    ## that we know we are just using the probability
    ## information. However, I haven't quite got estimating the map
    ## from the probabilities working, so I'm making naive thresholded
    ## calls as follows. These hard genotypes are used for the map
    ## only and then discarded at the end of this function.

    ## genos <- array(dim=c(length(indivs), length(markers)))
    ## if(backcross)
    ## Help read.cross recognise that we have backcross data; we
    ## will set these to missing later.
    ## genos[1:2] <- genotypes[1:2]

    genos <- prob.AA
    genos[] <- ifelse(prob.AA > 0.5, "A", "H")

    males <- which(sex == 1)
    X <- which(contigs.f == "X")
    isMaleX <- row(genos) %in% males & col(genos) %in% X
    genos[isMaleX][genos[isMaleX] == "H"] <- "B"
    
    ## Get phenotypes
    if(is.null(phenos)) phenos <- list(rep(0, length(indivs)))
    else {
        if(!is.list(phenos)) phenos <- list(phenos)
        stopifnot(sapply(phenos, length) == length(indivs))
    }
    phenotype.columns <-
        lapply(seq_along(phenos), function(i) c(sprintf("pheno%d", i), "", phenos[[i]]))
    phenotype.columns <- do.call("cbind", phenotype.columns)

    cross.data <- cbind(phenotype.columns, sex.column, rbind(bp, contigs, genos))

    ## Write this to disk as a csv file
    temp.file <- tempfile()
    on.exit(file.remove(temp.file))
    write.table(cross.data, file=temp.file, sep=",",
                row.names=FALSE, col.names=FALSE, quote=FALSE)

    ## Read it is as an R/qtl cross object

    ## TODO: We should estimate the map from the ancestry
    ## probabilities, but for now we allow R/qtl to do so using the
    ## simplistic thresholded genotype calls made above.

    cross <- read.cross(format="csv", file=temp.file,
                        ##  estimate.map=FALSE,
                        genotypes=genotypes, alleles=alleles)
    
    ## Construct arrays of probabilities, split up by contigs
    ## (Could use abind package rather than bindfuns below)
    split.probs <- function(p) {
        p <- split(as.data.frame(t(p)), f=contigs.f)
        lapply(lapply(p, as.matrix), t)
    }
    prob.AA <- split.probs(prob.AA)
    if(backcross) {
        prob.AB <- lapply(prob.AA, function(p) 1-p)
        bindfun <- function(pAA, pAB, markers) {
            p <- array(dim=c(dim(pAA), 2),
                       dimnames=list(indivs, markers, genotypes[1:2]))
            p[,,1] <- pAA
            p[,,2] <- pAB
            p
        }
        probs <- mapply(bindfun, prob.AA, prob.AB, split(markers, f=contigs.f), SIMPLIFY=FALSE)
    }
    else {
        prob.AB <- split.probs(prob.AB)
        prob.BB <- mapply(function(pAA, pBB) 1 - pAA - pBB, prob.AA, prob.AB, SIMPLIFY=FALSE)
        bindfun <- function(pAA, pAB, pBB, markers) {
            p <- array(dim=c(dim(pAA), 3),
                       dimnames=list(indivs, markers, genotypes))
            p[,,1] <- pAA
            p[,,2] <- pAB
            p[,,3] <- pBB
            p
        }
        probs <- mapply(bindfun, prob.AA, prob.AB, prob.BB, split(markers, f=contigs.f), SIMPLIFY=FALSE)
    }

    ## Set the genotype probability slots, in the same way as
    ## calc.genoprob does. We may need to work out suitable values for
    ## some of those attributes.

    step <- 0 ; off.end <- 0 ; stepwidth <- "fixed"
    
    for(contig in levels(contigs.f)) {
        cross$geno[[contig]]$data[] <- NA
        cross$geno[[contig]]$prob <- probs[[contig]]
        attr(cross$geno[[contig]]$prob, "map") <-
            create.map(cross$geno[[contig]]$map, step, off.end, stepwidth)
        attr(cross$geno[[contig]]$prob, "error.prob") <- NA
        attr(cross$geno[[contig]]$prob, "step") <- step
        attr(cross$geno[[contig]]$prob, "off.end") <- off.end
        attr(cross$geno[[contig]]$prob, "map.function") <- "MSG"
        attr(cross$geno[[contig]]$prob, "stepwidth") <- stepwidth
    }
    cross
}

removeNs <- function(bases,quals) {
	bases <- toupper(bases)
	paste(unlist(strsplit(quals,split=""))[which(unlist(strsplit(bases,split="")) %in% c("A","T","C","G"))],sep="")
}

### taken from write-hmm-data.R
cleanupReadPileup <- function(x, ref) {
    x <- toupper(x)

    ## x is a vector of pileup base codes. Examples are
    ## c..c   .+2cc   ^fa^f,^f,
    ## ref is a vector containing the corresponding reference alleles
    ## http://samtools.sourceforge.net/pileup.shtml
    x <- gsub("\\^.", "", x) ## begin contiguous something...
    x <- gsub("\\$", "", x)  ## end contiguous something...

	 ## get rid of indels for now
    x <- gsub("^[\\+-][ACGTNXMRWSYKVHDBacgtnxmrwsykvhdb]+", "", x)
    x <- gsub(".*\\*.*", "", x)

    len <- 1
    while(length(grep("[\\+-]", x)) > 0) {
        ## re <- paste("[\\+-]", len, paste(rep('[ACGTNacgtn]', len), collapse=""), sep="")
        re <- paste("[\\+-]", len, paste(rep('[ACGTNXMRWSYKVHDBacgtnxmrwsykvhdb]', len), collapse=""), sep="")
        x <- gsub(re, "", x)
        len <- len + 1
        ## We used to die if there were any indels > length 10 but we decided there's no harm in allowing them through.
        #if(len > 10) {
        #    print("hmmlib.R: replacing more than 10")
        #    print(x)
        #    stop()
        #}
    }

    x <- mapply(gsub, pattern=list("[\\.,]"), replacement=ref, x=x)
	 return(x);
}
