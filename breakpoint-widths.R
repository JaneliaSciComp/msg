#!/usr/bin/env Rscript

reportQuants <- function(data) {
	return(do.call("sprintf",c(fmt="quantiles \n.05=%.2f \n.50=%.2f \n.95=%.2f ",
		as.list(quantile(data,probs=c(.05,.5,.95),na.rm=T)))));
}


plotHists <- function(filename,widths,bin_breaks,bin_names,xlab_name) {
	pdf(sprintf("%s/%s",figdir,filename),width=5,height=3);
	par(bg="transparent",mar=c(3,2.5,1,.1),cex.lab=.68,cex.axis=.38,cex.main=.68,col="gray28",mgp=c(1.5,.5,0));

	a <- hist(widths,breaks=bin_breaks,plot=F)$counts;
	mp<-barplot(a/sum(a),names.arg=NULL,las=1,col="gray68",bor="transparent",axes=F,xpd=T);
	median_width <- median(widths,na.rm=T);

	axis(1,at=mp,bin_names);
	axis(2,line=0);
	mtext(side=1,xlab_name,line=1.8,font=2,col="black")
	mtext(side=2,"frequency",line=1.2,font=2,col="black");
	median_index <- which(bin_breaks[-1]>median_width)[1];
	points(mp[median_index],a[median_index]/sum(a),pch="*",cex=1.8,col="darkred",xpd=T)
	dev.off()

	cat(sprintf("%s: median %d",filename,median(widths)));
}


plotBarplot <- function(filename,data,data.names,xlab_name) {
	pdf(sprintf("%s/%s",figdir,filename),width=5,height=3)
	par(bg="transparent",mar=c(3,3,1,.5),cex.lab=.68,cex.axis=.38,cex.main=.68,col="gray28",mgp=c(1.5,.5,0));
	barplot(data,names.arg=data.names,las=2,cex.axis=.5,cex.lab=.8,cex.names=.28,main="",col="gray28",bor="transparent",ylab="");
	mtext(side=2,line=1.8,xlab_name,col="black")
	mtext(side=3,line=-2,reportQuants(data),col="black",adj=1,cex=.38)
	box(col="gray68")
	dev.off()
}


figdir <- "figure_summary"
if(!file.exists(figdir)) dir.create(figdir)


### plot the histogram of breakpoint widths
####################################################################################################
width_data <- retrieve.data(dir,"breakpoints.csv");
widths <- abs(as.vector(unlist(lapply(width_data,"[[","V9"))) - as.vector(unlist(lapply(width_data,"[[","V8"))));
widths <- widths[is.na(widths)==F]

bin_breaks <- c(c(0:20)*10^4,10^9);
bin_names <- (bin_breaks[-length(bin_breaks)]+(bin_breaks[-1]-bin_breaks[-length(bin_breaks)])/2) /10^3;
bin_names[length(bin_names)] <- ">200"
plotHists("hist_bp_widths.pdf",widths,bin_breaks,bin_names,"breakpoint width (kb)");


####################################################################################################
### plot the distribution of informative markers
numMarkers <- list();
globType <- "hmmprob.RData";
ff <- list.files(dir, paste("\\-",globType,"$",sep=""), recursive=TRUE, full.names=TRUE)
inames <- gsub(paste("-",globType,"$",sep=""), "", basename(ff))
for (indiv in inames) {
	hmmdata.file <- file.path(dir, indiv, paste(indiv, "hmmprob.RData", sep="-"))
	dataa <- read.object(hmmdata.file)
	numMarkers[[indiv]] = 0;

   for(contig in contigs) {
		if (sum(names(dataa) %in% contig)!=0) {
			contig_data <- dataa[[contig]];
			numMarkers[[indiv]] = numMarkers[[indiv]] + nrow(contig_data)
		}
	}
}

write.table(t(as.data.frame(numMarkers)),file="numInformativeMarkers.tsv",quote=F,col.names=F,row.names=T)
markerData <- unlist(numMarkers);
markerData <- markerData[order(markerData)]
plotBarplot("hist_informative_markers.pdf",markerData,names(markerData),"# informative markers");

### plot the distribution of informative reads
####################################################################################################
system(sprintf("python %s/countReads.py -d hmm_data -b %s -o filtered_read_stats.txt", dirname(dollar0), bc))
data <- read.csv("filtered_read_stats.txt",header=T,as.is=T)
filtered_sorted_data <- data[order(data$num_reads,decreasing=T),];
plotBarplot("hist_filtered_reads.pdf",filtered_sorted_data$num_reads,as.vector(filtered_sorted_data$indiv),"# filtered reads");


### plot the distribution of unfiltered reads
####################################################################################################
data <- read.csv(sprintf("%s_stats.txt",read),header=F,sep="\t"); 
names(data) <- c("indiv","num_reads");
data <- data[data$indiv %in% plate_order,];
unfiltered_sorted_data <- data[order(data$num_reads,decreasing=T),];
plotBarplot("hist_unfiltered_reads.pdf",unfiltered_sorted_data$num_reads,as.vector(unfiltered_sorted_data$indiv),"# reads");
