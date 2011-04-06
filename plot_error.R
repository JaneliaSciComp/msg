pdf("hmm_fit_images/error_gamma.pdf", width=6, height=10, pointsize=14)
par(mfcol= c(2,1))
err<-read.table("error_gamma", header=T)

if (sum(is.na(err$gam_par1)==F) > 0) {
   mean_gam1<-mean(err$gam_par1, na.rm=T)
   hist(err$gam_par1, main=mean_gam1, xlab="error (gamma) par1", xlim=c(0,0.2))
   abline(v=mean_gam1, col="blue", lwd=2, lty="dotted")
}
if (sum(is.na(err$gam_par2)==F) > 0) {
   mean_gam2<-mean(err$gam_par2, na.rm=T)
   hist(err$gam_par2, main=mean_gam2, xlab="error (gamma) par2", xlim=c(0,0.2))
   abline(v=mean_gam2, col="blue", lwd=2, lty="dotted")
}
dev.off()
