

SFS <- scan("HD_allsites.sfs")
sumSFS <- sum(SFS)
pctPoly = 100*(1-(SFS[1]/sumSFS))
plotSFS_HD <- SFS[2:length(SFS)]
barplot(plotSFS_HD, names.arg = seq(1,length(SFS)-1), xlab="Derived allele frequency", ylab="Number of sites", col=rgb(1, 0, 0, 0.5))

SFS <- scan("CW_allsites.sfs")
sumSFS <- sum(SFS)
pctPoly = 100*(1-(SFS[1]/sumSFS))
plotSFS_CW <- SFS[2:length(SFS)]
barplot(plotSFS_CW, col=rgb(0, 0, 1, .5),add=TRUE)

legend('topright', bty = 'n', title = 'Source climate',
       legend = c('HotDry', 'ColdWet'), fill = c("tomato3", "dodgerblue3"))


div <- read.table("HD_allsites.thetas.idx.pestPG")
colnames(div)=c("window","chrome","wincenter","tW","tP","tF","tH","tL","tajD","fulif","fuliD","fayH","zengsE","numSites")
div$tWpersite = div$tW/div$numSites
div$tPpersite = div$tP/div$numSites
pdf("XWS_diversity_stats2.pdf")
par(mfrow=c(2,2))

hist(div$tWpersite,col="darkorchid",xlab="Theta_W",main="")
hist(div$tPpersite,col="darkorchid", xlab="Theta-Pi",main="")
hist(div$tajD,col=rgb(1, 0, 0, 0.5),xlab="Tajima's D",main="")

summary(div)

barplot(plotSFS, main="SFS", xlab= "Derived allele frequency", ylab="Number of sites")
dev.off()

plot(div$tajD)

div <- read.table("CW_allsites.thetas.idx.pestPG")
colnames(div)=c("window","chrome","wincenter","tW","tP","tF","tH","tL","tajD","fulif","fuliD","fayH","zengsE","numSites")
div$tWpersite = div$tW/div$numSites
div$tPpersite = div$tP/div$numSites
pdf("XWS_diversity_stats2.pdf")
par(mfrow=c(2,2))

hist(div$tWpersite,col="darkorchid",xlab="Theta_W",main="")
hist(div$tPpersite,col="darkorchid", xlab="Theta-Pi",main="")
hist(div$tajD,col=rgb(0, 0, 1, .5),xlab="Tajima's D",main="", add=TRUE)
abline(v = 0, col="red", lwd=3, lty=2)
summary(div)

