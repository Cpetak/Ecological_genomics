data <- read.delim("XCV.thetas", header = TRUE, sep = "\t" )


all_pstW <- c(0)
for (i in seq(length(data$tW))){
  all_pstW[i] <- data$tW[i]/data$nSites[i]
}

all_pstP <- c(0)
for (i in seq(length(data$tP))){
  all_pstP[i] <- data$tP[i]/data$nSites[i]
}

hist(all_pstW, main="Histogram of per-site tW")
hist(all_pstP, main="Histogram of per-site tP")
hist(data$Tajima, main="Tajima's D")


