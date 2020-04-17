library(methylKit)
library(tidyverse)
library(ggplot2)
library(pheatmap)

# first, we want to read in the raw methylation calls with methylkit

# set directory with absolute path (why is this necessary? I have no idea, but gz files wont work with relative paths)
dir <- "/Users/csengepetak/Desktop/Ecological_genomics/epigenetics"

# read in the sample ids
samples <- read.table("~/Desktop/Ecological_genomics/epigenetics/sample_id.txt", header=FALSE)

# now point to coverage files
files <- file.path(dir,samples$V1)
all(file.exists(files)) #if true we have all files there

# convert to list
file.list <- as.list(files)

# get the names only for naming our samples
nmlist <- as.list(gsub("_1_bismark_bt2_pe.bismark.cov.gz","",samples$V1))

# use methRead to read in the coverage files
myobj <- methRead(location= file.list,
                  sample.id =   nmlist,
                  assembly = "atonsa", # this is just a string. no actual database
                  dbtype = "tabix", #so that not everything has to be in memory at the same time
                  context = "CpG",
                  resolution = "base", #we have SNP data
                  mincov = 20, #minimum coverage per base because it is data from pooled copepods
                  treatment = 
                    c(0,0,0,0,
                      1,1,1,1,
                      2,2,2,2,
                      3,3,3,3,
                      4,4,4,4), # eg we have 4 HA_F25 in "samples"
                  pipeline = "bismarkCoverage",
                  dbdir = "/Users/csengepetak/Desktop/Ecological_genomics/epigenetics")

######
# visualize coverage and filter
######

# We can look at the coverage for individual samples with getCoverageStats()
getCoverageStats(myobj[[1]],plot=TRUE) # we don't want to consider the very high coverage ones because that migth be due to eg pcr duplicates

# and can plot all of our samples at once to compare.

# filter samples by depth with filterByCoverage()
filtered.myobj <- filterByCoverage(myobj, lo.count = 20, lo.perc = NULL, hi.count = NULL, hi.perc = 97.5, db.dir = "/Users/csengepetak/Desktop/Ecological_genomics/epigenetics")


######
# merge samples
######

#Note! This takes a while and we're skipping it

# use unite() to merge all the samples. We will require sites to be present in each sample or else will drop it

meth <- unite(filtered.myobj, mc.cores=3, suffix="united") # we are going to analyse only samples that we have all data of

#the output of unite() is loaded in below from a file Reid created, unite() takes hours
meth <- methylKit:::readMethylBaseDB(
  dbpath = "/Users/csengepetak/Desktop/Ecological_genomics/epigenetics/methylBase_united.txt.bgz",
  dbtype = "tabix",
  sample.id =   unlist(nmlist),
  assembly = "atonsa", # this is just a string. no actual database
  context = "CpG",
  resolution = "base",
  treatment = c(0,0,0,0,
                1,1,1,1,
                2,2,2,2,
                3,3,3,3,
                4,4,4,4),
  destrand = FALSE)


# percMethylation() calculates the percent methylation for each site and sample
pm <- percMethylation(meth)

#plot methylation histograms
# might need to run dev.off() first
ggplot(gather(as.data.frame(pm)), aes(value)) + 
         geom_histogram(bins = 10, color="black", fill="grey") + 
         facet_wrap(~key)

# calculate and plot mean methylation
sp.means <- colMeans(pm)

p.df <- data.frame(sample=names(sp.means),
                   group = substr(names(sp.means), 1,6),
                   methylation = sp.means)

ggplot(p.df, aes(x=group, y=methylation, color=group)) + 
  stat_summary(color="black") + geom_jitter(width=0.1, size=3) 
#AA_F00 was also the one that was more mapped to reference -> difficult to tell we see real diff in average methlyation rate per site (SNP) here

# sample clustering
clusterSamples(meth, dist="correlation", method = "ward.D", plot=TRUE)

# subset with reorganize(), splitting data to a lot of pairwise alignments, no better way yet

meth_sub <- reorganize(meth,
                       sample.ids =c("AA_F25_1","AA_F25_2","AA_F25_3", "AA_F25_4",
                                     "HA_F25_1","HA_F25_2","HA_F25_3","HA_F25_4"),
                       treatment = c(0,0,0,0,1,1,1,1),
                       save.db=FALSE)

# calculate differential methylation between two groups specified above

myDiff <- calculateDiffMeth(meth_sub, overdispersion = "MN", mc.cores = 1, suffix = "AA_HH", adjust = "qvalue", test = "Chisq")

# get all differentially methylated bases, the ones that are significant from above
myDiff <- getMethylDiff(myDiff, qvalue = 0.05, difference = 10) #need to be at least 10 differences, kind of arbitrary
#myDiff4 <- getMethylDiff(myDiff2, qvalue = 0.05, difference = 0) #need to be at least 10 differences, kind of arbitrary

# we can visualize the changes in methylation frequencies quickly.
hist(getData(myDiff)$meth.diff) #-> less methylation in HH, more in -ve change

#hist(getData(myDiff3)$meth.diff, breaks = 10, col=rgb(1,0,0,0.5),xlim=c(-40,40), ylim=c(0,45), main="Frequency distribution of SNPs significantly differentially methylated in different treatment groups", xlab="% change in methylation rate")
#hist(getData(myDiff4)$meth.diff, col=rgb(0,0,1,0.5), add=T)
#legend("topright", c("Carbon dioxide response", "Heat response"), col=c("blue", "red"), lwd=10)

# get hyper methylated bases
hyper=getMethylDiff(myDiff,difference=10,qvalue=0.05,type="hyper")
#
# get hypo methylated bases
hypo=getMethylDiff(myDiff,difference=10,qvalue=0.05,type="hypo")

# heatmap
pm <- percMethylation(meth_sub)
# make a dataframe with snp id's, methylation, etc.
sig.in <- as.numeric(row.names(myDiff))
pm.sig <- pm[sig.in,]

# add snp, chr, start, stop

din <- getData(myDiff)[,1:3]
df.out <- cbind(paste(getData(myDiff)$chr, getData(myDiff)$start, sep=":"), din, pm.sig)
colnames(df.out) <- c("snp", colnames(din), colnames(df.out[5:ncol(df.out)]))
df.out <- (cbind(df.out,getData(myDiff)[,5:7]))

my_heatmap <- pheatmap(pm.sig,show_rownames = FALSE)

#normalising
ctrmean <- rowMeans(pm.sig[,1:4])
h.norm<- (pm.sig-ctrmean)
my_heatmap <- pheatmap(h.norm,show_rownames = FALSE)

#####
#let's look at methylation of specific gene or snp
####

df.out
df.plot <- df.out[,c(1,5:12)] %>% pivot_longer(-snp, values_to = "methylation")
df.plot$group <- substr(df.plot$name,1,2)
head(df.plot)

# looking at snp LS051659.1:1214
# if you choose a different snp, you can create different plots.

df.plot %>% filter(snp=="LS051734.1:5459") %>% 
  ggplot(., aes(x=group, y=methylation, color=group, fill=group)) +
  stat_summary(fun.data = "mean_se", size = 2) +
  geom_jitter(width = 0.1, size=3, pch=21, color="black")

## write bed file for intersection with genome annotation

write.table(file = "~/Desktop/diffmeth_AA_HA.bed",
            data.frame(chr= df.out$chr, start = df.out$start, end = df.out$end),
            row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

# to get overlap between SNP datasets
#j <- merge(x = HA_data, y = AH_data, by = c("chr", "start"))
