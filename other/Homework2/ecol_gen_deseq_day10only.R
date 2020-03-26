## Set your working directory
setwd("~/Desktop/Ecological_genomics/other/Homework2")

## Import the libraries that we're likely to need in this session
library(DESeq2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(ggpubr)
library(wesanderson)
library(vsn)  ### First: BiocManager::install("vsn") AND BiocManager::install("hexbin")

## Import the counts matrix
countsTable <- read.table("RS_cds2kb_countsMatrix.txt", header=TRUE, row.names=1)
head(countsTable)
dim(countsTable)
countsTableRound <- round(countsTable) # Need to round because DESeq wants only integers
head(countsTableRound)

## Import the samples description table - links each sample to factors of the experimental design.
# Need the colClasses otherwise imports "day" as numeric which DESeq doesn't like, coula altneratively change to d0, d5, d10
conds <- read.delim("RS_samples.txt", header=TRUE, stringsAsFactors = TRUE, row.names=1, colClasses=c('factor', 'factor', 'factor', 'factor'))
head(conds)
dim(conds)

############ Try with only Day 10 data
day10countstable <- countsTableRound %>% select(contains("10"))
dim(day10countstable)

conds10<- subset(conds, day=="10")
dim(conds10)
head(conds10)

## Let's see how many reads we have from each sample:
colSums(day10countstable)
mean(colSums(day10countstable))
barplot(colSums(day10countstable), las=3, cex.names=0.5,names.arg = substring(colnames(day10countstable),1,13))
abline(h=mean(colSums(day10countstable)), col="blue", lwd =2)

## What is the average number of counts per gene?
summary(rowSums(day10countstable))
mean(rowSums(day10countstable))
median(rowSums(day10countstable)) # median is much lower -> not normally distributed gene experssion, some very low/high expression

## What is the average number of counts per gene per sample?

summary(apply(day10countstable,2,mean))

## Create a DESeq object and define the experimental design here with the tilde

dds <- DESeqDataSetFromMatrix(countData = day10countstable, colData = conds10, design = ~ climate + treatment + climate:treatment)
dim(dds)

# Filter out genes with few reads, if we choose 76 as minimum that is on average 1 read per sample!
# this reduces the number of tests we have to from 66408 to 23887 -> have to correct less later because of the mulitple stat tests
dds <- dds[rowSums(counts(dds)) > 76]


## Run the DESeq model to test for differential gene expression: 1) estimate size factors (per sample), 2) estimate dispersion (per gene), 3) run negative binomial glm

dds <- DESeq(dds)

# List the results you've generated

resultsNames(dds) #factors it compared, based on line where we specified design (DESeqDataSetFromMatrix). here first factor was pop. since we had day0 first, it compared day5 to it and day10 to it. see below
#[1] "Intercept"            "pop_BRU_05_vs_ASC_06" "pop_CAM_02_vs_ASC_06"
#[4] "pop_ESC_01_vs_ASC_06" "pop_JAY_02_vs_ASC_06" "pop_KAN_04_vs_ASC_06"
#[7] "pop_LOL_02_vs_ASC_06" "pop_MMF_13_vs_ASC_06" "pop_NOR_02_vs_ASC_06"
#[10] "pop_XBM_07_vs_ASC_06" "day_10_vs_0"          "day_5_vs_0"          
#[13] "treatment_D_vs_C"     "treatment_H_vs_C"

#if we  have climate instead of pop in the design (DESeqDataSetFromMatrix):
resultsNames(dds)
#[1] "Intercept"        "climate_HD_vs_CW" "day_10_vs_0"      "day_5_vs_0"      
#[5] "treatment_D_vs_C" "treatment_H_vs_C"

-------------------------------------------------------------------------------
  
  # Order and list and summarize results from specific contrasts (i.e. the 6 above)
  # Here you set your adjusted p-value cutoff, can make summary tables of the number of genes differentially expressed (up- or down-regulated) for each contrast
  
res <- results(dds, alpha = 0.05) #automatically gives you last comparision
res <- res[order(res$padj),]
#-ve log change = more highly expressed in control! i.e. downregulated in treated

summary(res) #we get, "treatment H vs C"
#LFC > 0 (up)       : 16, 0.067%
#LFC < 0 (down)     : 3, 0.013% -> 3 were higher in control, down regulated in hot

#to pull out results from different contrasts
res_clim <- results(dds, name ="treatment_D_vs_C", alpha = 0.05)
res_clim <- res_clim[order(res_clim$padj),]

summary(res_clim) #we get, "treatment_D_vs_C"
#LFC > 0 (up)       : 678, 2.8% #up regulated in drought
#LFC < 0 (down)     : 424, 1.8% #up in control

-------------------------------------------------------------------------------
  
  ##### Data visualization #####
# MA plot
plotMA(res_clim, ylim=c(-3,3))
#below 0 = higher expression in control, red=significant

# PCA
vsd <- vst(dds, blind=FALSE)
#you can play with "nsub" operation of this function
data <- plotPCA(vsd, intgroup=c("climate", "treatment"), returnData =TRUE) #can play with ntop, default 500 genes!
percentVar <- round(100*attr(data, "percentVar"))

data$treatment <- factor(data$treatment, levels=c("C","H","D"), labels = c("C","H","D"))
data$climate <- factor(data$climate, levels=c("CW","HD"), labels = c("CW","HD"))


ggplot(data, aes(PC1, PC2, color=climate, shape=treatment)) +
  geom_point(size=4, alpha=0.85) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_minimal()

# Counts of specific top gene! (important validatition that the normalization, model is working)
#look at a few top genes, if they look similar to what you expect based on model significant result, result it is not driven by one spec gene
d <-plotCounts(dds, gene="MA_10425837g0010", intgroup = (c("treatment","climate")), returnData=TRUE)
d

p <-ggplot(d, aes(x=treatment, y=count, color=climate)) + 
  theme_minimal() + theme(text = element_text(size=20), panel.grid.major=element_line(colour="grey"))
p <- p + geom_point(position=position_jitter(w=0.3,h=0), size=3) +
  scale_x_discrete(limits=c("C","H","D"))
p

# Heatmap of top 20 genes sorted by pvalue
library(pheatmap)
topgenes <- head(rownames(res_clim),20) #we ordered results by pvalue alread
mat <- assay(vsd)[topgenes,]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(dds)[,c("treatment","climate")])
pheatmap(mat, annotation_col=df)
