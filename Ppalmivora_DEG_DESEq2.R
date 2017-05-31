# analyse ARI data the final figure with DESeq2
# test for differential expression by use of negative binomial generalized linear models

#0.-----
#configure the environment
source("http://bioconductor.org/biocLite.R")
biocLite()
library(ggplot2)
library(RColorBrewer)
library(DESeq2)
library(genefilter)
library(gplots)
library(pheatmap)

#1.-----
#Preapare input datasets, we need raw counts and sample table only
edata <- read.csv("/home/anna/anna/Labjournal/MarcAri/ARI_LILI_expression_table_and_annotation_tab.csv", header = FALSE, row.names = 1, sep = '\t') #TODO: remove paths!
edata <- edata[-c(8:11)] #remove redundant columns
names(edata) <- c('scaffold_id', 'start', 'stop', 'strand', 'annotation', 'AG_pred', 'ALI_pred', 'len', 'MZ', "NBIn6h","NBIn18h", "NBIn24h",
                  "NBIn30h", "NBIn48h", "NBIn72h", "MPIn1dA","MPIn1dB", "MPIn1dC", "MPIn2dA","MPIn2dB", "MPIn2dC", "MPIn3dA","MPIn3dB", "MPIn3dC",
                  "MPIn4dA","MPIn4dB", "MPIn4dC", "Myc1", "Myc2","Myc3")
lean_expression <- edata[9:30] #raw counts only for all the samples

#count matrix
cts <- lean_expression[c(14:22)] #ARI counts at 3 dai, 4 dai and myc samples only, too few reads in 1 dai and 2 dai samples, so we drop them

#sample table
sample_table <- read.csv("/home/anna/anna/Labjournal/MarcAri/REP/turboARI/data/ARI_LILIsample_table.csv", header = T, row.names = 1) #TODO: remove paths!
coldata <- sample_table[c(14:22),] #ARI sample table, 3-4 dai and myc samples only
#put the variable of interest at the end of the formula, the control level is the first level.
coldata$experiment <-factor(coldata$experiment, levels = c("myc", "in_planta")) 

#2.------
#Create DESeqDataSet object to store the read counts and the intermediate estimated quantities
#design formula is used to estimate the dispersions and to estimate the log2 fold changes of the model.
#We want to measure effect of the in planta vs myc state mainly time differences are npt that crucial

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ experiment)
dds #check what is there

#3-------
#Prefiltering, i.e removing rows with noreads or nearly no reads to reduce the memory size of the dds data object and increase the speed of the transformation and testing functions within DESeq2.
#remove rows that have only 10 reads in total, more strict filtering is applied via independent filtering on the mean of normalised counts within results function
dds <- dds[ rowSums(counts(dds)) > 10, ]
#we have only biological reps, no technical reps to collapse, so we proceed further

#4-------
#Diff expression analysis, the standard differential expression analysis steps are wrapped into a single function, DESeq
dds <- DESeq(dds)

#extract a results table with log2 fold changes, p values and adjusted p values. 
res <- results(dds, alpha = 0.001)
resultsNames(dds)

#order the results table by the smallest adjusted p value:
res_ordered <- res[order(res$padj),]

#summarize some basics using the summary function.
summary(res)

#How many adjusted p-values were less than 0.001?
sum(res$padj < 0.001, na.rm=TRUE)

# Export significant results: |LFC| >= 2 and adjusted p-value <= 10^-3
resSig <- as.data.frame(subset(res_ordered, padj < 0.001 & abs(log2FoldChange) >= 2))
# Merge with annotations
annotation <- edata[c(1:8)]
resSig_annotated <- merge(resSig, annotation, by = 'row.names') # suplementary table
write.table(resSig_annotated, row.names = FALSE, "/home/anna/anna/Labjournal/Manuscripts/Marchantia_ARI_small_paper/particles/SuppTable_DESeq2_ARI_results.csv", quote = F, sep = '\t')

#5.-----
# Plots for the figure

#Extract transformed values
#Regularized log transformation, inherently accounts for differences in sequencing depth and shrinks low counts
rld <- rlog(dds, blind=FALSE)
mat <- assay(rld)
#extract the rows corresponding to our set of significant DEGS
idx <- rownames(res)[which(res$padj < 0.001 & abs(res$log2FoldChange) >= 2)]
#index mat:
my_mat <- mat[idx,]  
my_mat <- my_mat - rowMeans(my_mat)
df <- as.data.frame(colData(rld)[,c('experiment', 'time')])

df
#plot the heatmap
#specify colors for the heatmap
RowSideColors<-colorRampPalette(c("navy", "white", "chartreuse4"))(20)

tiff("/home/anna/anna/Labjournal/Manuscripts/Marchantia_ARI_small_paper/PP_DEG_pheatmap.tiff", height=7, width=4, units='in', res=600)
pheatmap(my_mat, annotation_col=df,
         show_rownames = FALSE,
         cluster_cols = FALSE,
         color = RowSideColors,
         scale = 'none')

dev.off()

###
#6.------
#plot table with curated categories for the up-regulated secretome
library(gridExtra)
curated <- read.table("/home/anna/anna/Labjournal/Manuscripts/Marchantia_ARI_small_paper/particles/curated_up_secretome.csv", header = TRUE, sep = ',', check.names=FALSE)
tiff("/home/anna/anna/Labjournal/Manuscripts/Marchantia_ARI_small_paper/PP_SECRETOME_CURATION.tiff", height=5, width=3, units='in', res=600)
grid.table(curated, rows = NULL)
dev.off()

