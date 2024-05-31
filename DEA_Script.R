#This script contains a full workflow where raw RNA-Seq count data is taken, normalized and examined via Differential Gene Expression analysis.

#First we setup our working space for the analysis by installing and loading necessary packages and libraries.

setwd("C:/Users/Vesa/Desktop/RNA-seq_DEA") #This should be altered to match your own working directory.

if (!require("BiocManager", quietly = TRUE)) #Contains essential libraries needed for the analysis.
  install.packages("BiocManager")

BiocManager::install("DESeq2") #Will be used for managing and manipulating RNA count data, ie. filtering normalized counts. 

BiocManager::install("edgeR") #Used for normalizing the counts data with the TMM (trimmed mean of M values)

install.packages("pheatmap") #Used in creating the heatmap for sample clustering.

install.packages("ggplot2") #Used for creating plots based on the differential gene expression.

library(edgeR)

library(DESeq2)

library(pheatmap)

library(ggplot2)

#Start with reading our data into necessary data structures than can then be manipulated efficiently.

#library(edgeR)

Coinfection.targets<-read.delim("./data/fileDesc.txt") #Read and store the information about the location of our data and its contents.

rownames(Coinfection.targets)<-c("Ha1","Ha2","Ha3","Ctr1","Ctr2","Ctr3") #Add descriptive name for each of the samples.

Coinfection.orig <- edgeR::readDGE(Coinfection.targets, header=F) #Create a Large DGEList from our sample data.

dim(Coinfection.orig) #Check the dimensions of the resulting data frame.

head(Coinfection.orig) #Check the first 6 rows of our data.

Coinfection.rawCount <- Coinfection.orig$count #Extract the raw counts from our data stored in the DGEList we just created into its own matrix.

dim(Coinfection.rawCount) #Check the dimensions again...

head(Coinfection.rawCount) #...And the contents.




#Let's then prepare our metadata.
sampletype <- factor(c(rep("Ha",3), rep("Ctr", 3))) #Create a factor describing sample types for metadata, used later in PCA-plotting.

meta <- data.frame(sampletype, row.names = colnames(Coinfection.orig$count)) #Create a data frame using the sampletype factor and row names extracted from 

#Check that the column names in our raw count and meta matrices match.
#colnames(Coinfection.orig$count)

colnames(Coinfection.rawCount)

rownames(meta)

#all(colnames(Coinfection.orig$count) %in% rownames(meta))

all(colnames(Coinfection.rawCount) %in% rownames(meta))




#Move on to creating DESeqDataset object that is then normalized using median of ratios method in DESeq2 library.

dds <- DESeq2::DESeqDataSetFromMatrix(Coinfection.orig, colData = meta, design = ~ sampletype) #Set up the DESeqDataSet-object.

head(DESeq2::counts(dds)) #Check the contents of that object.

dds <- DESeq2::estimateSizeFactors(dds) #Generate the size factors based on our data.

DESeq2::sizeFactors(dds) #Check these generated size factors.

normalized_counts <- DESeq2::counts(dds, normalized=TRUE) #Apply median of ratios normalization on our count data.

write.csv(normalized_counts, file="./results/coinfection_normalized_counts_DESeq2.csv") #Extract the normalized read counts into its own .csv file needed for the submission.




#Now that we have normalized our count data, we can perform sample-level QC on it with a help of PCA-plot.
rld <- DESeq2::rlog(dds, blind=TRUE) #Apply a rlog transformation on our normalized count to improve the clustering of our data.

DESeq2::plotPCA(rld, intgroup="sampletype") #Use the resulting DESeqTransform-object to draw the PCA-plot displaying the clustering in our data.

pdf("./results/PlotPCA_dds.pdf") #Save the plot in a separate PDF file for the submission.
DESeq2::plotPCA(rld, intgroup="sampletype")
dev.off()

# 1.According to the plot, the infected samples and the control samples tend to form their own cluster, indicating that they have similar RNA read counts for different genes. 
#   However, one control sample seems to be an outlier with notably skewed location on the Y-axis. This might be due to plethora of reasons such as technical or experimental errors or other reasons.

# 2.Overall, the PCA plot is in agreement with our experimental design, with control and infected samples clustering as was expected.

# 3.The variance information of the PC1 and PC2 in the plot tells us how much each of these principle components account for the total observed variance in our data. 
#   Together they capture 75% of the total variance. However, out of these two, the PC1 is more important as it is twice as big as PC2 and accounts for half of the total variance in the data by itself.

#Continuing with the sample-level QC for our normalized count data, let's draw a hierarchical clustering heatmap.
rld_mat <- SummarizedExperiment::assay(rld) #First, read our transformed count data into a matrix.

rld_cor <- cor(rld_mat) #Then, compute pairwise correlation values for the samples.

head(rld_cor) #Check the resulting pairwise correlation matrix.

head(meta) #Check meta-object created earlier, as this will be used for annotating the heatmap.

#Display the hierarchical clustering heatmap. 
heat.colors <- RColorBrewer::brewer.pal(6, "Oranges") 
pheatmap::pheatmap(rld_cor, annotation = meta, color = heat.colors, border_color=NA, fontsize = 10, 
                   fontsize_row = 10, height=20)

#Save the resulting heatmap in a separate PDF file for the submission.
pdf("./results/PlotHeatmap_dds.pdf")
heat.colors <- RColorBrewer::brewer.pal(6, "Oranges")
pheatmap(rld_cor, annotation = meta, color = heat.colors, border_color=NA, fontsize = 10, 
         fontsize_row = 10, height=20)




#Moving on to Differential Expression Analysis using edgeR.
library(edgeR)
options(digits=3)

#Start by prepping the DGEList-object in the same manner as earlier with Coinfection.targets.

infection.targets<-read.delim("./data/fileDesc.txt") #Read and store the information about the location of our data and its contents.

#infection.targets #Check.

rownames(infection.targets)<-c("Ha1","Ha2","Ha3","Ctr1","Ctr2","Ctr3") #Add the necessary row names aswell.

infection.targets #And check that the information is all there.

infection <- edgeR::readDGE(infection.targets, header=F) #Create a Large DGEList from our sample data.

dim(infection) #Check its dimensions...

head(infection) #...And first six rows.

infection.rawCount <- infection$count #Extract raw counts from the DGEList in to its own matrix.

head(infection.rawCount) #Check the resulting matrix.

write.csv(infection.rawCount, file="./results/infection.rawCounts.csv") #Save the raw count data for the submission.

#library(ggplot2)
#Next, let's use the raw count data to create a histogram depicting the number of genes and their raw expression counts. 
#Draw the histogram using ggplot2-library.
ggplot(infection.rawCount) +
  geom_histogram(aes(x = Ha1), stat = "bin", bins = 200) +
  xlab("Raw expression counts") +
  ylab("Number of genes")

#Save a picture of the resulting histogram for the submission.
png("./results/count distribution.png", res=300, height=1800, width=1800)
ggplot(infection.rawCount) +
  geom_histogram(aes(x = Ha1), stat = "bin", bins = 200) +
  xlab("Raw expression counts") +
  ylab("Number of genes")
dev.off()




#Let's move on to calculating to CPM for our data before and after filtering and normalization.
#First we will normalize our data and calculate the CPM before filtering using the TMM method implemented by edgeR. 
infection.normCPM <- edgeR::cpm(edgeR::calcNormFactors(infection))

dim(infection.normCPM) #Check the dimensions of the resulting matrix.

head(infection.normCPM) #And first 6 rows, as you can see, there are some rows with only 0.0, which will be removed in the filtering.

write.csv(infection.normCPM, file="./results/infection.normCPM.csv") #Save the unfiltered CPM data in a separate file.

#Let's filter the data.

infection.filtered <- rowSums(edgeR::cpm(infection)>1) >=3 #Remove genes that are not expressed in 1 CPM in at least 3 samples.

table(infection.filtered) #Show how many genes are filtered away.

infection$samples$lib.size #Check the library size before filtering.

Infection <- infection[infection.filtered,] #Create a new DGEList with the filtered data.

colSums(Infection$counts) #Check the library size after filtering.

dim(Infection) #Check the dimensions of the filtered data.

Infection$samples$lib.size <- colSums(Infection$counts) #Update the filtered library size.

Infection$samples #And check it.

Infection = edgeR::calcNormFactors(Infection) #Calculate the normalization factors for the TMM normalization of filtered data.

Infection$samples #Check the normalization factors.

Infection.filtered.normCPM <-edgeR::cpm(edgeR::calcNormFactors(Infection)) #Normalize the filtered data.

write.csv(Infection.filtered.normCPM, file="./results/Infection.filtered.normCPM.csv") #And write it into a separate .csv file.




#Let's move on then to the final part, Differential Gene Expression Analysis.
group<-factor(c('Ha','Ha','Ha',"Ctr","Ctr","Ctr")) #Define a treatment factor.

#Describe the experimental design
Infection.design <- model.matrix(~group)   
rownames(Infection.design)<-colnames(Infection$counts)
Infection.design

#Identify possible outliers from our samples using a MDS plot. 
limma::plotMDS(Infection, main="MDS plot of RNA-Seq", labels=colnames(Infection$counts))

#Save the MDS plot in a separate file for the submission
png("./results/plotMDS.Infection.png", res=300, height=1800, width=1800)
limma::plotMDS(Infection, main="MDS plot of Infection RNA-Seq", labels=colnames(Infection$counts))
dev.off()



#library(edgeR)
#Based on the experimental design, estimate common, trended and tagwise dispersion for the filtered data.
Infection <- estimateGLMCommonDisp(Infection, Infection.design)

Infection <- estimateGLMTrendedDisp(Infection, Infection.design)

Infection <- estimateGLMTagwiseDisp(Infection, Infection.design)

#library(edgeR)
#Plot the means and variances of our data.
plotMeanVar(Infection, show.tagwise.vars=T,NBline=T)

#Plot the Common, Trended and Tagwise dispersions for the filtered data.
plotBCV(Infection)

#Fit Generalized Linear Models to our data.
Infection.fit <- glmFit(Infection, Infection.design)
colnames(Infection.fit)

#Conduct a Likelihood Ratio Test to assess the effects of treatment vs control.
lrt.Ha_vs_Ctr <- glmLRT(Infection.fit, coef=2)

#Get a table containing top differentially expressed genes.
t1<-topTags(lrt.Ha_vs_Ctr, n=nrow(Infection))
head(t1$table)

#See how many of the genes are up and down regulated.
summary(decideTests(lrt.Ha_vs_Ctr, adjust.method="BH", p.value=0.05))

nrow(subset(topTags(lrt.Ha_vs_Ctr, n=586)$table,  logFC > 0))

nrow(subset(topTags(lrt.Ha_vs_Ctr, n=586)$table,  logFC < 0))

#Output their own tables for up and downregulated genes.
lrt.Ha_vs_Ctr_UP <- subset(topTags(lrt.Ha_vs_Ctr, n=586)$table, logFC > 0)
lrt.Ha_vs_Ctr_DW <- subset(topTags(lrt.Ha_vs_Ctr, n=586)$table, logFC < 0)

#And write them into a separate .csv file.
write.csv(lrt.Ha_vs_Ctr_UP, file="./results/lrt.Ha_vs_Ctr_UP.csv")
write.csv(lrt.Ha_vs_Ctr_DW, file="./results/lrt.Ha_vs_Ctr_DW.csv")

#Specify DEGs for plotting.
DEtags.lrt.Ha_vs_Ctr <- rownames(Infection)[as.logical(decideTests(lrt.Ha_vs_Ctr, adjust.method="BH", p.value=0.05))]

#For plotting DEGs, label all genes grey.
Infection.colHavsCtr = rep('grey55', nrow(Infection))

#And then label upregulated genes red and downregulated blue.
Infection.colHavsCtr[lrt.Ha_vs_Ctr$table$PValue < 0.05 & lrt.Ha_vs_Ctr$table$logFC >0 ] <- "red"
Infection.colHavsCtr[lrt.Ha_vs_Ctr$table$PValue < 0.05 & lrt.Ha_vs_Ctr$table$logFC <0 ] <- "blue"

#Draw a plot visualizing DEGs, changes in their expression levels and log-counts per million.
par(omi=c(0.1,0.1,0.1,0.1), las=1, cex=0.5, mgp=c(3,1,0), cex.main=1.8, cex.lab=1.4, cex.axis=1.4)
plotSmear(lrt.Ha_vs_Ctr, xlab="log-counts per million (logCPM)", ylab="log2-fold change (log2FC)", main="Ha infection compared to Control", smearWidth=0.5, pch=21, cex=0.4, deCol="red", col=Infection.colHavsCtr, ylim=c(-7,7), yaxs="i")
abline(h=c(-1,1),col="dodgerblue")

#Save the said plot for submission.
png("./results/plotSmear.InfectionRNAseq.png", res=300, height=1800, width=1800)
par(omi=c(0.1,0.1,0.1,0.1), las=1, cex=0.5, mgp=c(3,1,0), cex.main=1.8, cex.lab=1.4, cex.axis=1.4)
plotSmear(lrt.Ha_vs_Ctr, xlab="log-counts per million (logCPM)", ylab="log2-fold change (log2FC)", main="Ha infection compared to Control", smearWidth=0.5, pch=21, cex=0.4, deCol="red", col=Infection.colHavsCtr, ylim=c(-7,7), yaxs="i")
abline(h=c(-1,1),col="dodgerblue")
dev.off()



