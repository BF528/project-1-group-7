#install.packages("BiocManager")
#BiocManager::install("affy")
library(affy)
path = "/Users/nitsuehkebere/Desktop/BF528/project-1/GSE39582/CEL_files"
Data <-ReadAffy(celfile.path = path)
#BiocManager::install(c("affyPLM","sva","AnnotationDbi","hgu133plus2.db"))

rma_data <- rma(Data) #rma normalization 

#browseVignettes("affyPLM")
library("affyPLM")
#fitPLM converts an AffyBatch into a PLMset by fitting a specified robust linear model to the probe level data
#PLMset = A class representation for Probe level Linear Models fitted to Affymetrix GeneChip probe level data
fit_data <- fitPLM(Data, normalize=TRUE, background=TRUE)

#outputs RLE median and IQR
rle_data <- RLE(fit_data, type = 'stats')

#outputs NUSE median and IQR
nuse_data <- NUSE(fit_data, type = 'stats')

#rle data frame + histogram 
rle_df <- data.frame(rle_data)
RLE_Median <- as.numeric(rle_df[1, ])
rle_hist <- hist(RLE_Median, main = "RLE")

#nuse data frame + histogram
nuse_df <- data.frame(nuse_data)
NUSE_Median <- as.numeric(nuse_df[1,])
nuse_hist <- hist(NUSE_Median, main = "NUSE")

library(readr)
library(sva)

#exprs = Retrieve Expression Data From ESets
#ESet= Class To Contain High-Throughput Assays And Experimental Metadata
rma_table <- exprs(rma_data)
proj_metadata <- read_csv("/Users/nitsuehkebere/Desktop/BF528/project-1/proj_metadata.csv")
batch = proj_metadata$normalizationcombatbatch

#model.matrix creates a design (or model) matrix, e.g., by expanding factors to a set of dummy variables (depending on the contrasts) 
mod = model.matrix(~normalizationcombatmod, data=proj_metadata) 

quality_data <-ComBat(dat = rma_table, batch=batch, mod=mod)
write.csv(quality_data, file = "Final_Data_1")

#getwd()

#scaling the transpose of expression data 
scale_quality_data <- scale(t(quality_data), center = TRUE, scale = TRUE) 
#retransposing 
#t_scale_rma_table <- t(scale_rma_table)

pca <-prcomp(scale_quality_data, scale = FALSE, center = FALSE)
#pca[["rotation"]]

library(ggplot2)


#In the pca analysis x contains the PCs for drawing a graph 
pca[["x"]]

#extracting PC1 and PC2
PC1 <- pca$x[,1]
PC2 <- pca$x[,2]
p <-plot(PC1,PC2)

#checking variation in the original data that each principal component accounts for to see how meaningful the clusters are 
pca.var <- pca$sdev^2 #outputs varience 
pca.var.per <-round(pca.var/sum(pca.var)*100, 1) #outputs percent varience 

#PCA Varience plot
pca.var.per.plot <-barplot(pca.var.per, main= "PCA VARIENCE PLOT", xlab= "Principal Component", ylab= "Percent Variation")

#pca dataframe to use for ggplot
pca.data <-data.frame(Sample=rownames(pca$x), PC1, PC2)

#PC1 vs PC2
#https://digibio.blogspot.com/2017/12/pca-plot-with-rnaseq-counts-using.html
pca_plot <-ggplot(data= pca.data, aes(x=PC1, y=PC2, label=Sample)) +
          geom_point() +
          xlab(paste("PC1-", pca.var.per[1], "%", sep = "")) +
          ylab(paste("PC2-", pca.var.per[2], "%", sep = "")) + ggtitle('PC1 vs PC2') + 
          theme(plot.title = element_text(hjust = 0.5))


                  



