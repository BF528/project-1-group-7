##Biologist Role

#Install and load Packages
library(BiocManager)
BiocManager::install(c("affy", "affyPLM", "sva", "AnnotationDbi", "hgu133plus2.db", "GSEABase"))
library(affy)
library(affyPLM)
library(sva)
library(AnnotationDbi)
library('hgu133plus2.db')
library(tidyverse)
library(GSEABase)

#Establish Keys
x <- hgu133plus2.db
keys(x)
columns(x)
keytypes(x)

#Map Probe IDs to Gene Symbols
cols <- c("SYMBOL", "GENENAME")
Sym_ID2 <- AnnotationDbi::select(x, keys=as.character(rownames(project_data2)), columns='SYMBOL')

#Download data and add new column
project_data2 <- read.csv('de_genes_4.4_p0.05_ordered_adj_p.csv', header=TRUE, row.names = 1)
ind <- match(rownames(project_data2), Sym_ID$PROBEID)
project_data2$Gene_Symbol <- Sym_ID$SYMBOL[ind]

#Find top up and down differentially expressed genes
project_data2_bot <- project_data2[order(project_data2$ttest),][1:1000,]
project_data2_top <- project_data2[order(-project_data2$ttest),][1:1000,]
ten_up <- project_data2_top[1:10,]
ten_down <- project_data2_bot[1:10,]
write.csv(ten_up, file='up_differentially_expressed.csv')
write.csv(ten_down, file = 'down_differentially_expressed.csv')



#Forming Geneset Collections
project_data3 <- read.csv('genes_all_4.2_notfiltered.csv', header=TRUE, row.names = 1)
colnames(project_data3) <- c('ttest', 'p.value', 'adjusted.value')
project_data3$Gene_Symbol <- Sym_ID$SYMBOL[match(rownames(project_data3), Sym_ID$PROBEID)]
project_data3 <- project_data3[order(project_data3$Gene_Symbol, project_data3$ttest),]
project_data3 <- project_data3[!duplicated(project_data3$Gene_Symbol),] #select the Probeset ID for unique gene symbol where the differential expression is largest
go <- getGmt('c5.all.v7.0.symbols.gmt') #9996 gene sets in collection
kegg <- getGmt('c2.cp.kegg.v7.0.symbols.gmt') #186 gene sets in collection
hall <- getGmt('h.all.v7.0.symbols.gmt') #50 gene sets in collection

#Fisher Test GO
estimates_go <- c()
p_values_go <- c()
names_go <- c()
for (i in 1:length(go)){
  name <- go[[i]]@setName
  ind_f <- match(go[[i]]@geneIds, project_data3$Gene_Symbol) #get indices for matching genes
  ind_f <- ind_f[!is.na(ind_f)] #remove nas
  in_dif <- length(which(project_data3$adjusted.value[ind_f] < 0.05)) #genes that are in gene set and differentially expressed 
  in_not <- length(project_data3$Gene_Symbol[ind_f]) - in_dif #genes that are in gene set but not differentially expressed
  prac <- project_data3[-(ind_f),]
  not_dif <- length(which(prac$adjusted.value < 0.05)) #genes that are not in gene set but are differentially expressed
  not_not <- length(prac$adjusted.value) - not_dif #gene that are not in gene set and not differentially expressed
  temp <- fisher.test(matrix(c(in_dif, in_not, not_dif, not_not), nrow = 2)) #create contigency table for gene set
  e <- temp[["estimate"]][["odds ratio"]]
  p <- temp$p.value #values from fisher test
  estimates_go[i] <- e
  p_values_go[i] <- p #place in vectors to concatentate later
  names_go[i] <- name
}
matrix_go <- data.frame(matrix(c(names_go, estimates_go, p_values_go), ncol = 3)) #make table of values
colnames(matrix_go) <- c('Pathway Name', 'Statistic Estimate', 'p_value')
matrix_go$padj <- p.adjust(matrix_go$p_value, method='fdr') #calculate padj

#Fisher Test Kegg (same code as above)
estimates_kegg <- c()
p_values_kegg <- c()
names_kegg <- c()
for (i in 1:length(kegg)){
  name <- kegg[[i]]@setName
  ind_f <- match(kegg[[i]]@geneIds, project_data3$Gene_Symbol)
  ind_f <- ind_f[!is.na(ind_f)]
  in_dif <- length(which(project_data3$adjusted.value[ind_f] < 0.05))
  in_not <- length(project_data3$Gene_Symbol[ind_f]) - in_dif
  prac <- project_data3[-(ind_f),]
  not_dif <- length(which(prac$adjusted.value < 0.05))
  not_not <- length(prac$adjusted.value) - not_dif
  temp <- fisher.test(matrix(c(in_dif, in_not, not_dif, not_not), nrow = 2))
  e <- temp[["estimate"]][["odds ratio"]]
  p <- temp$p.value
  estimates_kegg[i] <- e
  p_values_kegg[i] <- p
  names_kegg[i] <- name
}
matrix_kegg <- data.frame(matrix(c(names_kegg, estimates_kegg, p_values_kegg), ncol = 3))
colnames(matrix_kegg) <- c('Pathway Name', 'Statistic Estimate', 'p_value')
matrix_kegg$padj <- p.adjust(matrix_kegg$p_value, method='fdr')

#Fisher Test Hall (same code as above)
estimates_hall <- c()
p_values_hall <- c()
names_hall <- c()
for (i in 1:length(hall)){
  name <- hall[[i]]@setName
  ind_f <- match(hall[[i]]@geneIds, project_data3$Gene_Symbol)
  ind_f <- ind_f[!is.na(ind_f)]
  in_dif <- length(which(project_data3$adjusted.value[ind_f] < 0.05))
  in_not <- length(project_data3$Gene_Symbol[ind_f]) - in_dif
  prac <- project_data3[-(ind_f),]
  not_dif <- length(which(prac$adjusted.value < 0.05))
  not_not <- length(prac$adjusted.value) - not_dif
  temp <- fisher.test(matrix(c(in_dif, in_not, not_dif, not_not), nrow = 2))
  e <- temp[["estimate"]][["odds ratio"]]
  p <- temp$p.value
  estimates_hall[i] <- e
  p_values_hall[i] <- p
  names_hall[i] <- name
}
matrix_hall <- data.frame(matrix(c(names_hall, estimates_hall, p_values_hall), ncol = 3))
colnames(matrix_hall) <- c('Pathway Name', 'Statistic Estimate', 'p_value')
matrix_hall$padj <- p.adjust(matrix_hall$p_value, method='fdr')
#there are no gene sets significantly enriched according to an adjusted p-value cutoff of 0.05

#Prepping tables for paper
go_table <- matrix_go[order(matrix_go$p_value),][1:3,]
kegg_table <- matrix_kegg[order(matrix_kegg$p_value),][1:3,]
hall_table <- matrix_hall[order(matrix_hall$p_value),][1:3,]
write.csv(go_table, file="go_enriched.csv")
write.csv(hall_table, file="hall_enriched.csv")
write.csv(kegg_table, file="kegg_enriched.csv")