#import data
mydat <- read.table("Final_Data_1.csv",sep=",",header = TRUE)

#preprocess the data so that the row-name becomes probeset id
col_name <- mydat[,1]
mydat <- mydat[,-1]
rownames(mydat) <-col_name
subtype_matrix <- read.csv("proj_metadata.csv")

#4.1 select genes which over 20% of total samples have expression > log2(15)
boo_data <- mydat>log2(15)

boo_sum <- rowSums (boo_data, na.rm = FALSE, dims = 1)

pre1_data <- mydat[boo_sum>(0.2*ncol(mydat)), ]

#4.2 chi-square test, select genes that are significantly different from the median 
# of total genes

var_by_row <- apply(pre1_data, MARGIN = 1, var) # calculate variance for each probeset

#mean_by_row <- apply(pre1_data, MARGIN = 1, mean)

median_var <- median(var_by_row) # calculate the median for variance by row

degree_free <- ncol(mydat)-1 # degree of freeedom 

qt_.995 <- qchisq(0.99, df=degree_free) # find the 0.99 quartile, when df is 133

#qt_.005 <- qchisq(0.01, df=degree_free)

# http://www.itl.nist.gov/div898/handbook/eda/section3/eda358.htm
test_stat <- degree_free*(var_by_row/median_var) # formula that gives test stats for chi-square test
#4.2
pre2_data <- pre1_data[test_stat>qt_.995,]
write.csv(pre2_data,"processed_4.2")

#4.3

# select out genes that have coeffient of variance >0.186
pre2_var_row <- apply(pre2_data, MARGIN = 1, var)

pre2_mean_row <- apply(pre2_data, MARGIN = 1, mean)

#formula of coeffient of variance
coef_var <- sqrt(pre2_var_row)/pre2_mean_row

pre3_data<- pre2_data[coef_var>0.186,]
write.csv(pre3_data,"processed_4.4")

# Section 5

#use hclust to cluster matrix by *samples*
clusters <- hclust(as.dist(1-cor(pre3_data, method="spearman")), method="complete") 

#cut the hcluster into two groups
clst_2g <- cutree(clusters,k=2)
table(clst_2g)

#import subtype matrix that categarize either C3 or C4 to all samples
subtype_vec <- subtype_matrix$cit.coloncancermolecularsubtype

#create a new vector that has red or blue in it, the vector will be used to color samples
#in heatmap
library(plyr)
subtype_vec_color <- revalue(subtype_vec,c("C3"="red", "C4"="blue"))
#the above vector is a factor, which does not work in heatmap function, so 
#I convert it to char string
color <- as.character(subtype_vec_color)
#Again, imported data does not have the correct format, so I convert it to matrix 
data_mat <- as.matrix(pre3_data)
hm <- heatmap(data_mat, ColSideColors = color)

#subset the original table to two groups, based on the boolean vector I created before
clst_group1 <- pre3_data[,clst_2g==1]
clst_group2 <- pre3_data[,clst_2g==2]

#perform t test for fun, not related to project
ttest <- t.test(clst_group1, clst_group2, paired=FALSE)
print(ttest)

#now perform ttest on each row of the table, 
# in each row, two groups'(c3 and c4) are compared, using welch t test
p_value_all_genes <- apply(pre3_data, 1, 
                           function(x){ t.test(x[clst_2g==1], 
                                               x[clst_2g==2]) $p.value } )

#repeat the same procedure to get t stats 
tstats_all_genes <- apply(pre3_data, 1, function(x){ t.test(x[clst_2g==1], x[clst_2g==2]) $statistic } )
#adjust p-value
p_adj_all_genes <- p.adjust(p_value_all_genes,method = "fdr")
print(p_adj_all_genes)

#select p-value <0.05
de_genes <- p_adj_all_genes[p_adj_all_genes <0.05]
print(de_genes)

#organize my results and output a table
all_filter_de_genes <- cbind(tstats_all_genes,p_value_all_genes,p_adj_all_genes)
all_filter_de_genes <- all_filter_de_genes[p_adj_all_genes <0.05,]
colnames(all_filter_de_genes) <- c("ttest","p-value","adjusted-value")
order_de_genes <- all_filter_de_genes[order(all_filter_de_genes[,3]),]
write.csv(order_de_genes, "de_genes_4.4_p0.05_ordered_adj_p")


#4.2 de gene 
#the code down below is eaxctly the same as the above section, 
#except naming differently with "part_" as preffix
part_clusters <- hclust(as.dist(1-cor(pre2_data, method="spearman")), method="complete") 
part_clst_2g <- cutree(part_clusters,k=2)
table(part_clst_2g)

subtype_vec <- subtype_matrix$cit.coloncancermolecularsubtype

library(plyr)
subtype_vec_color <- revalue(subtype_vec,c("C3"="red", "C4"="blue"))
color <- as.character(subtype_vec_color)
part_data_mat <- as.matrix(pre2_data)
part_hm <- heatmap(part_data_mat, ColSideColors = color)

part_clst_group1 <- pre2_data[,part_clst_2g==1]
part_clst_group2 <- pre2_data[,part_clst_2g==2]

part_ttest <- t.test(part_clst_group1, part_clst_group2, paired=FALSE)

part_p_value_all_genes <- apply(pre2_data, 1, 
                           function(x){ t.test(x[part_clst_2g==1], 
                                               x[part_clst_2g==2]) $p.value } )
part_tstats_all_genes <- apply(pre2_data, 1, function(x){ t.test(x[part_clst_2g==1], x[part_clst_2g==2]) $statistic } )
part_p_adj_all_genes <- p.adjust(part_p_value_all_genes,method = "fdr")
print(p_adj_all_genes)
part_de_genes <- part_p_adj_all_genes[part_p_adj_all_genes <0.05]
print(part_de_genes)
part_all_filter_de_genes <- cbind(part_tstats_all_genes,part_p_value_all_genes,part_p_adj_all_genes)
part_all_filter_de_genes <- part_all_filter_de_genes[part_p_adj_all_genes <0.05,]
colnames(part_all_filter_de_genes) <- c("ttest","p-value","adjusted-value")
part_order_de_genes <- part_all_filter_de_genes[order(part_all_filter_de_genes[,3]),]
write.csv(part_order_de_genes, "de_genes_4.2_p0.05_ordered_adj_p")
