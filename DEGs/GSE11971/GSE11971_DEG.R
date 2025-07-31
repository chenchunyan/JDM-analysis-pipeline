setwd("D:/BaiduSyncdisk/GSE11971")
###############################1.Read the CEL data in the folder########################################
library(affy)
Data1 <- ReadAffy()
pData(Data1)

#############################2.Change the sample name####################################################
library(stringr)
raw.names<-sampleNames(Data1)
new<-str_split_fixed(raw.names, ".", 1)
new.names<-new[,1]
sampleNames(Data1)<-new.names
pData(Data1)

##############################3.The GCRMA algorithm can also be used for preprocessing.######################################### The same data can be used to obtain this method

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("gcrma") 
library(gcrma)
tum.gcrma<-gcrma(Data1) 
##########Repeat the above analysis and draw the graphs for quality control####
boxplot(tum.gcrma,col=cols) 
#Create a density curve graph
hist(tum.gcrma,col=cols)
####5.Filter probe##############################
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("genefilter")
library(genefilter)
gcrma_filter<-nsFilter(tum.gcrma,require.entrez = FALSE,remove.dupEntrez = FALSE)
gcrma_filter$filter.log
$numLowVar
[1] 27307

$feature.exclude

##############################4.Extract the expression matrix################################################
#tum_exp <- exprs(tum.rma)
gcrma_exp <- exprs(gcrma_filter$eset) 
pData(gcrma_filter$eset)
# mean-variance trend 
probes_expr <- na.omit(gcrma_exp)#eliminate rows with NAs
plotSA(lmFit(probes_expr), main="Mean variance trend, GSE11971")
##############################5.Comment on the probe of the matrix########################################
library(tinyarray)
library(GEOquery)
GPL96 <- getGEO("GPL96",destdir = ".")
ids <- GPL96@dataTable@table
length(unique(ids$`Gene Symbol`)) 
table(rownames(probes_expr) %in% ids$ID) # table(rownames(probes_expr) %in% ids$probe_id)

data<-probes_expr[rownames(probes_expr) %in% ids$ID,]
dim(data)
colnames(data) <- gsub("\\.CEL\\.gz", "", colnames(data)) 
colnames(data)
ids<-ids[match(rownames(data),ids$ID),]
tmp = by(data,ids$`Gene Symbol`,function(x) rownames(x)[which.max(rowMeans(x))] )
probes<-as.character(tmp)
data<-data[rownames(data) %in% probes,]
ids<-ids[ids$ID %in% probes,]
row.names(data)<-ids$`Gene Symbol`
max(data)

pd <- read.table("pd.txt",sep="\t",header=T,row.names = 1)
id <- pd
library(edgeR)
library(ggplot2)
library(dplyr)
library(cluster)
library(factoextra)
#Remove abnormal samples
dist_data <- dist(data[, -ncol(data)], method = "euclidean")
hclust_data <- hclust(dist_data, method = "ward.D2")
clusters <- cutree(hclust_data, k = 3) 
cluster_heights <- aggregate(data, by = list(clusters), FUN = mean)
threshold <- 1.5 * IQR(cluster_heights$mean_height)
outliers <- cluster_heights$clusters[cluster_heights$mean_height > threshold]
write.table(data,file="GSE11971_express.txt",sep="\t",quote = F,col.names = T)
data_1 <- t(data)
write.table(data_1,file="GSE11971_express_2.txt",sep="\t",quote = F,col.names = T)
merged_data <- merge(data_1, pd, by.x = "row.names", by.y ="row.names", all.x = TRUE)
#colnames(pd)[which(names(pd) == "source_name_ch1")] <- "tissue"
subset_pd <- subset(pd, group %in% c("control", "JDM"))
#Neutrophil_pd<-subset_pd[subset_pd$tissue=="neutrophil",]
#PBMC_pd<-subset_pd[subset_pd$tissue=="PBMC",]
#table(Neutrophil_pd$group)
#table(PBMC_pd$group)
#Neutrophil_data <- data[,rownames(Neutrophil_pd)]
subset_data <- data[,rownames(subset_pd)]
#merged_data <- merge(t(Neutrophil_data), Neutrophil_pd, by.x = "row.names", by.y = "geo_accession", all.x = TRUE)
merged_data$group <- ifelse(merged_data$group == "control", 1, 2)
write.table(merged_data,file="GSE11971_express_2.txt",sep="\t",quote = F,col.names = T)

########6.Perform differential expression analysis with Limma################
filtered_pd <- pd
filtered_data <- subset_data
group <- as.factor(filtered_pd$group)
samplename<-rownames(filtered_pd)
filtered_data <- data[, samplename]
design=model.matrix(~factor(group)) 
fit=lmFit(filtered_data,design)
fit=eBayes(fit)
DEG=topTable(fit,coef=2,n=Inf)
View(DEG)
#####火山图
df=DEG
attach(df)
df$v= -log10(q.Value)
df$g=ifelse(df$q.Value>0.01,'No_sig',
            ifelse( df$logFC >0,'up',
                    ifelse( df$logFC < 0,'down','No_sig') )
) 
table(df$g)
#Indicate the type of gene, only for protein-coding genes
genetype <- read.table("./human_geneID_symbol.txt",header=TRUE,sep="\t")
df$name=rownames(df)
df$type <- genetype[match(df$name, genetype$Symbol), "Type"]
library(dplyr)
newdata <- df %>%
  filter(type == "protein-coding")
newdf <- newdata[,1:8]
table(newdf$g)
write.table(newdf,file="GSE11971_deg.txt",sep="\t",quote = F,col.names = T)

library(pheatmap)
n=t(scale(t(filtered_data[cg,])))
n[n>2]=2
n[n< -2]= -2
n[1:4,1:4]
ac=data.frame(groupList=group)
rownames(ac)=colnames(n)  
pheatmap(n,show_colnames =F,show_rownames = F,
         annotation_col=ac)
