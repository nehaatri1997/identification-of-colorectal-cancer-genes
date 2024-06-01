#install.packages("BiocManager")
#install.packages("forcats")
#install.packages("stringr")
#install.packages("ggplot2")
#install.packages("ggrepel")
#install.packages("readr")
#install.packages("tidyr")
#install.packages("survminer")
#BiocManager::install("GEOquery")
#BiocManager::install("limma")
#BiocManager::install("pheatmap")
#BiocManager::install("org.Hs.eg.db")


library(GEOquery)
library(limma)
install.packages("umap")

library(umap)

# load series and platform data from GEO
gse<- getGEO("GSE110223")

## check how many platforms used
length(gse)

gse <- gse[[1]]
gse

pData(gse) ## print the sample information
fData(gse) ## print the gene annotation
exprs(gse) ## print the expression data

## exprs get the expression levels as a data frame and get the distribution
summary(exprs(gse))

exprs(gse) <- log2(exprs(gse))
boxplot(exprs(gse),outline=FALSE)

library(dplyr)

sampleInfo <- pData(gse)

sampleInfo

## source_name_ch1 and characteristics_ch1.1 seem to contain factors we might need for the analysis. Let's pick just those columns

#sampleInfo <- select(sampleInfo, characteristics_ch1,characteristics_ch1.1)

## Optionally, rename to more convenient column names
sampleInfo <- rename(sampleInfo,group = characteristics_ch1, patient=characteristics_ch1.1)

library(pheatmap)
## argument use="c" stops an error if there are any missing data points

corMatrix <- cor(exprs(gse),use="c")
pheatmap(corMatrix)


## Print the rownames of the sample information and check it matches the correlation matrix
rownames(sampleInfo)

colnames(corMatrix)


## If not, force the rownames to match the columns

rownames(sampleInfo) <- colnames(corMatrix)
pheatmap(corMatrix,
         annotation_col=sampleInfo)

library(ggplot2)
## MAKE SURE TO TRANSPOSE THE EXPRESSION MATRIX

pca <- prcomp(t(exprs(gse)))


## Join the PCs to the sample information
library(ggrepel)
cbind(sampleInfo, pca$x) %>%
  ggplot(aes(x = PC1, y=PC2, col=group,label=paste("patient", patient))) + geom_point() + geom_text_repel()



library(readr)
full_output <- cbind(fData(gse),exprs(gse))
write_csv(full_output, file ="gse_full_output.csv")

# Print example datA))
features <- fData(gse)
View(features)
str(features)
head(features)
features_clean_columns<- which(complete.cases(features))
features_clean_columns



features_clean<-na.omit(features)
View(features_clean)


# change column name for Gene Symbol



colnames(features_clean)[colnames(features_clean)== "Gene Symbol"]<-"Gene.Symbol"
View(features_clean)



#clean the feature file



#step 1: remove space on Column
#step 2: change GeneID to new column and run the code
#step 3: change GeneSymboL to Gene.Symbol ( how to rname column name)



### Look at the features data frame and decide the names of the columns you want to keep
features <- select(features_clean, ENTREZ_GENE_ID,Gene.Symbol)
write_csv(full_output, path="gse_full_output.csv")

library(limma)
design <- model.matrix(~0+sampleInfo$group)
design

## the column names are a bit ugly, so we will rename
colnames(design) <- c("Normal","Tumour")


summary(exprs(gse))

## calculate median expression level
cutoff <- median(exprs(gse))

## TRUE or FALSE for whether each gene is "expressed" in each sample
is_expressed <- exprs(gse) > cutoff

## Identify genes expressed in more than 2 samples

keep <- rowSums(is_expressed) > 1

## check how many genes are removed / retained.
table(keep)

## subset to just those expressed genes
gse <- gse[keep,]

fit <- lmFit(exprs(gse), design)
head(fit$coefficients)


contrasts <- makeContrasts(Tumour - Normal, levels=design)

## can define multiple contrasts
## e.g. makeContrasts(Group1 - Group2, Group2 - Group3,....levels=design)

fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)

topTable(fit2, coef=1)
### to see the results of the second contrast (if it exists)
## topTable(fit2, coef=2)

decideTests(fit2)
table(decideTests(fit2))



## calculate relative array weights
aw <- arrayWeights(exprs(gse),design)
aw


fit <- lmFit(exprs(gse), design,
             weights = aw)
contrasts <- makeContrasts(Tumour - Normal, levels=design)
fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)


anno <- features_clean
View(anno)

anno <- select(anno,Gene.Symbol,ENTREZ_GENE_ID)
fit2$genes <- anno
topTable(fit2)

full_results <- topTable(fit2, number=Inf)
colnames(gse)

## Make sure you have ggplot2 loaded
library(ggplot2)
ggplot(full_results,aes(x = logFC, y=B)) + geom_point()


## change according to your needs
p_cutoff <- 0.05
fc_cutoff <- 2

full_results <- full_results %>%
  mutate(Significant = adj.P.Val < p_cutoff & abs(logFC) > fc_cutoff)

# Scatter plot with colored points for significant results
ggplot(full_results, aes(x = logFC, y = B, col = Significant)) +
  geom_point() +
  scale_color_manual(values = c("gray", "red"), labels = c("Not Significant", "Significant"))

# View the updated full_results data frame
View(full_results)

# Add labels for the top N significant results using ggrepel
library(ggrepel)
topN <- 1000
top_results <- full_results %>% filter(Significant) %>% top_n(topN, wt = -adj.P.Val)
ggplot(full_results, aes(x = logFC, y = B, col = Significant)) +
  geom_point() +
  geom_label_repel(data = top_results, aes(label = Gene.Symbol))

full_results %>%
  mutate(Significant = adj.P.Val < p_cutoff, abs(logFC) > fc_cutoff ) %>%
  mutate(Rank = 1:n(), Label = ifelse(Rank < topN, Gene.Symbol,"")) %>%
  ggplot(aes(x = logFC, y = B, col=Significant,label=Label)) + geom_point() + geom_text_repel(col="black")

#Filtering and exporting the results table

library(readr)
filter(full_results, adj.P.Val < 0.01, abs(logFC) >2) %>%
  write_csv(path="filtered_de_results.csv")
View(full_results)

biomarkers <- filter(full_results, adj.P.Val < p_cutoff, abs(logFC) > fc_cutoff)



if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("STRINGdb")

library(STRINGdb)


# Establish a connection to the STRING database with custom settings
string_db <- STRINGdb$new(species = 9606, score_threshold = 500, version = "11.0")


# Retrieve PPI data for the gene IDs in numericIDs
mapped_data <- string_db$map(full_results, c("Gene.Symbol","logFC","P.Value"), removeUnmappedRows = TRUE)
# Explore the structure and preview of mapped_data
str(mapped_data)
head(mapped_data)

# Check the contents of the mapped_data
print(mapped_data)


# Extract the STRING IDs from the 'mapped_data' dataframe
string_ids <- mapped_data$STRING_id


# Analyze interaction scores
interaction_scores <- string_db$get_interactions(string_ids)
summary(interaction_scores)


# Create a data frame with the interaction scores
interaction_data <- data.frame(Interaction_Scores = interaction_scores)


# Load the required package
library(igraph)

# Create a graph object from the PPI data
ppi_graph <- graph.data.frame(interaction_data, directed = FALSE)

# Compute node degrees
node_degrees <- degree(ppi_graph)

# Sort nodes by degree in descending order
sorted_nodes <- sort(node_degrees, decreasing = TRUE)

# Select the top N hub genes (e.g., top 10)
N <- 10
hub_genes <- names(sorted_nodes)[1:N]

# Print the hub genes
print(hub_genes)

