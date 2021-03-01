library(GEOquery)
library(oligo)
library(oligoClasses)
library(tidyverse)
library(limma)
library(AnnotationDbi)
library(hugene20sttranscriptcluster.db)
library(RColorBrewer)

######STEP 1: IMPORTING GSE DATA FROM GEO########################################
#Retrieve GSE from GEO, and downloaded cel data from local folder.
gse <- getGEO('GSE48973')
gse48973 <- gse[[1]]

######STEP 2: PROCESSING DATA########################################
#Ensure order of cel and gse data is aligned.
phenodata <- pData(gse48973)
phenodata[,'cel_file'] <- str_split(phenodata$supplementary_file,"/") %>% map_chr(tail,1)
cel_path <- '/Users/a.h./Documents/University/Comp Bio/LSM3241/Assignments/GSE48973_RAW'
cel_data <- read.celfiles(paste0(cel_path,'/',phenodata$cel_file),phenoData=phenoData(gse48973))

#Overview of the experimental design. 
data.overview <- pData(cel_data)[,c("geo_accession","cell type:ch1","telomere length:ch1")]

######STEP 3: DIFFERENTIAL GENE EXPRESSION########################################
#Perform RMA: A 3 step process consisting of 1) Background correction 2) Normalisation 3) Summarisation.
gse48973_eset <- rma(cel_data)

#Generate a model matrix that specifies the experimental design. 
#Note that the resultant model matrix takes longer telomere length (15 kb) as reference for differential analysis.
expt_design <- model.matrix(~ gse48973_eset[["telomere length:ch1"]])
colnames(expt_design)[2] <- "Telomere Length"

#Fit a linear model for ever row of the expression matrix. 
#For a two-class comparison problem (short telomere length VS long telomere length), it is equivalent to performing a t-test.
model_fit <- lmFit(gse48973_eset,expt_design)

#Carry out differential analysis by performing Emprirical Bayes on the fitted model.
ebayes <- eBayes(model_fit)
ebayes_table <- topTable(ebayes, n = Inf)

######STEP 4: ANALYSIS OF RESULTS########################################
#We are only interested in genes that show significant differential gene expression.
#A threshold of p-adj val <= 0.05 is used to define signficant differential gene expression.
ebayes_sig <- ebayes_table[ebayes_table$adj.P.Val <= 0.05,]

#Map the probe IDs to genes using Affymetrix Human Gene 2.0 ST Array. 
#Probes that map to 'NA' are outdated IDs and thus removed.
probeIDs <- rownames(ebayes_sig)
genes <- mapIds(hugene20sttranscriptcluster.db, probeIDs, "SYMBOL", keytype = "PROBEID")
probe2genes <- as.data.frame(cbind(probeIDs, genes), row.names = FALSE)
probe2genes <- na.omit(probe2genes)
annotated_results <- ebayes_sig[intersect(rownames(ebayes_sig),probe2genes$probeIDs),]
annotated_results <- as.data.frame(cbind(annotated_results, probe2genes$genes))
colnames(annotated_results)[7] <- "Gene"

#Plot heatmap to show that differential gene expression cluster according to the length of telomere.
#To tidy up the heatmap visually, correlation instead of the default euclidean distance is used for clustering.
sig_eset <- gse48973_eset[rownames(annotated_results),]
heatmap(exprs(sig_eset), 
        labCol = sig_eset[["telomere length:ch1"]] ,labRow=NA,
        col = rev(brewer.pal(10, "RdBu")),
        distfun = function(x) as.dist(1-cor(t(x))))

#To find more meaningful results, we further subset the signficant genes to those that also display substantial levels of differential gene expression.
#We define substantial differential gene analysis as having an abs(logFC) >= 1.
#Genes that display both signficant and substantial differential gene expression is highlighted in red in a volcano plot.
interesting_genes <- annotated_results[abs(annotated_results$logFC) >= 1,]
volcanoplot(ebayes, coef=2)
points(interesting_genes[['logFC']],-log10(interesting_genes[['P.Value']]),col='red')
