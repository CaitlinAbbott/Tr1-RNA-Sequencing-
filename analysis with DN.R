#install packages
install.packages("BiocManager")
install.packages("devtools")
install.packages("tidyverse")
install.packages("magrittr")
BiocManager::install("limma")
BiocManager::install("edgeR")
BiocManager::install("ggrepel")
BiocManager::install("AnnotationHub")
BiocManager::install("ensembldb")
BiocManager::install("RColorBrewer")
install.packages("ggplot2")
install.packages("pheatmap")
install.packages("plotly")
install.packages("RColorBrewer")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Mm.eg.db")
BiocManager::install("variancePartition")
BiocManager::install("KEGG.db")
install.packages("statmod")


setwd("C:/Users/caitlin/Desktop/Tr1 RNA seq")
library(tidyverse)
library(magrittr)
library(limma)
library(edgeR)
library(ggrepel)
library(ggplot2)
library(pheatmap)
library(AnnotationHub)
library(ensembldb)
library(plotly)
library(RColorBrewer)
library(clusterProfiler)
library(org.Mm.eg.db)
library(variancePartition)
library(KEGG.db)
library(fgsea)
library(data.table)

file <- "Caitlin_genes.out"
counts <- read.delim(file, comment = "#")
head(counts)

#Converting counts data into a readable list called dgelist
dgeList <- counts %>%
  dplyr::select(Geneid, ends_with("bam")) %>%
  as.data.frame() %>%
  column_to_rownames("Geneid") %>%
  as.matrix() %>%
  DGEList()
dgeList$counts

dgeList
#change sample names
colnames(dgeList) <- basename(gsub("X.data.biohub.20190819_McColl_RNASeq.Caitlin.2_alignedData.bams.", "", colnames(dgeList)))
colnames(dgeList) <- basename(gsub("Aligned.sortedByCoord.out.bam", "", colnames(dgeList)))

dgeList$counts[,c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)] -> dgeList$counts
dgeList$samples
dgeList$samples <- dgeList$samples[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16),]


#Normfactors based on library size condition is from sample names
dgeList <- calcNormFactors(dgeList, method = c("TMM"))
dgeList$samples$condition <- gsub("M2_|M3_|M5_|M6_", "", rownames(dgeList$samples)) #add column to show populations
dgeList$samples$condition <- factor(dgeList$samples$condition, 
                                    levels= c("DN", "DP", "LAG_3", "CD49b"))
dgeList$samples$mouse <- str_extract(colnames(dgeList), "M[0-9]")
dgeList$samples

#Library size for each sample
dgeList$samples %>%
  rownames_to_column("Sample") %>%
  mutate(CellType = dgeList$samples$condition) %>%
  ggplot(aes(x = Sample, y = lib.size / 1e6, fill = CellType)) +
  geom_bar(stat = "identity") +
  labs(y = "Library Size (millions)") +
  theme_bw()

#MDS plot
mds <- plotMDS(dgeList[,], col = as.integer(dgeList$samples$condition))

#Filter out unnecessary genes
logcount <- cpm(dgeList, log = TRUE)
plotDensities(logcount)
genes2keep <- rowSums(logcount > 1) >4 #>1million reads in at least 3 samples
summary(genes2keep)
plotDensities(logcount[genes2keep, ])
dgeFilt <- dgeList[genes2keep, , keep.lib.sizes = FALSE]
dgeFilt

# PCA and filter out DN samples
pca <- dgeFilt %>%
  # .[,!grepl("DN", colnames(.))] %>%
  cpm(log = TRUE) %>%
  t() %>%
  prcomp()

x <- "PC1"
y <- "PC2"
dgeFilt$samples %>%
  rownames_to_column("Sample") %>%
  #dplyr::filter(!str_detect(Sample, "DN")) %>%
  mutate(Mouse = str_extract(Sample, "M[0-9]")) %>%
  as_tibble() %>%
  #dplyr::filter(condition != "DN") %>%
  cbind(pca$x[.$Sample,]) %>%
  ggplot(aes_string(x, y, shape = "Mouse", colour = "condition")) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = Sample), show.legend =FALSE) + 
  labs(
    x = paste0(x, " (", scales::percent(summary(pca)$importance[2, x]), ")"),
    y = paste0(y, " (", scales::percent(summary(pca)$importance[2, y]), ")")
  )

#Annotating
Ah <- AnnotationHub()
#Find which genome did we use
unique(Ah$dataprovider)
subset(Ah, dataprovider == "Ensembl") %>%  #In Ensembl databse
  subset(species == "Mus musculus") %>%  #under Mouse
  subset(rdataclass == "EnsDb") 
ensDb <- Ah[["AH69210"]] #This is the genome I used for the sequencing
genes <- genes(ensDb) %>% #extract the genes
  subset(seqnames %in% c(1:19, "MT", "X", "Y")) %>%
  keepStandardChromosomes()
seqlevels(genes) <- str_sort(seqlevels(genes), numeric = TRUE) #order it by numeric
genes
dgeFilt$genes <- genes[rownames(dgeFilt),]  
dgeFilt$genes

#vooming
design <- model.matrix(~0 + condition, data = dgeFilt$samples)
colnames(design) <- str_remove(colnames(design), "condition")
v = voomWithQualityWeights(dgeFilt, design = matrix(1, nrow = ncol(dgeFilt)), plot = TRUE)

mds$cmdscale.out %>%
  as.data.frame() %>%
  set_colnames(c("Dim1", "Dim2")) %>%
  rownames_to_column("sample") %>%
  cbind(v$targets[.$sample,]) %>%
  dplyr::rename(w = sample.weights) %>%
  ggplot(aes(Dim1, Dim2, colour = condition, label = sample)) +
  geom_point(aes(size = w)) +
  geom_text_repel() +
  labs(x = "Leading logFC dim 1",
       y = "Leading logFC dim 2",
       size = "Sample\nweight",
       colour = "Cell Type") +
  theme_bw(pca <- v$E) %>%
  #.[,!grepl("DN", colnames(.))] %>%
  cpm(log = TRUE) %>%
  t() %>%
  prcomp()

x <- "PC1"
y <- "PC2"
v$targets %>%
  rownames_to_column("Sample") %>%
  #dplyr::filter(!str_detect(Sample, "DN")) %>%
  mutate(Mouse = str_extract(Sample, "M[0-12]")) %>%
  as_tibble() %>%
  #dplyr::filter(condition != "DN") %>%
  cbind(pca$x[.$Sample,]) %>%
  ggplot(aes_string(x, y, shape = "Mouse", colour = "condition", size = "sample.weights")) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = Sample), show.legend =FALSE) +
  labs(
    x = paste0(x, " (", scales::percent(summary(pca)$importance[2, x]), ")"),
    y = paste0(y, " (", scales::percent(summary(pca)$importance[2, y]), ")")
  )

v$targets %>%
  ggplot(aes(condition, sample.weights, fill = mouse)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 1, linetype = 2) +
  facet_wrap(~mouse)

## Calculate correlations
# apply duplicateCorrelation is two rounds
tmp <- voom(dgeFilt, design, plot=FALSE)
dupcor <- duplicateCorrelation(tmp, design, block=dgeFilt$samples$mouse)

# run voom considering the duplicateCorrelation results
# in order to compute more accurate precision weights
v = voomWithQualityWeights(
  counts = dgeFilt, 
  design = matrix(1, nrow = ncol (dgeFilt)), 
  plot=FALSE, 
  block=dgeFilt$samples$mouse, 
  correlation=dupcor$consensus
)

# Estimate linear mixed model with a single variance component
# Fit the model for each gene, 
dupcor <- duplicateCorrelation(v, design, block=dgeFilt$samples$mouse)

# But this step uses only the genome-wide average for the random effect
#fitDupCor <- lmFit(vobj, design, block=metadata$Individual, correlation=dupcor$consensus)

cont.matrix <- makeContrasts(
  DNvLAG3 = DN - LAG_3,
  DNvDP = DN - DP,
  DNvCD49b = DN - CD49b,
  LAG3vCD49b = LAG_3 - CD49b,
  CD49bvDP = CD49b - DP,
  LAG3vDP = LAG_3 - DP,
  levels = design)

#limma
fit <- lmFit(v, design = design, block = v$targets$mouse, correlation = dupcor$consensus)
fit.cont <-  contrasts.fit(fit, cont.matrix) %>%
  eBayes()

summa.fit <- decideTests(fit.cont, lfc = 1)
summary(summa.fit) 
plotSA(fit.cont)


#####################################################fit cont 244 genes###########################

###################################################################################################

#Heatmap
mat_col <- data.frame(group = v$design)
mat_colours <- list(group = brewer.pal(4, "Set1"))
names(mat_colours$group) <- unique(mat_col)
valueRange <- fit$coefficients[HighConfGenes_CPM, ]%>%
  subtract(rowMeans(.)) %>%
  range()
myPalette <- colorRampPalette(c("blue", "white", "red"))(101)
myBreaks <- c(seq(min(valueRange),0, length.out = 50),
              seq(max(valueRange) / 101, max(valueRange), length.out = 50))

#interestedGenes = c("Eomes", "P2rx7", "Cd226")
#merged_DE_mt_geneIDs = merged_DE_mt$ID.gene_name

pheatmap(mat = fit.cont,
         color = myPalette,
         breaks = myBreaks,
         border_color = NA,
         # show_colnames = FALSE,
         show_rownames = FALSE,
         cluster_cols = FALSE,
         scale = "none",
         #labels_col = hello,
         #labels_row = merged_DE_mt$ID.gene_name,
         #cutree_rows = 3,
         drop_levels = TRUE,
         fontsize = 11)

###########################################################################################################

#MAke TAbles sorted by p value and lfc
topTable(fit.cont, coef = "DNvDP", number = Inf, sort.by = "p", p.value = 0.05, lfc = 1) -> DNvDP
DNvDP %>% select(ID.gene_id, ID.gene_name, P.Value) -> FiltDNDP

topTable(fit.cont, coef = "DNvLAG3", number = Inf, sort.by = "p", p.value = 0.05, lfc = 1) -> DNvLAG3
DNvLAG3 %>% select(ID.gene_id, ID.gene_name, P.Value) -> FiltDNLAG3

topTable(fit.cont, coef = "DNvCD49b", number = Inf, sort.by = "p", p.value = 0.05, lfc = 1) -> DNvCD49b
DNvCD49b %>% select(ID.gene_id, ID.gene_name, P.Value) -> FiltDNCD49b  

topTable(fit.cont, coef = "LAG3vDP", number = Inf, sort.by = "p", p.value = 0.05, lfc = 1) -> LAG3vDP
LAG3vDP %>% select(ID.gene_id, ID.gene_name, P.Value) -> FiltLAG3DP

topTable(fit.cont, coef = "CD49bvDP", number = Inf, sort.by = "p", p.value = 0.05, lfc = 1) -> CD49bvDP
CD49bvDP %>% select(ID.gene_id, ID.gene_name, P.Value) -> FiltCD49bDP

topTable(fit.cont, coef = "LAG3vCD49b", number = Inf, sort.by = "p", p.value = 0.05, lfc = 1) -> LAG3vCD49b
LAG3vCD49b %>% select(ID.gene_id, ID.gene_name, P.Value) -> FiltLAG3CD49b


###################### union matix #####################################################################

merged_DE_mt = rbind(FiltCD49bDP, FiltDNCD49b, FiltDNDP, FiltDNLAG3, FiltLAG3CD49b, FiltLAG3DP)
merged_DE_mt = merged_DE_mt[!duplicated(merged_DE_mt$ID.gene_id), ]

merged_T100_LFC %>% 
  select(ID.gene_id, ID.gene_name, AveExpr, logFC, adj.P.Val)  %>%
  arrange(logFC)-> T100LFCarr
T100genes <- head(T100LFCarr, 100) %>% arrange(logFC)

merged_T100_LFC %>% 
  select(ID.gene_id, ID.gene_name, AveExpr, logFC, adj.P.Val)  %>%
  arrange(adj.P.Val)-> T100LFCarrPV
T100genesPV <- head(T100LFCarr, 100) %>% arrange(adj.P.Val)


#Heatmap
mat_col <- data.frame(group = v$design)
mat_colours <- list(group = brewer.pal(4, "Set1"))
names(mat_colours$group) <- unique(mat_col)
valueRange <- fit.cont$coefficients %>%
  subtract(rowMeans(.)) %>%
  range()
myPalette <- colorRampPalette(c("blue", "white", "red"))(101)
myBreaks <- c(seq(min(merged_T100_LFC$logFC), 0, length.out = 35),
              seq(max(merged_T100_LFC$logFC) / 101, max(merged_T100_LFC$logFC), length.out = 40))

pheatmap(mat = fit$coefficients[T100genesPV$ID.gene_id, ],
         color = myPalette,
         breaks = myBreaks,
         border_color = NA,
         # show_colnames = FALSE,
         show_rownames = TRUE,
         cluster_cols = FALSE,
         scale = "none",
         #labels_col = hello,
         labels_row = T100genesPV$ID.gene_name,
         #cutree_rows = 3,
         drop_levels = TRUE,
         fontsize = 11)


###################################################################################################################
#Heatmap
mat_col <- data.frame(group = v$design)
mat_colours <- list(group = brewer.pal(4, "Set1"))
names(mat_colours$group) <- unique(mat_col)
valueRange <- fit.cont$coefficients %>%
  subtract(rowMeans(.)) %>%
  range()
myPalette <- colorRampPalette(c("blue", "white", "red"))(101)
myBreaks <- c(seq(min(valueRange), 0, length.out =40),
              seq(max(valueRange) / 101, max(valueRange), length.out = 50))

#interestedGenes = c("Eomes", "P2rx7", "Cd226")
#merged_DE_mt_geneIDs = merged_DE_mt$ID.gene_name

pheatmap(mat = fit$coefficients[merged_DE_mt$ID.gene_id, ],
         color = myPalette,
         breaks = myBreaks,
         border_color = NA,
         # show_colnames = FALSE,
         show_rownames = TRUE,
         cluster_cols = FALSE,
         scale = "none",
         #labels_col = hello,
         labels_row = merged_DE_mt$ID.gene_name,
         #cutree_rows = 3,
         drop_levels = TRUE,
         fontsize = 11)

#fit$coefficients[merged_DE_mt$ID.gene_id, ]

###############################Excel Spreadsheet##############################################################
# make one of these for each of the Filt contrasts 
merged_DE_df = data.frame(gene_id = merged_DE_mt$ID.gene_id, gene_name = merged_DE_mt$ID.gene_name, fit$coefficients[merged_DE_mt$ID.gene_id, ])
write.csv(merged_DE_df, file = "C:/Users/caitlin/Desktop/Tr1 RNA seq/merged_DE_mt.csv", col.names = T, row.names = F)

#Filt contrast spreadsheets
FiltCD49bDP = data.frame(gene_id = FiltCD49bDP$ID.gene_id, gene_name = FiltCD49bDP$ID.gene_name, fit$coefficients[FiltCD49bDP$ID.gene_id, ])
write.csv(FiltCD49bDP, file = "C:/Users/caitlin/Desktop/Tr1 RNA seq/FiltCD49bDP.csv", col.names = T, row.names = F)

##Tibble with all four comparisons
allCont <- colnames(fit.cont)
Comparisons <- lapply(allCont, function(x){
  topTable(fit.cont, coef = x, number = nrow(fit)) %>%
    rownames_to_column("Gene") %>%
    as_tibble()
}) %>% 
  set_names(allCont)

######################BIological replicates#####################################################################

###BR_CPMlog
dgeFilt.cpm = cpm(dgeFilt)
BR_CPM = dgeFilt.cpm[as.character(rownames(merged_T100_LFC)), ]
BR_CPM = dgeFilt.cpm[as.character(rownames(dgeFilt.cpm)) %in% as.character(rownames(merged_T100_LFC)), ]
dim(BR_CPM)
table(as.character(rownames(BR_CPM))==as.character(rownames(merged_T100_LFC)))
BR_CPMlog = log10(BR_CPM+1) 

#BR_CPMlog_minor = BR_CPMlog[, c("M2_DN", "M5_DN")] run this before heatmap to re-order replicates
#Heatmap BR
myPalette <- colorRampPalette(c("white", "red"))(101)
myBreaks <- c(seq(min(BR_CPMlog), 1, length.out = 40),
              seq(max(BR_CPMlog) / 10, max(BR_CPMlog), length.out = 70))

pheatmap(mat = BR_CPMlog, 
         color = myPalette,
         breaks = myBreaks,
         border_color = "black",
         # show_colnames = FALSE,
         show_rownames = TRUE,
         cluster_cols = FALSE,
         scale = "none",
         #labels_col = hello,
         labels_row = merged_T100_LFC$ID.gene_name,
         #cutree_rows = 3,
         drop_levels = TRUE,
         fontsize = 11)

###BR_LCPM###
dgeFilt.logcpm = cpm(dgeFilt, log = TRUE)
BR_LCPM = dgeFilt.logcpm[as.character(rownames(merged_T100_LFC)), ] ##for using in built log function
BR_CPM = dgeFilt.logcpm[as.character(rownames(dgeFilt.logcpm)) %in% as.character(rownames(merged_T100_LFC)), ]
dim(BR_LCPM)
table(as.character(rownames(BR_LCPM))==as.character(rownames(merged_T100_LFC)))


#Heatmap BR
myPalette <- colorRampPalette(c("blue","white", "red"))(101)
myBreaks <- c(seq(min(valueRange), 1, length.out = 51),
              seq(max(valueRange) / 101, max(valueRange), length.out = 50))
pheatmap(mat = BR_LCPM, 
         color = myPalette,
         breaks = myBreaks,
         border_color = "black",
         # show_colnames = FALSE,
         show_rownames = TRUE,
         cluster_cols = FALSE,
         scale = "none",
         #labels_col = hello,
         labels_row = merged_T100_LFC$ID.gene_name,
         #cutree_rows = 3,
         drop_levels = TRUE,
         fontsize = 11)


################ get averaged BR_CPM###############################################################
BR_LCPM_DN = BR_LCPM[, c("M2_DN", "M3_DN", "M5_DN", "M6_DN")]
BR_LCPM_DN_mean = rowMeans(BR_LCPM_DN)

BR_LCPM_DP = BR_LCPM[, c("M2_DP", "M3_DP", "M5_DP", "M6_DP")]
BR_LCPM_DP_mean = rowMeans(BR_LCPM_DP)


BR_LCPM_CD49b = BR_LCPM[, c("M2_CD49b", "M3_CD49b", "M5_CD49b", "M6_CD49b")]
BR_LCPM_CD49b_mean = rowMeans(BR_LCPM_CD49b)

BR_LCPM_LAG_3 = BR_LCPM[, c("M2_LAG_3", "M3_LAG_3", "M5_LAG_3", "M6_LAG_3")]
BR_LCPM_LAG_3_mean = rowMeans(BR_LCPM_LAG_3)

BR_LCPM_mean = cbind(BR_LCPM_DN_mean, 
                                  BR_LCPM_DP_mean,
                                  BR_LCPM_LAG_3_mean,
                                  BR_LCPM_CD49b_mean)

colnames(BR_LCPM_mean) = c("DN", "DP", "LAG_3", "CD49b")
BR_CPM_mean_log = log10(BR_CPM_mean + 1)

newgenenames <- lapply(
 rownames(merged_T100_LFC),
 function(x) bquote(italic(.(x))))

# make gene names in italics
#gene_list <- as.character(df$Gene)

#make_italics <- function(x) {
#  as.expression(lapply(x, function(y) bquote(italic(.(y)))))
#}

#Heatmaps
mat_col <- data.frame(group = v$design)
mat_colours <- list(group = brewer.pal(4, "Set1"))
names(mat_colours$group) <- unique(mat_col)
valueRange <- BR_LCPM_mean %>%
subtract(rowMeans(.)) %>%
range()
#myPalette <- colorRampPalette(c("blue","white", "red"))(101)
#myBreaks <- seq(min(BR_CPMlog), max(BR_CPMlog), 101)

myPalette <- colorRampPalette(c("blue", "white", "red"))(101)
myBreaks <- c(seq(min(BR_LCPM_mean), 0, length.out = 50),
              seq(max(BR_LCPM_mean) / 101, max(BR_LCPM_mean), length.out = 50))

pheatmap(mat = BR_LCPM_mean,
         color = myPalette,
         breaks = myBreaks,
         border_color = "black",
         #show_colnames = FALSE,
         show_rownames = TRUE, 
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         scale = "none",
         #labels_col = hello,
         labels_row = merged_T100_LFC$ID.gene_name,
         #cutree_rows = 3,
         drop_levels = TRUE,
         fontsize = 18)

###############################################################################################

#########################code for volcano plot comparisons######################################
options(ggrepel.max.overlaps = Inf)

currentComp <- "DNvDP"
Comparisons[[currentComp]] %>%
  dplyr::rename(fdr = adj.P.Val,
                Name = ID.gene_name) %>%
  mutate(DE = abs(logFC) >1 & fdr < 0.01,
         DE = sign(logFC) * DE,
         DE = as.factor(DE)) %>%
  arrange(desc(P.Value)) %>%
  ggplot(aes(logFC, -log(P.Value), label = Name, fdr = fdr)) +
  geom_point(aes(colour = DE)) +
  ggtitle(currentComp) +
  geom_text_repel(data = . %>% tail(25), size =4, aes(fontface=3)) +
  xlim(-5, 5) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  theme(text = element_text(size = 24)) 

currentComp <- "DNvLAG3"
Comparisons[[currentComp]] %>%
  dplyr::rename(fdr = adj.P.Val,
                Name = ID.gene_name) %>%
  mutate(DE = abs(logFC) >1 & fdr < 0.01,
         DE = sign(logFC) * DE,
         DE = as.factor(DE)) %>%
  arrange(desc(P.Value)) %>%
  ggplot(aes(logFC, -log(P.Value), label = Name, fdr = fdr)) +
  geom_point(aes(colour = DE)) +
  ggtitle(currentComp) +
  geom_text_repel(data = . %>% tail(25), size =4, aes(fontface=3)) +
  xlim(-5, 5) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  theme(text = element_text(size = 24)) 

currentComp <- "DNvCD49b"
Comparisons[[currentComp]] %>%
  dplyr::rename(fdr = adj.P.Val,
                Name = ID.gene_name) %>%
  mutate(DE = abs(logFC) >1 & fdr < 0.01,
         DE = sign(logFC) * DE,
         DE = as.factor(DE)) %>%
  arrange(desc(P.Value)) %>%
  ggplot(aes(logFC, -log(P.Value), label = Name, fdr = fdr)) +
  geom_point(aes(colour = DE)) +
  ggtitle(currentComp) +
  geom_text_repel(data = . %>% tail(25), size =4, aes(fontface=3)) +
  xlim(-5, 5) +
  scale_color_manual(values = c("grey", "blue")) +
  theme(text = element_text(size = 24)) 

currentComp <- "LAG3vDP"
Comparisons[[currentComp]] %>%
  dplyr::rename(fdr = adj.P.Val,
                Name = ID.gene_name) %>%
  mutate(DE = abs(logFC) >1 & fdr < 0.01,
         DE = sign(logFC) * DE,
         DE = as.factor(DE)) %>%
  arrange(desc(P.Value)) %>%
  ggplot(aes(logFC, -log(P.Value), label = Name, fdr = fdr)) +
  geom_point(aes(colour = DE)) +
  ggtitle(currentComp) +
  geom_text_repel(data = . %>% tail(25), size =4, aes(fontface=3)) +
  xlim(-5, 5) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  theme(text = element_text(size = 24)) 

currentComp <- "LAG3vCD49b"
Comparisons[[currentComp]] %>%
  dplyr::rename(fdr = adj.P.Val,
                Name = ID.gene_name) %>%
  mutate(DE = abs(logFC) >1 & fdr < 0.01,
         DE = sign(logFC) * DE,
         DE = as.factor(DE)) %>%
  arrange(desc(P.Value)) %>%
  ggplot(aes(logFC, -log(P.Value), label = Name, fdr = fdr)) +
  geom_point(aes(colour = DE)) +
  ggtitle(currentComp) +
  geom_text_repel(data = . %>% tail(25), size =4, aes(fontface=3)) +
  xlim(-5, 5) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  theme(text = element_text(size = 24))

currentComp <- "CD49bvDP"
Comparisons[[currentComp]] %>%
  dplyr::rename(fdr = adj.P.Val,
                Name = ID.gene_name) %>%
  mutate(DE = abs(logFC) >1 & fdr < 0.05,
         DE = sign(logFC) * DE,
         DE = as.factor(DE)) %>%
  arrange(desc(P.Value)) %>%
  ggplot(aes(logFC, -log(P.Value), label = Name, fdr = fdr)) +
  geom_point(aes(colour = DE)) +
  ggtitle(currentComp) +
  geom_text_repel(data = . %>% tail(25), size =4, aes(fontface=3)) +
  xlim(-5, 5) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  theme(text = element_text(size = 24))
  
##############updated code for heatmaps#############################

options(ggrepel.max.overlaps = Inf)

##change current comp to plot each comparison###
  
currentComp <- "CD49bvDP"
Comparisons[[currentComp]] %>%
dplyr::rename(fdr = adj.P.Val,
                Name = ID.gene_name) %>%
  mutate(DE = abs(logFC) >1 & fdr < 0.05,
         DE = sign(logFC) * DE,
         DE = as.factor(DE)) %>%
  arrange(desc(P.Value)) %>%
  ggplot(aes(logFC, -log(P.Value), label = Name, fdr = fdr)) +
  geom_point(aes(colour = DE,)) +
  ggtitle(currentComp) +
  geom_text_repel(data = . %>% tail(25), aes(size = 20)) +
  xlim(-5, 5) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  theme(text = element_text(size = 24)) 

##################################################################

currentComp <- "DNvLAG3"
Comparisons[[currentComp]] %>%
  dplyr::rename(fdr = adj.P.Val,
                Name = ID.gene_name) %>%
  mutate(DE = abs(logFC) >1 & fdr < 0.01,
         DE = sign(logFC) * DE,
         DE = as.factor(DE)) %>%
  arrange(desc(P.Value)) %>%
  ggplot(aes(logFC, -log(P.Value), label = Name, fdr = fdr)) +
  geom_point(aes(colour = DE)) +
  ggtitle(currentComp) +
  geom_text_repel(data = . %>% tail(30)) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  theme_bw()

currentComp <- "LAG3vDP"
Comparisons[[currentComp]] %>%
  dplyr::rename(fdr = adj.P.Val,
                Name = ID.gene_name) %>%
  mutate(DE = abs(logFC) >1 & fdr < 0.01,
         DE = sign(logFC) * DE,
         DE = as.factor(DE)) %>%
  arrange(desc(P.Value)) %>%
  ggplot(aes(logFC, -log(P.Value), label = Name, fdr = fdr)) +
  geom_point(aes(colour = DE)) +
  ggtitle(currentComp) +
  geom_text_repel(data = . %>% tail(30)) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  theme_bw()

currentComp <- "LAG3vCD49b"
Comparisons[[currentComp]] %>%
  dplyr::rename(fdr = adj.P.Val,
                Name = ID.gene_name) %>%
  mutate(DE = abs(logFC) >1 & fdr < 0.01,
         DE = sign(logFC) * DE,
         DE = as.factor(DE)) %>%
  arrange(desc(P.Value)) %>%
  ggplot(aes(logFC, -log(P.Value), label = Name, fdr = fdr)) +
  geom_point(aes(colour = DE)) +
  ggtitle(currentComp) +
  geom_text_repel(data = . %>% tail(30)) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  theme_bw()

currentComp <- "LAG3vCD49b"
Comparisons[[currentComp]] %>%
  dplyr::rename(fdr = adj.P.Val,
                Name = ID.gene_name) %>%
  mutate(DE = abs(logFC) >1 & fdr < 0.01,
         DE = sign(logFC) * DE,
         DE = as.factor(DE)) %>%
  arrange(desc(P.Value)) %>%
  ggplot(aes(logFC, -log(P.Value), label = Name, fdr = fdr)) +
  geom_point(aes(colour = DE)) +
  ggtitle(currentComp) +
  geom_text_repel(data = . %>% tail(30)) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  theme_bw()

currentComp <- "CD49bvDP"
Comparisons[[currentComp]] %>%
  dplyr::rename(fdr = adj.P.Val,
                Name = ID.gene_name) %>%
  mutate(DE = abs(logFC) >1 & fdr < 0.01,
         DE = sign(logFC) * DE,
         DE = as.factor(DE)) %>%
  arrange(desc(P.Value)) %>%
  ggplot(aes(logFC, -log(P.Value), label = Name, fdr = fdr)) +
  geom_point(aes(colour = DE)) +
  ggtitle(currentComp) +
  geom_text_repel(data = . %>% tail(30)) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  theme_bw()

##Tibble with all three comparisons
allCont <- colnames(fit.cont)
Comparisons <- lapply(allCont, function(x){
  topTable(fit.cont, coef = x, number = nrow(fit)) %>%
    rownames_to_column("Gene") %>%
    as_tibble()
}) %>% 
  set_names(allCont)

#High confidence table
highConfTable <- lapply(allCont, function(x){
  topTable(fit.cont, coef = x, number = nrow(fit), sort.by = "p", p.value = 0.05, lfc = 1) %>%
    rownames_to_column("Gene") %>%
    as_tibble()
}) %>% 
  set_names(allCont)

HighConfGenes <- highConfTable%>%
  lapply(extract2, "Gene") %>%   #isolate ensembl IDs
  unlist() %>%  #merge everything into one list
  unique() #retrieve only unique ensembl IDs
length(HighConfGenes)

###Heatmap genes for single comparison##########################################
## remember the name for LAG3 is LAG_3####

#Heatmap
mat_col <- data.frame(group = v$design)
mat_colours <- list(group = brewer.pal(4, "Set1"))
names(mat_colours$group) <- unique(mat_col)
valueRange <- fit.cont$coefficients[HighConfGenes,] %>%
 subtract(rowMeans(.)) %>%
  range()
myPalette <- colorRampPalette(c("blue", "white", "red"))(101)
myBreaks <- c(seq(min(valueRange), 0, length.out = 40),
              seq(max(valueRange) / 101, max(valueRange), length.out = 50))


pheatmap(mat = fit$coefficients[FiltLAG3CD49b$ID.gene_id, c("LAG_3", "CD49b")],
         color = myPalette,
         breaks = myBreaks,
         border_color = "black",
         # show_colnames = FALSE,
         show_rownames = TRUE,
         cluster_cols = FALSE,
         scale = "none",
         #labels_col = hello,
         labels_row = FiltLAG3CD49b$ID.gene_name,
         #cutree_rows = 3,
         drop_levels = TRUE,
         fontsize = 11)


##############################################################

#Heatmap
mat_col <- data.frame(group = v$design)
mat_colours <- list(group = brewer.pal(4, "Set1"))
names(mat_colours$group) <- unique(mat_col)
valueRange <- fit.cont$coefficients[HighConfGenes,] %>%
  subtract(rowMeans(.)) %>%
  range()
myPalette <- colorRampPalette(c("blue", "white", "red"))(101)
myBreaks <- c(seq(min(valueRange), 0, length.out = 50),
              seq(max(valueRange) / 101, max(valueRange), length.out = 50))


pheatmap(mat = fit$coefficients[FiltLAG3CD49b$ID.gene_id, c("LAG3", "CD49b")],
         color = myPalette,
         breaks = myBreaks,
         border_color = NA,
         # show_colnames = FALSE,
         show_rownames = TRUE,
         cluster_cols = FALSE,
         scale = "none",
         #labels_col = hello,
         labels_row = FiltLAG3CD49b$ID.gene_name,
         #cutree_rows = 3,
         drop_levels = TRUE,
         fontsize = 11)

####################################### CPM plots   ##################################
#cpm plots with pretty colours
#to get mean CPM for each population remove #before summarise line and change Sample on ggplot line with condition
cpm(dgeList, log = TRUE)["ENSMUSG00000029378",] %>%
  as.data.frame() %>%
  set_colnames("CPM") %>%
  rownames_to_column("Sample") %>%
  left_join(rownames_to_column(dgeFilt$samples, "Sample")) %>%
  #filter(condition != "DN") %>%
  droplevels() %>%
  arrange(Sample) %>%
  group_by(condition) %>%
  #summarise(CPM = mean(CPM)) %>%
  ggplot(aes(Sample, CPM, fill = condition)) +
  scale_fill_manual(values = c("#000066","#228B22", "#FFA500", "#B22222") ) +  
  ylim(-1, 12) +
  geom_bar(stat = "identity", color="black", position = position_dodge()) + 
  ##geom_text to label the value for each column CPM remove to visualise the pops change colour to "#0099FF"
  #geom_text(aes(label = stat(y), angle = 90)) +
  #scale_y_continuous(breaks = seq(0, 12, by = 2,)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

###HOW to get values from all samples for one gene####
cpm(dgeList, log = TRUE) ["ENSMUSG00000030342", ]


#cpm plots with pretty colours
#to get mean CPM for each population remove #before summarise line and change Sample on ggplot line with condition
cpm(dgeList, log = TRUE)["ENSMUSG00000030342",] %>%
  as.data.frame() %>%
  set_colnames("CPM") %>%
  rownames_to_column("Sample") %>%
  left_join(rownames_to_column(dgeFilt$samples, "Sample")) %>%
  #filter(condition != "DN") %>%
  droplevels() %>%
  arrange(Sample) %>%
  group_by(condition) %>%
  summarise(CPM = mean(CPM)) %>%
  ggplot(aes(condition, CPM, fill = condition)) +
  scale_fill_manual(values = c("#000066", "#228B22", "#FFA500", "#B22222") ) +  
  ylim(-1, 12) +
  geom_bar(stat = "identity") +
  #scale_y_continuous(breaks = seq(0, 12, by = 2,)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust =1, vjust = 0.5))

###### HOW TO rank genes by LFC #####
# to see which genes are upregulated and downregulated
DNvDPLFC <- head(DNvDP, 100) %>% arrange(logFC)
DNvLAG3LFC <- head(DNvLAG3, 100) %>% arrange(logFC)
DNvCD49bLFC <- head(DNvCD49b, 100) %>% arrange(logFC)
LAG3vDPLFC <- head(LAG3vDP, 100) %>% arrange(logFC)
LAG3vCD49bLFC <- head(LAG3vCD49b, 100) %>% arrange(logFC)
CD49bvDPLFC <- head(CD49bvDP, 100) %>% arrange(logFC)


#####need to work out how to plot DNup as a heatmap need to add gene id? different approach below
######UP on DN

#DN pop
UponDN = subset(DNvDPLFC, logFC >=0.1) %>% select(ID.gene_id, ID.gene_name, AveExpr)
UponDN = c(subset(DNvLAG3LFC, logFC>=0.1) %>% select(ID.gene_id, ID.gene_name, AveExpr), UponDN) 
UponDN = c(subset(DNvCD49bLFC, logFC>=0.1) %>% select(ID.gene_id, ID.gene_name, AveExpr), UponDN)
#UponDN %>% unlist() %>% unique() -> UponDN

#Heatmap UponDN
mat_col <- data.frame(group = v$design)
mat_colours <- list(group = brewer.pal(4, "Set1"))
names(mat_colours$group) <- unique(mat_col)
valueRange <- fit.cont$coefficients[HighConfGenes,] %>%
  subtract(rowMeans(.)) %>%
  range()
myPalette <- colorRampPalette(c("blue", "white", "red"))(101)
myBreaks <- c(seq(min(valueRange), 0, length.out = 51),
              seq(max(valueRange) / 101, max(valueRange), length.out = 35))

pheatmap(mat = fit$coefficients[UponDN$ID.gene_id, ] %>%
           subtract(rowMeans(.)),
         color = myPalette,
         breaks = myBreaks,
         border_color = "black", 
         # show_colnames = FALSE, 
         show_rownames = TRUE, 
         cluster_cols = FALSE,
         scale = "none",
         labels_row = UponDN$ID.gene_name,
         #cutree_rows = 3,
         drop_levels = TRUE, 
         fontsize = 11)  
 #########################UP on DP###################################
#DP pop
UponDP = subset(DNvDPLFC, logFC<=0.1) %>% select(ID.gene_id, ID.gene_name, AveExpr)
UponDP = c(subset(LAG3vDPLFC, logFC<=0.1) %>% select(ID.gene_id, ID.gene_name, AveExpr), UponDP) 
UponDP = c(subset(CD49bvDPLFC, logFC<=0.1) %>% select(ID.gene_id, ID.gene_name, AveExpr), UponDP)
#UponDN %>% unlist() %>% unique() -> UponDN

#Heatmap UponDP
mat_col <- data.frame(group = v$design)
mat_colours <- list(group = brewer.pal(4, "Set1"))
names(mat_colours$group) <- unique(mat_col)
valueRange <- fit.cont$coefficients[HighConfGenes,] %>%
  subtract(rowMeans(.)) %>%
  range()
myPalette <- colorRampPalette(c("blue", "white", "red"))(101)
myBreaks <- c(seq(min(valueRange), 0, length.out = 51),
              seq(max(valueRange) / 101, max(valueRange), length.out = 35))

pheatmap(mat = fit$coefficients[UponDP$ID.gene_id, ] %>%
           subtract(rowMeans(.)),
         color = myPalette,
         breaks = myBreaks,
         border_color = "black", 
         # show_colnames = FALSE, 
         show_rownames = TRUE, 
         cluster_cols = FALSE,
         scale = "none",
         labels_row = UponDP$ID.gene_name,
         #cutree_rows = 3,
         drop_levels = TRUE, 
         fontsize = 11)  

#########################UP on LAG3###################################
#LAG3 pop
UponLAG3 = subset(DNvLAG3LFC, logFC<=0.1) %>% select(ID.gene_id, ID.gene_name, AveExpr)
UponLAG3 = c(subset(LAG3vDPLFC, logFC>=0.1) %>% select(ID.gene_id, ID.gene_name, AveExpr), UponLAG3) 
UponLAG3 = c(subset(LAG3vCD49bLFC, logFC>=0.1) %>% select(ID.gene_id, ID.gene_name, AveExpr), UponLAG3)
#UponDN %>% unlist() %>% unique() -> UponDN

#Heatmap UponLAG3
mat_col <- data.frame(group = v$design)
mat_colours <- list(group = brewer.pal(4, "Set1"))
names(mat_colours$group) <- unique(mat_col)
valueRange <- fit.cont$coefficients[HighConfGenes,] %>%
  subtract(rowMeans(.)) %>%
  range()
myPalette <- colorRampPalette(c("blue", "white", "red"))(101)
myBreaks <- c(seq(min(valueRange), 0, length.out = 51),
              seq(max(valueRange) / 101, max(valueRange), length.out = 35))

pheatmap(mat = fit$coefficients[UponLAG3$ID.gene_id, ] %>%
           subtract(rowMeans(.)),
         color = myPalette,
         breaks = myBreaks,
         border_color = "black", 
         # show_colnames = FALSE, 
         show_rownames = TRUE, 
         cluster_cols = FALSE,
         scale = "none",
         labels_row = UponLAG3$ID.gene_name,
         #cutree_rows = 3,
         drop_levels = TRUE, 
         fontsize = 11)  

#########################UP on CD49b###################################
#CD49b pop
UponCD49b = subset(DNvCD49bLFC, logFC<=0.1) %>% select(ID.gene_id, ID.gene_name, AveExpr)
UponCD49b = c(subset(CD49bvDPLFC, logFC>=0.1) %>% select(ID.gene_id, ID.gene_name, AveExpr), UponCD49b) 
UponCD49b = c(subset(LAG3vCD49bLFC, logFC<=0.1) %>% select(ID.gene_id, ID.gene_name, AveExpr), UponCD49b)
#UponDN %>% unlist() %>% unique() -> UponDN

#Heatmap UponCD49b
mat_col <- data.frame(group = v$design)
mat_colours <- list(group = brewer.pal(4, "Set1"))
names(mat_colours$group) <- unique(mat_col)
valueRange <- fit.cont$coefficients[HighConfGenes,] %>%
  subtract(rowMeans(.)) %>%
  range()
myPalette <- colorRampPalette(c("blue", "white", "red"))(101)
myBreaks <- c(seq(min(valueRange), 0, length.out = 51),
              seq(max(valueRange) / 101, max(valueRange), length.out = 35))

pheatmap(mat = fit$coefficients[UponCD49b$ID.gene_id, ] %>%
           subtract(rowMeans(.)),
         color = myPalette,
         breaks = myBreaks,
         border_color = "black", 
         # show_colnames = FALSE, 
         show_rownames = TRUE, 
         cluster_cols = FALSE,
         scale = "none",
         labels_row = UponCD49b$ID.gene_name,
         #cutree_rows = 3,
         drop_levels = TRUE, 
         fontsize = 11)  

###################### union matix LFC #####################################################################

merged_T100_LFC = rbind(DNvDPLFC, DNvLAG3LFC, DNvCD49bLFC, LAG3vDPLFC, LAG3vCD49bLFC, CD49bvDPLFC)
merged_T100_LFC = merged_T100_LFC[!duplicated(merged_T100_LFC$ID.gene_id), ]


# make one of these for each of the Filt contrasts 
write.csv(merged_DE_df, file = "C:/Users/caitlin/Desktop/Tr1 RNA seq/merged_DE_df.csv", col.names = T, row.names = F)

write.csv(FiltCD49bDP, file = "C:/Users/caitlin/Desktop/Tr1 RNA seq/FiltCD49bDP.csv", col.names = T, row.names = F)

write.csv(DNvLAG3, file = "C:/Users/caitlin/Desktop/Tr1 RNA seq/DNvLAG3.csv", col.names = T, row.names = F)



###############################Excel Spreadsheet##############################################################
# make one of these for each of the Filt contrasts 
merged_DE_df = data.frame(gene_id = merged_DE_mt$ID.gene_id, gene_name = merged_DE_mt$ID.gene_name, fit$coefficients[merged_DE_mt$ID.gene_id, ])
write.csv(merged_DE_df, file = "C:/Users/caitlin/Desktop/Tr1 RNA seq/merged_DE_mt.csv", col.names = T, row.names = F)

#Filt contrast spreadsheets
FiltCD49bDP = data.frame(gene_id = FiltCD49bDP$ID.gene_id, gene_name = FiltCD49bDP$ID.gene_name, fit$coefficients[FiltCD49bDP$ID.gene_id, ])
write.csv(FiltCD49bDP, file = "C:/Users/caitlin/Desktop/Tr1 RNA seq/FiltCD49bDP.csv", col.names = T, row.names = F)

#All genes table
Allgenestable <- lapply(allCont, function(x){
  topTable(fit.cont, coef = x, number = nrow(fit), sort.by = "p", p.value = 0.05, lfc = 1) %>%
    rownames_to_column("Gene") %>%
    as_tibble()
}) %>% 
  set_names(allCont)


#Heatmap
mat_col <- data.frame(group = v$design)
mat_colours <- list(group = brewer.pal(4, "Set1"))
names(mat_colours$group) <- unique(mat_col)
valueRange <- fit.cont$coefficients[merged_DE_df$CD49b, ] %>%
 subtract(rowMeans(.)) %>%
  range()
myPalette <- colorRampPalette(c("blue", "white", "red"))(101)
myBreaks <- c(seq(min(valueRange), 0, length.out = 50),
              seq(max(valueRange) / 101, max(valueRange), length.out = 50))

pheatmap(mat = merged_T100_LFC,
         color = myPalette,
         breaks = myBreaks,
         border_color = NA,
         #show_colnames = FALSE,
         show_rownames = TRUE,
         cluster_cols = FALSE,
         scale = "none",
         #labels_col = hello,
         labels_row = merged_T100_LFC$ID.gene_name,
         #cutree_rows = 3,
         drop_levels = TRUE,
         fontsize = 11)


ImpGenes <- Allgenestable%>%
  lapply(extract2, "Gene") %>%   #isolate ensembl IDs
  unlist() %>%  #merge everything into one list
  unique() #retrieve only unique ensembl IDs
length(ImpGenes)

#only extract the confident genes from the highConfGenes list without the LFC
Allgenestable <- sapply(allCont, function(x){
  topTable(fit.cont, coef = x, number = Inf, p.value = 0.05, lfc = 1) %>%
    rownames_to_column("Gene") %>%
    as_tibble() %>%
    dplyr::filter(Gene %in% ImpGenes)
}, simplify = FALSE)
Allgenestable$genes <- genes[ImpGenes,] #table containing our high confidence genes

alltable <- topTable(fit.cont, coef = "DNvDP", number = Inf) %>% 
  rownames_to_column("Gene") %>%
  dplyr::filter(Gene %in% ImpGenes) %>%
  select(Gene, ID.gene_name, logFC, adj.P.Val, ID.entrezid)

#Try making a new table with all the highConfGenes in it regardless of P or logFC
#alltable <- topTable(fit.cont, coef = "LAG3vCD49b", number = Inf) %>% 
#  rownames_to_column("Gene") %>%
# dplyr::filter(Gene %in% HighConfGenes) %>%
#  select(Gene, ID.gene_name, logFC, adj.P.Val, ID.entrezid)

#Heatmap
mat_col <- data.frame(group = v$design)
mat_colours <- list(group = brewer.pal(4, "Set1"))
names(mat_colours$group) <- unique(mat_col)
valueRange <- fit.cont$coefficients[Heatmapgenes,] %>%
  subtract(rowMeans(.)) %>%
  range()
myPalette <- colorRampPalette(c("blue", "white", "red"))(101)
myBreaks <- c(seq(min(valueRange), 0, length.out = 51),
              seq(max(valueRange) / 101, max(valueRange), length.out = 50))

pheatmap(mat = fit$coefficients[HeatmapGenes$gene_id,] %>% subtract(rowMeans(.)),
         color = myPalette,
         breaks = myBreaks,
         border_color = NA,
         #show_colnames = FALSE,
         show_rownames = TRUE,
         cluster_cols = FALSE,
         scale = "none",
         #labels_col = hello,
         labels_row = HeatmapGenes$gene_name,
         #cutree_rows = 3,
         drop_levels = TRUE,
         fontsize = 11)

#list for heatmaps with selected genes*******************make a table containing all genes and gene IDs
highConfTable <- sapply(allCont, function(x){
  topTable(fit.cont, coef = x, number = Inf) %>%
    rownames_to_column("Gene") %>%
    as_tibble() %>%
    dplyr::filter(Gene %in% HighConfGenes)
}, simplify = FALSE) 

pheatmap(avg.mat)

########################################################### Making selected genes heatmap###########


HighConfGenes <- c("ENSMUSG00000039521", "ENSMUSG00000026770", "ENSMUSG00000068227", "ENSMUSG00000002578", 
                   "ENSMUSG00000026069", "ENSMUSG00000015619", "ENSMUSG00000000869", "ENSMUSG00000036117", 
                   "ENSMUSG00000020383", "ENSMUSG00000028150", "ENSMUSG00000025929", "ENSMUSG00000041872", 
                   "ENSMUSG00000040899", "ENSMUSG00000049103",  
                   "ENSMUSG00000050232", "ENSMUSG00000048521", "ENSMUSG00000024401", "ENSMUSG00000055170", 
                   "ENSMUSG00000022508", "ENSMUSG00000047880", "ENSMUSG00000000782", "ENSMUSG00000074607", 
                   "ENSMUSG00000079227", "ENSMUSG00000016529", "ENSMUSG00000027718", "ENSMUSG00000002603", 
                   "ENSMUSG00000015437", "ENSMUSG00000030124", "ENSMUSG00000015533", "ENSMUSG00000026009", 
                   "ENSMUSG00000026285", "ENSMUSG00000026011", "ENSMUSG00000020399", "ENSMUSG00000071552", 
                   "ENSMUSG00000034028", "ENSMUSG00000055435", "ENSMUSG00000038151", "ENSMUSG00000034266", 
                   "ENSMUSG00000019256", "ENSMUSG00000032446", "ENSMUSG00000018899", "ENSMUSG00000021356", 
                   "ENSMUSG00000026104", "ENSMUSG00000004040")

#logCPM of each population genes#

###############################################################################################################


HighConfGenes_logCPM_DN = HighConfGenes_logCPM[, c("M2_DN", "M3_DN", "M5_DN", "M6_DN")]
HighConfGenes_logCPM_DN_mean = rowMeans(HighConfGenes_logCPM_DN)

HighConfGenes_logCPM_DP = HighConfGenes_logCPM[, c("M2_DP", "M3_DP", "M5_DP", "M6_DP")]
HighConfGenes_logCPM_DP_mean = rowMeans(HighConfGenes_logCPM_DP)


HighConfGenes_logCPM_CD49b = HighConfGenes_logCPM[, c("M2_CD49b", "M3_CD49b", "M5_CD49b", "M6_CD49b")]
HighConfGenes_logCPM_CD49b_mean = rowMeans(HighConfGenes_logCPM_CD49b)

HighConfGenes_logCPM_LAG_3 = HighConfGenes_logCPM[, c("M2_LAG_3", "M3_LAG_3", "M5_LAG_3", "M6_LAG_3")]
HighConfGenes_logCPM_LAG_3_mean = rowMeans(HighConfGenes_logCPM_LAG_3)

HighConfGenes_logCPM_mean = cbind(HighConfGenes_logCPM_DN_mean, 
                                  HighConfGenes_logCPM_DP_mean,
                                  HighConfGenes_logCPM_LAG_3_mean,
                                  HighConfGenes_logCPM_CD49b_mean)


colnames(HighConfGenes_logCPM_mean) = c("DN", "DP", "LAG_3", "CD49b")
rownames(HighConfGenes_logCPM_mean) = c("Foxp3", "Il2ra", "Il2rb", "Ikzf4", "Il1rl1", 
                                        "Gata3", "Il4", "Il5", "Il13",  
                                        "Rorc","Il17a", "Il17f", "Ccr6", 
                                        "Tbx21", "Cxcr3", "Cxcr6", "Tnfa", "Ifng",
                                        "Bcl6", "Cxcr5","Tcf7","Tox2", 
                                        "Ccr5","Il10", "Il21" , "Tgfb1", "Gzmb", "Lag3", "Itga2", "Icos", "Pdcd1", "Ctla4", "Havcr2", "Tigit", "Cd226", "Maf", "Prdm1", "Batf", "Ahr", "Eomes", "Irf1", "Irf4", "Stat1", "Stat3")

newnames <- lapply(
  rownames(HighConfGenes_logCPM_mean),
  function(x) bquote(italic(.(x))))

myPalette <- colorRampPalette(c("blue", "white", "red"))(101)
myBreaks <- c(seq(min(HighConfGenes_logCPM_mean), 0, length.out = 45),
              seq(max(HighConfGenes_logCPM_mean) / 101, max(HighConfGenes_logCPM_mean), length.out = 50)) #table containing our high confidence genes

pheatmap(mat = HighConfGenes_logCPM_mean,
         color = myPalette,
         breaks = myBreaks,
         border_color = "black",
         #show_colnames = FALSE,
         show_rownames = TRUE, 
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         scale = "none",
         #labels_col = hello,
         labels_row = as.expression(newnames),
         #cutree_rows = 3,
         drop_levels = TRUE,
         fontsize = 18)



####################################################################################################

### temp codes for heatmap##############################################################################
filtLogCPM <- cpm(dgeList$counts, log = TRUE)
HighConfGenes_logCPM <- filtLogCPM[HighConfGenes, ]

myPalette <- colorRampPalette(c("blue", "white", "red"))(101)
myBreaks <- c(seq(min(HighConfGenes_logCPM), 0, length.out = 50),
              seq(max(HighConfGenes_logCPM) / 101, max(HighConfGenes_logCPM), length.out = 40)) #table containing our high confidence genes

pheatmap(mat = HighConfGenes_logCPM,
         color = myPalette,
         breaks = myBreaks,
         border_color = "black",
         #show_colnames = FALSE,
         show_rownames = TRUE,
         cluster_cols = FALSE,
         scale = "none",
         #labels_col = hello,
         labels_row = merged_DE_mt$ID.gene_name,
         #cutree_rows = 3,
         drop_levels = FALSE,
         fontsize = 11)



###############################################################################################################

############# GO ANALYSIS 4 populations #######################################################################

#FROM 3 pops seq

#GO ANALYSIS ##################################################################################################
#Define universe for go terms (Use dgefilt for all gene ID)
universe <- dgeFilt$genes$entrezid %>%
  unlist() %>%
  extract(!is.na(.)) %>% #keep non NA gene
  unique() %>%
  as.character()

#get the EnterezID from each contrast
DELAG3vCD49b <- highConfTable$LAG3vCD49b %>%
  dplyr::filter(adj.P.Val < 0.05) %>%
  dplyr::select(Gene, ID.gene_name, ID.entrezid)

DE49vDP <- highConfTable$CD49bvDP %>%
  dplyr::filter(adj.P.Val < 0.05) %>%
  dplyr::select(Gene, ID.gene_name, ID.entrezid)

DELag3vDP <- highConfTable$LAG3vDP %>%
  dplyr::filter(adj.P.Val < 0.05) %>%
  dplyr::select(Gene, ID.gene_name, ID.entrezid)

DEDNvDP <- highConfTable$DNvDP %>%
  dplyr::filter(adj.P.Val < 0.05) %>%
  dplyr::select(Gene, ID.gene_name, ID.entrezid)

DEDNvLag3 <- highConfTable$DNvLAG3 %>%
  dplyr::filter(adj.P.Val < 0.05) %>%
  dplyr::select(Gene, ID.gene_name, ID.entrezid)

DEDNvCD49b <- highConfTable$DNvCD49b %>%
  dplyr::filter(adj.P.Val < 0.05) %>%
  dplyr::select(Gene, ID.gene_name, ID.entrezid)



#ClusteredprofileR LAG-3vDP
clusterLag3vDP <- enrichGO(gene = DELag3vDP$ID.entrezid,
                           universe = universe,
                           OrgDb = org.Mm.eg.db,
                           ont = "ALL",
                           keyType = "ENTREZID",
                           pAdjustMethod = "bonferroni",
                           pvalueCutoff = 0.05,
                           qvalueCutoff = 0.05,
                           readable = TRUE)

FCLAG3DP <- structure(highConfTable$LAG3vDP$logFC, names = highConfTable$LAG3vDP$ID.entrezid %>%
                        vapply(function(x){as.character(x)[[1]]}, character(1))) #grab the fold change and make it into a character
FCLAG3DP <- FCLAG3DP[!is.na(names(FCLAG3DP))]

#make the cluster into a dataframe for easy readability
dLAG3vDP <- as.data.frame(clusterLag3vDP)
head(summary(LAG3vDP))

dotplot(clusterLag3vDP, showCategory = 10, font.size = 18) +
theme(title = element_text(size = 18),legend.text = element_text(size = 14))
  
cnetplot(clusterLag3vDP, node_label = "all", foldChange = FCLAG3DP) 
emapplot(clusterLag3vDP)

id <- clusterLag3vDP$ID[1:10]
id

clusterLag3vDP[[id[1]]]

geneInCategory(clusterLag3vDP)[id] 

GO_Lag3vDP = geneInCategory(clusterLag3vDP)[id] 


#ClusteredprofileR LAG-3 v CD49b
clusterLag3vCD49b <- enrichGO(gene = DELAG3vCD49b$ID.entrezid,
                              universe = universe,
                              OrgDb = org.Mm.eg.db,
                              ont = "ALL",
                              keyType = "ENTREZID",
                              pAdjustMethod = "bonferroni",
                              pvalueCutoff = 0.05,
                              qvalueCutoff = 0.05,
                              readable = TRUE)

FCLAG3CD49b <- structure(highConfTable$LAG3vCD49b$logFC, names = highConfTable$LAG3vCD49b$ID.entrezid %>%
                           vapply(function(x){as.character(x)[[1]]}, character(1))) #grab the fold change and make it into a character
FCLAG3CD49b <- FCLAG3CD49b[!is.na(names(FCLAG3CD49b))]



#make the cluster into a dataframe for easy readability
dCD49vLag3 <- as.data.frame(cluster49vLag3)
head(summary(dCP49vLag3))

dotplot(clusterLag3vCD49b, showCategory = 10, font.size = 18) +
theme(title = element_text(size = 18),legend.text = element_text(size = 14))

cnetplot(clusterLag3vCD49b, node_label = "all", foldChange = FCLAG3CD49b) 
emapplot(clusterLag3vDP)

id <- clusterLag3vCD49b$ID[1:10]
id

clusterLag3vCD49b[[id[1]]]

geneInCategory(clusterLag3vCD49b)[id] 

GO_Lag3vCD49b = geneInCategory(clusterLag3vCD49b)[id] 

#CLusterProfileR CD49b v DP

clusterCD49bvDP <- enrichGO(gene = DE49vDP$ID.entrezid,
                            universe = universe,
                            OrgDb = org.Mm.eg.db,
                            ont = "ALL",
                            keyType = "ENTREZID",
                            pAdjustMethod = "bonferroni",
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.05,
                            readable = TRUE)

FC49bvDP <- structure(highConfTable$CD49bvDP$logFC, names = highConfTable$CD49bvDP$ID.entrezid %>%
                        vapply(function(x){as.character(x)[[1]]}, character(1))) #grab the fold change and make it into a character
FC49bvDP <- FC49bvDP[!is.na(names(FC49bvDP))]


#make the cluster into a dataframe for easy readability
dCD49bvDP <- as.data.frame(clusterCD49bvDP)
head(summary(dCD49bvDP))

dotplot(clusterCD49bvDP, showCategory = 10, font.size = 18) +
  theme(title = element_text(size = 18),legend.text = element_text(size = 14))

cnetplot(clusterCD49bvDP, node_label = "all", foldChange = FC49bvDP) 

id <- clusterCD49bvDP$ID[1:10]
id

clusterCD49bvDP[[id[1]]]

geneInCategory(clusterCD49bvDP)[id] 

GO_CD49bvDP = geneInCategory(clusterCD49bvDP)[id] 


#CLusterProfileR DN v DP

clusterDNvDP <- enrichGO(gene = DEDNvDP$ID.entrezid,
                            universe = universe,
                            OrgDb = org.Mm.eg.db,
                            ont = "ALL",
                            keyType = "ENTREZID",
                            pAdjustMethod = "bonferroni",
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.05,
                            readable = TRUE)

FCDNvDP <- structure(highConfTable$DNvDP$logFC, names = highConfTable$DNvDP$ID.entrezid %>%
                        vapply(function(x){as.character(x)[[1]]}, character(1))) #grab the fold change and make it into a character
FCDNvDP <- FC49bvDP[!is.na(names(FCDNvDP))]


#make the cluster into a dataframe for easy readability
dDNvDP <- as.data.frame(clusterDNvDP)
head(summary(dDNvDP))
dotplot(clusterDNvDP, showCategory = 10, font.size = 18) +
  theme(title = element_text(size = 18),legend.text = element_text(size = 14))
cnetplot(clusterCD49bvDP, node_label = "all", foldChange = FC49bvDP) 

id <- clusterDNvDP$ID[1:10]
id

clusterDNvDP[[id[1]]]

geneInCategory(clusterDNvDP)[id] 

GO_DNvDP = geneInCategory(clusterDNvDP)[id] 


#CLusterProfileR DN v LAG3

clusterDNvLAG3 <- enrichGO(gene = DEDNvLag3$ID.entrezid,
                         universe = universe,
                         OrgDb = org.Mm.eg.db,
                         ont = "ALL",
                         keyType = "ENTREZID",
                         pAdjustMethod = "bonferroni",
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.05,
                         readable = TRUE)

FCDNvLAG3 <- structure(highConfTable$DNvLAG3$logFC, names = highConfTable$DNvLAG3$ID.entrezid %>%
                       vapply(function(x){as.character(x)[[1]]}, character(1))) #grab the fold change and make it into a character
FCDNvLAG3 <- FCDNvLAG3[!is.na(names(FCDNvLAG3))]


#make the cluster into a dataframe for easy readability
dDNvLAG3 <- as.data.frame(clusterDNvLAG3)
head(summary(dDNvLAG3))
dotplot(clusterDNvLAG3, showCategory = 10, font.size = 18) +
theme(title = element_text(size = 18),legend.text = element_text(size = 14))

cnetplot(clusterDNvLAG3, node_label = "all", foldChange = FCDNvLAG3) 

id <- clusterDNvLAG3$ID[1:10]
id

clusterDNvLAG3[[id[1]]]

geneInCategory(clusterDNvLAG3)[id] 

GO_DNvLAG3 = geneInCategory(clusterDNvLAG3)[id] 


#CLusterProfileR CD49b v DN

clusterCD49bvDN <- enrichGO(gene = DEDNvCD49b$ID.entrezid,
                            universe = universe,
                            OrgDb = org.Mm.eg.db,
                            ont = "ALL",
                            keyType = "ENTREZID",
                            pAdjustMethod = "bonferroni",
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.05,
                            readable = TRUE)

FC49bvDN <- structure(highConfTable$CD49bvDN$logFC, names = highConfTable$CD49bvDN$ID.entrezid %>%
                        vapply(function(x){as.character(x)[[1]]}, character(1))) #grab the fold change and make it into a character
FC49bvDN <- FC49bvDN[!is.na(names(FC49bvDN))]


#make the cluster into a dataframe for easy readability
dCD49bvDN <- as.data.frame(clusterCD49bvDN)
head(summary(dCD49bvDN))
dotplot(clusterCD49bvDN, showCategory = 10, font.size = 18) +
  theme(title = element_text(size = 18),legend.text = element_text(size = 14))

cnetplot(clusterCD49bvDN, node_label = "all", foldChange = FC49bvDP) 

id <- clusterCD49bvDN$ID[1:10]
id

clusterCD49bvDN[[id[1]]]

geneInCategory(clusterCD49bvDN)[id]

GO_CD49bvDN = geneInCategory(clusterCD49bvDN)[id] 




##################################################################END GO TERMS#####################################


mat_col <- data.frame(group = v$design)
mat_colours <- list(group = brewer.pal(4, "Set1"))
names(mat_colours$group) <- unique(mat_col)
valueRange <- fit.cont$coefficients[HighConfGenes, ] %>%
  subtract(rowMeans(.)) %>%
  range()
myPalette <- colorRampPalette(c("blue", "white", "red"))(101)
myBreaks <- c(seq(min(valueRange), 0, length.out = 51),
              seq(max(valueRange) / 101, max(valueRange), length.out = 50)) #table containing our high confidence genes


pheatmap(mat = fit$coefficients[HighConfGenes] %>% subtract(rowMeans(.)),
         color = myPalette,
         breaks = myBreaks,
         border_color = NA,
         #show_colnames = FALSE,
         show_rownames = TRUE,
         cluster_cols = FALSE,
         scale = "none",
         #labels_col = hello,
         labels_row = highConfTable$gene_name,
         #cutree_rows = 3,
         drop_levels = TRUE,
         fontsize = 11)

###Attempt to make heatmaps of shared genes #######################################################################
filt.mat = logcount[]
dat.split = lapply(split(1:ncol(filt.mat), idx), function(x) filt.mat[, x])

count = 0
split.mat = list()
for (i in unique(idx)) {
  count = count + 1
  split.mat[[count]] = filt.mat[, grep(i, colnames(filt.mat))]
}

unique(idx)

#### how to get pretty heatmaps ####

filt.mat = logcount[]
idx = sapply(stringr::str_split(colnames(filt.mat), "_"), "[[", 2)


dat.split = lapply(split(1:ncol(filt.mat), idx), function(x) filt.mat[, x])

####split.mat = lapply(split(filt.mat, names), rowMeans)####

count = 0
split.mat = list()
for (i in unique(idx)) {
  count = count + 1
  split.mat[[count]] = filt.mat[, grep(i, colnames(filt.mat))]
}

unique(idx)

count = 0
split.mat = list()
for (i in unique(idx)) {
 count = count + 1
 split.mat[[count]] = filt.mat[, grep(i, colnames(filt.mat))]
 }
 length(split.mat)
 
names(split.mat)

 avg = lapply(split.mat, rowMeans)
 names(avg) = unique(idx)
 avg.mat = do.call(cbind, split.mat)

avg.mat = do.call(cbind, avg)
head(avg.mat)

pheatmap(avg.mat)

###########################################################################################################


heatmapgenes <- c("ENSMUSG00000016529", "ENSMUSG00000030124")
mat_col <- data.frame(group = v$design)
mat_colours <- list(group = brewer.pal(4, "Set1"))
names(mat_colours$group) <- unique(mat_col)
valueRange <- fit.cont$coefficients[att2Get,] %>%
  subtract(rowMeans(.)) %>%
  range()
myPalette <- colorRampPalette(c("blue", "white", "red"))(101)
myBreaks <- c(seq(min(valueRange), 0, length.out = 51),
              seq(max(valueRange) / 101, max(valueRange), length.out = 50))


<- sapply(allCont, function(x){
  topTable(fit.cont, coef = x, number = Inf) %>%
    rownames_to_column("Gene") %>%
    as_tibble() %>%
    dplyr::filter(Gene %in% HighConfGenes)
}, simplify = FALSE)
highConfTable$genes <- genes[HighConfGenes,]

pheatmap(mat = fit$coefficients[heatmapgenes,] %>% subtract(rowMeans(.)),
         color = myPalette,
         breaks = myBreaks,
         border_color = NA,
         #show_colnames = FALSE,
         show_rownames = TRUE,
         cluster_cols = FALSE,
         scale = "none",
         #labels_col = hello,
         labels_row = heatmapgenes$gene_name,
         #cutree_rows = 3,
         drop_levels = TRUE,
         fontsize = 11)


######################TROUBLESHOOTING###############################################################


# Making isolated heatmaps for genes upregulated in the four populations ###LFC comparisons
# arrange top 100 genes  by log fold change
# to see which genes are upregulated and downregulated
DNvDPLFC <- head(DNvDP, 100) %>% arrange(logFC)
DNvLAG3LFC <- head(DNvLAG3, 100) %>% arrange(logFC)
DNvCD49bLFC <- head(DNvCD49b, 100) %>% arrange(logFC)
LAG3vDPLFC <- head(LAG3vDP, 100) %>% arrange(logFC)
LAG3vCD49bLFC <- head(LAG3vCD49b, 100) %>% arrange(logFC)
CD49bvDPLFC <- head(CD49bvDP, 100) %>% arrange(logFC)

DNup <- c(DNvDPLFC$ID.gene_name[35:61], DNvLAG3LFC$ID.gene_name[32:78], DNvCD49bLFC$ID.gene_name[3:18])
DPup <- c(DNvDPLFC$ID.gene_name[1:33], LAG3vDPLFC$ID.gene_name[1:18], CD49bvDPLFC$ID.gene_name[1:32])
LAG3up <- c(DNvLAG3LFC$ID.gene_name[1:31], LAG3vDPLFC$ID.gene_name[19:34], LAG3vCD49bLFC$ID.gene_name[30:100])
CD49bup <- c(DNvCD49bLFC$ID.gene_name[1:2], LAG3vCD49bLFC$ID.gene_name[1:29], CD49bvDPLFC$ID.gene_name[33:38])


#confirm that the genes are properly regulated
#DN pop
ConfirmDNvDP <- subset(DNvDPLFC, DNvDPLFC$ID.gene_name %in% DNup) %>% 
  select(ID.gene_name, logFC, adj.P.Val) %>%
  arrange(logFC)
ConfirmDNvLAG3 <- subset(DNvLAG3LFC, DNvLAG3LFC$ID.gene_name %in% DNup) %>% 
  select(ID.gene_name, logFC, adj.P.Val) %>%
  arrange(logFC)
ConfirmDNvCD49b <- subset(DNvCD49bLFC, DNvCD49bLFC$ID.gene_name %in% DNup) %>% 
  select(ID.gene_name, logFC, adj.P.Val) %>%
  arrange(logFC)

#LAG3 pop
ConfirmLAG3vDP <- subset(LAG3vDPLFC, LAG3vDPLFC$ID.gene_name %in% LAG3up) %>% 
  select(ID.gene_name, logFC, adj.P.Val) %>%
  arrange(logFC)
ConfirmLAG3vDN <- subset(DNvLAG3LFC, DNvLAG3LFC$ID.gene_name %in% LAG3up) %>% 
  select(ID.gene_name, logFC, adj.P.Val) %>%
  arrange(logFC)
ConfirmLAG3vCD49b <- subset(LAG3vCD49bLFC, LAG3vCD49bLFC$ID.gene_name %in% LAG3up) %>% 
  select(ID.gene_name, logFC, adj.P.Val) %>%
  arrange(logFC)


#DP pop
ConfirmDPvLAG3 <- subset(LAG3vDPLFC, LAG3vDPLFC$ID.gene_name %in% DPup) %>% 
  select(ID.gene_name, logFC, adj.P.Val) %>%
  arrange(logFC)
ConfirmDPvDN <- subset(DNvDPLFC, DNvDPLFC$ID.gene_name %in% DPup) %>% 
  select(ID.gene_name, logFC, adj.P.Val) %>%
  arrange(logFC)
ConfirmDPvCD49b <- subset(CD49bvDPLFC, CD49bvDPLFC$ID.gene_name %in% DPup) %>% 
  select(ID.gene_name, logFC, adj.P.Val) %>%
  arrange(logFC)

#CD49b pop
ConfirmCD49bvLAG3 <- subset(LAG3vCD49bLFC, LAG3vCD49bLFC$ID.gene_name %in% CD49bup) %>% 
  select(ID.gene_name, logFC, adj.P.Val) %>%
  arrange(logFC)
ConfirmCD49bvDN <- subset(DNvCD49bLFC, DNvCD49bLFC$ID.gene_name %in% CD49bup) %>% 
  select(ID.gene_name, logFC, adj.P.Val) %>%
  arrange(logFC)
ConfirmCD49bvDP <- subset(CD49bvDPLFC, CD49bvDPLFC$ID.gene_name %in% CD49bup) %>% 
  select(ID.gene_name, logFC, adj.P.Val) %>%
  arrange(logFC)


#Heatmap
mat_col <- data.frame(group = v$design)
mat_colours <- list(group = brewer.pal(4, "Set1"))
names(mat_colours$group) <- unique(mat_col)
valueRange <- fit.cont$coefficients[HighConfGenes,] %>%
  subtract(rowMeans(.)) %>%
  range()
myPalette <- colorRampPalette(c("blue", "white", "red"))(101)
myBreaks <- c(seq(min(valueRange), 0, length.out = 50),
              seq(max(valueRange) / 101, max(valueRange), length.out = 50))

pheatmap(mat = fit$coefficients[FinalDNup$ID.gene_id, ] %>%
           subtract(rowMeans(.)),
         color = myPalette,
         breaks = myBreaks,
         border_color = NA, 
         # show_colnames = FALSE, 
         show_rownames = TRUE, 
         cluster_cols = FALSE,
         scale = "none",
         labels_row = FinalDNup$ID.gene_name,
         #cutree_rows = 3,
         drop_levels = TRUE, 
         fontsize = 11)  

#Heatmap
mat_col <- data.frame(group = v$design)
mat_colours <- list(group = brewer.pal(4, "Set1"))
names(mat_colours$group) <- unique(mat_col)
#valueRange <- fit.cont$coefficients[HighConfGenes,] %>%
#subtract(rowMeans(.)) %>%
#range()
myPalette <- colorRampPalette(c("blue", "white", "red"))(101)
myBreaks <- c(seq(min(valueRange), 0, length.out = 51),
              seq(max(valueRange) / 101, max(valueRange), length.out = 35))

#interestedGenes = c("Eomes", "P2rx7", "Cd226")
#merged_DE_mt_geneIDs = merged_DE_mt$ID.gene_name

pheatmap(mat = merged_DE_mt_geneIDs,
         color = myPalette,
         breaks = myBreaks,
         border_color = NA,
         # show_colnames = FALSE,
         show_rownames = TRUE,
         cluster_cols = FALSE,
         scale = "none",
         #labels_col = hello,
         labels_row = merged_DE_mt$ID.gene_name,
         #cutree_rows = 3,
         drop_levels = TRUE,
         fontsize = 11)

####################################################################################################################
#MORE TROUBLESHOOTING#

#######################################################################################################

#Heatmap
mat_col <- data.frame(group = v$design)
mat_colours <- list(group = brewer.pal(4, "Set1"))
names(mat_colours$group) <- unique(mat_col)
valueRange <- fit.cont$coefficients[HighConfGenes,] %>%
  subtract(rowMeans(.)) %>%
  range()
myPalette <- colorRampPalette(c("blue", "white", "red"))(101)
myBreaks <- c(seq(min(valueRange), 0, length.out = 50),
              seq(max(valueRange) / 101, max(valueRange), length.out = 50))


pheatmap(mat = fit$coefficients[FCLAG3CD49b$ID.gene_id, c("LAG-3", "CD49b")],
         color = myPalette,
         breaks = myBreaks,
         border_color = NA,
         # show_colnames = FALSE,
         show_rownames = TRUE,
         cluster_cols = FALSE,
         scale = "none",
         #labels_col = hello,
         labels_row = FCLAG3CD49b$ID.gene_name,
         #cutree_rows = 3,
         drop_levels = TRUE,
         fontsize = 11)


##########  Heatmap of biological replicates- 244 genes of interest #########################

##LogCPM##

filtLogCPM <- cpm(dgeList$counts, log = TRUE)
HighConfGenes_logCPM <- filtLogCPM[HighConfGenes, ]

HighConfGenes_logCPM_DN = HighConfGenes_logCPM[, c("M2_DN", "M3_DN", "M5_DN", "M6_DN")]

HighConfGenes_logCPM_DP = HighConfGenes_logCPM[, c("M2_DP", "M3_DP", "M5_DP", "M6_DP")]

HighConfGenes_logCPM_CD49b = HighConfGenes_logCPM[, c("M2_CD49b", "M3_CD49b", "M5_CD49b", "M6_CD49b")]

HighConfGenes_logCPM_LAG_3 = HighConfGenes_logCPM[, c("M2_LAG_3", "M3_LAG_3", "M5_LAG_3", "M6_LAG_3")]

HighConfGenes_logCPM = cbind(HighConfGenes_logCPM_DN, 
                             HighConfGenes_logCPM_DP,
                             HighConfGenes_logCPM_LAG_3,
                             HighConfGenes_logCPM_CD49b)

myPalette <- colorRampPalette(c("blue", "white", "red"))(101)
myBreaks <- c(seq(min(HighConfGenes_logCPM), 0, length.out = 48),
              seq(max(HighConfGenes_logCPM) / 101, max(HighConfGenes_logCPM), length.out = 40)) #table containing our high confidence genes

pheatmap(mat = HighConfGenes_logCPM,
         color = myPalette,
         breaks = myBreaks,
         border_color = "black",
         #show_colnames = FALSE,
         show_rownames = TRUE, 
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         scale = "none",
         #labels_col = hello,
         labels_row = merged_DE_mt$ID.gene_name,
         #cutree_rows = 3,
         drop_levels = FALSE,
         fontsize = 24)




