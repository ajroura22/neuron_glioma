## ---------------------------

## Purpose of script: normalizing non-tumour and tumour cells from https://doi.org/10.1093/noajnl/vdae005 and data integration per patient
##
## Author: Adria-Jaume Roura
##
## Date Created: 18-02-2025
## 

####################### 0. Libraries and directories ############################################
library(Seurat)
library(future)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
#devtools::install_github("diegoalexespi/pochi", force = T)
library(pochi)

options(width = 300)
options(future.globals.maxSize = 4000 * 1024^2)
options(future.fork.enable = TRUE)
supportsMulticore()
plan("multisession", workers = 20)

# Paths
setwd("/data/ebatlle/aroura/projects/21.GBM_review/")
mainDir <- "/data/ebatlle/aroura/projects/21.GBM_review/"
dir.create(dataDir <- "/data/ebatlle/aroura/projects/21.GBM_review/data/")
dir.create(figDir <- paste0(mainDir, "figs_GBM/"))


custom_colors <- list()
colors1 <- c('#FFC312','#C4E538','#12CBC4','#FDA7DF','#ED4C67','#F79F1F','#A3CB38','#1289A7','#D980FA',
             '#B53471','#EE5A24','#009432','#0652DD','#9980FA','#833471',#EA2027','#006266','#1B1464',
             '#5758BB','#6F1E51')

colors2 <- c("#2c2c54", "#706fd3", "#40407a", "#474787", "#aaa69d", "#ffda79", "#f7f1e3")

set.seed(639245)
colors3 <- brewer.pal(9, "Set1") 
colors4 <- brewer.pal(12, "Set3")
custom_colors$discrete <- c(colors1, colors2, colors3)
custom_colors$cell_cycle <- setNames(c('mediumseagreen', 'gold2', 'darkorchid1'),c('G1',      'S',       'G2M'))

col_vector<-c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', 'darkblue', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080')
final_colors <- c(col_vector, colors4[-c(2,4)])
final_colors2 <- c(final_colors, brewer.pal(8, "Dark2"))
final_colors2 <- c(final_colors2, "red", "black", "blue", "yellow", "grey")

median.stat <- function(x){
  out <- quantile(x, probs = c(.5))
  names(out) <- c("ymed")
  return(out) 
}


####################### 1. Merge TME and tumour cells ############################################
seu_tumour <- DietSeurat(object = readRDS(paste0(dataDir, "tumorcells.integrated.rds")), assays="RNA")
seu_TME <- DietSeurat(object = readRDS(paste0(dataDir, "nontumorcells.integrated.rds")), assays="RNA")
merged <- merge(seu_tumour, y = seu_TME)

# Complete the metadata based on author's separated annotations
merged$annotation <- paste(merged$celltype, merged$Suva_A, sep = "_")
merged$annotation <- gsub("_NA|NA_", "", merged$annotation)

# Only GBM samples
merged <- subset(merged, subset = patient %in% c("Gr4-GBM-1", "Gr4-GBM-2", "Gr4-GBM-3", "Gr4-GBM-4"))

# Filtering genes which are expressed in fewer than X cells
f_genes <- rownames(merged)[Matrix::rowSums(merged) > 5]
merged <- subset(merged, features = f_genes)

# Run UMAP
#merged <- RunUMAP(merged, reduction.name = 'UMAP',reduction = 'pca', dims = 1:15, n.components = 2)
#Idents(merged) <- "annotation"
#DimPlot(merged)
#saveRDS(merged, paste0(dataDir, "merged_tumor_nontumor_SCT.rds"))


####################### Data integration per patient on merged objects ############################################
merged[["RNA"]] <- split(merged[["RNA"]], f = merged$patient)
merged <- NormalizeData(merged)
merged <- FindVariableFeatures(merged)
merged <- ScaleData(merged)
merged <- RunPCA(merged)

# Check the number of PC to retain
p <- tibble(PC = 1:50, stdev = merged@reductions$pca@stdev) %>%
  ggplot(aes(PC, stdev)) +
  geom_point() +
  geom_vline(xintercept = 10, color = 'darkgreen',linetype = "dashed" ) +
  geom_vline(xintercept = 15, color = 'darkgreen',linetype = "dashed") +
  theme_bw() +
  labs(x = 'PCs', y = 'Standard deviation') +
  theme(axis.text = element_text(size = 12)) # 15 looks fine
p

# KNN, clustering and dimensionality reduction
merged <- FindNeighbors(merged, reduction = "pca",
                        dims = 1:15, #dimensions of reduction to use as input
                        k.param = 20, #defines k for the k-nearest neighbor (lower K par, more clusters)
                        prune.SNN = 1/15) #defines a cutoff for acceptable Jaccard index when computing neighborhood overlap for the SNN construction

merged <- FindClusters(merged,verbose = FALSE, resolution = c(0.6, 0.8, 1, 1.2), )
merged <- RunUMAP(merged, dims = 1:15, reduction = "pca", reduction.name = "umap.unintegrated")


pdf(paste0(figDir, "UMAP_unintegrated.pdf"), height = 6, width = 14, onefile = T)
DimPlot(merged, reduction = "umap.unintegrated", group.by = c("patient", "annotation"), cols = final_colors2)
dev.off()

png(paste0(figDir, "UMAP_unintegrated.png"), height = 600, width = 1400, res = 150)
DimPlot(merged, reduction = "umap.unintegrated", group.by = c("patient", "annotation"), cols = final_colors2)
dev.off()

# Anchor-based CCA integration
obj <- IntegrateLayers(object = merged, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca", verbose = FALSE)
obj <- FindNeighbors(obj, reduction = "integrated.cca", dims = 1:15)
obj <- FindClusters(obj, resolution = c(1), cluster.name = "cca_clusters")
obj <- RunUMAP(obj, reduction = "integrated.cca", dims = 1:15, reduction.name = "umap.cca")

pdf(paste0(figDir, "UMAP_patient_integrated.pdf"), height = 6, width = 14, onefile = T)
DimPlot(obj, reduction = "umap.cca", group.by = c("patient", "annotation"), cols = final_colors2)
dev.off()

png(paste0(figDir, "UMAP_patient_integrated.png"), height = 600, width = 1400, res = 150)
DimPlot(obj, reduction = "umap.cca", group.by = c("patient", "annotation"), cols = final_colors2)
dev.off()

pdf(paste0(figDir, "UMAP_patient_integrated_clusters.pdf"), height = 6, width = 7, onefile = T)
DimPlot(obj, reduction = "umap.cca", group.by = c("cca_clusters"), cols = final_colors2)
dev.off()

# Join layers, which were splitting by patient (counts and data matrices)
obj <- JoinLayers(obj)


####################### Data imputation (MAGIC) ############################################
library(reticulate)
use_virtualenv("/home/aroura/R/x86_64-pc-linux-gnu-library/4.3/envMagic")
keep_rows <- rowSums(obj) > 1
tmp <- obj[keep_rows,]
tmp <- Rmagic::magic(tmp, n.jobs = 20) # imputing all genes
obj@assays$MAGIC_RNA <- tmp@assays$MAGIC_RNA
rm(tmp)
DefaultAssay(obj) <- "MAGIC_RNA"

saveRDS(obj, paste0(dataDir, "obj_cca_integrated_patient_GBM.rds"))
#obj <- readRDS(paste0(dataDir, "obj_cca_integrated_patient_GBM.rds"))


####################### Scores from senders and receivers from CellChat ############################################
signaling_pathway <- list(
  TULP = list(sender = "TUB", receiver = "MERTK"),
  THBS = list(sender = "THBS2", receiver = "CD47"),
  TENASCIN = list(sender = "ITGA9", receiver = "ITGB1"),
  SLITRK = list(sender = "SLITRK5", receiver = "PTPRS"),
  SEMA6 = list(sender = "SEMA6A", receiver = "PLXNA4, PLXNA2"),
  SEMA4 = list(sender = "SEMA4D", receiver = "PLXNB1"),
  SEMA3 = list(sender = c("SEMA3A", "SEMA3D"), receiver = c("NRP1, PLXNA2, NRP2, PLXNA4")),
  PROS = list(sender = "PROS1", receiver = "TYRO3"),
  NRG = list(sender = c("NRG3", "NRG1", "NRG2"), receiver = "ERBB4"),
  NG4 = list(sender = "LRRC4B", receiver = "PTPRF"),
  NEGR = list(sender = "NEGR1", receiver = "NEGR1"),
  ncWNT = list(sender = c("WNT5A", "WNT5B"), receiver = "FZD3"),
  IGF = list(sender = "IGF1", receiver = "IGF1R"),
  GABA_B = list(sender = c("GAD1","SLC6A1", "GAD2", "SLC6A6"), receiver = c("GABBR2", "GABBR1")),
  FLRT = list(sender = "FLRT2", receiver = "FLRT2"),
  FGF = list(sender = "FGF9", receiver = c("FGFR1, FGFR2, FGFR3")),
  EPHB = list(sender = "EFNB2", receiver = c("EPHA4, EPHB1")),
  EPHA = list(sender = "EFNA2", receiver = c("EPHA3, EPHA4, EPHA7, EPHA5, EPHB2")),
  DHEAS = list(sender = "SULT2B1", receiver = c("PPARA, PPARD, PPARG")),
  DHEA = list(sender = "STS", receiver = "PPARA"),
  Cholesterol = list(sender = "LIPA", receiver = "RORA"),
  ACTIVIN = list(sender = "INHBA", receiver = c("ACVR1B","ACVR2A")),
  Testosterone = list(sender = "HSD17B12", receiver = "AR"))


# Senders
genesets_senders <- lapply(signaling_pathway, function(pathway) {
  sender_genes <- unlist(strsplit(pathway$sender, ", "))
  receiver_genes <- unlist(strsplit(pathway$receiver, ", "))
  sender_genes
})
names(genesets_senders) <- paste0(names(genesets_senders), "_S")
genesets_senders_f <- Filter(function(genes) length(genes) > 1, genesets_senders) # only pathways to calculate a signature score
genesets_senders_f2 <- Filter(function(genes) length(genes) == 1, genesets_senders)

# Receivers
genesets_receivers <- lapply(signaling_pathway, function(pathway) {
  sender_genes <- unlist(strsplit(pathway$sender, ", "))
  receiver_genes <- unlist(strsplit(pathway$receiver, ", "))
  receiver_genes
})
names(genesets_receivers) <- paste0(names(genesets_receivers), "_R")
genesets_receivers_f <- Filter(function(genes) length(genes) > 1, genesets_receivers) # only pathways to calculate a signature score
genesets_receivers_f2 <- Filter(function(genes) length(genes) == 1, genesets_receivers) 

obj$annotation <- factor(obj$annotation, levels = names(table(obj$annotation)))
obj <- RunAUCell(obj,
                 assay = "MAGIC_RNA",
                 slot = "data",
                 genesets = c(genesets_senders_f, genesets_receivers_f),
                 ranking.save = TRUE,
                 ranking.key = NULL,
                 normAUC = TRUE,
                 aucMaxRank = 0.05,
                 verbose = TRUE,
                 auc_assay_name = "AUC")

AUC_df <- as.data.frame(t(as.matrix(obj@assays$AUC@data)))
colnames(AUC_df) <- names(c(genesets_senders_f, genesets_receivers_f))
obj <- AddMetaData(obj, metadata = AUC_df)

# 1) UMAP-violin
plots <- FeaturePlot(obj, 
                     reduction = "umap.cca", 
                     dims = 1:2, 
                     pt.size = 0.1, 
                     features = names(c(genesets_senders_f, genesets_receivers_f)),
                     order = TRUE, 
                     combine = F) 

#plots <- lapply(plots, function(p) {
#  p + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())
#})

Idents(obj) <- "annotation"
plots2 <- VlnPlot(obj, assay = "MAGIC_RNA",
                  pt.size = 0, 
                  features = names(c(genesets_senders_f, genesets_receivers_f)),
                  combine = FALSE, 
                  cols = final_colors2)

plots2 <- lapply(plots2, function(p) {
  p + stat_summary(fun = median.stat, geom='point', size = 2, colour = c("black")) + theme(legend.position = "none", axis.title = element_blank(), axis.text = element_text(size=9))
})

combined_vector <- c(rbind(plots, plots2))

pdf(paste0(figDir, "UMAP_Violins_senders_receivers.pdf"), width = 8, height = 40)
CombinePlots(plots = combined_vector, ncol = 2) + theme(legend.position = "none", axis.title.x = element_text(face = "plain", colour = "black", size = 0))
dev.off()

png(paste0(figDir, "UMAP_Violins_senders_receivers.png"), width = 800, height = 4000)
CombinePlots(plots = combined_vector, ncol = 2) + theme(legend.position = "none", axis.title.x = element_text(face = "plain", colour = "black", size = 0))
dev.off()

# 2) UMAP-violin, rest of candidates
plots <- FeaturePlot(obj, 
                     reduction = "umap.cca", 
                     dims = 1:2, 
                     pt.size = 0.1, 
                     features = as.character(unlist(c(genesets_senders_f2, genesets_receivers_f2))),
                     order = TRUE, 
                     combine = F) 

#plots <- lapply(plots, function(p) {
#  p + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())
#})

Idents(obj) <- "annotation"
plots2 <- VlnPlot(obj, assay = "MAGIC_RNA",
                  pt.size = 0, 
                  features = as.character(unlist(c(genesets_senders_f2, genesets_receivers_f2))),
                  combine = FALSE, 
                  cols = final_colors2)

plots2 <- lapply(plots2, function(p) {
  p +     stat_summary(fun = median.stat, geom='point', size = 2, colour = c("black")) + theme(legend.position = "none", axis.title = element_blank(), axis.text = element_text(size=9))
})

combined_vector <- c(rbind(plots, plots2))

pdf(paste0(figDir, "UMAP_Violins_senders_receivers2.pdf"), width = 8, height = 110)
CombinePlots(plots = combined_vector, ncol = 2) + theme(legend.position = "none", axis.title.x = element_text(face = "plain", colour = "black", size = 0))
dev.off()

png(paste0(figDir, "UMAP_Violins_senders_receivers2.png"), width = 800, height = 8000)
CombinePlots(plots = combined_vector, ncol = 2) + theme(legend.position = "none", axis.title.x = element_text(face = "plain", colour = "black", size = 0))
dev.off()


# Stacked SENDERS
pdf(paste0(figDir, "UMAP_Violins_senders_stacked.pdf"), width = 6, height = 15)
VlnPlot(obj, pt.size = 0, 
        features = c("TUB", "THBS2", "SLITRK5", "SEMA6A", "SEMA4D", "SEMA3_S", "PROS1", "NRG3", "LRRC4B", "NEGR1", "ncWNT_S", "IGF1",
                     "GABA_B_S", "FLRT2", "FGF9", "EFNB2", "EFNA2", "SULT2B1", "STS", "LIPA", "INHBA", "HSD17B12"),
        split.by = "annotation", group.by = "annotation",
        cols = final_colors2, 
        stack = TRUE, 
        sort = F, 
        flip = T) +
  stat_summary(fun = median.stat, geom='point', size = 2, colour = c("black")) +
  theme(legend.position = "none", plot.title = element_text(size = 10))
dev.off()


# Stacked SENDERS filtered
pdf(paste0(figDir, "UMAP_Violins_senders_stacked_filtered.pdf"), width = 5, height = 10)
VlnPlot(obj, pt.size = 0, 
        features = c("TUB", "SLITRK5", "SEMA4D", "SEMA3_S", "IGF1", "FGF9", "EFNB2", "SULT2B1", "STS", "LIPA", "HSD17B12"),
        split.by = "annotation", group.by = "annotation",
        cols = final_colors2, 
        stack = TRUE, 
        sort = F, 
        flip = T) +
  stat_summary(fun = median.stat, geom='point', size = 2, colour = c("black")) +
  theme(legend.position = "none", plot.title = element_text(size = 10), axis.title.x = element_blank())
dev.off()

png(paste0(figDir, "UMAP_Violins_senders_stacked_filtered.png"), width = 500, height = 1000)
VlnPlot(obj, pt.size = 0, 
        features = c("TUB", "SLITRK5", "SEMA4D", "SEMA3_S", "IGF1", "FGF9", "EFNB2", "SULT2B1", "STS", "LIPA", "HSD17B12"),
        split.by = "annotation", group.by = "annotation",
        cols = final_colors2[c(4,5,1,7,9,11,2,3,6,8,10,12,13)], 
        stack = TRUE, 
        sort = F, 
        flip = T) +
  stat_summary(fun = median.stat, geom='point', size = 2, colour = c("black")) +
  theme(axis.text.y = element_text(size=14),legend.position = "none", plot.title = element_text(size = 10), axis.title.x = element_blank()) +
  annotate("rect", xmin = 0.5, xmax = 2.5, ymin = -Inf, ymax = Inf, fill = "orange1", alpha = 0.2) + 
  annotate("rect", xmin = 2.5, xmax = 6.5, ymin = -Inf, ymax = Inf, fill = "darkorchid1", alpha = 0.2) +  
  annotate("rect", xmin = 6.5, xmax = 13.5, ymin = -Inf, ymax = Inf, fill = "darkgreen", alpha = 0.2)
dev.off()



signaling_pathway <- list(
  TULP = list(sender = "TUB", receiver = "MERTK"),
  THBS = list(sender = "THBS2", receiver = "CD47"),
  TENASCIN = list(sender = "THBS2", receiver = "ITGB1"),
  SLITRK = list(sender = "SLITRK5", receiver = "PTPRS"),
  SEMA6 = list(sender = "SEMA6A", receiver = "PLXNA4, PLXNA2"),
  SEMA4 = list(sender = "SEMA4D", receiver = "PLXNB1"),
  SEMA3 = list(sender = c("SEMA3A", "SEMA3D"), receiver = c("NRP1, PLXNA2, NRP2, PLXNA4")),
  PROS = list(sender = "PROS1", receiver = "TYRO3"),
  NRG = list(sender = c("NRG3", "NRG1", "NRG2"), receiver = "ERBB4"),
  NG4 = list(sender = "LRRC4B", receiver = "PTPRF"),
  NEGR = list(sender = "NEGR1", receiver = "NEGR1"),
  ncWNT = list(sender = c("WNT5A", "WNT5B"), receiver = "FZD3"),
  IGF = list(sender = "IGF1", receiver = "IGF1R"),
  GABA_B = list(sender = c("GAD1","SLC6A1", "GAD2", "SLC6A6"), receiver = c("GABBR2", "GABBR1")),
  FLRT = list(sender = "FLRT2", receiver = "FLRT2"),
  FGF = list(sender = "FGF9", receiver = c("FGFR1, FGFR2, FGFR3")),
  EPHB = list(sender = "EFNB2", receiver = c("EPHA4, EPHB1")),
  EPHA = list(sender = "EFNA2", receiver = c("EPHA3, EPHA4, EPHA7, EPHA5, EPHB2")),
  DHEAS = list(sender = "SULT2B1", receiver = c("PPARA, PPARD, PPARG")),
  DHEA = list(sender = "STS", receiver = "PPARA"),
  Cholesterol = list(sender = "LIPA", receiver = "RORA"),
  ACTIVIN = list(sender = "INHBA", receiver = c("ACVR1B","ACVR2A")),
  Testosterone = list(sender = "HSD17B12", receiver = "AR"))


pdf(paste0(figDir, "UMAP_Violins_receivers_stacked.pdf"), width = 6, height = 15)
VlnPlot(obj, pt.size = 0, 
        features = c("MERTK", "CD47", "ITGB1", "PTPRS", "SEMA6_R", "PLXNB1", "SEMA3_R", "TYRO3", "ERBB4", "PTPRF", "NEGR1", "FZD3",
                     "IGF1R", "GABA_B_R", "FLRT2", "FGF_R", "EPHB_R", "EPHA_R", "DHEAS_R", "PPARA", "RORA", "ACTIVIN_R", "AR"),
        split.by = "annotation", group.by = "annotation",
        cols = final_colors2, 
        stack = TRUE, 
        sort = F, 
        flip = T) +
  stat_summary(fun = median.stat, geom='point', size = 2, colour = c("black")) +
  theme(legend.position = "none", plot.title = element_text(size = 10))
dev.off()



pdf(paste0(figDir, "UMAP_Violins_receivers_stacked_filtered.pdf"), width = 5, height = 10)
VlnPlot(obj, pt.size = 0, 
        features = c("MERTK", "PTPRS", "PLXNB1", "SEMA3_R", "IGF1R", "FGF_R", "EPHB_R", "DHEAS_R", "PPARA", "RORA", "AR"),
        split.by = "annotation", group.by = "annotation",
        cols = final_colors2, 
        stack = TRUE, 
        sort = F, 
        flip = T) +
  stat_summary(fun = median.stat, geom='point', size = 2, colour = c("black")) +
  theme(legend.position = "none", plot.title = element_text(size = 9), axis.title.x = element_blank())
dev.off()


png(paste0(figDir, "UMAP_Violins_receivers_stacked_filtered.png"), width = 500, height = 1000)
VlnPlot(obj, pt.size = 0, 
        features = c("MERTK", "PTPRS", "PLXNB1", "SEMA3_R", "IGF1R", "FGF_R", "EPHB_R", "DHEAS_R", "PPARA", "RORA", "AR"),
        split.by = "annotation", group.by = "annotation",
        cols = final_colors2, 
        stack = TRUE, 
        sort = F, 
        flip = T) +
  stat_summary(fun = median.stat, geom='point', size = 2, colour = c("black")) +
  theme(axis.text.y = element_text(size=14),legend.position = "none", plot.title = element_text(size = 10), axis.title.x = element_blank())
dev.off()

#Order
obj$annotation <- factor(obj$annotation, levels = c("ExN", "InN", "AC", "MES", "NPC", "OPC", "Astro", "Endo", "Mac", "MG", "Oligo", "T_cell", "unclassified"))
png(paste0(figDir, "UMAP_Violins_receivers_stacked_filtered.png"), width = 500, height = 1000)
VlnPlot(obj, pt.size = 0, 
        features = c("MERTK", "PTPRS", "PLXNB1", "SEMA3_R", "IGF1R", "FGF_R", "EPHB_R", "DHEAS_R", "PPARA", "RORA", "AR"),
        split.by = "annotation", group.by = "annotation",
        cols = final_colors2[c(4,5,1,7,9,11,2,3,6,8,10,12,13)], 
        stack = TRUE, 
        sort = F, 
        flip = T) +
  stat_summary(fun = median.stat, geom='point', size = 2, colour = c("black")) +
  theme(axis.text.y = element_text(size=14),legend.position = "none", plot.title = element_text(size = 10), axis.title.x = element_blank()) +
  annotate("rect", xmin = 0.5, xmax = 2.5, ymin = -Inf, ymax = Inf, fill = "orange1", alpha = 0.2) + 
  annotate("rect", xmin = 2.5, xmax = 6.5, ymin = -Inf, ymax = Inf, fill = "darkorchid1", alpha = 0.2) +  
  annotate("rect", xmin = 6.5, xmax = 13.5, ymin = -Inf, ymax = Inf, fill = "darkgreen", alpha = 0.2)
dev.off()



# Final UMAPs
plots <- FeaturePlot(obj, 
                     reduction = "umap.cca", 
                     dims = 1:2, 
                     pt.size = 0.3, 
                     features = c("TUB", "SLITRK5", "SEMA4D", "SEMA3_S", "IGF1", "FGF9", "EFNB2", "SULT2B1", "STS", "LIPA", "HSD17B12",
                                  "MERTK", "PTPRS", "PLXNB1", "SEMA3_R", "IGF1R", "FGF_R", "EPHB_R", "DHEAS_R", "PPARA", "RORA", "AR"),
                     order = TRUE, combine = F) 

plots <- lapply(plots, function(p) {
  p + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank()) + NoAxes() + labs(title="")
})

genes <- c("TUB", "SLITRK5", "SEMA4D", "SEMA3_S", "IGF1", "FGF9", "EFNB2", "SULT2B1", "STS", "LIPA", "HSD17B12",
           "MERTK", "PTPRS", "PLXNB1", "SEMA3_R", "IGF1R", "FGF_R", "EPHB_R", "DHEAS_R", "PPARA", "RORA", "AR")
for (i in seq_along(plots)) {
  # Generate a unique file name for each plot
  file_name <- paste0(figDir, "UMAP_plot_", genes[i], ".png")
  
  # Save the current plot as a PNG
  png(file_name, width = 600, height = 600)
  print(plots[[i]])  # Render the plot
  dev.off()  # Close the device
}






saveRDS(obj, paste0(dataDir, "obj_cca_integrated_patient_051224.rds"))
obj <- readRDS(paste0(dataDir, "obj_cca_integrated_patient_051224.rds"))






#### END ####



