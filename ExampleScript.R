#this is the seurat cluster script that will be basically manual


.libPaths("/global/scratch/johnblair/rlibs/3.5")
setwd("/global/scratch/johnblair/Seurat")

library(Seurat)
library(dplyr)
library(cowplot)

# Load the datasets
d120_etoh_CTL_1.data <- Read10X(data.dir = "/global/scratch/johnblair/Seurat/filtered_gene_matrices/120_etoh_CTL_1/filtered_gene_bc_matrices/newref_Ai9")
d120_etoh_TdTom_1.data <- Read10X(data.dir = 
"/global/scratch/johnblair/Seurat/filtered_gene_matrices/120_etoh_tdTom_1/filtered_gene_bc_matrices/newref_Ai9")
d120_etoh_CTL_2.data <- Read10X(data.dir = "/global/scratch/johnblair/Seurat/filtered_gene_matrices/120_etoh_CTL_2/filtered_gene_bc_matrices/newref_Ai9")
d120_etoh_TdTom_2.data <- Read10X(data.dir = 
"/global/scratch/johnblair/Seurat/filtered_gene_matrices/120_etoh_tdTom_2/filtered_gene_bc_matrices/newref_Ai9")
d120_etoh_CTL_3.data <- Read10X(data.dir = "/global/scratch/johnblair/Seurat/filtered_gene_matrices/120_etoh_CTL_3/filtered_gene_bc_matrices/newref_Ai9")
d120_etoh_TdTom_3.data <- Read10X(data.dir = 
"/global/scratch/johnblair/Seurat/filtered_gene_matrices/120_etoh_tdTom_3/filtered_gene_bc_matrices/newref_Ai9")
d50_etoh_CTL_1.data <- Read10X(data.dir = "/global/scratch/johnblair/Seurat/filtered_gene_matrices/50_etoh_CTL_1/filtered_gene_bc_matrices/newref_Ai9")
d50_etoh_TdTom_1.data <- Read10X(data.dir = "/global/scratch/johnblair/Seurat/filtered_gene_matrices/50_etoh_tdTom_1/filtered_gene_bc_matrices/newref_Ai9")
d50_etoh_CTL_2.data <- Read10X(data.dir = "/global/scratch/johnblair/Seurat/filtered_gene_matrices/50_etoh_CTL_2/filtered_gene_bc_matrices/newref_Ai9")
d50_etoh_TdTom_2.data <- Read10X(data.dir = "/global/scratch/johnblair/Seurat/filtered_gene_matrices/50_etoh_tdTom_2/filtered_gene_bc_matrices/newref_Ai9")
d50_etoh_CTL_3.data <- Read10X(data.dir = "/global/scratch/johnblair/Seurat/filtered_gene_matrices/50_etoh_CTL_3/filtered_gene_bc_matrices/newref_Ai9")
d50_etoh_TdTom_3.data <- Read10X(data.dir = "/global/scratch/johnblair/Seurat/filtered_gene_matrices/50_etoh_tdTom_3/filtered_gene_bc_matrices/newref_Ai9")

#create Seurat objects
d120_etoh_CTL_1 <- CreateSeuratObject(raw.data = d120_etoh_CTL_1.data, min.cells = 1, min.genes = 200, project = "d120_all")
d120_etoh_CTL_1@meta.data$geno <- "CTRL"
d120_etoh_CTL_1@meta.data$treatment <- "ETOH"
d120_etoh_CTL_1@meta.data$batch <- "1"
d120_etoh_CTL_1@meta.data$timepoint <- "120"
d120_etoh_CTL_1@meta.data$geno.time <- "CTRL_120"
d120_etoh_CTL_1@meta.data$name <- "CTRL_120_1"

d120_etoh_TdTom_1 <- CreateSeuratObject(raw.data = d120_etoh_TdTom_1.data, min.cells = 1, min.genes = 200, project = "d120_all")
d120_etoh_TdTom_1@meta.data$geno <- "TSC2KO"
d120_etoh_TdTom_1@meta.data$treatment <- "ETOH"
d120_etoh_TdTom_1@meta.data$batch <- "1"
d120_etoh_TdTom_1@meta.data$timepoint <- "120"
d120_etoh_TdTom_1@meta.data$geno.time <- "TSC2KO_120"
d120_etoh_TdTom_1@meta.data$name <- "TSC2KO_120_1"

d120_etoh_CTL_2 <- CreateSeuratObject(raw.data = d120_etoh_CTL_2.data, min.cells = 1, min.genes = 200, project = "d120_all")
d120_etoh_CTL_2@meta.data$geno <- "CTRL"
d120_etoh_CTL_2@meta.data$treatment <- "ETOH"
d120_etoh_CTL_2@meta.data$batch <- "2"
d120_etoh_CTL_2@meta.data$timepoint <- "120"
d120_etoh_CTL_2@meta.data$geno.time <- "CTRL_120"
d120_etoh_CTL_2@meta.data$name <- "CTRL_120_2"

d120_etoh_TdTom_2 <- CreateSeuratObject(raw.data = d120_etoh_TdTom_2.data, min.cells = 1, min.genes = 200, project = "d120_all")
d120_etoh_TdTom_2@meta.data$geno <- "TSC2KO"
d120_etoh_TdTom_2@meta.data$treatment <- "ETOH"
d120_etoh_TdTom_2@meta.data$batch <- "2"
d120_etoh_TdTom_2@meta.data$timepoint <- "120"
d120_etoh_TdTom_2@meta.data$geno.time <- "TSC2KO_120"
d120_etoh_TdTom_2@meta.data$name <- "TSC2KO_120_2"

d120_etoh_CTL_3 <- CreateSeuratObject(raw.data = d120_etoh_CTL_3.data, min.cells = 1, min.genes = 200, project = "d120_all")
d120_etoh_CTL_3@meta.data$geno <- "CTRL"
d120_etoh_CTL_3@meta.data$treatment <- "ETOH"
d120_etoh_CTL_3@meta.data$batch <- "3"
d120_etoh_CTL_3@meta.data$timepoint <- "120"
d120_etoh_CTL_3@meta.data$geno.time <- "CTRL_120"
d120_etoh_CTL_3@meta.data$name <- "CTRL_120_3"

d120_etoh_TdTom_3 <- CreateSeuratObject(raw.data = d120_etoh_TdTom_3.data, min.cells = 1, min.genes = 200, project = "d120_all")
d120_etoh_TdTom_3@meta.data$geno <- "TSC2KO"
d120_etoh_TdTom_3@meta.data$treatment <- "ETOH"
d120_etoh_TdTom_3@meta.data$batch <- "3"
d120_etoh_TdTom_3@meta.data$timepoint <- "120"
d120_etoh_TdTom_3@meta.data$geno.time <- "TSC2KO_120"
d120_etoh_TdTom_3@meta.data$name <- "TSC2KO_120_3"

d50_etoh_CTL_1 <- CreateSeuratObject(raw.data = d50_etoh_CTL_1.data, min.cells = 1, min.genes = 200, project = "d50_all")
d50_etoh_CTL_1@meta.data$geno <- "CTRL"
d50_etoh_CTL_1@meta.data$treatment <- "ETOH"
d50_etoh_CTL_1@meta.data$batch <- "1"
d50_etoh_CTL_1@meta.data$timepoint <- "50"
d50_etoh_CTL_1@meta.data$geno.time <- "CTRL_50"
d50_etoh_CTL_1@meta.data$name <- "CTRL_50_1"

d50_etoh_TdTom_1 <- CreateSeuratObject(raw.data = d50_etoh_TdTom_1.data, min.cells = 1, min.genes = 200, project = "d50_all")
d50_etoh_TdTom_1@meta.data$geno <- "TSC2KO"
d50_etoh_TdTom_1@meta.data$treatment <- "ETOH"
d50_etoh_TdTom_1@meta.data$batch <- "1"
d50_etoh_TdTom_1@meta.data$timepoint <- "50"
d50_etoh_TdTom_1@meta.data$geno.time <- "TSC2KO_50"
d50_etoh_TdTom_1@meta.data$name <- "TSC2KO_50_1"

d50_etoh_CTL_2 <- CreateSeuratObject(raw.data = d50_etoh_CTL_2.data, min.cells = 1, min.genes = 200, project = "d50_all")
d50_etoh_CTL_2@meta.data$geno <- "CTRL"
d50_etoh_CTL_2@meta.data$treatment <- "ETOH"
d50_etoh_CTL_2@meta.data$batch <- "2"
d50_etoh_CTL_2@meta.data$timepoint <- "50"
d50_etoh_CTL_2@meta.data$geno.time <- "CTRL_50"
d50_etoh_CTL_2@meta.data$name <- "CTRL_50_2"

d50_etoh_TdTom_2 <- CreateSeuratObject(raw.data = d50_etoh_TdTom_2.data, min.cells = 1, min.genes = 200, project = "d50_all")
d50_etoh_TdTom_2@meta.data$geno <- "TSC2KO"
d50_etoh_TdTom_2@meta.data$treatment <- "ETOH"
d50_etoh_TdTom_2@meta.data$batch <- "2"
d50_etoh_TdTom_2@meta.data$timepoint <- "50"
d50_etoh_TdTom_2@meta.data$geno.time <- "TSC2KO_50"
d50_etoh_TdTom_2@meta.data$name <- "TSC2KO_50_2"

d50_etoh_CTL_3 <- CreateSeuratObject(raw.data = d50_etoh_CTL_3.data, min.cells = 1, min.genes = 200, project = "d50_all")
d50_etoh_CTL_3@meta.data$geno <- "CTRL"
d50_etoh_CTL_3@meta.data$treatment <- "ETOH"
d50_etoh_CTL_3@meta.data$batch <- "3"
d50_etoh_CTL_3@meta.data$timepoint <- "50"
d50_etoh_CTL_3@meta.data$geno.time <- "CTRL_50"
d50_etoh_CTL_3@meta.data$name <- "CTRL_50_3"

d50_etoh_TdTom_3 <- CreateSeuratObject(raw.data = d50_etoh_TdTom_3.data, min.cells = 1, min.genes = 200, project = "d50_all")
d50_etoh_TdTom_3@meta.data$geno <- "TSC2KO"
d50_etoh_TdTom_3@meta.data$treatment <- "ETOH"
d50_etoh_TdTom_3@meta.data$batch <- "3"
d50_etoh_TdTom_3@meta.data$timepoint <- "50"
d50_etoh_TdTom_3@meta.data$geno.time <- "TSC2KO_50"
d50_etoh_TdTom_3@meta.data$name <- "TSC2KO_50_3"

#alternate combining to make the weird things go away

combctrl_etoh <- MergeSeurat(object1 = d120_etoh_CTL_1, object2 = d120_etoh_CTL_2, add.cell.id1 = "CTRL1", add.cell.id2 = "CTRL2",do.normalize = FALSE, 
project = "2samples")
d120_etoh_CTL <- MergeSeurat(object1 = combctrl_etoh, object2 = d120_etoh_CTL_3, add.cell.id2 = "CTRL3", do.normalize = FALSE, project ="d120_etOH")
d120_etoh_CTL <- NormalizeData(d120_etoh_CTL)
d120_etoh_CTL <- FilterCells(d120_etoh_CTL, subset.names = "nGene", low.thresholds = 500, high.thresholds = Inf)
d120_etoh_CTL <- ScaleData(d120_etoh_CTL, display.progress = F)
d120_etoh_CTL <- FindVariableGenes(d120_etoh_CTL, do.plot = F)

combTdTom_etoh <- MergeSeurat(object1 = d120_etoh_TdTom_1, object2 = d120_etoh_TdTom_2, add.cell.id1 = "TSCKO1", add.cell.id2 = "TSCKO2",do.normalize = 
FALSE, project = "2samples")
d120_etoh_TdTom <- MergeSeurat(object1 = combTdTom_etoh, object2 = d120_etoh_TdTom_3, add.cell.id2 = "TSCKO3", do.normalize = FALSE, project ="d120_etOH")
d120_etoh_TdTom <- NormalizeData(d120_etoh_TdTom)
d120_etoh_TdTom <- FilterCells(d120_etoh_TdTom, subset.names = "nGene", low.thresholds = 500, high.thresholds = Inf)
d120_etoh_TdTom <- ScaleData(d120_etoh_TdTom, display.progress = F)
d120_etoh_TdTom <- FindVariableGenes(d120_etoh_TdTom, do.plot = F)

combctrl_rap <- MergeSeurat(object1 = d50_etoh_CTL_1, object2 = d50_etoh_CTL_2, add.cell.id1 = "CTRL1", add.cell.id2 = "CTRL2",do.normalize = FALSE, project 
= "2samples")
d50_etoh_CTL <- MergeSeurat(object1 = combctrl_rap, object2 = d50_etoh_CTL_3, add.cell.id2 = "CTRL3", do.normalize = FALSE, project ="d50_etoh")
d50_etoh_CTL <- NormalizeData(d50_etoh_CTL)
d50_etoh_CTL <- FilterCells(d50_etoh_CTL, subset.names = "nGene", low.thresholds = 500, high.thresholds = Inf)
d50_etoh_CTL <- ScaleData(d50_etoh_CTL, display.progress = F)
d50_etoh_CTL <- FindVariableGenes(d50_etoh_CTL, do.plot = F)

combTdTom_rap <- MergeSeurat(object1 = d50_etoh_TdTom_1, object2 = d50_etoh_TdTom_2, add.cell.id1 = "TSCKO1", add.cell.id2 = "TSCKO2",do.normalize = FALSE, 
project = "2samples")
d50_etoh_TdTom <- MergeSeurat(object1 = combTdTom_rap, object2 = d50_etoh_TdTom_3, add.cell.id2 = "TSCKO3", do.normalize = FALSE, project ="d50_etoh")
d50_etoh_TdTom <- NormalizeData(d50_etoh_TdTom)
d50_etoh_TdTom <- FilterCells(d50_etoh_TdTom, subset.names = "nGene", low.thresholds = 500, high.thresholds = Inf)
d50_etoh_TdTom <- ScaleData(d50_etoh_TdTom, display.progress = F)
d50_etoh_TdTom <- FindVariableGenes(d50_etoh_TdTom, do.plot = F)

# Gene selection for input to CCA
g.1 <- head(rownames(d120_etoh_CTL@hvg.info), 1000)
g.2 <- head(rownames(d120_etoh_TdTom@hvg.info), 1000)
g.3 <- head(rownames(d50_etoh_CTL@hvg.info), 1000)
g.4 <- head(rownames(d50_etoh_TdTom@hvg.info), 1000)

genes.use <- unique(c(g.1, g.2,g.3,g.4))
genes.use <- intersect(genes.use, rownames(d120_etoh_CTL@scale.data))
genes.use <- intersect(genes.use, rownames(d120_etoh_TdTom@scale.data))
genes.use <- intersect(genes.use, rownames(d50_etoh_CTL@scale.data))
genes.use <- intersect(genes.use, rownames(d50_etoh_TdTom@scale.data))

etoh.list <- list(d120_etoh_CTL, d120_etoh_TdTom, d50_etoh_CTL,d50_etoh_TdTom)
id.list <- list("CTRL_ETOH","TSC2KO_ETOH","CTRL_RAP","TSC2KO_RAP")

#lets run the CCA
etoh.combined <- RunMultiCCA(object.list = etoh.list, add.cell.ids = id.list, genes.use = genes.use, num.cc = 30)

#now we align the subspaces
etoh.combined <- AlignSubspace(etoh.combined, reduction.type = "cca", grouping.var = "geno.time", dims.align = 1:20)

#now we can perform an integrated analysis
# t-SNE and Clustering
etoh.combined <- RunTSNE(etoh.combined, reduction.use = "cca.aligned", dims.use = 1:20, do.fast = T,check_duplicates = FALSE)
etoh.combined <- FindClusters(etoh.combined, reduction.type = "cca.aligned", resolution = 0.6, dims.use = 1:20)

#find all markers for the cluster
etoh.markers <- FindAllMarkers(object = etoh.combined, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
clusters <- etoh.markers %>% group_by(cluster) %>% top_n(50, avg_logFC)
write.csv(clusters, file = "clusters_etoh.csv")


saveRDS(etoh.combined, file = "etoh_combined.rds")
