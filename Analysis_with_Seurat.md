# Summary

This pipeline is based on the R package Seurat and R version 4.3. You need to start an Rstudio session to run the code.

# Launch Rstudio

To launch Rstudio, you need to use the scRNA-seq.sif in the cluster. (Slide 37-41)

    cp /project/def-sponsor00/scRNA-seq/test_cluster.sh .
    sbatch test_cluster.sh

##### open job_xxx.log and find the line with "ssh -N -L ...". Open another terminal or windows powershell, copy the "ssh" line and run. Please note you need to enter the password twice. 

##### Then find the link with localhost, open this link in your browser to launch rstudio.
    

# R packages needed 

    require("DropletUtils")
    require("Seurat")
    require("Matrix")
    library(SingleCellExperiment)
    library(celda)
    library(scuttle)
    library(scater)
    library(DoubletFinder)
    library(SingleR)
    library(celldex)
    library(tidyr)
    library(dplyr)

# Step 1: Load raw matrix and understand the input (Slide 42)

    require("DropletUtils")
    require("Seurat")
    require("Matrix")
    raw_dat <- Seurat::Read10X(data.dir = "/usr/local/10x_data/sample_raw_feature_bc_matrix")

> change data.dir to where your cellranger outputs locate

##### Row is the gene ID, column is the barcode of each droplet, the numbers are UMI counts

    dim(raw_dat)
    raw_dat[1:5, 1:5]

##### Easy filtration: remove those droplet with 0 or 1 UMI

    raw_dat <- raw_dat[,colSums(raw_dat)>1]
    dim(raw_dat)

# Step 2: Identify the empty droplets (Slide 43-47)

    br.out <- barcodeRanks(raw_dat)

##### barcode rank plot:

    plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total")
    abline(h=metadata(br.out)$knee, col="dodgerblue", lty=2)

##### Identify empty drops with Poisson method, and fill them out from the matrix:

    e.out <- emptyDrops(raw_dat, lower=100, niters=10000, ignore=NULL, retain=2*br.out$knee)
    
> niters=10000: change it based on your needs, usually 100000

    e.out <- e.out[!is.na(e.out$PValue),]   # filter out those not-caculated droplets
    is.cell <- e.out$FDR <= 0.01            # Among all barcodes we call "cells", at most ~1% are expected to be false positives.
    e.cells <- rownames(e.out)[is.cell]
    head(e.cells)
    filtered_counts <- raw_dat[, e.cells]

# Step 3: Remove ambient RNA (Slide 48-50)

##### Load the R packages:

    library(SingleCellExperiment)
    library(celda)
    library(scuttle)

##### Convert the matrix into a single cell experiment (SCE) object

    sce <- SingleCellExperiment(assays = list(counts = filtered_counts))

##### Run DecontX

    sce <- decontX(sce)
    sce <- logNormCounts(sce)

##### Convert SCE object to Seurat object (with "RNA" assay)

    seur_obj <- as.Seurat(sce, counts = "counts", data = "logcounts")
    seur_obj[["RNA"]] <- seur_obj[[DefaultAssay(seur_obj)]]
    DefaultAssay(seur_obj) <- "RNA"

# Step 4: Identify dead/broken cells (Slide 51-53)

    library(SingleCellExperiment)
    library(DropletUtils)
    library(scater)

##### Convert Seurat → SCE

    sce <- as.SingleCellExperiment(seur_obj)

##### Generate mitochondrial gene list

    rowData(sce)$Symbol <- rownames(sce)
    is.mt <- grepl("^MT-", rowData(sce)$Symbol)

##### run QC metrics per cell

    sce <- addPerCellQC(sce, subsets = list(Mito = is.mt))

##### Define Filtering Criteria. Find the threshold through two histograms.

##### Histogram 1: Number of Detected Genes per Cell

    hist(sce$detected, breaks = 200, main = "Number of Genes Detected per Cell", xlab = "Detected Genes", col = "lightgray", border = "white")
    abline(v = 300, col = "red", lwd = 2)
    abline(v = 100, col = "green", lwd = 2)

##### Histogram 2: Percentage of Mitochondrial UMIs per Cell

    hist(sce$subsets_Mito_percent, breaks = 100, main = "Mitochondrial % per Cell", xlab = "Mitochondrial Percent", col = "lightblue")
    abline(v = 30, col = "red", lwd = 2)
    abline(v = 50, col = "green", lwd = 2)

##### Mark dead/broken cells

    cell_filter_detect <- sce$detected < 100
    cell_filter_MT <- sce$subsets_Mito_percent > 30
    sce$discard <- (cell_filter_detect | cell_filter_MT)

> Change 100 and 30 based on the two histograms above
    
##### visualize and inspect

    plotColData(sce, x = "sum", y = "subsets_Mito_percent", colour_by = "discard")

##### filter out dead/broken cells

    sce_filtered <- sce[, !sce$discard]
    sce_filtered <- logNormCounts(sce_filtered)
    
## convert back to seurat object

    seur_filtered <- as.Seurat(sce_filtered, counts = "counts", data = "logcounts")

# Step 5 Expression Normalization (Slide 54-55)

    seur_filtered <- NormalizeData(seur_filtered, normalization.method = "LogNormalize", scale.factor = 10000)

# Step 6: Feature Selection - Identify highly variable genes (Slide 56-59)

    summary(seur_filtered$nFeature_originalexp)
    seur_filtered <- FindVariableFeatures(seur_filtered, selection.method = "vst", nfeatures = 500)

> nfeatures = 2000 in default. Change it based on your sample size and the detected gene numbers per sample.

# Step 7 Dimentionality Reduction (Slide 60-61)

##### Scale the data

    seur_filtered[["percent.mt"]] <- PercentageFeatureSet(seur_filtered, pattern = "^MT-")    # Calculates the percentage of mitochondrial gene expression for each cell
    seur_filtered <- ScaleData(seur_filtered, vars.to.regress = c("percent.mt"))        # center and scale all genes but regressing out the expression of MT genes

##### Run PCA 

    seur_filtered <- RunPCA(seur_filtered, features = VariableFeatures(seur_filtered))
    ElbowPlot(seur_filtered)
    PCs <- 7

> Determine PCs based on ElbowPlot diagram.

# Step 8. Cell clustering (Slide 62-63)

    seur_filtered <- FindNeighbors(seur_filtered, dims = 1:PCs)
    seur_filtered <- FindClusters(seur_filtered, resolution = 0.3)

> resolution=0.3 can be set between 0 and 1. Need fine tuning to find the best number for each project.

# Step 9. Filter out doublets (Slide 64-71)

    library(DoubletFinder)

##### Estimate expected number of detectable doublets

    type.freq <- table(seur_filtered@meta.data$seurat_clusters) / ncol(seur_filtered)        # No. of cells in each cluster
    homotypic.prop <- sum(type.freq^2)            # homotypic proportion — the probability that two randomly selected cells belong to the same cluster
    nEXP <- 0.009 * (ncol(seur_filtered) / 1000) * (1 - homotypic.prop) * ncol(seur_filtered) 
    nEXP

> 0.009: expected number of doublets per 1000 cells (recommended by 10x Genomics)
> ncol(seur_filtered): total number of cells
> (1 - homotypic.prop): adjusts for undetectable homotypic doublets

##### Load Custom DoubletFinder patch for Seurat v5

    source("/usr/local/10x_data/DoubletFinder_custom.R")

##### Parameter sweep to find best pK: neighborhood size - how many nearest neighbors

    sweep.out <- paramSweep_v3(seur_filtered, PCs = 1:PCs)
    sweep.stats <- summarizeSweep(sweep.out)
    plot(sweep.stats[, 2:3])

##### Automatically select best pK

    best.pK <- as.numeric(as.character(sweep.stats[which.max(sweep.stats$BCreal), "pK"]))
    best.pK

##### Inspect sweep results

    sweep.stats[as.numeric(as.character(sweep.stats[, 2])) < 0.075 & sweep.stats[, 3] > 0.8, ]

##### Run DoubletFinder

    seur_filtered <- doubletFinder_v3_custom(
      seur_filtered,
      PCs = 1:PCs,
      pN = 0.05,
      pK = best.pK,
      nExp = round(nEXP)
    )

> Change pN and pK based on sweep results.

##### Visualize doublets on UMAP 

    df_col <- grep("DF.classifications", colnames(seur_filtered@meta.data), value = TRUE)    # Find the column name that contains predicted doublet status
    table(seur_filtered@meta.data$seurat_clusters, seur_filtered@meta.data[[df_col]])
    seur_filtered <- RunUMAP(seur_filtered, dims = 1:PCs)
    DimPlot(seur_filtered, reduction = "umap", group.by = df_col)
    prop.table(table(seur_filtered@meta.data[[df_col]])) * 100

##### Filter out detected doublets (keep only singlets)

    seur_filtered <- subset(seur_filtered, subset = !!as.name(df_col[1]) == "Singlet")
    seur_filtered <- FindNeighbors(seur_filtered, dims = 1:PCs)
    seur_filtered <- FindClusters(seur_filtered, resolution = 0.3)

# Step 10. Run UMAP for visualization (Slide 72-73)

    seur_filtered <- RunUMAP(seur_filtered, dims = 1:PCs)
    umap_cluster <- DimPlot(seur_filtered, reduction = "umap", group.by = "seurat_clusters", label = TRUE, label.size = 5)
    ggsave("umap_cluster_plot.png", plot = umap_cluster, width = 8, height = 6, dpi = 300)

# Step 11. Identify cluster-specific markers (Slide 74-75)

    markers <- FindAllMarkers(seur_filtered, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
    head(markers)
    library(dplyr)
    top3 <- markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
    dotplot <- DotPlot(seur_filtered, features = unique(top3$gene)) + RotatedAxis()
    ggsave("dotplot_top3_markers.png", plot = dotplot, width = 8, height = 6, dpi = 300)

# Step 12. Cell Type Annotation (Slide 76)

    library(SingleR)
    library(celldex)

##### Assign cluster ID to each cell

    Idents(seur_filtered) <- "seurat_clusters"

##### Download the reference
    
    ref <- celldex::HumanPrimaryCellAtlasData()

##### Retrieve the normalized expression profile and run SingleR
    
    expr <- GetAssayData(seur_filtered, slot = "data")
    cluster_labels <- SingleR(test = expr, ref = ref, labels = ref$label.main, clusters = Idents(seur_filtered))

##### Get the cluster ID, cell type label, and pass it back to Seurat object
    
    celltype_labels <- setNames(cluster_labels$labels, rownames(cluster_labels))
    seur_filtered <- RenameIdents(seur_filtered, celltype_labels)
    seur_filtered$predicted_celltype <- as.character(Idents(seur_filtered))

##### Plot the results
    
    umap_annoated <- DimPlot(seur_filtered, group.by = "predicted_celltype", reduction = "umap", label = TRUE, repel = TRUE)
    ggsave("umap_annoated.png", plot = umap_annoated, width = 8, height = 6, dpi = 300)
