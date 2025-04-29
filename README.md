# Summary
This workshop offers a balanced blend of theoretical foundation and hands-on experience in single-cell RNA sequencing (scRNA-seq) data analysis using the Seurat R package. It is designed for researchers new to single-cell data or looking to deepen their practical skills in cell clustering, visualization, and interpretation.
# 
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

# Step 1: Load raw matrix and understand the input

    require("DropletUtils")
    require("Seurat")
    require("Matrix")
    raw_dat <- Seurat::Read10X(data.dir = "/usr/local/10x_data/sample_raw_feature_bc_matrix")

> change data.dir to where your cellranger outputs locate

##### Row is gene, column is droplet, the number is how many UMI

    dim(raw_dat)

##### Easy filtration: remove those droplet with 0 or 1 UMI

    raw_dat <- raw_dat[,colSums(raw_dat)>1]
    dim(raw_dat)

# Step 2: Identify the empty droplets

    br.out <- barcodeRanks(raw_dat)

##### barcode rank plot:

    plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total")
    abline(h=metadata(br.out)$knee, col="dodgerblue", lty=2)
    abline(h=metadata(br.out)$inflection, col="forestgreen", lty=2)
    legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"), legend=c("knee", "inflection"))

##### Identify empty drops with Poisson method, and fill them out from the matrix:

    e.out <- emptyDrops(raw_dat, lower=100, niters=10000, ignore=NULL, retain=2*br.out$knee)
    
> niters=10000: change it based on your needs, usually 100000

    e.out <- e.out[!is.na(e.out$PValue),] 
    is.cell <- e.out$FDR <= 0.01
    e.cells <- rownames(e.out)[is.cell]
    head(e.cells)
    filtered_counts <- raw_dat[, e.cells]

# Step 3: Remove ambient RNA

##### Load the R packages:

    library(SingleCellExperiment)
    library(celda)
    library(scuttle)

##### Convert the seurat object into a single cell experiment (SCE) object

    sce <- SingleCellExperiment(assays = list(counts = filtered_counts))

##### Run DecontX

    sce <- decontX(sce)
    sce <- logNormCounts(sce)

##### Convert SCE object back to Seurat object

    seur_obj <- as.Seurat(sce, counts = "counts", data = "logcounts")
    seur_obj[["RNA"]] <- seur_obj[[DefaultAssay(seur_obj)]]
    DefaultAssay(seur_obj) <- "RNA"



