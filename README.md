# Spline-DV
A spline-based scRNA-seq method for identifying differentially variable (DV) genes across two experimental conditions.

[Read Preprint Here!](https://doi.org/10.1101/2024.08.08.607086)

# Why to use Spline-DV?
One of the most intuitive ways to evaluate a gene expression change is using Differential Expression (DE) analysis. Traditionally, DE analysis focuses on identifying genes that are up- or down-regulated (increased or decreased expression) between conditions, typically employing a basic mean-difference approach. We propose a paradigm shift that acknowledges the central role of gene expression variability in cellular function and challenges the current dominance of mean-based DE analysis in single-cell studies. We suggest that scRNA-seq data analysis should embrace the role of inherent gene expression variability in defining cellular function and move beyond mean-based approaches. 

# Installation 
## Stable Bioconductor release 
```R
BiocManager::install("Xenon8778/SplineDV")
```
## Developer version 
```R
if (!require("devtools")) install.packages("remotes")
remotes::install_github("Xenon8778/SplineDV")
```

# Loading Input
The expression matrix should be formatted with the genes/features as rows and cells as columns.
The two input matrices for DV_splinefit must be stored as two separate count expression matrices or SingleCellExperiment objects.

## Loading .h5 files
Reading 10X .h5 files as SingleCellExperiment objects using the [DropletUtils](https://bioconductor.org/packages/release/bioc/html/DropletUtils.html) package
```{r loadingSingleExprMat2, eval=FALSE}
# Reading a 10X .h5 file
exprMatrix <- DropletUtils::read10xCounts(samples = 'mySingleCellExperiment.h5', type = 'HDF5')

# Reading two 10X .h5 files 
exprMatrix_1 <- DropletUtils::read10xCounts(samples = 'mySingleCellExperiment_1.h5', type = 'HDF5')
exprMatrix_2 <- DropletUtils::read10xCounts(samples = 'mySingleCellExperiment_2.h5', type = 'HDF5')
```

## Loading Seurat files
Reading .rds files generated using the [Seurat](https://satijalab.org/seurat/) package
```R
# Reading a Seurat RDS file
exprMatrix <- readRDS('mySingleCellExperiment.rds')
exprMatrix <- Seurat::GetAssayData(exprMatrix, layers = 'count') # Extract counts

# Reading two Seurat RDS file
exprMatrix_1 <- readRDS('mySingleCellExperiment_1.rds')
exprMatrix_1 <- Seurat::GetAssayData(exprMatrix_1, layers = 'count') # Extract counts
exprMatrix_2 <- readRDS('mySingleCellExperiment_2.rds')
exprMatrix_2 <- Seurat::GetAssayData(exprMatrix_2, layers = 'count') # Extract counts
```
# Tutorial - Spline-DV
## Loading scRNAseq count example data
The example data is borrowed from an experimental *Nkx2-1* Gene knockout scRNA-seq study by Liebler *et al.* [1]
```R
# Load Data
library(SplineDV)
load(system.file("extdata", "WT_count.rda", package = "SplineDV")) # WT Sample
load(system.file("extdata", "KO_count.rda", package = "SplineDV")) # KO Sample
```

## Running Spline-DV
For the analysis, the test data (X) is always use in contrast with the control data (Y).
```R
DV_res <- DV_splinefit(X = KO_count, Y = WT_count, ncells = 3, ncounts = 200)
head(DV_res)
```
## Visualize Gene Expression statistics
```R
DV_plot(DV_res)
```
![image](https://github.com/user-attachments/assets/4d17a58b-5ce0-4ad4-a65f-f7b7b01f3ebf)

# Tutorial - Spline-HVG
```R
## Loading Data
load(system.file("extdata", "WT_count.rda", package = "SplineDV")) # WT Sample

## Running Spline-HVG
HVG_res <- HVG_splinefit(WT_count, nHVGs = 20, ncells = 3, ncounts = 200)
head(HVG_res)
```
## Visualize Highly Variable Genes
```R
HVG_plot(HVG_res)
```
![image](https://github.com/user-attachments/assets/5942ef6e-cdd8-496c-a316-b3cfa60826e7)

```R
# Plot individual genes
HVG_plot(HVG_res,'Eln') 
```
![image](https://github.com/user-attachments/assets/7862108c-06cf-4769-8b32-b3248a1ce464)

# References
1. Liebler JM, Marconett CN, Juul N, et al. Combinations of differentiation markers distinguish subpopulations of alveolar epithelial cells in adult lung. Am J Physiol Lung Cell Mol Physiol. 2016;310(2):L114-L120. doi:10.1152/ajplung.00337.2015
