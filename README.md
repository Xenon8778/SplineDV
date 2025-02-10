# Spline-DV
A spline-based scRNA-seq method for identifying differentially variable (DV) genes across two experimental conditions.

[Read Preprint Here!](https://doi.org/10.1101/2024.08.08.607086)

# Why to use Spline-DV?
One of the most intuitive ways to evaluate a gene expression change is using Differential Expression (DE) analysis. Traditionally, DE analysis focuses on identifying genes that are up- or down-regulated (increased or decreased expression) between conditions, typically employing a basic mean-difference approach. We propose a paradigm shift that acknowledges the central role of gene expression variability in cellular function and challenges the current dominance of mean-based DE analysis in single-cell studies. We suggest that scRNA-seq data analysis should embrace the role of inherent gene expression variability in defining cellular function and move beyond mean-based approaches. 

# Citation 
Gatlin, V., Gupta, S., Romero, S., Chapkin, R., & Cai, J. J. (2024). Beyond Differential Expression: Embracing Cell-to-Cell Variability in Single-Cell Gene Expression Data Analysis. bioRxiv, 2024-08. doi: [doi.org/10.1101/2024.08.08.607086](doi.org/10.1101/2024.08.08.607086)

# Installation 
### Stable Bioconductor release (currently only in the devel branch)
```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# The following initializes usage of Bioc devel
BiocManager::install(version='devel')

BiocManager::install("Xenon8778/SplineDV")
```
### Developer version 
```R
if (!require("remotes")) install.packages("remotes")
remotes::install_github("Xenon8778/SplineDV")
```

# Loading Input
The expression matrix should be formatted with the genes/features as rows and cells as columns.
The two input matrices for DV_splinefit must be stored as two separate count expression matrices or SingleCellExperiment objects.

# Tutorial
### Input - scRNA-seq expression matrix from two condtions

Spline-DV is designed to compare expression variability between two experimental conditions. To illustrate this, we'll utilize a publicly available mouse lung single-cell RNA-seq dataset from the scRNAseq Bioconductor package by Zilionis et al. [1]. This dataset includes tumor-infiltrating myeloid cells from both tumor and healthy conditions.

```R
sce <- fetchDataset('zilionis-lung-2019','2023-12-20', path="mouse")
sce
```

The expression matrix should be formatted with the genes/features as rows and cells as columns. The two input matrices for splineDV must be stored as two separate count expression matrices or SingleCellExperiment objects. We can isolate the two experimental conditions. To focus solely on changes within a specific cell type and eliminate potential variability arising from shifts in cell type distribution, we select only Neutrophils for further analysis.

```R
# Extract Healthy data
healthyCount <- sce[, sce$Tissue == "healthy"]

# Extract Tumor data
tumorCount <- sce[, sce$Tissue == "tumor"]

# Select a specific cell type - Neutrophils
healthyCount <- healthyCount[, which(healthyCount$`Major cell type` == "Neutrophils")]
print(healthyCount)

tumorCount <- tumorCount[, which(tumorCount$`Major cell type` == "Neutrophils")]
print(tumorCount)
```

Given our isolation of a single cell type, any observed changes in gene expression variability should be attributed to experimental factors. To mitigate the potential impact of background distribution shifts caused by varying cell type and state proportions, we recommend conducting Spline-DV analysis in a cell state/type-specific manner.

### Running Spline-DV

The main function of SplineDV ios the splineDV function. For the analysis, the test data (X) is always use in contrast with the control data (Y). We use smaller QC parameters for the small example data sets to preserve enough cells and genes. **We recommend using the default QC parameters for large data sets.**

```R
dvRes <- splineDV(X = tumorCount, Y = healthyCount)
head(dvRes)
```
The output is a DataFrame containing statistical measures for each gene in the differential variability (DV) analysis. Genes with large vectorDist values exhibit substantial changes in expression variability between the two experimental conditions. Genes with FDR values below 0.05 are considered significantly differentially variable. The Direction column indicates whether the variability increased or decreased. By leveraging this information, we can identify biologically relevant genes that not only demonstrate shifts in mean expression but also alterations in expression variability across conditions.

### Visualize Gene Expression statistics

To visualize the differential variability of genes, we can employ a 3D scatter plot. The `DVPlot` function allows us to visualize the computed spline for each condition, represented by a distinct color, along with the distance vectors of genes from their respective spline. This plot provides a clear representation of the shift in expression variability between the two conditions.

```R
fig <- DVPlot(dvRes, targetgene='Spp1')
fig 
```
## Highly Variable Genes (HVGs) using Spline-HVG

Next, we'll introduce the Spline HVG algorithm, a feature selection technique designed to identify Highly Variable Genes (HVGs) in scRNA-seq data. This algorithm computes the distance from the spline to identify genes with significant variability. To demonstrate this, we'll utilize the same sample dataset of healthy neutrophils from mice lungs used in the previous section.


### Input - scRNA-seq Expression matrix

```R
# Healthy neutrophils data
print(healthyCount)
```

### Running Spline-HVG

The expression matrix should be formatted with the genes/features as rows and cells as columns. The input matrix for `splineHVG` must be stored as a count expression matrices or SingleCellExperiment objects. To focus solely on changes within a specific cell type and eliminate potential variability arising from shifts in cell type distribution, we select only Neutrophils for further analysis. Smaller QC parameters can be used for the small data sets (e.g., cells < 500 and genes < 5000) to preserve enough cells and genes for splineHVG analysis. However, **we recommend using the default QC parameters for large data sets**, if not more stringent.
```R
HVGRes <- splineHVG(healthyCount, nHVGs = 200)
head(HVGRes)
```
The output is a DataFrame containing statistical measures for each gene in the Spline HVG analysis. Genes with large Distance values exhibit significant deviation from expected behavior, as represented by the 3D spline. The top n HVGs can be identified by sorting the DataFrame by Distance in descending order. These highly variable genes, which are biologically relevant, can be utilized in various downstream analyses, such as for dimensional reduction and PCA.

```R
# Extracting HVG Gene list
HVGList <- rownames(HVGRes)[HVGRes$HVG == TRUE]
head(HVGList)
```

### Visualize Highly Variable Genes

To visualize the highly variable genes, we can employ a 3D scatter plot. The `HVGPlot` function allows us to visualize the computed spline, which represents the expected gene behavior. Each dot on the plot represents a single gene. This plot provides a clear representation of HVGs (highlighted in red), which deviate significantly from the 3D spline, indicating a substantial shift in their gene expression.

```R
fig <- HVGPlot(HVGRes)
fig
```

# Loading .h5 files
Reading 10X .h5 files as SingleCellExperiment objects using the [DropletUtils](https://bioconductor.org/packages/release/bioc/html/DropletUtils.html) package
```{r loadingSingleExprMat2, eval=FALSE}
# Reading a 10X .h5 file
exprMatrix <- DropletUtils::read10xCounts(samples = 'mySingleCellExperiment.h5', type = 'HDF5')

# Reading two 10X .h5 files 
exprMatrix_1 <- DropletUtils::read10xCounts(samples = 'mySingleCellExperiment_1.h5', type = 'HDF5')
exprMatrix_2 <- DropletUtils::read10xCounts(samples = 'mySingleCellExperiment_2.h5', type = 'HDF5')
```

# Working with Seurat files
Read files generated using the [Seurat](https://satijalab.org/seurat/) package. Seuratobject is not compatible with Spline-DV as it needs matrix or sparse matrix objects. Hence we extract the count matrix stored within the Seurat object. **The two experimental conditions to be tested should be input as two separate files.** 
```R
# Extract dgCMatrix from two Seurat RDS files
exprMatrix_1 <- Seurat::GetAssayData(seuratObj_1, layers = 'count') # Extract counts
exprMatrix_2 <- Seurat::GetAssayData(seuratObj_1, layers = 'count') # Extract counts
```
### Running Spline-DV
For the analysis, the test data (X) is always used in contrast with the control data (Y).
```R
DV_res <- splineDV(X = exprMatrix_1, Y = exprMatrix_2)
head(DV_res)
```
### Visualize Gene Expression statistics
```R
DVPlot(DV_res, targetgene='Ager')
```
![image](https://github.com/user-attachments/assets/4d17a58b-5ce0-4ad4-a65f-f7b7b01f3ebf)


# Conclusion
SplineDV is a valuable tool for single-cell RNA sequencing (scRNA-seq) analysis, specifically designed to identify differential gene expression variability between conditions. By examining changes in expression dispersion, SplineDV complements traditional differential expression analysis methods. We recommend using SplineDV in conjunction with other tools from Seurat and Bioconductor to gain a more comprehensive understanding of biological processes.

