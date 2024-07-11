# SplineDV
A spline-based scRNA-seq method for identifying differentially variable (DV) genes across two experimental conditions.

# Installation 
```R
devtools::install_github("Xenon8778/SplineDV")
```
# Tutorial - Spline-DV
## Loading scRNAseq count example data
The example data is borrowed from an experimental *Nkx2-1* Gene knockout scRNA-seq study by Liebler *et al.* ^1^
```R
# Load Data
library(SplineDV)
exprMatrix_WT = read.csv('../data/WTdata.csv', row.names = 1) # WT Sample
exprMatrix_KO = read.csv('../data/KOdata.csv', row.names = 1) # KO Sample
```
## Running Spline-DV
For the analysis, the test data (X) is always use in contrast with the control data (Y).
```R
DV_res = DV_splinefit(X = exprMatrix_KO, Y = exprMatrix_WT)
head(DV_res)
```

# Tutorial - Spline-HVG

```R
## Loading Data
exprMatrix = read.csv('../data/WTdata.csv', row.names = 1) # WT Sample

## Running Spline-HVG
HVG_res = HVG_splinefit(exprMatrix, nHVGs = 100)
head(HVG_res)
```
