# SplineDV
A spline-based scRNA-seq method for identifying differentially variable (DV) genes across two experimental conditions.

# Installation 
```R
devtools::install_github("Xenon8778/SplineDV")
```
# Example
```R
# Load Data
library(SplineDV)
exprMatrix_WT = read.csv('../data/WTdata.csv', row.names = 1) # WT Sample
exprMatrix_KO = read.csv('../data/KOdata.csv', row.names = 1) # KO Sample
DV_res = DV_splinefit(X = exprMatrix_KO, Y = exprMatrix_WT)
head(DV_res)
```
