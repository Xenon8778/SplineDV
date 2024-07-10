# SplineDV
A spline-based scRNA-seq method for identifying differentially variable (DV) genes across two experimental conditions.

# Example
```R
so = readRDS('data/Nkx2-1_ENDO.rds')
#' so_WT = subset(so, subset = Batch == 'WT')
#' so_KO = subset(so, subset = Batch == 'KO')
#' DV_res = DV_splinefit(so_KO,so_WT)
```
