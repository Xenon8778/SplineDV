## ----loadingSingleExprMat, eval=FALSE-----------------------------------------
#  # i.e. Reading a Seurat RDS file
#  exprMatrix <- readRDS('mySingleCellExperiment.rds')
#  
#  # or .h5 file experiment
#  exprMatrix <- Seurat::Read10X_h5('mySingleCellExperiment.h5')

## ----loadingSeparateExprMat, eval=FALSE---------------------------------------
#  # i.e. Reading a Seurat RDS file
#  exprMatrix_1 <- readRDS('mySingleCellExperiment_1.rds')
#  exprMatrix_2 <- readRDS('mySingleCellExperiment_2.rds')
#  
#  # or .h5 file experiment
#  exprMatrix_1 <- Seurat::Read10X_h5('mySingleCellExperiment_1.h5')
#  exprMatrix_2 <- Seurat::Read10X_h5('mySingleCellExperiment_2.h5')

## ----loadingExampledata-------------------------------------------------------
# Load Data
exprMatrix_WT = read.csv('../data/WTdata.csv', row.names = 1) # WT Sample
exprMatrix_KO = read.csv('../data/KOdata.csv', row.names = 1) # KO Sample

## ----runSplineDV, results='hide'----------------------------------------------
library(SplineDV)
DV_res = DV_splinefit(X = exprMatrix_KO, Y = exprMatrix_WT)

## ----showResults, echo=FALSE--------------------------------------------------
head(DV_res)

## -----------------------------------------------------------------------------
library(MASS)
library(scales)
library(ggplot2)
library(ggpubr)

getdensity <- function(data, gene, col = 'firebrick3',
                       ident = NULL, plot.mu = F){

  dat = log(as.numeric(data[gene,])+1) # Log1p transforming
  dat.line = density(dat)

  df = data.frame(x = dat.line$x, y = dat.line$y, Genotype = ident)
  p1 = ggplot(df, aes(x = x, y = y))+ theme_classic() +
    geom_area(aes(x = x, y = y), color = col,
              fill = alpha(col,0.3), linewidth = 1) +
    labs(x = 'Expression',y = 'Density') +
    ggtitle(paste0(ident,' ',gene, ' Density'))
  if (is.null(legend) == F){
    p1 = p1
    p1
  }
  if (plot.mu == T){
    p1 = p1 + geom_vline(xintercept = log(mean(dat)+1),
                         col = col, lty = 'dashed', linewidth = 1)
  }
  return(p1)
}

gene = 'Ager'
p1 = getdensity(exprMatrix_KO, gene = gene, ident = 'KO')
p2 = getdensity(exprMatrix_WT, gene = gene, col = 'steelblue', ident = 'WT')
ggarrange(p1,p2,ncol = 1, align = 'hv')

## -----------------------------------------------------------------------------
# Loading Data
exprMatrix = read.csv('../data/WTdata.csv', row.names = 1) # WT Sample

## ----runSplineHVG, results='hide', fig.show='animate'-------------------------
HVG_res = HVG_splinefit(exprMatrix, nHVGs = 100)

## ----showHVGResults-----------------------------------------------------------
head(HVG_res)

## ----showHVGList--------------------------------------------------------------
# Extracting HVG Gene list
HVG_list = rownames(HVG_res)[HVG_res$HVG == TRUE]
HVG_list

## -----------------------------------------------------------------------------
date()
sessionInfo()

