#' @title Spline-HVG
#' @description Compute Highly Variable Genes from an scRNAseq expression data.
#' @export HVG_splinefit
#' @author Shreyan Gupta <xenon8778@tamu.edu>
#' @import dplyr
#' @import Seurat
#' @import plotly
#' @importFrom sparseMatrixStats rowSds
#' @importFrom Matrix rowSums rowMeans
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom Biobase matchpt
#' @importFrom stats smooth.spline
#' @importFrom methods as is
#' @return A DataFrame with HVG selection Statistics
#' @param X Count matrix or SeuratObject
#' @param QC A Boolean value (TRUE/FALSE), if TRUE, a quality control is applied over the data.
#' @param nfeatures An integer value. Defines the minimum features (genes) required for a cell to be included in the analysis.
#' @param ncells An integer value. Defines the minimum cells required for a gene to be included in the analysis.
#' @param mt.perc A double value. Defines the minimum percent mitochondrial genes expression required for a cell to be excluded from the analysis.
#' @param degf An integer value. Degrees of freedom for the Spline.
#' @param spar A double value. Smoothing parameter for Spline
#' @param nHVGs An integer value. Number of top Highly Variable Genes (HVGs) to select.
#' @param show.spline A Boolean value (TRUE/FALSE), if TRUE, shows the 3D Spline as Plotly figure.
#' @param use.ndist A Boolean value (TRUE/FALSE), if TRUE, computes the nearest point on spline by nearest-neighbor search (TRUE Recommended). Else, uses the position of the corresponding gene on the spline for distance computation.

HVG_splinefit <- function(X = x, QC = T,
                          nfeatures = 500, ncells = 15, mt.perc = 15,
                          degf = 15, spar = 0.75, nHVGs = 2000, show.spline = FALSE,
                          use.ndist = T){

  # Check if object is SeuratObject
  if (class(X) == "Seurat"){
    adata = X
  } else {
    adata = CreateSeuratObject(CreateAssayObject(X))
  }

  # Checking if object has enough cells
  if(ncol(X) < 50){
    print('Sample has too few cells')
    stop()
  }

  # QC Filtering
  if (QC == T){
    adata = HVG_QC(adata, nfeatures = nfeatures, ncells = ncells,
                   mt.perc = mt.perc)
    print('QC Done')
  }

  # Normalizing
  adata = NormalizeData(adata, normalization.method = 'RC',
                        scale.factor = mean(adata$nCount_RNA, na.rm = T),
                        verbose = F)

  ## Highly Variable Genes
  Dropout = rowSums(GetAssayData(adata, layer = 'counts') == 0)
  Dropout = Dropout/ncol(adata)
  Means = rowMeans(GetAssayData(adata, layer = 'data'))
  SDs = rowSds(GetAssayData(adata, layer = 'data'))
  CV = SDs/Means
  splinefit_df = as.data.frame(Means)
  splinefit_df$CV = CV
  splinefit_df$Dropout = Dropout
  splinefit_df$logMean = log(Means+1)
  splinefit_df$logCV = log(CV+1)
  splinefit_df = splinefit_df %>% arrange(Dropout)
  # splinefit_df = splinefit_df %>% filter(Dropout>0.01 & Dropout<0.99) # Trimming ends

  # Calculate the differences and the squared differences between consecutive elements
  diff_lgu <- diff(splinefit_df$logMean)
  diff_lgcv <- diff(splinefit_df$logCV)
  diff_dropr <- diff(splinefit_df$Dropout)
  diff_squared_sum <- diff_lgu^2 + diff_lgcv^2 + diff_dropr^2

  # Calculate the cumulative sum
  s <- c(0, cumsum(sqrt(diff_squared_sum)))
  xyz = splinefit_df %>% select(logMean,logCV,Dropout)

  fitx <- smooth.spline(s, splinefit_df$logMean, df = degf, spar = spar)
  fity <- smooth.spline(s, splinefit_df$logCV, df = degf, spar = spar)
  fitz <- smooth.spline(s, splinefit_df$Dropout, df = degf, spar = spar)

  # Computing Distances
  xyz1 = cbind(predict(fitx,s)$y,predict(fity,s)$y,predict(fitz,s)$y)
  xyz1 = as.data.frame(xyz1)
  colnames(xyz1) = c('logMean','logCV','Dropout')
  rownames(xyz1) = rownames(splinefit_df)
  euclidean <- function(a, b) sqrt(sum((a - b)^2))
  if (use.ndist == F){
    Dist_HVG = c()
    pb = txtProgressBar(min = 0, max = nrow(xyz), initial = 0, style = 3)
    for (i in 1:nrow(xyz)){
      Dist_HVG = append(Dist_HVG,euclidean(xyz[i,],xyz1[i,]))
      setTxtProgressBar(pb,i)
    }
    close(pb)
  } else {
    Dist_HVG = c()
    dvecx = c()
    dvecy = c()
    dvecz = c()
    nearidx = c()
    df = as.matrix(xyz1)
    pb = txtProgressBar(min = 0, max = nrow(xyz), initial = 0, style = 3)
    for (i in 1:nrow(xyz)){
      p1 = as.matrix(xyz[i,])
      near_point = matchpt(p1, df)
      d1 = p1-df[near_point$index,]
      Dist_HVG = append(Dist_HVG,near_point$distance)
      dvecx = append(dvecx, d1[1])
      dvecy = append(dvecy, d1[2])
      dvecz = append(dvecz, d1[3])
      nearidx = append(nearidx,near_point$index)
      setTxtProgressBar(pb,i)
    }
    close(pb)
  }
  splinefit_df$Distance = Dist_HVG
  splinefit_df$nearidx = nearidx
  splinefit_df$dvecx = dvecx
  splinefit_df$dvecy = dvecy
  splinefit_df$dvecz = dvecz
  splinefit_df$splinex = xyz1$logMean
  splinefit_df$spliney = xyz1$logCV
  splinefit_df$splinez = xyz1$Dropout
  splinefit_df$Distance[splinefit_df$Dropout > 0.95 | splinefit_df$Dropout < 0.01] = 0
  top_HVG = splinefit_df %>% top_n(n = nHVGs, wt = Distance)
  mask = rownames(splinefit_df) %in% rownames(top_HVG)
  splinefit_df$HVG = mask
  splinefit_df = splinefit_df %>% arrange(Dropout)

  # Plotting in 3D
  if(show.spline == TRUE){
    fig = plot_ly(splinefit_df, x = ~logMean, y = ~logCV, z = ~Dropout, mode = 'lines+marker',
                  color = ~HVG, colors = c('grey','red'))
    fig <- fig %>% add_markers(mode = 'marker', type = "scatter3d", opacity = 0.5,
                               marker = list(size = 2))
    fig <- fig %>% add_trace(splinefit_df, x = ~splinex, y = ~spliney, z = ~splinez, type="scatter3d", mode = 'lines+marker',
                             opacity = 1, line = list(width = 5, color = 'black', reverscale = FALSE))
    fig <- fig %>% layout(scene = list(xaxis = list(title = 'log(Mean)'),
                                       yaxis = list(title = 'log(CV)'),
                                       zaxis = list(title = 'Dropout rate')))
    print(fig)
  }

  splinefit_df = splinefit_df %>% filter(nearidx != 1 & nearidx != max(nearidx))
  splinefit_df = splinefit_df %>% arrange(-Distance)

  return(splinefit_df)
}


#' @export HVG_QC
#' @title Quality control
#' @importFrom Matrix rowSums
#' @description QC filter scRNA-seq expression data
#' @author Shreyan Gupta <xenon8778@tamu.edu>
#' @import Seurat
#' @return QC filtered SeuratObject
#' @param nfeatures An integer value. Defines the minimum features (genes) required for a cell to be included in the analysis.
#' @param ncells An integer value. Defines the minimum cells required for a gene to be included in the analysis.
#' @param mt.perc A double value. Defines the minimum percent mitochondrial genes expression required for a cell to be excluded from the analysis.

HVG_QC <- function(X = x, nfeatures = 500, ncells = 15, mt.perc = 15){

  # Check if object is SeuratObject
  if (class(X) == "Seurat"){
    adata = X
  } else {
    adata = CreateSeuratObject(CreateAssayObject(X))
  }

  ## Calculating mitochondrial expression
  if(dim(table(startsWith(rownames(adata),'mt-'))) == 2){
    adata[["percent.mt"]] <- PercentageFeatureSet(adata, pattern = "^mt-")
  }else {
    adata[["percent.mt"]] <- PercentageFeatureSet(adata, pattern = "^MT-")
  }
  adata <- subset(adata, subset = nFeature_RNA > nfeatures & percent.mt < mt.perc)
  selected_f <- rownames(adata)[rowSums(adata) > ncells]
  adata = subset(adata, features = selected_f)

  return(adata)
}



#' @export DV_splinefit
#' @title Spline-DV
#' @description Differential Variability analysis
#' @author Shreyan Gupta <xenon8778@tamu.edu>
#' @import dplyr
#' @importFrom stats pnorm
#' @return A DataFrame with DV Statistics
#' @param X Count matrix or SeuratObject (Test sample)
#' @param Y Count matrix or SeuratObject (Control sample)
#' @param nfeatures An integer value. Defines the minimum features (genes) required for a cell to be included in the analysis.
#' @param ncells An integer value. Defines the minimum cells required for a gene to be included in the analysis.
#' @param mt.perc A double value. Defines the minimum percent mitochondrial genes expression required for a cell to be excluded from the analysis.

DV_splinefit <- function(X = x, Y = y,  nfeatures = 500, ncells = 15,
                         mt.perc = 15) {

  # QC Filtering
  X = HVG_QC(X, nfeatures = nfeatures, ncells = ncells,
             mt.perc = mt.perc)
  Y = HVG_QC(Y, nfeatures = nfeatures, ncells = ncells,
             mt.perc = mt.perc)

  # Intersect gene space of the two data sets
  feat = intersect(rownames(X),rownames(Y))
  X = X[feat,]
  Y = Y[feat,]

  # Computing distances
  print('Computing distances for Data 1')
  res_X = HVG_splinefit(X, QC = F)
  print('Computing distances for Data 2')
  res_Y = HVG_splinefit(Y, QC = F)

  # Intersect gene space of the two data sets
  feat = intersect(rownames(res_X),rownames(res_Y))
  X_filt = res_X[feat,]
  Y_filt = res_Y[feat,]

  # Creating summary table
  df = as.data.frame(list(genes = feat,
                          dist1 = X_filt$Distance,
                          dist2 = Y_filt$Distance,
                          mu1 = X_filt$logMean,
                          mu2 = Y_filt$logMean,
                          CV1 = X_filt$logCV,
                          CV2 = Y_filt$logCV,
                          drop1 = X_filt$Dropout,
                          drop2 = Y_filt$Dropout,
                          X_splinex = X_filt$splinex,
                          X_spliney = X_filt$spliney,
                          X_splinez = X_filt$splinez,
                          Y_splinex = Y_filt$splinex,
                          Y_spliney = Y_filt$spliney,
                          Y_splinez = Y_filt$splinez,
                          X_dvecx = X_filt$dvecx,
                          X_dvecy = X_filt$dvecy,
                          X_dvecz = X_filt$dvecz,
                          Y_dvecx = Y_filt$dvecx,
                          Y_dvecy = Y_filt$dvecy,
                          Y_dvecz = Y_filt$dvecz),
                     row.names = feat)

  # Computing Diff Distance
  DV_res = df %>% as_tibble %>%
    mutate(dist_diff = dist1 - dist2) %>%
    mutate(Vector_Dist = sqrt((X_dvecx-Y_dvecx)^2+(X_dvecy-Y_dvecy)^2+(X_dvecz-Y_dvecz)^2)+1)  %>%
    mutate(Direction = dist_diff/abs(dist_diff)) %>%
    mutate(Z= as.numeric(scale(Spline_Dist))) %>%
    mutate(Pval = 2*pnorm(abs(Z),lower.tail = F)) %>%
    arrange(Pval)

  res_out = DV_res %>% select(genes,mu1,mu2,CV1,CV2,drop1,drop2,dist1,dist2,
                              Vector_Dist,Direction,Pval) %>% as.data.frame()
  output = list(HVG_X = res_X, HVG_Y = res_Y, DV = res_out)
  return(res_out)
}



#' # Test Code
#' so = readRDS('data/Nkx2-1_ENDO.rds')
#' so_WT = subset(so, subset = Batch == 'WT')
#' so_KO = subset(so, subset = Batch == 'KO')
#' DV_res = DV_splinefit(so_KO,so_WT)
