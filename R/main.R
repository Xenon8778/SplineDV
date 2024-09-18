#' @title Spline-HVG
#' @description Compute Highly Variable Genes from an scRNAseq expression data.
#' @export HVG_splinefit
#' @author Shreyan Gupta <xenon8778@tamu.edu>
#' @import dplyr
#' @import SingleCellExperiment
#' @import SummarizedExperiment
#' @importFrom sparseMatrixStats rowSds
#' @importFrom Matrix rowSums rowMeans
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom Biobase matchpt
#' @importFrom stats smooth.spline predict
#' @importFrom methods as is
#' @return A data.frame with HVG selection Statistics
#' @param X Count matrix or SingleCellExperiment
#' @param QC A Boolean value (TRUE/FALSE), if TRUE, a quality control is applied over the data.
#' @param ncounts An integer value. Defines the minimum reads required for a cell to be included in the analysis.
#' @param ncells An integer value. Defines the minimum cells required for a gene to be included in the analysis.
#' @param mtPerc A double value. Defines the minimum percent mitochondrial genes expression required for a cell to be excluded from the analysis.
#' @param degf An integer value. Degrees of freedom for the Spline.
#' @param spar A double value. Smoothing parameter for Spline
#' @param nHVGs An integer value. Number of top Highly Variable Genes (HVGs) to select.
#' @param use.ndist A Boolean value (TRUE/FALSE), if TRUE, computes the nearest point on spline by nearest-neighbor search (TRUE Recommended). Else, uses the position of the corresponding gene on the spline for distance computation.
#' @examples
#' # example code
#' ## Load Data
#' load(system.file("extdata", "WT_count.rda", package = "SplineDV")) # WT Sample
#' HVG_res <- HVG_splinefit(WT_count, nHVGs = 100)

HVG_splinefit <- function(X, QC=TRUE,
                          ncounts=500, ncells=15, mtPerc=15,
                          degf=15, spar=0.75, nHVGs=2000, use.ndist=TRUE){

  # Check if object is SingleCellExperiment
  if (is(X,"SingleCellExperiment")){
    adata <- X
  } else {
    adata <- SingleCellExperiment(list(counts=as.matrix(X)))
  }

  # Checking if object has enough cells
  if(ncol(X) < 50){
    stop('Sample has too few cells')
  }

  # QC Filtering
  if (QC == TRUE){
    adata <- HVG_QC(adata, ncounts=ncounts, ncells=ncells,
                   mtPerc=mtPerc)
    message('QC Done')
  }

  # Checking if object has enough genes
  if(nrow(adata) < 50){
    stop('Sample has too few genes')
  }

  # Normalizing
  SummarizedExperiment::assay(adata, 'data') <- t(t(SummarizedExperiment::assay(adata, 'counts')) / colSums(SummarizedExperiment::assay(adata, 'counts'))) * mean(colSums(SummarizedExperiment::assay(adata, 'counts')), na.rm=TRUE)

  ## Highly Variable Genes
  Dropout <- rowSums(SummarizedExperiment::assay(adata, 'counts') == 0)
  Dropout <- Dropout / ncol(adata)
  Means <- rowMeans(SummarizedExperiment::assay(adata, 'data'))
  SDs <- rowSds(SummarizedExperiment::assay(adata, 'data'))
  CV <- SDs / Means
  splinefit_df <- as.data.frame(Means)
  splinefit_df$CV <- CV
  splinefit_df$Dropout <- Dropout
  splinefit_df$logMean <- log(Means+1)
  splinefit_df$logCV <- log(CV+1)
  splinefit_df <- splinefit_df %>% arrange(Dropout)
  # splinefit_df <- splinefit_df %>% filter(Dropout>0.01 & Dropout<0.99) # Trimming ends

  # Calculate the differences and the squared differences between consecutive elements
  diff_lgu <- diff(splinefit_df$logMean)
  diff_lgcv <- diff(splinefit_df$logCV)
  diff_dropr <- diff(splinefit_df$Dropout)
  diff_squared_sum <- diff_lgu^2 + diff_lgcv^2 + diff_dropr^2

  # Calculate the cumulative sum
  s <- c(0, cumsum(sqrt(diff_squared_sum)))
  xyz <- splinefit_df %>% dplyr::select(logMean,logCV,Dropout)

  fitx <- smooth.spline(s, splinefit_df$logMean, df=degf, spar=spar)
  fity <- smooth.spline(s, splinefit_df$logCV, df=degf, spar=spar)
  fitz <- smooth.spline(s, splinefit_df$Dropout, df=degf, spar=spar)

  # Computing Distances
  xyz1 <- cbind(predict(fitx,s)$y, predict(fity,s)$y, predict(fitz,s)$y)
  xyz1 <- as.data.frame(xyz1)
  colnames(xyz1) <- c('logMean','logCV','Dropout')
  rownames(xyz1) <- rownames(splinefit_df)
  euclidean <- function(a, b) sqrt(sum((a - b) ^ 2))
  if (!use.ndist){
    df <- cbind(xyz, xyz1)
    p1 <- as.list(as.data.frame(t(df)))
    Dist_df <- as.data.frame(
      do.call('rbind', lapply(p1, FUN=function(p){
        d <- euclidean(c(p[1:3]), c(p[4:6]))# Computing nearest point on spline
        return(d)
        }))
    )
    names(Dist_df) <- c('Dist_HVG')
    splinefit_df$Distance <- Dist_df$Dist_HVG
  } else {
    df <- as.matrix(xyz1)
    p1 <- as.list(as.data.frame(t(xyz)))
    Dist_df <- as.data.frame(
      do.call('rbind', lapply(p1, FUN=function(p){
        p <- t(as.matrix(p))
        near_point <- matchpt(p, df) # Computing nearest point on spline
        d <- p - df[near_point$index,] # Computing distance vector
        return(c(near_point$distance, d[1], d[2], d[3], near_point$index))
        }))
    )
    colnames(Dist_df) <- c('Dist_HVG', 'dvecx', 'dvecy', 'dvecz', 'nearidx')

    splinefit_df$Distance <- Dist_df$Dist_HVG
    splinefit_df$nearidx <- Dist_df$nearidx
    splinefit_df$dvecx <- Dist_df$dvecx
    splinefit_df$dvecy <- Dist_df$dvecy
    splinefit_df$dvecz <- Dist_df$dvecz
  }
  splinefit_df$splinex <- xyz1$logMean
  splinefit_df$spliney <- xyz1$logCV
  splinefit_df$splinez <- xyz1$Dropout

  # Removing Genes with too high or low dropout
  splinefit_df$Distance[splinefit_df$Dropout > 0.95 | splinefit_df$Dropout < 0.01] <- 0
  top_HVG <- splinefit_df %>% top_n(n=nHVGs, wt=Distance)
  mask <- rownames(splinefit_df) %in% rownames(top_HVG)
  splinefit_df$HVG <- mask
  splinefit_df <- splinefit_df %>% arrange(Dropout)

  # Trimming genes at the ends of spline
  if (use.ndist){
    splinefit_df <- splinefit_df %>% filter(nearidx != 1 & nearidx != max(nearidx))
  }
  splinefit_df <- splinefit_df %>% arrange(-Distance)

  return(splinefit_df)
}


#' @export HVG_QC
#' @title Quality control
#' @importFrom Matrix rowSums
#' @description QC filter scRNA-seq expression data
#' @author Shreyan Gupta <xenon8778@tamu.edu>
#' @import scuttle
#' @import SummarizedExperiment
#' @import SingleCellExperiment
#' @return QC filtered SingleCellExperiment
#' @param X Count matrix or SingleCellExperiment
#' @param ncounts An integer value. Defines the minimum reads required for a cell to be included in the analysis.
#' @param ncells An integer value. Defines the minimum cells required for a gene to be included in the analysis.
#' @param mtPerc A double value. Defines the minimum percent mitochondrial genes expression required for a cell to be excluded from the analysis.
#' @examples
#' # example code
#' ## Load Data
#' load(system.file("extdata", "WT_count.rda", package="SplineDV")) # WT Sample
#' adata = HVG_QC(WT_count)

HVG_QC <- function(X, ncounts=1000, ncells=15, mtPerc=15){

  # Check if object is SingleCellExperiment
  if (is(X,"SingleCellExperiment")){
    sce <- X
  } else {
    sce <- SingleCellExperiment(list(counts=X))
  }

  ## Calculating mitochondrial expression
  if(dim(table(startsWith(rownames(sce), 'mt-'))) == 2){
    isMito <- any(rownames(sce) == "mt")
    qc_df <- perCellQCMetrics(sce, subsets=list(Mito=isMito))
  } else {
    isMito <- any(rownames(sce) == "MT")
    qc_df <- perCellQCMetrics(sce, subsets=list(Mito=isMito))
  }

  qc.lib <- qc_df$sum < ncounts
  qc.mito <- qc_df$subsets_Mito_percent > mtPerc
  discard_cells <- qc.lib | qc.mito
  sce <- sce[, !discard_cells]

  discard_genes <- rowSums(SummarizedExperiment::assay(sce,'counts')) < ncells
  sce <- sce[!discard_genes,]

  return(sce)
}


#' @export HVG_plot
#' @title Plot HVG 3D scatter plot
#' @description Plots 3D scatter plot with HVGs or target gene highlighted
#' @author Shreyan Gupta <xenon8778@tamu.edu>
#' @import plotly
#' @import dplyr
#' @return 3D plotly scatter plot
#' @param df Resultant data.frame after HVG_splinefit analysis
#' @param targetgene An integer value. Defines the minimum reads required for a cell to be included in the analysis.
#' @param lwd An integer value. Defines line width for spline.
#' @param ptSize An integer value. Defines point size for dots
#' @param dlwd An integer value. Defines Line width for target gene distance.
#' @examples
#' # example code
#' ## Load Data
#' load(system.file("extdata", "WT_count.rda", package="SplineDV")) # WT Sample
#' HVG_df <- HVG_splinefit(WT_count)
#' fig <- HVG_plot(HVG_df)
#' print(fig)

HVG_plot <- function(df, targetgene=NULL, ptSize=3, lwd=5, dlwd=7){
  df$genenames <- rownames(df)
  if(is.null(targetgene)){
    fig <- plot_ly(df, x=~logMean, y=~logCV, z=~Dropout,
                   mode='markers', color=~HVG, colors=c('grey','red'),
                   type="scatter3d", text=~genenames,
                   marker=list(size=ptSize, opacity=0.5)) %>%
      add_trace(data = df %>% arrange(Dropout),x=~splinex, y=~spliney, z=~splinez,
                type="scatter3d", mode='lines+markers', showlegend = FALSE,
                marker=list(size=1, color='black'),
                line=list(width=lwd, color='black', opacity=1, reverscale=FALSE)) %>%
      layout(legend=list(title=list(text='<b> HVG </b>')),
             scene = list(camera = list(eye = list(x = 2, y = 2, z = 2))))
  } else {
    df_sub <- df[targetgene,]
    df_sub <- rbind(end = df_sub[,c('logMean', 'logCV', 'Dropout')],
                     start = unlist(df_sub[,c('splinex', 'spliney', 'splinez')]))

    fig <- plot_ly(df, x=~logMean, y=~logCV, z=~Dropout,
                   mode='markers', type="scatter3d", text=~genenames,
                   marker=list(size=ptSize, color='grey', opacity=0.5)) %>%
      add_trace(data = df %>% arrange(Dropout),x=~splinex, y=~spliney, z=~splinez,
                mode='lines+markers',
                marker=list(size=1, color='black', opacity=1),
                line=list(width=lwd, color='black', opacity=1, reverscale=FALSE)) %>%
      add_trace(data = df_sub, x=~logMean, y=~logCV, z=~Dropout,
                mode='lines+markers',
                text=targetgene,
                marker=list(size=ptSize+3, color='red',opacity=1),
                line=list(width=dlwd, color='red')) %>%
      add_text(data = df_sub['end',], x=~logMean, y=~logCV, z=~Dropout,
               text=targetgene, mode='text', inherit=FALSE) %>%
      layout(showlegend = FALSE,
             scene = list(camera = list(eye = list(x = 2, y = 2, z = 2))))
  }
  return(fig)
}


#' @export DV_plot
#' @title Plot 3D scatter plot in two conditions
#' @description Plots 3D gene statistic scatter plot of target genes in twop conditions
#' @author Shreyan Gupta <xenon8778@tamu.edu>
#' @import plotly
#' @import dplyr
#' @return 3D plotly scatter plot
#' @param df Resultant data.frame after DV_splinefit analysis
#' @param targetgene An integer value. Defines the minimum reads required for a cell to be included in the analysis.
#' @param lwd An integer value. Defines line width for spline.
#' @param ptSize An integer value. Defines point size for dots
#' @param dlwd An integer value. Defines Line width for target gene distance.
#' @examples
#' # example code
#' ## Load Data
#' load(system.file("extdata", "WT_count.rda", package="SplineDV")) # WT Sample
#' load(system.file("extdata", "KO_count.rda", package="SplineDV")) # KO Sample
#' DV_res <- DV_splinefit(X=KO_count, Y=WT_count)
#' fig <- DV_plot(DV_res)
#' print(fig)

DV_plot <- function(df, targetgene=NULL, ptSize=3, lwd=5, dlwd=7){
  if(isFALSE(targetgene %in% df$gene)){
    stop('Gene not in analysis!')
  }
  if (is.null(targetgene)){
    df_sub <- df[1,]
  } else {
    df_sub <- df[df$gene == targetgene,]
  }

  df_subX <- rbind(end = df_sub[,c('mu1', 'CV1', 'drop1')],
                  start = unlist(df_sub[,c('X_splinex', 'X_spliney',
                                           'X_splinez')]))
  colnames(df_subX) <- c('logMean', 'logCV', 'Dropout')
  df_subY <- rbind(end = df_sub[,c('mu2', 'CV2', 'drop2')],
                   start = unlist(df_sub[,c('Y_splinex', 'Y_spliney',
                                            'Y_splinez')]))
  colnames(df_subY) <- c('logMean', 'logCV', 'Dropout')

  fig <- plot_ly(df %>% arrange(drop1), type="scatter3d", mode='markers') %>%
    add_trace(data = df %>% arrange(drop1),x=~X_splinex, y=~X_spliney, z=~X_splinez,
              mode='lines',
              line=list(width=lwd, color='firebrick', opacity=1, reverscale=FALSE)) %>%
    add_trace(data = df %>% arrange(drop2),x=~Y_splinex, y=~Y_spliney, z=~Y_splinez,
              mode='lines',
              line=list(width=lwd, color='steelblue', opacity=1, reverscale=FALSE)) %>%
    add_trace(data = df_subX, x=~logMean, y=~logCV, z=~Dropout,
              mode='lines+markers',
              opacity=1, text=df_sub$genes,
              marker=list(size=ptSize, color='firebrick', reverscale=FALSE),
              line=list(width=dlwd, color='firebrick')) %>%
    add_trace(data = df_subY, x=~logMean, y=~logCV, z=~Dropout,
              mode='lines+markers',
              opacity=1, text=df_sub$genes,
              marker=list(size=ptSize, color='steelblue', reverscale=FALSE),
              line=list(width=dlwd, color='steelblue')) %>%
    add_text(data = df_subX['end',], x=~logMean, y=~logCV, z=~Dropout,
             mode='text', opacity=1, text=df_sub$genes, inherit=FALSE) %>%
    add_text(data = df_subY['end',], x=~logMean, y=~logCV, z=~Dropout,
             mode='text', opacity=1, text=df_sub$genes, inherit=FALSE) %>%
    layout(showlegend = FALSE,
           scene=list(xaxis = list(title = 'logMean'),
                      yaxis=list(title = 'logCV'),
                      zaxis=list(title = 'Dropout'),
                      camera = list(eye = list(x = 2, y = 2, z = 2))
                      ))
  return(fig)
}


#' @export DV_splinefit
#' @title Spline-DV
#' @description Differential Variability analysis
#' @author Shreyan Gupta <xenon8778@tamu.edu>
#' @import dplyr
#' @importFrom stats pnorm
#' @return A data.frame with DV Statistics
#' @param X Count matrix or SingleCellExperiment (Test sample)
#' @param Y Count matrix or SingleCellExperiment (Control sample)
#' @param ncounts An integer value. Defines the minimum reads required for a cell to be included in the analysis.
#' @param ncells An integer value. Defines the minimum cells required for a gene to be included in the analysis.
#' @param mtPerc A double value. Defines the minimum percent mitochondrial genes expression required for a cell to be excluded from the analysis.
#' @param detailed A boolean value. Defines whether to add individual HVG_splinefit data.frame to the output.
#' @examples
#' # example code
#' # Load Data
#' load(system.file("extdata", "WT_count.rda", package="SplineDV")) # WT Sample
#' load(system.file("extdata", "KO_count.rda", package="SplineDV")) # KO Sample
#' DV_res <- DV_splinefit(X=KO_count, Y=WT_count)

DV_splinefit <- function(X, Y, ncounts=500, ncells=15,
                         mtPerc=15, detailed=FALSE) {

  # QC Filtering
  X <- HVG_QC(X, ncounts = ncounts, ncells = ncells,
             mtPerc = mtPerc)
  Y <- HVG_QC(Y, ncounts = ncounts, ncells = ncells,
             mtPerc = mtPerc)

  # Intersect gene space of the two data sets
  feat <- intersect(rownames(X),rownames(Y))
  X <- X[feat,]
  Y <- Y[feat,]

  # Computing distances
  message('Computing distances for Data 1')
  res_X <- HVG_splinefit(X, QC = FALSE)
  message('Computing distances for Data 2')
  res_Y <- HVG_splinefit(Y, QC = FALSE)

  # Intersect gene space of the two data sets
  feat <- intersect(rownames(res_X),rownames(res_Y))
  X_filt <- res_X[feat,]
  Y_filt <- res_Y[feat,]

  # Creating summary table
  df <- as.data.frame(list('genes' = feat,
                          'dist1' = X_filt$Distance,
                          'dist2' = Y_filt$Distance,
                          'mu1' = X_filt$logMean,
                          'mu2' = Y_filt$logMean,
                          'CV1' = X_filt$logCV,
                          'CV2' = Y_filt$logCV,
                          'drop1' = X_filt$Dropout,
                          'drop2' = Y_filt$Dropout,
                          'X_splinex' = X_filt$splinex,
                          'X_spliney' = X_filt$spliney,
                          'X_splinez' = X_filt$splinez,
                          'Y_splinex' = Y_filt$splinex,
                          'Y_spliney' = Y_filt$spliney,
                          'Y_splinez' = Y_filt$splinez,
                          'X_dvecx' = X_filt$dvecx,
                          'X_dvecy' = X_filt$dvecy,
                          'X_dvecz' = X_filt$dvecz,
                          'Y_dvecx' = Y_filt$dvecx,
                          'Y_dvecy' = Y_filt$dvecy,
                          'Y_dvecz' = Y_filt$dvecz),
                     row.names = feat)

  # Trimming zeros spline differences
  df <- df[!c(df$dist1 == 0 & df$dist2 == 0),]

  # Computing Diff Distance
  DV_res <- df %>% as_tibble %>%
    mutate('dist1'=sqrt(X_dvecx^2 + X_dvecy^2 + X_dvecz^2)) %>%
    mutate('dist2'=sqrt(Y_dvecx^2 + Y_dvecy^2 + Y_dvecz^2)) %>%
    mutate('Vector_Dist'=sqrt((X_dvecx - Y_dvecx)^2 +
                                (X_dvecy - Y_dvecy)^2 +
                                (X_dvecz - Y_dvecz)^2))  %>%
    mutate('Direction'=sign(dist1 - dist2)) %>%
    mutate('Z'= as.numeric(scale(Vector_Dist))) %>%
    mutate('Pval' = 2*pnorm(abs(Z), lower.tail = FALSE)) %>%
    arrange(Pval)

  res_out <- DV_res %>% select(genes, mu1, mu2, CV1, CV2, drop1, drop2,
                               dist1, dist2, X_splinex, X_spliney, X_splinez,
                               Y_splinex, Y_spliney, Y_splinez,
                              Vector_Dist, Direction, Pval) %>% as.data.frame()
  output <- list(HVG_X=res_X, HVG_Y=res_Y, DV=res_out)
  if(detailed){
    res_out <- output
  }
  return(res_out)
}
