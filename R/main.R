#' @title Spline-HVG
#' @description Compute Highly Variable Genes from an scRNAseq expression data.
#' @export splineHVG
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
#' HVG_res <- splineHVG(WT_count, nHVGs = 100)

splineHVG <- function(X, QC=TRUE,
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
    adata <- hvgQC(adata, ncounts=ncounts, ncells=ncells,
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
  splinefitDF <- as.data.frame(Means)
  splinefitDF$CV <- CV
  splinefitDF$Dropout <- Dropout
  splinefitDF$logMean <- log(Means+1)
  splinefitDF$logCV <- log(CV+1)
  splinefitDF <- splinefitDF %>% arrange(logMean) # arrange by Mean expression
  # splinefitDF <- splinefitDF %>% filter(Dropout>0.01 & Dropout<0.99) # Trimming ends

  # Calculate the differences and the squared differences between consecutive elements
  diffLgu <- diff(splinefitDF$logMean)
  diffLgcv <- diff(splinefitDF$logCV)
  diffDropr <- diff(splinefitDF$Dropout)
  diffSquaredSum <- diffLgu^2 + diffLgcv^2 + diffDropr^2

  # Calculate the cumulative sum
  s <- c(0, cumsum(sqrt(diffSquaredSum)))
  xyz <- splinefitDF %>% dplyr::select(logMean,logCV,Dropout)

  fitx <- smooth.spline(s, splinefitDF$logMean, df=degf, spar=spar)
  fity <- smooth.spline(s, splinefitDF$logCV, df=degf, spar=spar)
  fitz <- smooth.spline(s, splinefitDF$Dropout, df=degf, spar=spar)

  # Computing Distances
  xyz1 <- cbind(predict(fitx,s)$y, predict(fity,s)$y, predict(fitz,s)$y)
  xyz1 <- as.data.frame(xyz1)
  colnames(xyz1) <- c('logMean','logCV','Dropout')
  rownames(xyz1) <- rownames(splinefitDF)
  euclidean <- function(a, b) sqrt(sum((a - b) ^ 2))
  if (!use.ndist){
    df <- cbind(xyz, xyz1)
    p1 <- as.list(as.data.frame(t(df)))
    distDF <- as.data.frame(
      do.call('rbind', lapply(p1, FUN=function(p){
        d <- euclidean(c(p[1:3]), c(p[4:6]))# Computing nearest point on spline
        return(d)
        }))
    )
    names(distDF) <- c('distHVG')
    splinefitDF$Distance <- distDF$distHVG
  } else {
    df <- as.matrix(xyz1)
    p1 <- as.list(as.data.frame(t(xyz)))
    distDF <- as.data.frame(
      do.call('rbind', lapply(p1, FUN=function(p){
        p <- t(as.matrix(p))
        nearPoint <- matchpt(p, df) # Computing nearest point on spline
        d <- p - df[nearPoint$index,] # Computing distance vector
        return(c(nearPoint$distance, d[1], d[2], d[3], nearPoint$index))
        }))
    )
    colnames(distDF) <- c('distHVG', 'dvecx', 'dvecy', 'dvecz', 'nearidx')

    splinefitDF$Distance <- distDF$distHVG
    splinefitDF$nearidx <- distDF$nearidx
    splinefitDF$dvecx <- distDF$dvecx
    splinefitDF$dvecy <- distDF$dvecy
    splinefitDF$dvecz <- distDF$dvecz
  }
  splinefitDF$splinex <- xyz1$logMean
  splinefitDF$spliney <- xyz1$logCV
  splinefitDF$splinez <- xyz1$Dropout

  # Removing Genes with too high or low dropout
  splinefitDF$Distance[splinefitDF$Dropout > 0.95 | splinefitDF$Dropout < 0.01] <- 0
  topHVG <- splinefitDF %>% top_n(n=nHVGs, wt=Distance)
  mask <- rownames(splinefitDF) %in% rownames(topHVG)
  splinefitDF$HVG <- mask
  splinefitDF <- splinefitDF %>% arrange(Dropout)

  # Trimming genes at the ends of spline
  if (use.ndist){
    splinefitDF <- splinefitDF %>% filter(nearidx != 1 & nearidx != max(nearidx))
  }
  splinefitDF <- splinefitDF %>% arrange(-Distance)

  return(splinefitDF)
}


#' @export hvgQC
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
#' adata = hvgQC(WT_count)

hvgQC <- function(X, ncounts=1000, ncells=15, mtPerc=15){

  # Check if object is SingleCellExperiment
  if (is(X,"SingleCellExperiment")){
    sce <- X
  } else {
    sce <- SingleCellExperiment(list(counts=X))
  }

  ## Calculating mitochondrial expression
  if(dim(table(startsWith(rownames(sce), 'mt-'))) == 2){
    isMito <- any(rownames(sce) == "mt")
    qcDF <- perCellQCMetrics(sce, subsets=list(Mito=isMito))
  } else {
    isMito <- any(rownames(sce) == "MT")
    qcDF <- perCellQCMetrics(sce, subsets=list(Mito=isMito))
  }

  qc.lib <- qcDF$sum < ncounts
  qc.mito <- qcDF$subsets_Mito_percent > mtPerc
  discard_cells <- qc.lib | qc.mito
  sce <- sce[, !discard_cells]

  discardGenes <- rowSums(SummarizedExperiment::assay(sce,'counts')) < ncells
  sce <- sce[!discardGenes,]

  return(sce)
}


#' @export HVGPlot
#' @title Plot HVG 3D scatter plot
#' @description Plots 3D scatter plot with HVGs or target gene highlighted
#' @author Shreyan Gupta <xenon8778@tamu.edu>
#' @import plotly
#' @import dplyr
#' @return 3D plotly scatter plot
#' @param df Resultant data.frame after splineHVG analysis
#' @param targetgene An integer value. Defines the minimum reads required for a cell to be included in the analysis.
#' @param lwd An integer value. Defines line width for spline.
#' @param ptSize An integer value. Defines point size for dots
#' @param dlwd An integer value. Defines Line width for target gene distance.
#' @examples
#' # example code
#' ## Load Data
#' load(system.file("extdata", "WT_count.rda", package="SplineDV")) # WT Sample
#' hvgDF <- splineHVG(WT_count)
#' fig <- HVGPlot(hvgDF)
#' print(fig)

HVGPlot <- function(df, targetgene=NULL, ptSize=3, lwd=5, dlwd=7){
  df$genenames <- rownames(df)
  if(is.null(targetgene)){
    fig <- plot_ly(df, x=~logMean, y=~logCV, z=~Dropout,
                   mode='markers', color=~HVG, colors=c('grey','red'),
                   type="scatter3d", text=~genenames,
                   marker=list(size=ptSize, opacity=0.5)) %>%
      add_trace(data = df %>% arrange(logMean),x=~splinex, y=~spliney, z=~splinez,
                type="scatter3d", mode='lines+markers', showlegend = FALSE,
                marker=list(size=1, color='black'),
                line=list(width=lwd, color='black', opacity=1, reverscale=FALSE)) %>%
      layout(legend=list(title=list(text='<b> HVG </b>')),
             scene = list(camera = list(eye = list(x = 2, y = 2, z = 2))))
  } else {
    dfSub <- df[targetgene,]
    dfSub <- rbind(end = dfSub[,c('logMean', 'logCV', 'Dropout')],
                     start = unlist(dfSub[,c('splinex', 'spliney', 'splinez')]))

    fig <- plot_ly(df, x=~logMean, y=~logCV, z=~Dropout,
                   mode='markers', type="scatter3d", text=~genenames,
                   marker=list(size=ptSize, color='grey', opacity=0.5)) %>%
      add_trace(data = df %>% arrange(logMean),x=~splinex, y=~spliney, z=~splinez,
                mode='lines+markers',
                marker=list(size=1, color='black', opacity=1),
                line=list(width=lwd, color='black', opacity=1, reverscale=FALSE)) %>%
      add_trace(data = dfSub, x=~logMean, y=~logCV, z=~Dropout,
                mode='lines+markers',
                text=targetgene,
                marker=list(size=ptSize+3, color='red',opacity=1),
                line=list(width=dlwd, color='red')) %>%
      add_text(data = dfSub['end',], x=~logMean, y=~logCV, z=~Dropout,
               text=targetgene, mode='text', inherit=FALSE) %>%
      layout(showlegend = FALSE,
             scene = list(camera = list(eye = list(x = 2, y = 2, z = 2))))
  }
  return(fig)
}


#' @export DVPlot
#' @title Plot 3D scatter plot in two conditions
#' @description Plots 3D gene statistic scatter plot of target genes in twop conditions
#' @author Shreyan Gupta <xenon8778@tamu.edu>
#' @import plotly
#' @import dplyr
#' @return 3D plotly scatter plot
#' @param df Resultant data.frame after splineDV analysis
#' @param targetgene An integer value. Defines the minimum reads required for a cell to be included in the analysis.
#' @param lwd An integer value. Defines line width for spline.
#' @param ptSize An integer value. Defines point size for dots
#' @param dlwd An integer value. Defines Line width for target gene distance.
#' @examples
#' # example code
#' ## Load Data
#' load(system.file("extdata", "WT_count.rda", package="SplineDV")) # WT Sample
#' load(system.file("extdata", "KO_count.rda", package="SplineDV")) # KO Sample
#' DVres <- splineDV(X=KO_count, Y=WT_count)
#' fig <- DVPlot(DVres)
#' print(fig)

DVPlot <- function(df, targetgene=NULL, ptSize=3, lwd=5, dlwd=7){
  if(isFALSE(targetgene %in% df$gene)){
    stop('Gene not in analysis!')
  }
  if (is.null(targetgene)){
    dfSub <- df[1,]
  } else {
    dfSub <- df[df$gene == targetgene,]
  }

  dfSubX <- rbind(end = dfSub[,c('mu1', 'CV1', 'drop1')],
                  start = unlist(dfSub[,c('X_splinex', 'X_spliney',
                                           'X_splinez')]))
  colnames(dfSubX) <- c('logMean', 'logCV', 'Dropout')
  dfSubY <- rbind(end = dfSub[,c('mu2', 'CV2', 'drop2')],
                   start = unlist(dfSub[,c('Y_splinex', 'Y_spliney',
                                            'Y_splinez')]))
  colnames(dfSubY) <- c('logMean', 'logCV', 'Dropout')

  fig <- plot_ly(df %>% arrange(mu1), type="scatter3d", mode='markers') %>%
    add_trace(data = df %>% arrange(mu1),x=~X_splinex, y=~X_spliney, z=~X_splinez,
              mode='lines',
              line=list(width=lwd, color='firebrick', opacity=1, reverscale=FALSE)) %>%
    add_trace(data = df %>% arrange(mu2),x=~Y_splinex, y=~Y_spliney, z=~Y_splinez,
              mode='lines',
              line=list(width=lwd, color='steelblue', opacity=1, reverscale=FALSE)) %>%
    add_trace(data = dfSubX, x=~logMean, y=~logCV, z=~Dropout,
              mode='lines+markers',
              opacity=1, text=dfSub$genes,
              marker=list(size=ptSize, color='firebrick', reverscale=FALSE),
              line=list(width=dlwd, color='firebrick')) %>%
    add_trace(data = dfSubY, x=~logMean, y=~logCV, z=~Dropout,
              mode='lines+markers',
              opacity=1, text=dfSub$genes,
              marker=list(size=ptSize, color='steelblue', reverscale=FALSE),
              line=list(width=dlwd, color='steelblue')) %>%
    add_text(data = dfSubX['end',], x=~logMean, y=~logCV, z=~Dropout,
             mode='text', opacity=1, text=dfSub$genes, inherit=FALSE) %>%
    add_text(data = dfSubY['end',], x=~logMean, y=~logCV, z=~Dropout,
             mode='text', opacity=1, text=dfSub$genes, inherit=FALSE) %>%
    layout(showlegend = FALSE,
           scene=list(xaxis = list(title = 'logMean'),
                      yaxis=list(title = 'logCV'),
                      zaxis=list(title = 'Dropout'),
                      camera = list(eye = list(x = 2, y = 2, z = 2))
                      ))
  return(fig)
}


#' @export splineDV
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
#' @param detailed A boolean value. Defines whether to add individual splineHVG data.frame to the output.
#' @examples
#' # example code
#' # Load Data
#' load(system.file("extdata", "WT_count.rda", package="SplineDV")) # WT Sample
#' load(system.file("extdata", "KO_count.rda", package="SplineDV")) # KO Sample
#' DVres <- splineDV(X=KO_count, Y=WT_count)

splineDV <- function(X, Y, ncounts=500, ncells=15,
                         mtPerc=15, detailed=FALSE) {

  # QC Filtering
  X <- hvgQC(X, ncounts = ncounts, ncells = ncells,
             mtPerc = mtPerc)
  Y <- hvgQC(Y, ncounts = ncounts, ncells = ncells,
             mtPerc = mtPerc)

  # Intersect gene space of the two data sets
  feat <- intersect(rownames(X),rownames(Y))
  X <- X[feat,]
  Y <- Y[feat,]

  # Computing distances
  message('Computing distances for Data 1')
  res_X <- splineHVG(X, QC = FALSE)
  message('Computing distances for Data 2')
  res_Y <- splineHVG(Y, QC = FALSE)

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
  DVres <- df %>% as_tibble %>%
    mutate('dist1'=sqrt(X_dvecx^2 + X_dvecy^2 + X_dvecz^2)) %>%
    mutate('dist2'=sqrt(Y_dvecx^2 + Y_dvecy^2 + Y_dvecz^2)) %>%
    mutate('vectorDist'=sqrt((X_dvecx - Y_dvecx)^2 +
                                (X_dvecy - Y_dvecy)^2 +
                                (X_dvecz - Y_dvecz)^2))  %>%
    mutate('Direction'=sign(dist1 - dist2)) %>%
    mutate('Z'= as.numeric(scale(vectorDist))) %>%
    mutate('Pval' = 2*pnorm(abs(Z), lower.tail = FALSE)) %>%
    arrange(Pval)

  resOut <- DVres %>% select(genes, mu1, mu2, CV1, CV2, drop1, drop2,
                               dist1, dist2, X_splinex, X_spliney, X_splinez,
                               Y_splinex, Y_spliney, Y_splinez,
                              vectorDist, Direction, Pval) %>% as.data.frame()
  output <- list(HVG_X=res_X, HVG_Y=res_Y, DV=resOut)
  if(detailed){
    resOut <- output
  }
  return(resOut)
}
