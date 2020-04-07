#' Sample outlier tagging via PCA
#'
#' @param my.vector PC1=PC$rotation[,1] where PC <- prcomp(df, scale.=T, center = T)
#' @param sampleqc sample names fetched from rownames(PC$rotation) 
#'
#' @return data frame with sample Id's and PCA outlier tagging columns 
#' @export
#'
#' @examples samplepca.outlier <- tag.pca.outlier(my.vector = PC$PC1, sampleqc = rownames(PC$rotation))
tag.pca.outlier <- function(my.vector, sampleqc){
  # PCA outlier tagging 
  pop.sd <- sd(my.vector) * sqrt((length(my.vector) - 1)/(length(my.vector)))
  pop.mean <- mean(my.vector)
  pop.median <- median(my.vector)
  pop.mad <- mad(my.vector)
  
  # zscore of PC1s ----------------------------------------------
  z.score <- lapply(seq_along(my.vector), function(i) {
    data.point <- my.vector[i]
    z <- (data.point - pop.mean) / pop.sd
    if (is.finite(z)) {
      ifelse((abs(z) < 3), "non-outlier", "outlier")
    } else {
      "non-outlier"
    }
  })
  vec <- as.vector(do.call(cbind, z.score))
  
  samplepca.outliers <- as.data.frame(cbind(sampleqc, vec))
  names(samplepca.outliers) <- c("SampleID", "pca_outlier")
  samplepca.outliers[samplepca.outliers=="outlier"] <- 1
  samplepca.outliers[samplepca.outliers=="non-outlier"] <- 0
  
  return(samplepca.outliers)
}

