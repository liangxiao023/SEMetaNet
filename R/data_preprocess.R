


metaOverlap <- function(expr){
  check_exp(expr)

  K <- length(expr)
  for(k in 1:K){
    rownames(expr[[k]]) <- toupper(rownames(expr[[k]]))
  }
  interG <- rownames(expr[[1]])
  for(k in 2:K){
    interG <- intersect(interG,rownames(expr[[k]]))
  }
  if(length(interG)<=5){
    stop(paste("Error: Number of overlap genes: ",interG,", less than 5.",sep=""))
  }
  for(k in 1:K){
    expr[[k]] <- expr[[k]][interG,]
  }

  return(expr)
}






remove_batch_effects <- function(expression_matrices, batch_labels){
  
  library(sva)
  library(limma)
  
  # Combine expression matrices into a single data frame
  df <- do.call(cbind, expression_matrices)
  ncols <- sapply(expression_matrices, ncol)
  
  # Perform combat correction
  corrected_df <- ComBat(df, batch_labels)
  
  corrected_tmp <- corrected_df
  for (i in 1:length(expression_matrices)) {
    expression_matrices[[i]] <- corrected_tmp[,1:ncol(expression_matrices[[i]])]
    corrected_tmp <- corrected_tmp[,-(1:ncol(expression_matrices[[i]]))]
  }
  
  return(expression_matrices)
}



#'  data <- remove_batch_effects(expression_matrices, batch_labels)






#to-do: 把array和rnaseq的过滤函数合并

## filter-microarray
Filter.list <- function(datasets, data.type="microarray", del.perc=c(0.3,0.3)) {
  if (data.type == "microarray") {
    mean.rank <- sapply(datasets,function(z)rank(apply(z, 1, mean, na.rm=T)))
    mean.r.mv <- rowMeans(mean.rank, na.rm=T)
    mean.r.mv <- mean.r.mv[order(mean.r.mv, decreasing=T)]
    index <- which(mean.r.mv > quantile(mean.r.mv, del.perc[1]))
    gene.mv <- names(mean.r.mv)[index]
    sd.rank <- sapply(datasets,function(z)rank(apply(z[gene.mv,], 1, sd, na.rm=T)))
    mean.r.sd <- rowMeans(sd.rank, na.rm=T)
    mean.r.sd <- mean.r.sd[order(mean.r.sd, decreasing=T)]
    index <- which(mean.r.sd > quantile(mean.r.sd, del.perc[2]))
    final.genes <- names(mean.r.sd)[index]
    res <- lapply(datasets, function(x) x[final.genes, ])
  } 
  return(res)
}


## filter-rnaseq
Filter.rnaseq <- function(datasets, data.type="RNAseq", del.perc=c(0.3,0.3), threshold=5) {
  min.mean <- apply(sapply(datasets, function(z)rowMeans(z)), 1, min)
  index <- min.mean > threshold
  res <- lapply(datasets, function(x) x[index, ])
  return(res)
}
