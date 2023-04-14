

#-----------------------------------------------------------------------------#
# check expression matrix for meta-analysis                                   #
#-----------------------------------------------------------------------------#
check_exp <- function(expr){
  K <- length(expr)

  ## 检查数据类型
  if (!is.list(expr)){
    stop("Please enter a list containing data frames")
  }
  for (k in 1:K) {
    if (!is.data.frame(expr[[k]])){
      stop("Please enter a list containing data frames")
    }
  }

  ## 检查名称
  if (is.null(row.names(expr[[1]]))) {
    ng<-nrow(expr[[1]])
    for (k in 1:K) {
      row.names(x[[k]])<-paste("gene",1:ng)
    }
  }

  return(expr)
}


check_group <- function(group){

  ## 与expr是否匹配？
  K <- length(group)

  ## 是否是二因子型变量？
}


## 检查参数的维度是否相同？ eg: expr和group


#------------------------------------------------------------------------------#
# check dimensions and size of argument
#------------------------------------------------------------------------------#
check.dim <- function(x,y,ind.method,meta.method,paired){
  K<-length(x)
  nperstudy<-sapply(x,function(y) ncol(y))
  nlabels<-sapply(y,function(z) length(z))
  if(sum(nperstudy==nlabels)!=K) {
    stop(cat("The number of samples does not match with the dimension of labels in study(s)",paste((1:K)[nperstudy!=nlabels],"",collapse=","),"!"))
  }
  if(sum(meta.method%in%c("FEM","REM","minMCC","rankProd"))<1) {
    if(length(ind.method)!=K)stop(paste('Argument "ind.method" should be a character vector of size',K))
  }
  if(("REM"%in%meta.method|"FEM"%in%meta.method)&(length(paired)!=K)) {
    stop(paste('Argument "paired" should be a logical vecter of size',K))
  }
}
