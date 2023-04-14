

FEMR <- function(expr, group, cutoff_FEM = 0.05){
  expr <- check_exp(expr)
  check_group(group)
  if (!is.numeric(cutoff_FEM) | length(cutoff_FEM) != 1){
    stop("cutoff_FEM should be a single number")
  }
  if (cutoff_FEM <= 0 | cutoff_FEM >= 1){
    stop("cutoff_FEM should be between 0 and 1")
  }

  ## FEM
  ind.ES <- cal_ES(expr, group)
  expr <- delete_expr(expr, ind.ES)
  meta.FEM.ES <- MetaDE.ES(ind.ES)$pval %>% matrix()             ## %>%
  input_boot <- which(meta.FEM.ES < cutoff.FEM)

  ## bootMRMR
  dataExpr <- data.frame(do.call(cbind,expr))[input_boot,]
  Samples <- ncol(dataExpr)
  label <- vector("numeric",0)
  for (k in 1:K) {
    tmp <- as.numeric(as.factor(group[[k]]))
    label <- c(label,tmp)
  }
  label[which(label==2)] <- -1
  pval <- c(BootMRMR::pval.mbmr(dataExpr, group, round(0.9*Samples), 200, Q=0.5))
  FDR <- p.adjust(pval,method = "BH")
  res <- list(pval=pval, FDR=FDR)
  return(res)
}


cal_ES <- function(expr, group){
  expr <- check_exp(expr)
  check_group(group)
  K <- length(expr)
  dataES <- vector("list",K)
  for (k in 1:K) {
    dataES[[k]][[1]] <- as.matrix(expr[[k]])
    dataES[[k]][[2]] <- group[[k]]
  }
  ind.ES <- ind.cal.ES(dataES, paired=rep(FALSE,K))
  return(ind.ES)
}


delete_expr <- function(expr, ind.ES){
  temp <- NULL
  temp <- which(rowSums(is.na(ind.ES$ES))==0)                    ## temp要重新设置一下，考虑Var为0而ES不为0的情况
  if(!is.null(temp)){
    ind.ES$ES <- ind.ES$ES[temp,]
    ind.ES$Var <- ind.ES$Var[temp,]
  }
  for (k in 1:K) {
    expr[[k]] <- expr[[k]][temp, ]
  }
  return(expr)
}


ind.cal.ES<-function(x,paired,nperm=NULL){
  K<-length(x)
  res<-get.ES(x,paired=paired)
  if(!is.null(nperm)){
    perm.ES<-perm.Var<-NULL
    for(i in 1:nperm){
      for(k in 1:K){
        x[[k]][[2]]<-perm.lab(x[[k]][[2]],paired[k])
      }
      tempRes<-get.ES(x,paired=paired)
      perm.ES<-rbind(perm.ES,tempRes$ES)
      perm.Var<-rbind(perm.Var,tempRes$Var)
    }
  }else{
    perm.ES<-perm.Var<-NULL
  }
  if(is.null(names(x))){colnames(res$ES)<-colnames(res$Var)<-paste("dataset",
                                                                   1:K,sep="")
  }else{colnames(res$ES)<-colnames(res$Var)<-names(x)}
  result<-list(ES=res$ES,Var=res$Var,perm.ES=perm.ES,perm.Var=perm.Var)
  attr(result,"nperstudy")<-attr(res,"nperstudy")
  attr(result,"nperlabelperstudy")<-attr(res,"nperlabelperstudy")
  return(result)
}


get.ES<-function(x,paired){
  K<-length(x)
  ES.m<-Var.m<- N<-n<-NULL
  for (k in 1:K){
    y<-x[[k]][[1]]      ## counts matrix
    l<-x[[k]][[2]]      ## design matrix
    temp<-cal.ES(y,l,paired[k])
    ES.m<-cbind(ES.m,temp[,"dprime"])
    Var.m<-cbind(Var.m,temp[,"vardprime"])
    N<-c(N,length(l))
    n<-c(n,table(l))
  }
  rownames(ES.m)<-rownames(y)
  rownames(Var.m)<-rownames(y)
  colnames(ES.m)<-paste("study",1:K,sep="")
  colnames(Var.m)<-paste("study",1:K,sep="")
  res<-list(ES=ES.m,Var=Var.m)
  attr(res,"nperstudy")<-N
  attr(res,"nperlabelperstudy")<-n
  return(res)
}


cal.ES<-function(y,l,paired=FALSE){
  l<-unclass(factor(l))
  n<-table(factor(l))
  if(paired){
    if (n[1]!=n[2]) {
      stop("The study is not paired design")
    }
    DM<-y[,l==2] ## disease is label = 2
    CM<-y[,l==1] ## control is label = 1, originally ref level
    ydiff<-DM-CM
    den<-sqrt(1/(n[1]-1)*(rowSums(ydiff^2))-1/(n[1]^2-n[1])*(rowSums(ydiff))^2)
    t<-rowMeans(ydiff)/(den/sqrt(n[1]))
    rnum<-n[1]*rowSums(DM*CM)-rowSums(DM)*rowSums(CM)
    rden<-sqrt((n[1]*rowSums(DM^2)-(rowSums(DM))^2)*(n[1]*rowSums(CM^2)-
                                                       (rowSums(CM))^2))
    r<-rnum/rden
    d<-t*sqrt(2*(1-r)/n[1])
    m<-n[1]-1
    cm = gamma(min(m/2, 100))/(sqrt(m/2) * gamma(min((m - 1)/2,100) ))
    dprime=cm*d
    vard=(2*(1-r)/n[1])*((n[1]-1)/(n[1]-3))*
      (1+n[1]*dprime^2/(2*(1-r)))-dprime^2/cm^2
    vardprime=cm^2*vard
  }else{
    ind<-diag(rep(1,length(n)))[l,]                          ## diag()
    ym<-y%*%ind%*%diag(1/n)
    ntilde<-1/sum(1/n)
    m=sum(n)-2
    cm = gamma(min(m/2,100))/(sqrt(m/2) * gamma(min((m - 1)/2,100) ))
    s<-sqrt((1/(sum(n)-2)*((y^2%*%ind)%*%diag(1/(n-1))-
                             ym^2%*%diag(n/(n-1)))%*%(n-1)))
    d<-(ym[,2]-ym[,1])/s
    dprime=d-3*d/(4*(sum(n)-2)-1)
    terme1=1/ntilde
    vard = terme1 + d^2 * (terme1 * ntilde - 1/cm^2)
    vardprime=sum(1/n)+dprime^2/(2*sum(n))
  }
  result = cbind( dprime, vardprime)
  colnames(result) = c( "dprime", "vardprime")
  rownames(result)<-rownames(y)
  result
}


get.FEM<-function(em,vm,n, pe=NULL,pv=NULL){
  wt<-1/vm                                      ## ## if vm is NA, wt is 0
  mu.hat<-rowSums(wt*em)/rowSums(wt)
  mu.var<-1/rowSums(wt)
  z.score<-mu.hat/sqrt(mu.var)
  if(!is.null(pe)&!is.null(pv)){
    rnum<-which(apply(em,1,function(x) !any(is.na(x))))
    Z0<-matrix(get.REM2(pe,pv, n=n, REM.type="HO")$zval,nrow(em),
               nrow(pe)/nrow(em))
    z.p<-rep(NA,nrow(em))
    z.p[rnum]<-perm.p(z.score[rnum],Z0[rnum,],"abs")
  }else{
    z.p<-2*(1-pnorm(abs(z.score)))
  }
  qval<-p.adjust(z.p,method="BH")
  res<-list(mu.hat=mu.hat,mu.var=mu.var,zval=z.score,pval=z.p,FDR=qval)
  return(res)
}


MetaDE.ES<-function(x) {

  #meta.method<-match.arg(meta.method,c("FEM","REM"))
  K<-ncol(x$ES)
  n <- attr(x,"nperstudy")
  res<-get.FEM(x$ES,x$Var,n=n, x$perm.ES,x$perm.Var)
  tempFDR<-matrix(res$FDR,ncol=1)
  rownames(tempFDR)<-rownames(x$ES)
  colnames(tempFDR)<-"FEM"
  meta.res<-list(mu.hat=res$mu.hat,mu.var=res$mu.var,zval=res$zval,pval=res$pval,FDR=tempFDR)

  attr(meta.res,"nstudy")<-K

  return(meta.res)
}




