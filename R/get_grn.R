
## get network

## to-do: 此处data为DEG的expr，记得检查
get_grn <- function(data, label, species, database_type, database, method="BH", alpha=0.05){
    DEG_list <- rownames(data)
    prior_grn <- get_prior_grn(DEG_list, species, database)
    check_id(DEG_list)

    tmp <- induced_subgraph(prior_grn,vids = rownames(data)[colnames(t(data)) %in% V(reactome)$name])
    dest <- SEMgraph::SEMace(tmp,t(data),label,type = "parents",effect = "direct",method = method, alpha = alpha, boot = NULL)
    ftm <- data.frame(from = dest$source, to = dest$sink)
    gD = graph_from_data_frame(ftm)

    dest$color <- dest$d_z
    dest$color[which(dest$color>0)] <- "red"
    dest$color[which(dest$color<0)] <- "green"
    E(gD)$color <- dest$color
    V(gD)$degree <- degree(gD)

}

##database -2023/4/28 写完功能+配套check函数
## to-do: user-defined 接口功能+配套check函数

get_prior_grn <- function(data, species, database_type, database){
    DEG_list <- rownames(data)
    
    ## pathway database

    ## tf database

    ## protein database


    check_id(DEG_list)
}



## pathway database 
get_pathway_grn <- function(DEG_list,species, database){                                 ## to-do: 把DEG_list整合到函数里
  check_support(species, database)
  species_pathways <- graphtie::pathways(species, database)
  graph <- vector("list",length(species_pathways))
  for (i in 1:length(species_pathways)) {
    pathway <- species_pathways[[i]]
    g <- pathwayGraph(pathway)
    g <- graph_from_graphnel(g)
    graph[[i]] <- g
  }
  x <- as_data_frame(graph[[1]])
  for (i in 1:length(graph)-1) {
    x <- rbind(x,as_data_frame(graph[[i+1]]))
  }
  x <- x%>% unique()
  src <- destin <- vector("numeric",nrow(x))
  for (i in 1:nrow(x)) {
    src[i] <- unlist(strsplit(x[i,1],":"))[2]
    destin[i] <- unlist(strsplit(x[i,2],":"))[2]
  }
  x[,1] <- src
  x[,2] <- destin
  prior_grn <- igraph::graph_from_data_frame(x)
  check_id(prior_grn)
  return(prior_grn)
}


## tf database 
## to-do: 统一一下各个数据库中先验网络的变量名称
## to-do: 改check_id()函数，使其输入参数

get_tf_grn <- function(data, DEG_list, top=100, species, database){
  gene_list <- rownames(data)
  check_support(species, database)
  check_id()
  
  #POST to ChEA3 server
  payload <- list(query_name = "myQuery", gene_set = DEG_list)
  response <- httr::POST(url = "https://maayanlab.cloud/chea3/api/enrich/", body = payload, encode = "json")
  json <- httr::content(response, "text")
  results <- jsonlite::fromJSON(json)

  ## 1. 选择top100的tf（在这一步设置tf数据库的选择参数）
  switch

  ## 2. 在gene_list里搜寻在top100的tf
  gene_list <- rownames(data)
  gene_list[gene_list %in% results$`Integrated--meanRank`[1:top,3]]

  ## 3. 取调控关系，转换成图（以ChEA3为例）
  tf_net <- data.frame(TF=0,Target=0)
  for (i in 1:top) {
    tmp <- data.frame(TF=results$`Integrated--meanRank`[i,3],Target=strsplit(results$`Integrated--meanRank`[i,6],","))
    colnames(tmp) <- c("TF","Target")
    tf_net <- rbind(tf_net,tmp)
  }
  tf_net <- igraph::graph_from_data_frame(tf_net)
  return(tf_net)
  
}




## protein database STRING











## get_DEG_expr
get_DEG_expr <- function(expr, FEMR.res, cutoff = 0.05){
  if (!is.numeric(cutoff_FEM) | length(cutoff_FEM) != 1){
    stop("cutoff_FEM should be a single number")
  }
  if (cutoff_FEM <= 0 | cutoff_FEM >= 1){
    stop("cutoff_FEM should be between 0 and 1")
  }

  if(!is.null(rownames(expr))){
    index <- which(FEMR.res$FDR < cutoff)
    gene_id <- names(FEMR.res[index])
    DEG_expr <- expr[which(rownames(expr) %in% gene_id),]
  }else{
    stop("Please entering a matrix with rownames")
  }

  return(DEG_expr)
}



## get_TF
get_TF <- function(expr, species){
  if(!is.null(species)){
    TF_dir <- paste0("../FEMboot/",species,".txt")
    TF_list <- read.csv(TF_dir)
    if(!is.null(rownames(expr))){
      TFs <- rownames(expr)[which(rownames(expr) %in% TF_list)]
    }
    else{
      stop("Please entering a matrix with rownames")
    }
  }else{
    stop("Please entering your species")
  }

  return(TFs)
}

