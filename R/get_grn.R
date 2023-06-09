
## get network

## to-do: 此处data为DEG的expr，记得检查

#' Title
#'
#' @param data
#' @param label
#' @param species
#' @param database_type
#' @param database
#' @param method
#' @param alpha
#'
#' @return
#' @export
#'
#' @examples
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
    return(gD)
}



##database -2023/4/28 写完功能+配套check函数
## to-do: user-defined 接口功能+配套check函数


#' Title
#'
#' @param DEG_list
#' @param species
#' @param database_type
#' @param database
#'
#' @return
#' @export
#'
#' @examples
get_prior_grn <- function(DEG_list, species, database_type, database){
  check_support(species, database)
  check_id(DEG_list)
  if(database_type == "pathway"){
    prior_grn <- get_pathway_grn(species, database)
    prior_grn <- induced_graph(prior_grn, DEG_list)
  }
  else if(database_type == "tf"){
    prior_grn <- get_tf_grn(species, database,top=100)
    prior_grn <- induced_graph(prior_grn, DEG_list)
  }
  else if(database_type == "protein"){
    prior_grn <- get_protein_grn(DEG_list, species, database, score_threshold = 400)
  }
  else{
    print("Please enter the correct database type")
  }
  return(prior_grn)
}






## pathway database

#' Title
#'
#' @param DEG_list
#' @param species
#' @param database
#'
#' @return
#' @export
#'
#' @examples
get_pathway_grn <- function(DEG_list,species, database){
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
  pathway_net <- igraph::graph_from_data_frame(x)
  check_idmatch(DEG_list)

  return(pathway_net)
}


## tf database
## to-do: 统一一下各个数据库中先验网络的变量名称
## to-do: 改check_id()函数，使其输入参数

#' Title
#'
#' @param species
#' @param database
#' @param top
#'
#' @return
#' @export
#'
#' @examples
get_tf_grn <- function(species, database,top=100){
  check_support(species, database)

  payload <- list(query_name = "myQuery", gene_set = DEG_list)
  response <- httr::POST(url = "https://maayanlab.cloud/chea3/api/enrich/", body = payload, encode = "json")
  json <- httr::content(response, "text")
  results <- jsonlite::fromJSON(json)
  tf_net <- data.frame(TF=0,Target=0)
  for (i in 1:top) {
    tmp <- data.frame(TF=results$`Integrated--meanRank`[i,3],Target=strsplit(results$`Integrated--meanRank`[i,6],","))
    colnames(tmp) <- c("TF","Target")
    tf_net <- rbind(tf_net,tmp)
  }
  tf_net <- igraph::graph_from_data_frame(tf_net)
  check_idmatch(DEG_list)

  return(tf_net)
}




## protein database STRING

#' Title
#'
#' @param DEG_list
#' @param species
#' @param database
#'
#' @return
#' @export
#'
#' @examples
get_protein_grn <- function(DEG_list, species, database){
  check_support(species, database)

  string_db <- STRINGdb$new(version = "10", species = 10116, score_threshold = 400, input_directory = "")
  data_mapped <- string_db$map(data.frame(DEG_list), "geneid", removeUnmappedRows = TRUE)$STRING_id
  protein_interactions <- string_db$get_interactions(data_mapped)
  protein_net <- igraph::graph_from_data_frame(protein_interactions)

  return(protein_net)
}










#' Title
#'
#' @param expr
#' @param FEMR.res
#' @param cutoff
#'
#' @return
#' @export
#'
#' @examples
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




#' Title
#'
#' @param expr
#' @param species
#'
#' @return
#' @export
#'
#' @examples
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

