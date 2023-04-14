
## net visualization

# pathway network
net <- induced_subgraph(gD, vids = names(membership)[membership == cnames[1]])
#net <- gD
tmp_id <- as_ids(V(net))
#group_id <- rep("Target",length(tmp_id))
#group_id[which(as_ids(V(net)) %in% TFs)] <- "TF"

#meta.FEM.ES <- MetaDE.ES(ind.ES,meta.method = "FEM")
#table(meta.FEM.ES$zval[tmp_id]>0)
#nodes_color <- meta.FEM.ES$zval[tmp_id]
#nodes_color[which(nodes_color>0)] <- "pink"
#nodes_color[which(nodes_color<0)] <- "lightblue"

#edges_color <- E(net)$es
#edges_color[which(edges_color>0)] <- "red"
#edges_color[which(edges_color<0)] <- "blue"
#edges_size <- abs(E(net)$es)


vis_nodes <- data.frame(id=tmp_id,group=group_id,label=tmp_id, 
                        value = 2*degree(net))
vis_edges <- cbind(igraph::as_data_frame(net),
                   color = E(net)$color,
                   dashes = c(E(net)$color == "green")
                   )
visNetwork(vis_nodes, vis_edges) %>% 
  visEdges(arrows ="to")  %>%                           # arrow "to" for all edges
  visGroups(groupname = "Target", shape = "square") %>%      # square for group "A"
  visGroups(groupname = "TF", shape = "triangle")           


# tf network
##区分tf和target



## module analysis visualization
#gprofiler2

