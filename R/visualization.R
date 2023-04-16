
## net visualization

# pathway network
net <- induced_subgraph(gD, vids = names(membership)[membership == cnames[1]])
tmp_id <- as_ids(V(net))
group_id <- rep("Target",length(tmp_id))
group_id[which(as_ids(V(net)) %in% TFs)] <- "TF"

vis_nodes <- data.frame(id=tmp_id,group=group_id,label=tmp_id, 
                        value = 2*degree(net))
vis_edges <- cbind(igraph::as_data_frame(net),
                   color = E(net)$color)
visNetwork(vis_nodes, vis_edges) %>% 
  visEdges(arrows ="to")  %>%                      
  visGroups(groupname = "Target") %>%      
  visGroups(groupname = "TF", shape = "triangle")           




## module analysis visualization
#gprofiler2

