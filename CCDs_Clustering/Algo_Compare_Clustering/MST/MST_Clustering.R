library(igraph)

# use MST for clustering
# data: the input data set
# thresh: the ratio to cut a edge when comparing its adjacent edges
MST_Clustering = function(data, thresh=2){
  n = dim(data)[1]
  
  # construct a full graph
  g = graph_from_adjacency_matrix(as.matrix(dist(data)), mode = "undirected", weighted = TRUE)
  
  # Calculate MST
  mst = mst(g)
  
  # the end points of each edge
  end_pts = ends(mst, E(mst))
  
  # Calculate edge lengths
  edge_lengths = E(mst)$weight
  
  # Calculate average length of adjacent edges for each edge
  adj_avg = sapply(1:length(edge_lengths), function(i){
    end_temp = end_pts[i,]
    # the two end points of current edge
    end1 = end_temp[1]
    end2 = end_temp[2]
    adj_ind = apply(end_pts, 1, function(row){
      sum(row == end1 | row == end2) == 1})
    return(mean(E(mst)$weight[adj_ind]))
  })
  
  # Remove edges significantly longer than adjacent average
  remove_ind <- which(E(mst)$weight > thresh*adj_avg)
  mst_pruned <- delete_edges(mst, remove_ind)
  
  # Extract clusters
  cls = components(mst_pruned)$membership
  n_cls = unique(cls) # number of clusters
  size_cls = table(cls) # the sizes of each clusters
  return(labels = cls)
}