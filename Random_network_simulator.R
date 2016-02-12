library(igraph)

n=2194
e=44680
  
for (i in 1:1000) {
  g <- erdos.renyi.game(n, e,'gnm',weight=T,mode="undirected")
  
  # Global toplogical features
  c <- cluster_walktrap(g)
  md <- modularity(g, membership(c), weights = NULL)
  cc <- transitivity(g, vids = NULL,
                     weights = NULL)
  spl <- average.path.length(g, directed=FALSE, unconnected=TRUE)
  gd  <- graph.density(g, loops=FALSE)
  nd  <- diameter(g, directed = FALSE, unconnected = TRUE, weights = NULL)
  
  ND <- degree(g, v = V(g), mode="all")
  ad  <- mean(node.degree)
  
  global.topol <- data.frame(n,e,cc,spl,md,gd,nd,ad)
  
  write.table(global.topol, file = sprintf("N%dE%d.er.random.network.xls",n,e), 
              append = TRUE, sep = "\t",row.names = FALSE, col.names = TRUE) }

# print node distribution statistics

degree <- data.frame(table(degree=factor(ND, levels=seq_len(max(ND)))))
plot(degree)


