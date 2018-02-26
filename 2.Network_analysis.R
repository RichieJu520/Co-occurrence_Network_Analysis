###@author: Feng Ju
###@email: richieju520@gmail.com
###@cite Ju F, Xia Y, Guo F, Wang ZP, Zhang T. 2014. 
###@Taxonomic relatedness shapes bacterial assembly in activated sludge of globally distributed wastewater treatment plants.
###@Environmental Microbiology. 16(8):2421-2432

# Co-occurrence-network-analysis
################## OTU filtering, network generation, topological analysis and export OTU table ###############################
library(igraph)
library(Hmisc)

Abu=read.table('NW.txt',header=T)
Abu<-as.matrix(Abu)

###1. Filtering OTUs by occurrence frequency (i.e.,number of samples an OTU is Present)
table<-Abu
table[table>0]<-1
table.generalist<-Abu[which(rowSums(table)>=12),]
Abu<-table.generalist

###2. Creating gml files of network (to be visulized in Gephi or Cytoscape)
pattern<-co_occurrence_network(Abu,0.6,0.01)  ## cutoffs for correlation coefficient and P-value

write.graph(pattern$graph1,'Pos0.6-NW.gml',format='gml')    #network file for positive association
#write.graph(pattern$graph2,'Neg0.6-NW.gml',format='gml')   #network file for negative association (if any)
write.graph(pattern$graph3,'PosNeg0.6-NW.gml',format='gml') #network file for all association

###3. Calculating network topological properties
g<-pattern$graph1   ###positive network
#g<-pattern$graph1   ###negative network

c <- cluster_walktrap(g)
# Global toplogical features
modularity(c)
md <- modularity(g, membership(c), weights = NULL)
cc <- transitivity(g, vids = NULL,
                   weights = NULL)
spl <- average.path.length(g, directed=FALSE, unconnected=TRUE)
gd  <- graph.density(g, loops=FALSE)
nd  <- diameter(g, directed = FALSE, unconnected = TRUE, weights = NA)

node.degree <- degree(g, v = V(g), mode="all")
ad  <- mean(node.degree)

e <- ecount(g)
v <- vcount(g)
global.topology <- data.frame(e,v,cc,spl,md,gd,nd,ad)
write.csv(global.topology, file="Pos0.6-NW-global.topology.csv")

# Node toplogical features
betweenness.centrality <- betweenness(g, v=V(g), 
                                      directed = FALSE, weights = NA,
                                      nobigint = TRUE, normalized = FALSE)
closeness.centrality <- closeness(g, vids = V(g),
                                  weights = NA, normalized = FALSE)
node.transitivity <- transitivity(g, type = c("local"), vids = NULL,
                                  weights = NA)

node.topology <- data.frame(node.degree, betweenness.centrality, closeness.centrality, node.transitivity)
write.csv(node.topology, file="Pos0.6-NW-node.topology.csv")

# Ploting node degreee distribution in a log-log plot
degree.df <- data.frame(table(degree=factor(node.degree, levels=seq_len(max(node.degree)))))
degree.df$degree <- as.numeric(as.character(degree.df$degree))

#4. Creating an abundance table for OTUs present in the positive and negative network
my.list1 <- row.names(pattern$matrix.cor1)
###my.list2 <- row.names(pattern$matrix.cor2)

logical1 <- row.names(Abu)  %in% my.list1
###logical2 <- row.names(Abu)  %in% my.list2

tab.subset1 <- subset(Abu,logical1)
###tab.subset2 <- subset(Abu,logical2)

write.table(tab.subset1,'Pos0.6-NW.txt',sep="\t")
###write.table(tab.subset2,'Neg0.6-NW.txt',sep="\t")
