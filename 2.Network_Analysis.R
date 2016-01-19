# Co-occurrence-network-analysis
################## OTU filtering, network generation, topological analysis and export OTU table ###############################

library(vegan)
library(igraph)
library(Hmisc)

Abu=read.table('NW.xls',header=T)
Abu<-as.matrix(Abu)

#1. Filtering OTUs by occurrence frequency (i.e.,number of samples an OTU is deemed as Present)
table<-Abu
table[table>0]<-1
table.generalist<-Abu[which(rowSums(table)>=12),]
Abu<-table.generalist

#2. Creating gml files of network (to be visulized in Gephi or Cytoscape)
pattern<-co_occurrence_network(Abu,0.8,0.01)

write.graph(pattern$graph1,'Pos0.8-NW.gml',format='gml')
write.graph(pattern$graph2,'Neg0.8-NW.gml',format='gml')

#3. Calcuating partial topological properties for positive co-occurrence network
g<-pattern$graph1
ecount(g)
degree(g, v = V(g), mode="all")
betweenness.centrality<-betweenness(g, v=V(g), directed = FALSE, weights = NULL,
            nobigint = TRUE, normalized = FALSE)
closeness.centrality<-closeness(g, vids = V(g), mode = c("out", "in", "all", "total"),
          weights = NULL, normalized = FALSE)

#plot(degree_distribution(g, cumulative = FALSE))

#4. Creating an abundance table for OTUs present in the positive and negative network
my.list1 <- row.names(pattern$matrix.cor1)
my.list2 <- row.names(pattern$matrix.cor2)

logical1 <- row.names(Abu)  %in% my.list1
logical2 <- row.names(Abu)  %in% my.list2

tab.subset1 <- subset(Abu,logical1)
tab.subset2 <- subset(Abu,logical2)

write.table(tab.subset1,'Pos0.8-NW.txt',sep="\t")
write.table(tab.subset2,'Neg0.8-NW.txt',sep="\t")
