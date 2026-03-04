

library(igraph)
library(Hmisc)
# read into the OTU abundance table 
asv <- read.csv("Genus.csv",row.names = 1)
genus1 <-asv[,19:27]
genus2 <-asv[,19:27]
genus1[genus1>0] <- 1
rt1<- genus2[which(rowSums(genus1) >=2), ]    

# to calculate whether there is a correlation between the abundance changes between the two genera, taking the spearman correlation coefficient as an example.
net_corr <- rcorr(t(rt1), type = 'pearson')
r <- net_corr$r
r[abs(r) < 0.6] <- 0
r[is.na(r)] <- 0

p <- net_corr$P
p <- p.adjust(p, method = 'BH')    
p[is.na(p)] <- -1
p[p>=0.05] <- -1
p[p<0.05 & p>=0] <- 1
p[p==-1] <- 0


z <- r * p
diag(z) <- 0    
head(z)[1:6,1:6]


# Construct a weighted undirected network, and the weight represents the spearman correlation coefficient of microbial abundance.
g <- graph.adjacency(z, weighted = TRUE, mode = 'undirected',diag = FALSE)
g

g <- simplify(g)

#Delete isolated nodes ( delete nodes with a degree of 0 )
g <- delete.vertices(g, names(degree(g)[degree(g) == 0]))
g
plot(g)

E(g)$correlation <- E(g)$weight
E(g)$weight <- abs(E(g)$weight)

#vertex_phylum <- tax_network2$Phylum[match(tax_network2$Genus,V(g)$name)]
#V(g)$Phylum <- vertex_phylum
E(g)
V(g)$label <- V(g)$name
V(g)$degree <- degree(g)
#library(rgexf)
#write.graph(g,"Fun.gml",format = "gml")

# View network diagram
g
plot.igraph(g,
            vertex.size=4, 
            vertex.label=NA,
            edge.curved=T,
            edge.size=1.5,
            layout = layout.fruchterman.reingold)


#adj_matrix <- as.matrix(get.adjacency(g, attr = 'correlation'))
#write.table(data.frame(adj_matrix, check.names = FALSE), 'network.adj_matrix.txt', col.names = NA, sep = '\t', quote = FALSE)

# Edge list
edge <- data.frame(as_edgelist(g))    

edge_list <- data.frame(
  source = edge[[1]],
  target = edge[[2]],
  weight = E(g)$weight,
  correlation = E(g)$correlation
)
head(edge_list)
# Add a column to distinguish positive and negative correlations for Gephi beautification
a <- edge_list[,4]
a[a>0] <- 1
a[a<0] <- -1
edge_list$cor <- a
edge_list <- subset(edge_list,weight<1)

write.csv(edge_list, 'CK-B.csv',  row.names = FALSE)

#List of node attributes 
# Add attribute information for nodes ( microbial genera ) ( annotated at the level of family and genus of phylum ) 
# ' genus _ taxonomy.txt ' records the attributes of microorganisms. After reading the table, the corresponding rows are matched according to the known network nodes.
node_list <- data.frame(ID = names(V(g)))
tax <- read.table("Tax2.txt",sep="\t", check.names = FALSE, stringsAsFactors = FALSE)
colnames(tax)=tax[1,]
tax<-tax[-1,]
#tax_network <- tax[which(!is.na(tax$Genus)),]## Removal of genus level NA

#tax_network1 <- tax_network[!duplicated(tax_network$Genus),]
#tax_network2<-tax_network[,row.names(tax_network)]
tax_network2 <- tax[tax$ID %in% V(g)$name,]
library(dplyr)
# Sorting by genus level
#tax_network2 <- tax_network2[order(tax_network2$Genus),]

row.names(tax_network2)=tax_network2$ID
tax2<-tax_network2[as.character(V(g)$name), ]
#V(g)$Kingdom <- tax2$Kingdom
V(g)$Phylum <- tax2$Phylum
#V(g)$Class <- tax2$Class
#V(g)$Order <- tax2$Order
#V(g)$Family <- tax2$Family
#V(g)$Genus <- tax2$Genus

node_list <- merge(node_list,tax2,by.x="ID",by.y="row.names")

#write.csv(tax, 'asv_node_rt2.csv',  row.names = T)

node_list <- data.frame(
  label = names(V(g)),
  #kingdom = V(g)$Kingdom,
  phylum = V(g)$Phylum
  #class = V(g)$Class,
  #order = V(g)$Order,
  #family = V(g)$Family
  #genus = V(g)$Genus
)
head(node_list)
colnames(node_list)[1] <- "Id"
#colnames(node_list)[6]<-"Label"
#
write.csv(node_list, 'CK-D.csv',  row.names = FALSE, )

