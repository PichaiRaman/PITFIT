#############################################################
#Code to generate matrix of proximity of each gene
#to each of the genes in the cancer gene census. Matrix will be
#Gene on x-axis, oncogenes as columns, as distance as values
#
#This will be a .RData object to be loaded by application
#
#############################################################



#call libraries
library("stringr");
library("RNeo4j");
library("gdata");
library("snow");

graph = startGraph("http://localhost:7474/db/data/")


#Read & format cancer data#####################
cancGeneCens <- read.delim("../../data/CancerGeneCensus.tsv", stringsAsFactors=F);
cancGene <- cancGeneCens[,c("Gene.Symbol", "Tumour.Types.Somatic.")];
colnames(cancGene) <- c("cgc_Gene", "cgc_Tumour.Type")
cancGene <- cancGene[,1];

#only take genes that are in network
query <- "match (n) return n.name"
allNetGenes <- as.character(cypher(graph, query)[,1]);
cancGene <- intersect(cancGene, allNetGenes);
#####################################################


#Read in list of gene names
featureList <- read.delim("/srv/shiny-server/PITFIT/data/TCGA_Data/ov/tcga/data_RNA_Seq_v2_expression_median.txt", stringsAsFactors=F)[,1];
featureList <- intersect(featureList, allNetGenes);




###############################################################
#Oncogenic Proximity Section
###############################################################

distGenes <- function(tmpGeneA, tmpGeneB)
{
query = paste("MATCH (p:Gene) WHERE p.name ='",tmpGeneA,"' RETURN p", sep="");
tmpNodeA <-getSingleNode(graph, query);
query = paste("MATCH (p:Gene) WHERE p.name ='",tmpGeneB,"' RETURN p", sep="");
tmpNodeB <-getSingleNode(graph, query);
p = shortestPath(tmpNodeA, "ULINK", tmpNodeB, max_depth=5)
out <- c(p[[1]]);
if(is.null(out)) { out <- 100; }
return(out);
}

distOncogene <- function(x)
{
print(x);
output <- as.numeric(sapply(cancGene, FUN=distGenes, x));
}


clus <- makeCluster(20);
print("Started cluster");
clusterExport(clus, "featureList");
clusterExport(clus, "cancGene");
clusterExport(clus, "distGenes");
clusterExport(clus, "distOncogene");
clusterExport(clus, "featureList");
clusterExport(clus, "graph");
clusterExport(clus, "getSingleNode");
clusterExport(clus, "shortestPath");
myoutput <- parSapply(clus, featureList, FUN=distOncogene);
stopCluster(clus);

OncProxMat <- data.frame(t(myoutput));
colnames(OncProxMat) <- cancGene;
keep(OncProxMat, sure=T);
save.image("OncProx.RData");













