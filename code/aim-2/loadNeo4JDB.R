###########################################
#Code to load Neo4J
#
#
###########################################

#Install Neo4J first and make sure you have JDK8
#http://tecadmin.net/install-java-8-on-centos-rhel-and-fedora/

#using package RNeo4J
#https://github.com/nicolewhite/RNeo4j
#Please note when installing you may need to install openssldevel first  : yum install -y openssl-devel
#Also make sure to disable authentication in conf/neo4j-server : dbms.security.auth_enabled=false


#Call library
library(RNeo4j)

#Start a graph & if its full get rid of 
graph = startGraph("http://localhost:7474/db/data/")
clear(graph, F)


#First let's start parse genemania data, there are 553 files
setwd("/home/ramanp/pitfit/data/GeneMania/genemania.org/data/current/Homo_sapiens");
allFiles <-list.files();
interactionFilesType <- c("Physi", "Predi", "Genet", "Co-lo", "Co-ex", "Pathw");
interactionFiles <- allFiles[substring(allFiles, 1,5)%in%interactionFilesType];

#Need Cancer Gene Census
cancGeneCens <- read.delim("/home/ramanp/pitfit/data/CancerGeneCensus.tsv", stringsAsFactors=F);
cancGene <- cancGeneCens[,c("Gene.Symbol", "Tumour.Types.Somatic.")];
colnames(cancGene) <- c("cgc_Gene", "cgc_Tumour.Type")
cancGene <- cancGene[,1]

#Mapping file to convert ID's
mappingFile <- read.delim("identifier_mappings.txt");
mappingFile <- mappingFile[mappingFile[,3]=="Gene Name",1:2];






#Function to convert ID's
convertID <- function(x)
{
rownames(mappingFile) <- mappingFile[,1];
y <- mappingFile[x,2];
return(as.character(y));
}
#Function to pull out stem
pullOutStem <- function(x)
{
x <-substring(x, 1,5);
x <- gsub("Physi", "Physical_Interaction", x);
x <- gsub("Predi", "Predicted_Interaction", x);
x <- gsub("Genet", "Genetic_Interactions", x);
x <- gsub("Co-lo", "Co-localization", x);
x <- gsub("Co-ex", "Co-expression", x);
x <- gsub("Pathw", "Pathway", x);
return(x);
}

########################################################
#Load all the nodes
########################################################
print("started adding nodes");

#Find the unique nodes
uniqueGenes <- c();
for(i in 1:length(interactionFiles))
{
intDFTmp <- read.delim(paste("/home/ramanp/pitfit/data/GeneMania/genemania.org/data/current/Homo_sapiens/", interactionFiles[i], sep=""), stringsAsFactors=F);
uniqueGenes <- unique(c(uniqueGenes, intDFTmp[,1], intDFTmp[,2]));
print(paste(i, " ", length(uniqueGenes)));
}
uniqueGenes <- unique(convertID(uniqueGenes));
uniqueGenes <- data.frame(uniqueGenes);
rownames(uniqueGenes) <- uniqueGenes[,1];
uniqueGenes[,"isCancerGene"] <- uniqueGenes[,1]%in%cancGene;

for(i in 1:length(uniqueGenes[,1]))
{
createNode(graph, "Gene", name=as.character(uniqueGenes[i,1]), CancerGene=as.character(uniqueGenes[i,2]))
}

########################################################
#Load all the edges
########################################################
print("started adding edges");
#Main Function to add edges
addEdge <- function(x)
{
tmpGeneA <- as.character(x[1]);
tmpGeneB <- as.character(x[2]);
tmpWeight <- as.numeric(x[3]);
tmpSource <- as.character(x[4]);
tmpType <- as.character(x[5]);

print(paste("adding genes", tmpGeneA, tmpGeneB));
query = paste("MATCH (p:Gene) WHERE p.name ='",tmpGeneA,"' RETURN p", sep="");
tmpNodeA <-getSingleNode(graph, query);
query = paste("MATCH (p:Gene) WHERE p.name ='",tmpGeneB,"' RETURN p", sep="");
tmpNodeB <-getSingleNode(graph, query);
  if(!is.null(tmpNodeA)&&!is.null(tmpNodeB))
  {
  createRel(tmpNodeA, "ULINK", tmpNodeB, weight=tmpWeight, source=tmpSource, sourceType=tmpType)
  }
}

for(i in 1:length(interactionFiles))
{
intDFTmp <- data.frame(read.delim(paste("/home/ramanp/pitfit/data/GeneMania/genemania.org/data/current/Homo_sapiens/", interactionFiles[i], sep=""), stringsAsFactors=F), gsub(".txt", "", interactionFiles[i]));
colnames(intDFTmp)[4] <- "source"
intDFTmp[,"Type"] <- pullOutStem(intDFTmp[,"source"]);
intDFTmp[,1] <- convertID(intDFTmp[,1]);
intDFTmp[,2] <- convertID(intDFTmp[,2]);
apply(intDFTmp, FUN=addEdge, MARGIN=1);
print(paste("loaded file", i, "file size is", nrow(intDFTmp)));
}











