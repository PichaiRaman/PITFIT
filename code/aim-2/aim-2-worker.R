##############################
#Code to support the Aim 2 tool
#which is for target prioritization
#
#Pichai Raman
#11/28/2015
##############################

#call libraries
library("stringr");
library("RNeo4j")
graph = startGraph("http://localhost:7474/db/data/")


#Read in list of gene names

featureList <- read.delim("/srv/shiny-server/PITFIT/data/TCGA_Data/ov/tcga/data_RNA_Seq_v2_expression_median.txt", stringsAsFactors=F)[,1];


#Load gene mania data
load("../../data/GeneManiaROBJ.RData");

#Read & format druggability data#####################
drugInt <- read.delim("../../data/categories.tsv", stringsAsFactors=F);
druggable <- drugInt[drugInt[,"category"]=="DRUGGABLE GENOME",]
rownames(druggable) <- druggable[,1];
druggable[,"Drug_Score"] <- str_count(druggable[,"category_sources"], ",")+1
druggable <- druggable[,c("entrez_gene_symbol", "category_sources", "Drug_Score")]
colnames(druggable) <- c("Gene", "Sources_Druggability", "Score_Druggability");
druggable <- druggable[1:3];
#####################################################

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

#Read & format tm data#####################
tmData <- read.delim("../../data/tmList.txt");
rownames(tmData) <- tmData[,1];
#####################################################


#########Functions to add data and create score#############

###############################################################
#Oncogenic Proximity Section
###############################################################

distGenes <- function(tmpGeneA, tmpGeneB)
{
print(tmpGeneA);
query = paste("MATCH (p:Gene) WHERE p.name ='",tmpGeneA,"' RETURN p", sep="");
tmpNodeA <-getSingleNode(graph, query);
query = paste("MATCH (p:Gene) WHERE p.name ='",tmpGeneB,"' RETURN p", sep="");
tmpNodeB <-getSingleNode(graph, query);
p = shortestPath(tmpNodeA, "ULINK", tmpNodeB, max_depth=5)
out <- c(p[[1]]);
return(out);
}

distOncogene <- function(x)
{
output <- sapply(cancGene, FUN=distGenes, x);
}


minDistOncogene <- function(x)
{




}


###############################################################
#End Oncogenic Proximity Section
###############################################################




#Function to add druggability score
isDruggable <- function(output)
{
tmpDruggable <- merge(output, druggable, by.x="Gene", by.y="Gene", all.x=T);
tmpDruggable[is.na(tmpDruggable[,"Sources_Druggability"]),"Sources_Druggability"] <- "NONE";
tmpDruggable[is.na(tmpDruggable[,"Score_Druggability"]),"Score_Druggability"] <- 0;
return(tmpDruggable);
}

#Function to see if known to be regulated by a txnFactor
isRegByCancerTxnFactor <- function(output)
{
tmpGenes <- as.character(output[,1]);
myOut <- data.frame(sapply(tmpGenes, FUN=txnRegHelper), sapply(tmpGenes, FUN=txnRegHelperScore))
colnames(myOut) <- c("TF_Regulation", "Score_Regulation");
myOut <- data.frame(output, myOut);
}
txnRegHelper <- function(geneName)
{
tmpOut <- txnFactor_genes[[geneName]]
tmpOut <- paste(tmpOut, collapse=", ")
}
txnRegHelperScore <- function(geneName)
{
tmpOut <- txnFactor_genes[[geneName]]
tmpOut <- length(intersect(tmpOut, cancGene));
}

#Function to see if it is a TM Protein 
isTM <- function(output)
{
tmpGenes <- as.character(output[,1]);
myOut <- tmData[tmpGenes,6];
myOut[is.na(myOut)] <- 0;
ifelse(myOut==1, "Yes", "No");
myOut <- data.frame(myOut);
colnames(myOut) <- c("isTM");
myOut <- data.frame(output, myOut);
}







##############################
#Main Function to prioritize and rank
#Targets
#
##############################
pitfitAnalyzeAim2 <- function(myCancer, geneList)
{
output <- data.frame(myCancer, geneList);
colnames(output) <- c("Cancer", "Gene");

#Add Druggability score to genes
output <- isDruggable(output);

#Add tm piece
output <- isTM(output);

#Add regulation piece
output <- isRegByCancerTxnFactor(output);

#Add Oncogenic Proximity
#output <- minDistOncogene(output);


return(output);
}

