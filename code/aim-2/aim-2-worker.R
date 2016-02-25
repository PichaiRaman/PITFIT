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

#Read & format normal tissue data###################

normExp <- read.delim("../../data/normalExpressionScore.txt");
rownames(normExp) <- normExp[,1];
#####################################################

#Read & format normal cancer Relevence data###################

cancRel <- read.delim("../../data/CancerRelSurv.txt");
#####################################################


#########Functions to add data and create score#############

###############################################################
#1. Oncogenic Proximity Section
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

###############################################################
#2. Oncogenic Regulation Section
###############################################################

#Function to see if known to be regulated by a txnFactor
isRegByCancerTxnFactor <- function(output)
{
tmpGenes <- as.character(output[,"Gene"]);
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
ifelse(tmpOut>1, 1, 0);
}


###############################################################
#End Oncogenic Regulation Section
###############################################################



###############################################################
#3. Cancer Relevance Section
###############################################################

isCancRel <- function(output)
{
tmpGenes <- as.character(output[,"Gene"]);
survPvals <- cancRel[tmpGenes,unique(output[,"Cancer"])]
cancRelScore <- ifelse(survPvals<.05, 1, 0);
output <- data.frame(output, survPvals, cancRelScore);
}
###############################################################
#End Cancer Relevance Section
###############################################################



###############################################################
#4. Druggability Section
###############################################################

#Function to add druggability score
isDruggable <- function(output)
{
tmpDruggable <- merge(output, druggable, by.x="Gene", by.y="Gene", all.x=T);
tmpDruggable[is.na(tmpDruggable[,"Sources_Druggability"]),"Sources_Druggability"] <- "NONE";
tmpDruggable[is.na(tmpDruggable[,"Score_Druggability"]),"Score_Druggability"] <- 0;
return(tmpDruggable);
}

###############################################################
#End Druggability Section
###############################################################


###############################################################
#5. TM Section
###############################################################


#Function to see if it is a TM Protein 
isTM <- function(output)
{
tmpGenes <- as.character(output[,"Gene"]);
myOut <- tmData[tmpGenes,6];
myOut[is.na(myOut)] <- 0;
ifelse(myOut==1, "Yes", "No");
myOut <- data.frame(myOut);
colnames(myOut) <- c("isTM");
myOut <- data.frame(output, myOut);
}


###############################################################
#5. End TM Section
###############################################################



###############################################################
#6. Norm Tissue Section
###############################################################


#Function to see if it is a TM Protein 
normExpProf <- function(output)
{
tmpGenes <- as.character(output[,"Gene"]);
myOut <- normExp[tmpGenes,5];
myOut[is.na(myOut)] <- 0;
myOut <- data.frame(myOut);
colnames(myOut) <- c("NormExpProfile");
myOut <- data.frame(output, myOut);
}


###############################################################
#6. End Norm Tissue 
###############################################################





##############################
#Main Function to prioritize and rank
#Targets
#
##############################
pitfitAnalyzeAim2 <- function(myCancer, geneList)
{
output <- data.frame(myCancer, geneList);
colnames(output) <- c("Cancer", "Gene");

#Add 1. Oncogenic Proximity piece
#output <- minDistOncogene(output);

#Add 2. Oncogenic regulation piece
output <- isRegByCancerTxnFactor(output);

#Add 3. Cancer Relevance piece
output <- isCancRel(output);

#Add 4. Druggability piece
output <- isDruggable(output);

#Add 5. TM piece
output <- isTM(output);

#Add 6. Normal Expression piece
output <- normExpProf(output);

return(output);
}

