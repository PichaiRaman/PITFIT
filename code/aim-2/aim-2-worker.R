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
library("lsa")
library("gtools")
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

#Read Onc Proximity Matrix ##################

load("../../data/OncProx.RData");

#####################################################

#########Functions to add data and create score#############

###############################################################
#1. Oncogenic Proximity Section
###############################################################

getOncProx <- function(output)
{
tmpGenes <- as.character(output[,"Gene"]);
tmpOncProxMat <- OncProxMat[c(tmpGenes),]
minDist <- apply(tmpOncProxMat, FUN=min, MARGIN=1);
closeOncogenes <- apply(tmpOncProxMat, FUN=getMinOncogenes, MARGIN=1);
tmpDF <- data.frame(minDist, closeOncogenes);
colnames(tmpDF) <- c("MinDist", "CloseOncogenes");
output <- cbind(output, tmpDF);
return(output);
}

getMinOncogenes <- function(x)
{
minVal <- min(x);
out <- names(x)[x==minVal]
out <- paste(out, collapse=", ");
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

###############################################################
#7. Calculate Scores
###############################################################

#create all permutaiton vector
allComb <- data.frame(permutations(4, 6, c(0:3), T, T));
allComb <- allComb[allComb[,2]%in%c(1,0),];
allComb <- allComb[allComb[,3]%in%c(1,0),];
allComb <- allComb[allComb[,4]%in%c(2,1,0),];
allComb <- allComb[allComb[,5]%in%c(1,0),];

#Filter ADC
ADC <- allComb[allComb[,5]==1,];
ADC <- ADC[ADC[,6]%in%c(3,2),];
ADC <- as.matrix(ADC);

#Filter CAR
CAR <- allComb[allComb[,5]==1,];
CAR <- CAR[CAR[,6]==3,];
CAR <- as.matrix(CAR);

#Filter LMW
LMW <- allComb[allComb[,1]==3,];
LMW <- LMW[LMW[,2]==1,];
LMW <- LMW[LMW[,3]==1,];
LMW <- LMW[LMW[,4]==2,];
LMW <- as.matrix(LMW);


MAB <- c(3,1,1,0,1,2);

getADCScore <- function(outputTmp)
{
tmpDat <- apply(outputTmp, FUN=cosine, MARGIN=1, y=ADC[1,]);
  for(i in 2:(dim(ADC)[1]))
  {
   tmpDat <-  rbind(tmpDat, apply(outputTmp, FUN=cosine, MARGIN=1, y=ADC[i,]))
  }
  finOut <- apply(tmpDat, FUN=max, MARGIN=2);
}

getCARScore <- function(outputTmp)
{
tmpDat <- apply(outputTmp, FUN=cosine, MARGIN=1, y=CAR[1,]);
  for(i in 2:(dim(ADC)[1]))
  {
   tmpDat <-  rbind(tmpDat, apply(outputTmp, FUN=cosine, MARGIN=1, y=CAR[i,]))
  }
  finOut <- apply(tmpDat, FUN=max, MARGIN=2);
}


getLMWScore <- function(outputTmp)
{
tmpDat <- apply(outputTmp, FUN=cosine, MARGIN=1, y=LMW[1,]);
  for(i in 2:(dim(ADC)[1]))
  {
   tmpDat <-  rbind(tmpDat, apply(outputTmp, FUN=cosine, MARGIN=1, y=LMW[i,]))
  }
  finOut <- apply(tmpDat, FUN=max, MARGIN=2);
}


getMABScore <- function(outputTmp)
{
tmpDat <- apply(outputTmp, FUN=cosine, MARGIN=1, y=MAB[1,]);
  for(i in 2:(dim(ADC)[1]))
  {
   tmpDat <-  rbind(tmpDat, apply(outputTmp, FUN=cosine, MARGIN=1, y=MAB[i,]))
  }
  finOut <- apply(tmpDat, FUN=max, MARGIN=2);
}


getScore <- function(output)
{
outputTmp <- output[,c("MinDist", "Score_Regulation", "cancRelScore", "Score_Druggability", "isTM", "NormExpProfile")];
output[,"ADC_Score"] <- getADCScore(outputTmp);
output[,"LMW_Score"] <- getLMWScore(outputTmp);
output[,"MAB_Score"] <- getMABScore(outputTmp);
output[,"CAR_Score"] <- getCARScore(outputTmp);
output[,"MAX_Score"] <- apply(output[13:16], FUN=max, MARGIN=1);
return(output);
}




###############################################################
#7. End Calculate scores
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
output <- getOncProx(output);

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

#Add 7. Scores
output <- getScore(output);


return(output);
}

