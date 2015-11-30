##############################
#Code to support the Aim 2 tool
#which is for target prioritization
#
#Pichai Raman
#11/28/2015
##############################

#call libraries
library("stringr");


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
#####################################################



#########Functions to add data and create score#############

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
tmpOut <- length(intersect(tmpOut, cancGene[,1]));
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

#Add regulation piece
output <- isRegByCancerTxnFactor(output);

return(output);
}

