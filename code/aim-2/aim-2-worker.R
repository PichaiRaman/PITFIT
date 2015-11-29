##############################
#Code to support the Aim 2 tool
#which is for target prioritization
#
#Pichai Raman
#11/28/2015
##############################

#call libraries


#Read & format druggability data#####################
drugInt <- read.delim("../../data/categories.tsv", stringsAsFactors=F);
druggable <- drugInt[drugInt[,"category"]=="DRUGGABLE GENOME",]
rownames(druggable) <- druggable[,1];
druggable[,"Drug_Score"] <- str_count(druggable[,"category_sources"], ",")+1
druggable <- druggable[,c("entrez_gene_symbol", "category_sources", "Drug_Score")]
colnames(druggable) <- c("Gene", "Sources_Druggability", "Score_Druggability");
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
return(output);
}

