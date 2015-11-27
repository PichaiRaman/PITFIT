
#############################
#Code to parse TCGA Data
#Author : Pichai Raman
#Date : 11/8/2015
#############################

#Call libraries
library("gdata");
library("tidyr");

#Find name of top directory
BaseDir <- "/srv/shiny-server/PITFIT/data/TCGA_Data/"

#Vector to hold all Features
featureVector <- c();

#Vector of cancers & lists to store all data
cancVec <- c("ov", "paad", "prad")
clinDataList <- list();
mutDataList <- list();
cnaDataList <- list();
expDataList <- list();

for(i in 1:length(cancVec))
{
#First set basedir
setwd(paste(BaseDir, cancVec[i], "/", "tcga", sep=""));

#Get Expression data 
expDataTmp <- read.delim("data_RNA_Seq_v2_expression_median.txt", stringsAsFactors=F);
rownames(expDataTmp) <- expDataTmp[,1];
expDataTmp <- expDataTmp[-1:-2];
expDataList[[cancVec[i]]] <- expDataTmp;
featureVector <- c(featureVector, paste(rownames(expDataTmp), "Exp", sep="_"));

#Get Copy Number data 
cnaDataTmp <- read.delim("data_log2CNA.txt", stringsAsFactors=F);
rownames(cnaDataTmp) <- cnaDataTmp[,1];
cnaDataTmp <- cnaDataTmp[-1:-2];
cnaDataList[[cancVec[i]]] <- cnaDataTmp;
featureVector <- c(featureVector, paste(rownames(cnaDataTmp), "Cna", sep="_"));


#Get Mutation data 
mutDataTmp <- read.delim("data_mutations_extended.txt", header=T, skip=1);
mutDataTmp <- mutDataTmp[,c("Hugo_Symbol", "IMPACT", "Tumor_Sample_Barcode")];
mutDataTmp[,"tmpConc"] <- paste(mutDataTmp[,"Hugo_Symbol"], mutDataTmp[,"Tumor_Sample_Barcode"], sep="");
mutDataTmp <- mutDataTmp[mutDataTmp[,2]!="",];
mutDataTmp[,2] <- factor(mutDataTmp[,2], levels=c("HIGH", "MODERATE","MODIFIER", "LOW"));
mutDataTmp <- mutDataTmp[order(mutDataTmp[,2]),];
mutDataTmp <- mutDataTmp[!duplicated(mutDataTmp[,4]),];
mutDataTmp <- mutDataTmp[1:3];
mutDataTmp <- spread(mutDataTmp, Tumor_Sample_Barcode, IMPACT)
mutDataTmp <- as.matrix(mutDataTmp)
mutDataTmp[is.na(mutDataTmp)]<- "NONE";
mutDataTmp <- as.data.frame(mutDataTmp);
rownames(mutDataTmp) <- mutDataTmp[,1];
mutDataTmp <- mutDataTmp[-1];
mutDataList[[cancVec[i]]] <- mutDataTmp;
featureVector <- c(featureVector, paste(rownames(mutDataTmp), "Mut", sep="_"));

#Get Clinical data 
clinDataTmp <- read.delim("data_bcr_clinical_data.txt");
clinCols <- c("PATIENT_ID", "GENDER", "CLINICAL_STAGE", "DAYS_TO_BIRTH", "OS_STATUS", "OS_MONTHS", "DFS_STATUS", "DFS_MONTHS")
clinDataTmp <- clinDataTmp[,clinCols];  
clinDataList[[cancVec[i]]] <- clinDataTmp;
print(paste("Done ", cancVec[i], sep=""));
}
featureVector <- unique(featureVector, c("Overall-Survival_Clin", "Progression-Free-Survival_Clin", "Age_Clin", "Stage_Clin", "Gender_Clin");

keep(clinDataList, expDataList, mutDataList,cnaDataList, featureVector, sure=T);
save.image("/bigdata/PITFIT_Data/ParsedTCGA.RData");














