
#############################
#Code to parse TCGA Data
#Author : Pichai Raman
#Date : 11/8/2015
#############################

#Find name of top directory
BaseDir <- "/srv/shiny-server/PITFIT/data/"
topDir <- list.files(BaseDir)

#Vector of cancers & lists to store all data
cancVec <- c("OV", "PAAD", "PRAD", "READ")
clinDataList <- list();
mutDataList <- list();
cnaDataList <- list();
expDataList <- list();

#Function to collaspe mutation files into 1
collapseMut <- function()
{
	
myFiles <- list.files()[grep("TCGA", list.files())];

allDat <- data.frame(myFiles[i], read.delim(myFiles[i], skip=3))

for(i in 2:length(myFiles))
{	
	tmpDat <- data.frame(myFiles[i], read.delim(myFiles[i], skip=3))
	allDat <- rbind(allDat, tmpDat);
}

return(allDat);
}



for(i in 1:length(cancVec))
{
#Get Clinical data 
setwd(paste(BaseDir , topDir, "/", cancVec[i], sep=""));
midDir <- list.files()
setwd(paste(getwd(), "/", midDir, "/", sep=""));
clinDirN <- list.files()[grep("Clinical", list.files())];
setwd(paste(getwd(), "/", clinDirN, "/", sep=""));
clinDataTmp <- read.delim(paste(getwd(), "/", list.files()[grep('merged_only_clinical_clin_format.txt', list.files())], sep=""));
rownames(clinDataTmp) <- clinDataTmp[,1];
clinDataTmp <- clinDataTmp[-1];
clinDataTmp <- data.frame(t(clinDataTmp));
clinDataList[[cancVec[i]]] <- clinDataTmp;

#Get Expression data 
setwd(paste(BaseDir , topDir, "/", cancVec[i], sep=""));
midDir <- list.files()
setwd(paste(getwd(), "/", midDir, "/", sep=""));
expDirN <- list.files()[grep("rnaseqv2", list.files())];
setwd(paste(getwd(), "/", expDirN, "/", sep=""));
expDataTmp <- read.delim(paste(getwd(), "/", list.files()[grep('rnaseqv2__illuminahiseq', list.files())], sep=""));
expDataTmp <- expDataTmp[-1,];
expDataList[[cancVec[i]]] <- expDataTmp;

#Get Mutation data 
setwd(paste(BaseDir , topDir, "/", cancVec[i], sep=""));
midDir <- list.files()
setwd(paste(getwd(), "/", midDir, "/", sep=""));
mutDirN <- list.files()[grep("Mutation", list.files())];
setwd(paste(getwd(), "/", mutDirN, "/", sep=""));
mutDataTmp <- collapseMut();
mutDataList[[cancVec[i]]] <- mutDataTmp;

#Get Copy Number data 
setwd(paste(BaseDir , topDir, "/", cancVec[i], sep=""));
midDir <- list.files()
setwd(paste(getwd(), "/", midDir, "/", sep=""));
cnaDirN <- list.files()[grep("Merge_snp", list.files())];
setwd(paste(getwd(), "/", cnaDirN, "/", sep=""));
cnaDataTmp <- read.delim(paste(getwd(), "/", list.files()[grep('snp', list.files())], sep=""));
cnaDataTmp <- data.frame(cnaDataTmp);
cnaDataList[[cancVec[i]]] <- cnaDataTmp;
}

keep(clinDataList, expDataList, mutDataList,cnaDataList, sure=T);
save.image("/bigdata/PITFIT_Data/ParsedTCGA.RData");














