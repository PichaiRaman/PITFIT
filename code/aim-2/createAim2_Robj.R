##############################
#Code to create R object for Aim 2
#Pichai Raman
#11/28/2015
############################

#First let's start parse genemania data
setwd("/home/ramanp/pitfit/data/GeneMania/genemania.org/data/current/Homo_sapiens");
allFiles <-list.files();
mappingFile <- read.delim("identifier_mappings.txt");
mappingFile <- mappingFile[mappingFile[,3]=="Gene Name",1:2]

interactionFilesType <- c("Physi", "Predi", "Genet", "Co-lo", "Co-ex");
interactionFiles <- allFiles[substring(allFiles, 1,5)%in%interactionFilesType];

interactionDF <- data.frame(read.delim(interactionFiles[1]), gsub(".txt", "", interactionFiles[1]))
colnames(interactionDF)[4] <- "source"
#for(i in 2:length(interactionFiles))
for(i in 2:3)
{
intDFTmp <- data.frame(read.delim(interactionFiles[i], stringsAsFactors=F), gsub(".txt", "", interactionFiles[i]))
colnames(intDFTmp)[4] <- "source"
interactionDF <- rbind(interactionDF, intDFTmp);
print(i/length(interactionFiles));
}

#Map it to ensembl
interactionDF <- merge(interactionDF, mappingFile, by.x="Gene_A", by.y="Preferred_Name");
interactionDF <- merge(interactionDF, mappingFile, by.x="Gene_B", by.y="Preferred_Name");
interactionDF <- interactionDF[,c(5,6,3,4)];
colnames(interactionDF)[1:2] <- c("Gene_A", "Gene_B");
interactionDF[,"Type"] <- sapply(interactionDF[,"source"], FUN=pullOutStem);


pullOutStem <- function(x)
{
myEnd <- gregexpr("\\.", x)[[1]][1]-1; 
out <- substring(x, 1, myEnd);
}









#"networks.txt"
#"Pathway"
#".gmt"
