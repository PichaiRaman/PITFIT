##############################
#Code to create R object for Aim 2
#Pichai Raman
#11/28/2015
############################

#Call library
library("sqldf")
library("RSQLite")
library("GSEABase");
library("stringr");

#First let's start parse genemania data
setwd("/home/ramanp/pitfit/data/GeneMania/genemania.org/data/current/Homo_sapiens");
allFiles <-list.files();

#Network file to get info on netorks
networks <- read.delim("networks.txt"); #not used currently

#Mapping file to convert ID's
mappingFile <- read.delim("identifier_mappings.txt");
mappingFile <- mappingFile[mappingFile[,3]=="Gene Name",1:2];

#now let's read in gmt files
gmtFiles <- allFiles[grep(".gmt", allFiles)];
interpro <- geneIds(getGmt("Attributes.InterPro.gmt")); #not used currently
txnFactor <- geneIds(getGmt("Attributes.Transcriptional-factor-targets-2013.gmt"));
drugs <- geneIds(getGmt("Attributes.Drug-interactions-2013.gmt")); #not used currently
cPathways <- geneIds(getGmt("Attributes.Consolidated-Pathways-2013.gmt"));

#Function to pull out just gene symbol TF's
formatTxn <- function(x)
{
x <- as.character(sapply(x, FUN=pullGene));
x <- as.character(na.omit(x));
}
pullGene <- function(x)
{
y <- substring(x, (str_locate(x, "\\$")[1]+1), str_length(x));
y <- substring(y, 1, (str_locate(y, "_")[1]-1))
return(y);
}
txnFactor_genes <- lapply(txnFactor, FUN=formatTxn);

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

#Convert to gene symbol
names(drugs) <- convertID(names(drugs));
names(txnFactor) <- convertID(names(txnFactor));
names(txnFactor_genes) <- convertID(names(txnFactor_genes));
names(interpro) <- convertID(names(interpro));
names(cPathways) <- convertID(names(cPathways));


print("Finished writing R Object");
print("Started Reading files");
interactionFilesType <- c("Physi", "Predi", "Genet", "Co-lo", "Co-ex", "Pathw");
interactionFiles <- allFiles[substring(allFiles, 1,5)%in%interactionFilesType];
interactionDF <- data.frame(read.delim(interactionFiles[1]), gsub(".txt", "", interactionFiles[1]))
colnames(interactionDF)[4] <- "source"
interactionDF[,"Type"] <- pullOutStem(interactionDF[,"source"]);

#Map it to Gene Symbol
interactionDF[,1] <- convertID(interactionDF[,1]);
interactionDF[,2] <- convertID(interactionDF[,2]);

setwd("/bigdata/PITFIT_Data/");
db <- dbConnect(SQLite(), dbname="aim2.sqlite")
dbWriteTable(db, "GeneMania_Interactions", interactionDF)
for(i in 2:length(interactionFiles))
{
intDFTmp <- read.delim(paste("/home/ramanp/pitfit/data/GeneMania/genemania.org/data/current/Homo_sapiens/", interactionFiles[i], sep=""), stringsAsFactors=F);
intDFTmp[,"source"] <- gsub(".txt", "", interactionFiles[i])
colnames(intDFTmp)[4] <- "source"
intDFTmp[,"Type"] <- pullOutStem(intDFTmp[,"source"]);
intDFTmp[,1] <- convertID(intDFTmp[,1]);
intDFTmp[,2] <- convertID(intDFTmp[,2]);
dbWriteTable(db, "GeneMania_Interactions", intDFTmp, append=T)
print(i/length(interactionFiles));
}
dbDisconnect(db)
print("Finished reading file");

save.image("/home/ramanp/pitfit/data/GeneManiaROBJ.RData");



