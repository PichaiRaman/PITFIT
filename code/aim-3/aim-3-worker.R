##################################
#Either find best set of lines for a single gene
#or Greedy Algorithm to find the best set of lines
#
##################################

#Read in the data
load("../../data/CCLEData.RData");

#call libraries
library(ggplot2);
library(reshape2);

cancLookUp <- function(canc)
{
  out <- "";
  if(canc=="ov")
  {
    out <- "OVARY";
  }
  if(canc=="paad")
  {
    out <- "PANCREAS";
  }
  if(canc=="prad")
  {
    out <- "PROSTATE";
  }
return(out);
}

#Function to use if single gene
getSingleBestLine <- function(canc, gene)
{
canc <- cancLookUp(canc);
tmpDat <- data.frame(t(ccleData[gene,grep(canc, colnames(ccleData))]));
tmpDat[,"Cell_Line"] <- rownames(tmpDat);
colnames(tmpDat)[1]<- "Value";
tmpDat <- tmpDat[order(-tmpDat[,"Value"]),];
tmpDat[,"Cell_Line"] <- gsub(paste("_", canc, sep=""), "", tmpDat[,"Cell_Line"]);
tmpDat[,"Cell_Line"] <- factor(tmpDat[,"Cell_Line"], levels=tmpDat[,"Cell_Line"]);
tmpDat <- tmpDat[,c(2,1)];
p <- ggplot(tmpDat, aes(Cell_Line, 2^Value))+geom_bar(stat="identity")+theme_bw();
p <- p+xlab("Cell Line")+ylab("Expression Value")
p <- p+theme(axis.text.x=element_text(angle = -75, hjust = 0));
return(list(tmpDat,p));

}

##################################
#Function to binarize matrix
myZ <- function(x) 
{ 
x <- (x-mean(x))/sd(x);
x <- ifelse(x<1, 0, 1);
}
createBinMat <- function(x)
{
x <- data.frame(t(apply(x, FUN=myZ, MARGIN=1))); 
}
##################################


#Function to use if multiple gene
getBestSetOfLines <- function(canc, genes, numLines)
{
canc <- cancLookUp(canc);
tmpDat <- ccleData[genes,grep(canc, colnames(ccleData))];

#Turn to Binary
tmpDatBin <- createBinMat(tmpDat);
bestLines <- greedyCL(tmpDatBin, numLines);


binMat <- tmpDatBin[, bestLines];
contDat <- tmpDat[, bestLines];
contDat[,"geneName"] <- rownames(contDat);
contDat <- melt(contDat, id.var="geneName");
colnames(contDat) <- c("Gene", "Cell_Line", "Value");
contDat[,"Value"] <- 2^contDat[,"Value"];
return(list(contDat, binMat));
}





##################################
#Functions to run greedy algorithm
##################################

greedyHelper <- function(sol, matR, numLines)
{
nMatr <- matR[,!colnames(matR)%in%sol]
mySol <- matR[,colnames(matR)%in%sol]
if(class(mySol)=="data.frame")
{
mySol <- as.numeric(rowSums(mySol));
}
nMatrAdd <- as.matrix(nMatr+mySol);
nMatrAdd  <- data.frame(ifelse(nMatrAdd> numLines, numLines, nMatrAdd));
nfSum <- colSums(nMatrAdd);
nfLine <- names(sort(nfSum, T))[1]
tmpScore <- sort(nfSum, T)[1]
return(list(nfLine, tmpScore));
}

greedyCL <- function(x, numLines)
{
bestScore <- numLines*nrow(x);
#Let's choose the best line first
fSum <- colSums(x);
cnSolution <- names(sort(fSum, T))[1]
origScore <- sort(fSum, T)[1]
	for(i in 1:ncol(x))
	{
		
		out <- greedyHelper(cnSolution, x, numLines);
		cnSolution <- c(cnSolution, out[[1]]);
		nScore <- out[[2]];
		if(out[[2]]>=bestScore|out[[2]]==origScore) { break; } 
		origScore <- nScore
	}
return(cnSolution); 

}
