##################################
#Either find best set of lines for a single gene
#or Greedy Algorithm to find the best set of lines
#
##################################

#Read in the data
load("../../data/CCLEData.RData");

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

getSingleBestLine <- function(canc, gene)
{
canc <- cancLookUp(canc);
tmpDat <- ccleData[gene,grep(canc, colnames(ccleData))];
tmpDat <- tmpDat[-order(tmpDat)];

}










