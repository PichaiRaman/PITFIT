##############################
#Code to support the Aim 1 tool
#which is for target discovery.
#
#Pichai Raman
#5/31/2015
##############################

#call libraries
library("limma");
library("scales");
library("ggplot2");
library("pheatmap");
library("gplots");


#Read in data
expData <- read.delim("../data/ovarianData.txt");
annotData <- read.delim("../data/ovarianDataAnnot.txt");



SigLimmaTrain <- function(myData, myGene, thresh=.20, pvalThresh=.25, logFCThresh=1)
{

if(myData=="ov")
{
myData = expData;
}

#pull our vector of values for gene
myVect <- myData[myGene,];
myEcdf <- ecdf(myVect);
myVectECDF <- myEcdf(myVect);
names(myVectECDF) <- names(myVect);

low <- names(myVectECDF[myVectECDF<thresh]);
high <- names(myVectECDF[myVectECDF>(1-thresh)]); 


#Create targets to run limma
tmpLow <- data.frame(low, rep("LOW", length(low)))
colnames(tmpLow) <- c("SAMP", "CLASS");

tmpHigh <- data.frame(high, rep("HIGH", length(high)))
colnames(tmpHigh) <- c("SAMP", "CLASS");

targets <- rbind(tmpLow, tmpHigh);
targets[,2] <- as.factor(targets[,2]);
rownames(targets) <- targets[,1];
targets <- targets[-1];

fTarget <- factor(targets[,"CLASS"]);
design <- model.matrix(~0+fTarget);
colnames(design) <- gsub("fTarget", "", colnames(design));

#Run Limma
trainDatatmp <- myData[,rownames(targets)];
fit <- lmFit(trainDatatmp, design);
myConts <- c("HIGH-LOW");


contrast.matrix <- makeContrasts(contrasts=myConts, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

#Okay now combine all into one result data frame
result <- topTable(fit2, number=60000)[,c("logFC", "P.Value", "adj.P.Val")];
tmpLimmaOut <- result[abs(result[,"logFC"])>logFCThresh,];
tmpLimmaOut <- tmpLimmaOut[tmpLimmaOut[,"adj.P.Val"]<pvalThresh,];

output <- list();
output$all <- result;
output$hits <- tmpLimmaOut;
return(output);
}

#accessory for volcano plot
reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
  log_breaks(base = base), 
  domain = c(1e-100, Inf))
}


#volcano plot, takes in limma analysis
plotVolcanoTrain <- function(result, hitp=.05, hitlfc=1)
{
result[,"HIT"] <- result[,"adj.P.Val"]<hitp&abs(result[,"logFC"])>hitlfc;
p <- ggplot(result, aes(x=logFC, y= adj.P.Val, color=HIT))+geom_point()+ scale_y_continuous(trans=reverselog_trans(10))+theme_bw();
return(p);
}

plotHeatmap <- function(hits, myData)
{
myRows <- rownames(hits);
myHM <- myData[myRows,];
hmPlot <- heatmap.2(myHM);
}


result <- SigLimmaTrain(expData, "121_at", .05)
png("cha.png", width=1440, height=1440,  res=216);
plotVolcanoTrain(result[[1]])
dev.off();

png("cha2.png", width=1440, height=1440,  res=216);
plotHeatmap(result[[2]], as.matrix(expData))
dev.off();






