##############################
#Code to support the Aim 1 tool
#which is for target discovery.
#
#Pichai Raman
#5/31/2015
##############################

#Read in data
expData <- read.delim("../data/ovarianData.txt");
annotData <- read.delim("../data/ovarianDataAnnot.txt");



SigLimmaTrain <- function(myData, myGene, thresh=.20, pvalThresh=.25, logFCThresh=1)
{
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
tmpLimmaOut <- topTable(fit2, number=60000)[,c("logFC", "P.Value", "adj.P.Val")];
tmpLimmaOut <- tmpLimmaOut[abs(tmpLimmaOut[,"logFC"])>logFCThresh,];
tmpLimmaOut <- tmpLimmaOut[tmpLimmaOut[,"adj.P.Val"]<pvalThresh,];


}


result <- SigLimmaTrain(expData, "121_at", .20)








