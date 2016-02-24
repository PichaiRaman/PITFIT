##################################
#Code to generate normal tissue score
#Pichai Raman
#
##################################

#Read in the data
load("../../data/GTEXNormFPKM.RData");

rownames(normDataAnnot) <- normDataAnnot[,1];
rownames(normDataAnnot) <- gsub("-", ".", rownames(normDataAnnot));
normDataAnnot <- normDataAnnot[(colnames(normData)[3:1643]),];


#Split into 2 groups, brain and non-brain

#Get all brain samples
bSamps <- rownames(normDataAnnot[normDataAnnot[,"SMTS"]=="Brain",])
bSamps <- intersect(bSamps, colnames(normData)[3:1643]);
brainData <- normData[,c("Name", "Description", bSamps)];
maxBrainValues <- apply(brainData[3:315], FUN=max, MARGIN=1);

#Get all non-critical tissue normals
ncList <- c("Breast", "Fallopian Tube", "Testis", "Cervix Uteri", "Prostate", "Spleen", "Ovary", "Uterus");
ncSamps <- rownames(normDataAnnot[normDataAnnot[,"SMTS"]%in%ncList,])
ncSamps <- intersect(ncSamps, colnames(normData)[3:1643]);
ncData <- normData[,c("Name", "Description", ncSamps)];
maxNcValues <- apply(ncData[3:66], FUN=max, MARGIN=1);



nonBNCData <- normData[,setdiff(colnames(normData), c(bSamps,ncSamps))];

maxNonBNCValues <- apply(nonBNCData[3:1266], FUN=max, MARGIN=1);

normExpScore <- data.frame(normData[1:2], maxBrainValues, maxNcValues, maxNonBNCValues);
normExpScore <- normExpScore[order(normExpScore[,"maxNonBNCValues"]),];
normExpScore <- normExpScore[!duplicated(normExpScore[,2]),];
normExpScore <- normExpScore[-1];

getScore <- function(x)
{
x <- x[2:4];
score <- 0;

#No expression
if(max(x)<10){ score <- 3}

#Expression only in brain
if(x[1]>10&max(x[2:3])<10) { score <- 2}

#Expression only in non-critical tissues
if(x[2]>10&max(x[c(1,3)])<10) { score <- 1}
return(score);

}

normExpScore[,"Score"] <- apply(normExpScore, FUN=getScore, MARGIN=1);

write.table(normExpScore, "../../data/normalExpressionScore.txt", sep="\t", row.names=F);

