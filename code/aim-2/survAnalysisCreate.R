##############################################
#Code to run & create Survival analysis
#Matrix to demonstrate cancer relevance
#
##############################################

#Call Libraries
library("survival");

cleanFormat <- function(annot, exprs)
{
#Get appropriate columns for annot and recode
annot <- annot[,c("PATIENT_ID", "OS_MONTHS", "OS_STATUS")];
colnames(annot) <- c("Patient", "TimeVar", "EventVar");
annot <- na.omit(annot);
annot[,"EventVar"] <- ifelse(annot[,"EventVar"]=="DECEASED", 1, 0);
rownames(annot) <- annot[,1];
annot <- annot[-1];


#Format expression
rownames(exprs) <- exprs[,1] 
exprs <- exprs[-1:-2];
colnames(exprs) <- gsub("\\.","-", colnames(exprs));
colnames(exprs) <- substring(colnames(exprs), 1, 12);

#Now get intersection of id's and conserve order
intSamps <- intersect(colnames(exprs), rownames(annot));
output <- list(exprs[,intSamps], annot[intSamps,]);
return(output);

}


#ovarian
annot_ov <- read.delim("/bigdata/PITFIT_Data/TCGA_Data/ov/tcga/data_bcr_clinical_data.txt");
exprs_ov <- read.delim("/bigdata/PITFIT_Data/TCGA_Data/ov/tcga/data_RNA_Seq_v2_expression_median.txt")
ov <- cleanFormat(annot_ov, exprs_ov);


#prostate
annot_pr <- read.delim("../data/raw/prca/annot.txt");
exprs_pr <- read.delim("../data/raw/prca/exprs.txt")
pr <- list(exprs_pr, annot_pr);

#head and neck
annot_pa <- read.delim("../data/raw/paad/annot.txt");
exprs_pa <- read.delim("../data/raw/paad/exprs.txt")
hn <- list(exprs_pa, annot_pa);



coxReg <- function(genes, myData)
{
    #Get metadata
    tmpMeta <- myData[[2]];
    
    myGene <- myData[[1]][genes,];
    tmpMeta[,"Gene"] <- as.numeric(myGene);
    temp <- table(tmpMeta[,"Gene"])
    out <- c(genes, 1);
    if(names(temp)[temp == max(temp)]!=0)
    {
    tmpMeta[,"Gene"] <- log2(tmpMeta[,"Gene"]+1);
    coxExpAnalysis <- coxph(formula = Surv(TimeVar, EventVar) ~ Gene, data = tmpMeta)
    pVal <- summary(coxExpAnalysis)[7][[1]][5];
    out <- c(genes, pVal);
    }
    out;
    
}



