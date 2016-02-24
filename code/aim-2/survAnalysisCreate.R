##############################################
#Code to run & create Survival analysis
#Matrix to demonstrate cancer relevance
#
##############################################

#Call Libraries
library("survival");

cleanFormat <- function(annot, expr)
{
#Get appropriate columns for annot and recode
annot <- annot_ov[,c("PATIENT_ID", "OS_MONTHS", "OS_STATUS")];
annot <- na.omit(annot);



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
    coxExpAnalysis <- coxph(formula = Surv(TimeVar, eventVar) ~ Gene, data = tmpMeta)
    pVal <- summary(coxExpAnalysis)[7][[1]][5];
    out <- c(genes, pVal);
    }
    out;
    
}



