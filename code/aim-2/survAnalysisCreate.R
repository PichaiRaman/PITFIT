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
annot_pr <- read.delim("/bigdata/PITFIT_Data/TCGA_Data/prad/tcga/data_bcr_clinical_data.txt");
exprs_pr <- read.delim("/bigdata/PITFIT_Data/TCGA_Data/prad/tcga/data_RNA_Seq_v2_expression_median.txt")
pr <- cleanFormat(annot_pr, exprs_pr);

#head and neck
annot_pa <- read.delim("/bigdata/PITFIT_Data/TCGA_Data/paad/tcga/data_bcr_clinical_data.txt");
exprs_pa <- read.delim("/bigdata/PITFIT_Data/TCGA_Data/paad/tcga/data_RNA_Seq_v2_expression_median.txt")
pa <- cleanFormat(annot_pa, exprs_pa);


geneCV1 <- rownames(ov[[1]]);
geneCV2 <- rownames(pr[[1]]);
geneCV3 <- rownames(pa[[1]]);
geneCV <- intersect(geneCV1, geneCV2);
geneCV <- intersect(geneCV, geneCV3);


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
        if(mean(tmpMeta[tmpMeta[,2]==0,3])>mean(tmpMeta[tmpMeta[,2]==1,3]))
        {
            pVal=1;
        }
    out <- c(genes, pVal);
    }
    out;
    
}

coxReg_ov <- sapply(geneCV, FUN=coxReg, ov);
coxReg_ov <- data.frame(t(data.frame(coxReg_ov)));
colnames(coxReg_ov) <- c("Gene", "P.Value");

coxReg_pr <- sapply(geneCV, FUN=coxReg, pr);
coxReg_pr <- data.frame(t(data.frame(coxReg_pr)));
colnames(coxReg_pr) <- c("Gene", "P.Value");

coxReg_pa <- sapply(geneCV, FUN=coxReg, pa);
coxReg_pa <- data.frame(t(data.frame(coxReg_pa)));
colnames(coxReg_pa) <- c("Gene", "P.Value");

coxReg_mat <- data.frame(coxReg_ov[2], coxReg_pr[2], coxReg_pa[2]);
colnames(coxReg_mat) <- c("ov", "pr", "pa");




