##############################################
#Code to run & create Survival analysis
#Matrix to demonstrate cancer relevance
#
##############################################



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




#ovarian
annot_ov <- read.delim("../data/raw/ovca/annot.txt");
exprs_ov <- read.delim("../data/raw/ovca/exprs.txt")
ov <- list(exprs_ov, annot_ov);

#prostate
annot_pr <- read.delim("../data/raw/prca/annot.txt");
exprs_pr <- read.delim("../data/raw/prca/exprs.txt")
pr <- list(exprs_pr, annot_pr);

#head and neck
annot_pa <- read.delim("../data/raw/hnca/annot.txt");
exprs_pa <- read.delim("../data/raw/hnca/exprs.txt")
hn <- list(exprs_hn, annot_hn);
