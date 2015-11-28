##############################
#Code to support the Aim 2 tool
#which is for target prioritization
#
#Pichai Raman
#11/28/2015
##############################

#call libraries

##############################
#Main Function to prioritize and rank
#Targets
#
##############################
pitfitAnalyzeAim2 <- function(myCancer, geneList)
{

output <- data.frame(geneList, myCancer);

return(output);
}

