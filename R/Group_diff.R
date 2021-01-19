#' Calculates differences in fold induction for two groups of samples
#' This function takes as input a data frame where the first n columns contain
#' categorical variables describing each sample (treatment, experiment, cell,...)
#' and the remaining columns containing expression values (as fold over mean) 
#' for different genes (one gene per column). The function calculates the 4
#' difference in fold expression for each gene between two groups of samples 
#' defined by the values in one of the columns describing the samples
#'  @param df data matrix containing sample info in the first n columns and Ct
#'   for measured genes (one gene per column)
#'  @param nInfCol number of columns containing sample information. This columns
#'  must precede those with Ct. The data frame must contain at least one column 
#'  containing sample information to group data.
#'  @param colGroup number of column containing the variable to group data (only)
#'   two categories are allowed.
#'  @return a data frame containing the first n columns with sample information
#'   followed by the differences in fold values for each gene
#'   @export

Group_diff<-function(df,nInfCol,colGroup){
  tmpTxt<-factor(df[,colGroup])
  tmpList<-split(df,tmpTxt)
  GED1<-as.matrix(tmpList[[1]][,(nInfCol+1):dim(tmpList[[1]])[2]])
  GED2<-as.matrix(tmpList[[2]][,(nInfCol+1):dim(tmpList[[2]])[2]])
  Diff<-as.data.frame(GED1-GED2)
  infocol<-setdiff(1:nInfCol,colGroup)
  return(cbind(tmpList[[1]][,infocol],Diff))
}

