#' Transforms Ct values into fold over mean
#' This function takes as input a data frame where the first n columns contain
#' categorical variables describing each sample (treatment, experiment, cell,...)
#' and the remaining columns containing the Ct values for different genes (one gene
#'  per column). The function calculates the fold expression of each sample over
#'  the mean (taken over the samples in the column) as :  Eff^(Ct_mean-Ct); where
#'  Eff is the PCR efficiency for that gene (default 2)
#'  
#'  @param df data matrix containing sample info in the first n columns and Ct
#'   for measured genes (one gene per column)
#'  @param nInfCol number of columns containing sample information. This columns
#'  must precede those with Ct. If data frame contains no sample information,
#'  nInfCol must be 0 (default)
#'  @param Eff efficiency of the PCR. Note that it is applied to all columns.
#'  Default value 2.
#'  @return a data frame containing the first n columns with sample information
#'   followed by fold values for each gene
#'   @export
FOM<-function(df,nInfCol=0,Eff=2){
  #transforms Ct into folds over mean Ct
  GED<-as.matrix(df[,(nInfCol+1):dim(df)[2]])
  Ct_centered<-scale(GED,center=colMeans(GED,na.rm=T),scale=F)
  Fold<-Eff**(-Ct_centered)
  if (nInfCol==0){
    return(as.data.frame(Fold))
  }else{
    return(cbind(df[,1:nInfCol],as.data.frame(Fold))) 
  }
}
