#' @title Benjamini-Hochberg Correction in LAWS-HiC
#' @description Applies the BH procedure for input p-values
#' @param pv A numeric vector of p-values
#' @param q The false discovery rate (FDR) threshold.
#' @return A list with the number of rejected hypotheses, threshold, and decision rule.
#' @examples
#' p_vector <- c(0.01, 0.02, 0.03, 0.04, 0.00001, 0.002, 0.003, 0.6)
#' bh.func(p_vector, 0.005)
#' @export 
bh.func<-function(pv, q)
{ 
  # the input 
  # pv: the p-values
  # q: the FDR level
  # the output 
  # nr: the number of hypothesis to be rejected
  # th: the p-value threshold
  # de: the decision rule
  
  m=length(pv)
  st.pv<-sort(pv)   
  pvi<-st.pv/1:m
  de<-rep(0, m)
  if (sum(pvi<=q/m)==0)
  {
    k<-0
    pk<-1
  }
  else
  {
    k<-max(which(pvi<=(q/m)))
    pk<-st.pv[k]
    de[which(pv<=pk)]<-1
  }
  y<-list(nr=k, th=pk, de=de)
  return (y)
}