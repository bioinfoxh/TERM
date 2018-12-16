#' The main function performs the whole pipeline of TERM to identify TE-regulated modules.
#'
#' @param raw_rna raw RNA-seq read counts matrix
#' @param raw_ribo raw Ribo-seq read counts matrix
#' @param rna_label the experiment conditions corresponding to the samples in raw_rna
#' @param ribo_label the experiment conditions corresponding to the samples in raw_ribo
#' @param baseLevel the condition used as the base level for comparisons of differential expression and differential translation 
#' @param raw_ppi a matrix/dataframe containing gene-gene interactions in two columns
#' @param minCounts the number of minimum read count of a gene in a sample, with default value 10
#' @param minCountsProportion the proportion of samples where a gene exhibits read counts larger than the minCounts in a same condition, with default value 1 
#' @param num_seed the number of seed genes used for searching modules from multi-layer networks, with default value 100
#' @param minModSize the number of minimum genes within a module, with default value 10
#' @param permTimes permutation times for calculating modulatrity significance of modules, with default value 100
#' @param modularity_p threshold for modulatrity significance, with default value 0.05
#' @param maxModOvlp threshold for meetmin index defining the overlapping between modules, with default value 0.5
#' @return a list of results including a dataframe of mRNA and RPF differential expression and differential translation of genes, 
#' a dataframe of the multi-layer network with edge weights in each network as well as the average weight, 
#' a list of identified TE-regulated modules   
#' @export
#' @examples
#' term()
#' term(raw_rna, raw_ribo, rna_label, ribo_label, baseLevel, raw_ppi, minCounts=10, minCountsProportion=1,num_seed=100,minModSize=10,permTimes=100, modularity_p=0.05,maxModOvlp=0.5 )





term <- function(raw_rna, raw_ribo, rna_label, ribo_label, baseLevel, raw_ppi, minCounts=10, minCountsProportion=1,num_seed=100,minModSize=10,permTimes=100, modularity_p=0.05,maxModOvlp=0.5 ){
  
  
  result_TEgene = calculateDTE(raw_rna, raw_ribo, rna_label, ribo_label, baseLevel, minCounts, minCountsProportion )
  
  
  
  ##################################
  ########## network analysis
  ##################################
  
  
  
  DE = result_TEgene[ (!is.na(result_TEgene$stat_RNA) ) & (!is.na(result_TEgene$stat_Ribo)) ,]
  gene_nonNA = rownames(DE)
  
  
  temp1 = raw_ppi[ (!is.na(raw_ppi[,1])) & (!is.na(raw_ppi[,2])) & (raw_ppi[,1]!=raw_ppi[,2]) ,]
  temp2 = t( apply( temp1, 1, function(x) sort(x) ) )
  temp_ppi = unique(temp2)
  ppi = temp_ppi[ (temp_ppi[,1]%in%gene_nonNA) & (temp_ppi[,2]%in%gene_nonNA), ]
  
  rm(temp1,temp2,temp_ppi)
  
  
  gene_netwk = union(ppi[,1],ppi[,2])
  
  
  DE = DE[ match(gene_netwk, rownames(DE)) ,]
  
  temp1 = cbind(DE$stat_RNA,DE$pvalue_RNA)
  rownames(temp1) = gene_netwk
  temp2 = cbind(DE$stat_Ribo,DE$pvalue_Ribo)
  rownames(temp2) = gene_netwk
  
  DEstat = list( rna=temp1, ribo=temp2  )
  
  TEstat = DE$pvalue_TE
  names(TEstat)=gene_netwk
  
  
  result_TERM = detectTERM(ppi, DEstat, TEstat, num_seed, minModSize, permTimes, modularity_p, maxModOvlp)
  
  
  return( list( DTE=result_TEgene, multiNetwk=result_TERM$multiNetwk, TERM=result_TERM$TERM  ) )
  
  
  
}