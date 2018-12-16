

get_permutationFDR <- function( obsScore, permScore ){
  perm_mean = apply( permScore, 1, mean )
  perm_sd = apply( permScore, 1, sd )
  
  obsScoreNorm = (obsScore-perm_mean)/perm_sd
  permScoreNorm = apply( permScore, 2, function(x) (x-perm_mean)/perm_sd  )
  
  fdr=c()
  for(i in 1:length(obsScore)){
    ratio_perm = sum(rowSums( permScoreNorm<obsScoreNorm[i]) ) / (ncol(permScoreNorm)*nrow(permScoreNorm))
    ratio_obs = sum( obsScoreNorm<=obsScoreNorm[i] )/length(obsScoreNorm)
    tempFDR = ratio_perm/ratio_obs
    if(tempFDR>1){
      tempFDR=1
    }
    fdr = c(fdr, tempFDR)
  }
  return(fdr)
}





getModStat <- function(n_modules,multi_netwk,DEstat, ppi,permTimes=100){
  
  
  
  module_entropy = entropy(n_modules,multi_netwk)
  
  #permTimes=100
  entropy_mrna_perm = c()
  entropy_rpf_perm = c()
  entropy_all_perm = c()
  
  for( i in 1:permTimes){
    print( paste("permutation times ",i," / ",permTimes,sep="")  )
    #if(i%%10==0){print( paste("permutation times ",i," / ",permTimes,sep="")  )}
    
    perm_DEstat = list()
    for(j in 1:length(DEstat)){
      perm_DEstat[[j]] = DEstat[[j]][  sample( 1:nrow(DEstat[[j]]) ) , ]
      rownames(perm_DEstat[[j]]) = rownames(DEstat[[j]])
    }
    
    perm_multi_netwk_g = constructMultiNetwork(perm_DEstat, ppi)
    
    gene_netwk =  rownames( DEstat[[1]] )
    
    perm_multi_netwk = array(0, dim=c(length(gene_netwk),length(gene_netwk),length(perm_multi_netwk_g) )  )
    for(k in 1:length(perm_multi_netwk_g)){
      temp_adj =  as.matrix( as_adjacency_matrix(perm_multi_netwk_g[[k]],type="both",attr="weight",names=TRUE) )
      perm_multi_netwk[,,k] = temp_adj
    }
    rm(temp_adj)
    
    temp_entropy = entropy(n_modules, perm_multi_netwk)
    
    entropy_mrna_perm = cbind(entropy_mrna_perm, temp_entropy[,1])
    entropy_rpf_perm = cbind(entropy_rpf_perm, temp_entropy[,2])
    entropy_all_perm = cbind(entropy_all_perm, temp_entropy[,3])
    
    rm(perm_multi_netwk)
    
  }
  
  
  module_p_mrna=apply(cbind(module_entropy[,1],entropy_mrna_perm),1,function(x){sum( x[2:length(x)]<x[1] )/(length(x)-1)} )
  module_p_rpf=apply(cbind(module_entropy[,2],entropy_rpf_perm),1,function(x){sum( x[2:length(x)]<x[1] )/(length(x)-1)} )
  module_p_all=apply(cbind(module_entropy[,3],entropy_all_perm),1,function(x){sum( x[2:length(x)]<x[1] )/(length(x)-1)} )
  
  module_p = data.frame( mrna=module_p_mrna, ribo=module_p_rpf,  average=module_p_all)
  
  
  module_fdr_mrna = get_permutationFDR( module_entropy[,1], entropy_mrna_perm  )
  module_fdr_rpf = get_permutationFDR( module_entropy[,2], entropy_rpf_perm  )
  module_fdr_all = get_permutationFDR( module_entropy[,3], entropy_all_perm  )
  
  module_fdr = data.frame( mrna=module_fdr_mrna, ribo=module_fdr_rpf, average=module_fdr_all)
  
  modularity_significance = list( entropy = module_entropy,
                                  pvalue = module_p,
                                  FDR = module_fdr)
  return(modularity_significance)
}


