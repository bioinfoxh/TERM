#' Detecting TE-regulation modules from the constructed multi-layer differential expression network
#'
#' @param ppi a matrix/dataframe containing gene-gene interactions in two columns
#' @param DEstat a list including two dataframes of statistics and p-values of differential expression for mRNAs and RPFs respectively
#' @param TEstat a vector of differential translation significance of genes
#' @param num_seed the number of seed genes used for searching modules from multi-layer networks, with default value 100
#' @param minModSize the number of minimum genes within a module, with default value 10
#' @param permTimes permutation times for calculating modulatrity significance of modules, with default value 100
#' @param modularity_p threshold for modulatrity significance, with default value 0.05
#' @param maxModOvlp threshold for meetmin index defining the overlapping between modules, with default value 0.5
#' @return a list of results including a dataframe of the multi-layer network with edge weights in each network as well as the average weight and a list of identified TE-regulated modules   
#' @export
#' @examples
#' detectTERM()
#' detectTERM(ppi, DEstat, TEstat,num_seed=100,minModSize=10,permTimes=100,modularity_p=0.05,maxModOvlp=0.5)

detectTERM <- function(ppi, DEstat, TEstat,num_seed=100,minModSize=10,permTimes=100,modularity_p=0.05,maxModOvlp=0.5){
  
  
  gene_netwk =  rownames( DEstat[[1]] )

multi_netwk_g = constructMultiNetwork( DEstat, ppi  )

multi_netwk = array(0, dim=c(length(gene_netwk),length(gene_netwk),length(multi_netwk_g) )  )
for(i in 1:length(multi_netwk_g)){
  temp_adj =  as.matrix( as_adjacency_matrix(multi_netwk_g[[i]],type="both",attr="weight",names=TRUE) )
  multi_netwk[,,i] = temp_adj
}

rm(temp_adj)



##########

#num_seed = 100

generank = TEstat
seedgene = order(generank)
seeds = sort(seedgene[1:num_seed]) 

n_modules = nModule(multi_netwk,seeds)


temp_len = sapply( n_modules, function(x) length(x$members) )
n_modules = n_modules[ temp_len>=minModSize ]
seeds = seeds[ temp_len>=minModSize ]


######### calculate the p-values for each module

print("evaluating modularity significance for modules:")

#module_significance = getModStat(n_modules,multi_netwk,multi_network_g,DEstat, ppi,permTimes)

module_significance = getModStat(n_modules,multi_netwk,DEstat, ppi,permTimes)

module_entropy = module_significance$entropy
module_p = module_significance$pvalue
module_fdr = module_significance$FDR



###################### process modules

sigModuleIdx = which(rowSums(module_p<=modularity_p)>0  )
sigModule = n_modules[sigModuleIdx]
names(sigModule) = sigModuleIdx



sigModule_delOvlp = deleteModOvlp(sigModule,multi_netwk,maxModOvlp)

TEmod_idx = as.numeric( names(sigModule_delOvlp) )
TEmod_idx = TEmod_idx[order(TEmod_idx)]

TEmod = n_modules[TEmod_idx]
names(TEmod) = TEmod_idx

TEmod_seeds = seeds[TEmod_idx]

TEmod_p = module_p[TEmod_idx,]
TEmod_fdr = module_fdr[TEmod_idx,]

TEmod_entropy = data.frame(module_entropy[ TEmod_idx ,],stringsAsFactors = F)
rownames(TEmod_entropy) = TEmod_idx
colnames(TEmod_entropy) = colnames(TEmod_p)




temp_w = sapply( multi_netwk_g, function(x) E(x)$weight )
ppi_w = data.frame(ppi, temp_w, rowMeans(temp_w), stringsAsFactors = F)
colnames(ppi_w) = c( "Gene1","Gene2","weight_rna","weight_ribo","weight_average"   )

result_TEmodule = list()
for( i in 1:length(TEmod)){
  temp_list=list()
  
  temp_member_idx = TEmod[[i]]$member
  temp_member_name = gene_netwk[ TEmod[[i]]$members ]
  
  temp_list$seed = gene_netwk[ TEmod_seeds[i] ]
  
  temp_list$members = temp_member_name
  temp_list$entropy = TEmod_entropy[i,]
  temp_list$entropy_pvalue = TEmod_p[i,]
  temp_list$entropy_FDR = TEmod_fdr[i,]
  
  temp_list$gene_DEstat = lapply(DEstat, function(x) x[ temp_member_idx ,]   )
  temp_list$gene_TEstat = TEstat[ temp_member_idx ]
  
  temp_list$links = ppi_w[ (ppi_w[,1]%in%temp_member_name)&(ppi_w[,2]%in%temp_member_name) ,]
  
  
  result_TEmodule[[i]] = temp_list
}
names(result_TEmodule) = gene_netwk[ TEmod_seeds ]

return( list( multiNetwk = ppi_w,TERM = result_TEmodule  )  )


}

