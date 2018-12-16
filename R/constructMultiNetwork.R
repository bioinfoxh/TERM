



get_wEdge <- function( a,b  ){
  a=abs(a)
  b=abs(b)
  temp1 = 2*min(a,b)/(a+b) 
  temp2 = max(a,b)
  wEdge = sqrt( temp1*temp2 )
  return(wEdge)
}








constructMultiNetwork <- function( DEstat, ppi){
  
  gene_netwk = rownames( DEstat[[1]] )
  
  ppi_g = graph.data.frame( list( ppi[,1],ppi[,2] ), directed=F, vertices=gene_netwk   )
  
  multi_netwk_g = list()
  for (i in 1:length(DEstat)){
    
    temp_stat = DEstat[[i]][,1]
    names(temp_stat) = rownames(DEstat[[i]])
    
    temp_p = DEstat[[i]][,2]
    names(temp_p) = rownames(DEstat[[i]])
    
    temp_ppi_w = cbind( temp_stat[ match(ppi[,1],names(temp_stat) )], temp_stat[ match(ppi[,2],names(temp_stat) ) ] )
    temp_wEdge = apply( temp_ppi_w, 1, function(x) get_wEdge(x[1],x[2]) )
    
    temp_ppi_g = ppi_g
    temp_ppi_g = set_edge_attr(temp_ppi_g,"weight", value = temp_wEdge )
    temp_ppi_g = set_vertex_attr(temp_ppi_g,"statistics", value = temp_stat )
    temp_ppi_g = set_vertex_attr(temp_ppi_g,"pvalue", value = temp_p )
    
    multi_netwk_g[[i]] = temp_ppi_g
    
  }
  names(multi_netwk_g) = names(DEstat)
  
  return(multi_netwk_g)
  
}


