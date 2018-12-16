#' Plot a module identified by TERM.
#'
#' @param TEmodule a module identified by TERM.
#' @param layout_style the layout style of the graph for the module, with the default value "layout.fruchterman.reingold".
#' @param gene2symbol a file converting the original gene identifiers to alternative identifiers, with the default value NULL.
#' @return visualization for the network of the input module 
#' where the core and border colors of the nodes denoting differential expression of mRNA and RPF 
#' and the square shape of the nodes indicating the significantly translated genes. 
#' @export
#' @examples
#' plot_TERM()
#' plot_TERM(TEmodule, layout_style="layout.fruchterman.reingold", gene2symbol=NULL  )


plot_TERM <- function(TEmodule,layout_style="layout.fruchterman.reingold", gene2symbol=NULL){

  
  #library(igraph)
  #library(marray)
  #library(corrplot)
  
  map <- function(x, range = c(0,1), from.range=NA) {
    if(any(is.na(from.range))) from.range <- range(x, na.rm=TRUE)
    
    ## check if all values are the same
    if(!diff(from.range)) return(
      matrix(mean(range), ncol=ncol(x), nrow=nrow(x), 
             dimnames = dimnames(x)))
    
    ## map to [0,1]
    x <- (x-from.range[1])
    x <- x/diff(from.range)
    ## handle single values
    if(diff(from.range) == 0) x <- 0 
    
    ## map from [0,1] to [range]
    if (range[1]>range[2]) x <- 1-x
    x <- x*(abs(diff(range))) + min(range)
    
    x[x<min(range) | x>max(range)] <- NA
    
    x
  }

mycircle <- function(coords, v=NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size  <- 1/200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  vertex.frame.color <- params("vertex", "frame.color")
  if (length(vertex.frame.color) != 1 && !is.null(v)) {
    vertex.frame.color <- vertex.frame.color[v]
  }
  vertex.frame.width <- params("vertex", "frame.width")
  if (length(vertex.frame.width) != 1 && !is.null(v)) {
    vertex.frame.width <- vertex.frame.width[v]
  }
  
  mapply(coords[,1], coords[,2], vertex.color, vertex.frame.color,
         vertex.size, vertex.frame.width,
         FUN=function(x, y, bg, fg, size, lwd) {
           symbols(x=x, y=y, bg=bg, fg=fg, lwd=lwd,
                   circles=size, add=TRUE, inches=FALSE)
         })
}


mysquare <- function(coords, v=NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size <-  1/200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  vertex.frame.color <- params("vertex", "frame.color")
  if (length(vertex.frame.color) != 1 && !is.null(v)) {
    vertex.frame.color <- vertex.frame.color[v]
  }
  vertex.frame.width <- params("vertex", "frame.width")
  if (length(vertex.frame.width) != 1 && !is.null(v)) {
    vertex.frame.width <- vertex.frame.width[v]
  }
  
  mapply(coords[,1], coords[,2], vertex.color, vertex.frame.color,
         vertex.size, vertex.frame.width,
         FUN=function(x, y, bg, fg, size, lwd) {
           symbols(x=x, y=y, bg=bg, fg=fg, lwd=lwd,
                   squares=2*size, add=TRUE, inches=FALSE) }
  )
}



ppi=TEmodule$links

mod.graph = graph.data.frame( list( ppi[,1],ppi[,2] ), directed=F, vertices=TEmodule$members  )
mod.graph = set_edge_attr(mod.graph,"weight", value = ppi[,ncol(ppi)])

E(mod.graph)$edgewidth=E(mod.graph)$weight


V(mod.graph)$mtval=TEmodule$gene_DEstat$rna[,1]
V(mod.graph)$rtval=TEmodule$gene_DEstat$ribo[,1]

maxStat = 5

# # add the vm.color, vr.color

mrnaPalette = colorRampPalette(c("green","grey90","red"))
riboPalette = colorRampPalette(c("blue","grey90","orange"))


tmcolor.scheme = mrnaPalette(100)

#tmcolor.scheme[41:60]="#A9A9A9";
#tmcolor.scheme[1:40]=maPalette(low = "yellow", high="lightyellow", k =40);
#tmcolor.scheme[61:100]=maPalette(low = "lightblue", high="blue", k =40);

temp_mStat = data.frame( rownames(TEmodule$gene_DEstat$rna), TEmodule$gene_DEstat$rna,stringsAsFactors = F)
colnames(temp_mStat) = c("GeneID","RNA_DEG_stat","RNA_DEG_pvalue")
temp_mStat$RNA_DEG_stat[temp_mStat$RNA_DEG_stat>=maxStat] = maxStat
temp_mStat$RNA_DEG_stat[temp_mStat$RNA_DEG_stat<=-maxStat] = -maxStat


tmcolor.position = 100 - ceiling(   (maxStat-temp_mStat$RNA_DEG_stat)/((maxStat+maxStat)/100) )
tmcolor.position[tmcolor.position==0] = 1

tmcolor_attr = tmcolor.scheme[tmcolor.position]
#tmcolor_attr[ (temp_mStat$RNA_DEG_pvalue>0.05) | is.na(temp_mStat$RNA_DEG_pvalue) ] = "#A9A9A9"

V(mod.graph)$vmcolor <- tmcolor_attr




trcolor.scheme = riboPalette(100)

#trcolor.scheme[41:60]=  "#A9A9A9";  ### darkgrey;         ## grey is "#BEBEBE";
#trcolor.scheme[1:40]=maPalette(low = "green", high="lightgreen", k =40);
#trcolor.scheme[61:100]=maPalette(low = "#DC6868", high="red", k =40);

temp_rStat = data.frame( rownames(TEmodule$gene_DEstat$ribo), TEmodule$gene_DEstat$ribo, stringsAsFactors = F )
colnames(temp_rStat) = c("GeneID", "Ribo_DEG_stat", "Ribo_DEG_pvalue")
temp_rStat$Ribo_DEG_stat[temp_rStat$Ribo_DEG_stat>=maxStat] = maxStat
temp_rStat$Ribo_DEG_stat[temp_rStat$Ribo_DEG_stat<=-maxStat] = -maxStat

trcolor.position = 100 - ceiling(   (maxStat-temp_rStat$Ribo_DEG_stat)/((maxStat+maxStat)/100) )
trcolor.position[trcolor.position==0] = 1

trcolor_attr = trcolor.scheme[trcolor.position]
#trcolor_attr[ (temp_rStat$Ribo_DEG_pvalue>0.05) | is.na(temp_rStat$Ribo_DEG_pvalue) ] = "#A9A9A9"

V(mod.graph)$vrcolor <- trcolor_attr





#################################################################################
#add the mod label value

if(is.null(gene2symbol)){
  label.v =  TEmodule$members
}else{
  label.v =  gene2symbol[match(TEmodule$members,gene2symbol[,1]) ,2]
  label.v[ is.na(label.v) ] = TEmodule$members[ is.na(label.v) ]
}

 


#create subgraph label.cex value
V(mod.graph)$label.cex=rep(0.5,length(as.vector(V(mod.graph)))); #all the cex first set as 0.7
# V(mod.graph)$label.cex[which(as.vector(V(mod.graph)$name)==as.vector(mod[1,1]))]=0.8 #only the firt mod name was set as 1


#################################################################################
# generate the plot 
# when you want to plot the vertex shape, and its frame width first you should load the api script api bellow
add.vertex.shape("fcircle", clip=igraph.shape.noclip,plot=mycircle, parameters=list(vertex.frame.color=1))
add.vertex.shape("fsquare", clip=igraph.shape.noclip,plot=mysquare, parameters=list(vertex.frame.color=1))

v_shape_TE = rep("fcircle", length(TEmodule$members) )
v_shape_TE[TEmodule$gene_TEstat<=0.05] = "fsquare"
V(mod.graph)$vertex.shape = v_shape_TE



v_size_base = 10
if( length(V(mod.graph))<10 ){
  v_size_scale = 1
}else if( length(V(mod.graph))>1000 ){
  v_size_scale = 10
}else{
  v_size_scale = sqrt( floor(length(V(mod.graph))/10 ) )
}
v_size = v_size_base/v_size_scale


V(mod.graph)$vertex.size = rep( v_size, length( V(mod.graph)$name ) )
#V(mod.graph)$vertex.size[ V(mod.graph)$vertex.shape=="fcircle"]=v_size
#V(mod.graph)$vertex.size[ V(mod.graph)$vertex.shape=="fsquare"]=2*v_size

V(mod.graph)$vertex.frame.width = rep( 0.5*v_size, length( V(mod.graph)$name ) )

V(mod.graph)$vertex.size[ V(mod.graph)$name==TEmodule$seed  ] = 1.5*V(mod.graph)$vertex.size[ V(mod.graph)$name==TEmodule$seed  ]
V(mod.graph)$vertex.frame.width[ V(mod.graph)$name==TEmodule$seed  ] = 1.5*V(mod.graph)$vertex.frame.width[ V(mod.graph)$name==TEmodule$seed  ]



plot(mod.graph,
     layout=get(layout_style),
     vertex.shape=V(mod.graph)$vertex.shape ,
     vertex.size=V(mod.graph)$vertex.size,
     vertex.color=V(mod.graph)$vmcolor,
     vertex.frame.width=V(mod.graph)$vertex.frame.width,
     vertex.frame.color=V(mod.graph)$vrcolor,
     vertex.label=label.v,
     vertex.label.dist=0,
     vertex.label.cex=V(mod.graph)$label.cex,
     vertex.label.font=3,
     edge.color="grey90",
     edge.width=map(E(mod.graph)$weight,c(0.5,v_size) )
)

colorlegend(tmcolor.scheme,seq(-maxStat,maxStat,1),ratio.colbar=0.3 ,xlim=c(-1.55,-1.4),ylim=c(0.5,1),align="r",cex=0.5)
colorlegend(trcolor.scheme,seq(-maxStat,maxStat,1),ratio.colbar=0.3 ,xlim=c(-1.55,-1.4),ylim=c(-0.5,0),align="r",cex=0.5)

text(-1.47, 0.43, c("stat(mRNA)\nCore"),cex=0.5)
text(-1.47, -0.57,c("stat(Ribo)\nBorder"),cex=0.5)

}