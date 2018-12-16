#' Calculate the differential expression of mRNAs and RPFs, and the differential translation of genes.
#'
#' @param raw_rna raw RNA-seq read counts
#' @param raw_ribo raw Ribo-seq read counts
#' @param rna_label the experiment conditions corresponding to the samples in raw_rna
#' @param ribo_label the experiment conditions corresponding to the samples in raw_ribo
#' @param baseLevel the condition used as the base level for comparisons of differential expression and differential translation 
#' @param minCounts the number of minimum read count of a gene in a sample, with default value 10 
#' @param minCountsProportion the proportion of samples where a gene exhibits read counts larger than the minCounts in a same condition, with default value 1 
#' @return a dataframe consists of DESeq2 results of mRNA and RPF differential expression and Riborex results of differential translation 
#' @export
#' @examples
#' calculateDTE()
#' calculateDTE(raw_rna, raw_ribo, rna_label, ribo_label, baseLevel, minCounts=10, minCountsProportion=1 )



calculateDTE <- function(raw_rna, raw_ribo, rna_label, ribo_label, baseLevel, minCounts=10, minCountsProportion=1 ){


sampleClass_rna = unique(rna_label)
sampleClass_ribo = unique(ribo_label)

keepGene_rna = rownames(raw_rna)[ rowSums(raw_rna[,rna_label==sampleClass_rna[1] ]>=minCounts) >= (sum(rna_label==sampleClass_rna[1])*minCountsProportion) | 
                                    rowSums(raw_rna[,rna_label==sampleClass_rna[2] ]>=minCounts) >= (sum(rna_label==sampleClass_rna[2])*minCountsProportion)  ]

keepGene_ribo = rownames(raw_ribo)[ rowSums(raw_ribo[,ribo_label==sampleClass_ribo[1] ]>=minCounts) >= (sum(ribo_label==sampleClass_ribo[1])*minCountsProportion) | 
                                      rowSums(raw_ribo[,ribo_label==sampleClass_ribo[2] ]>=minCounts) >= (sum(ribo_label==sampleClass_ribo[2])*minCountsProportion)  ]


keepGene = intersect(keepGene_rna,keepGene_ribo)

rna = raw_rna[ match( keepGene, rownames(raw_rna)  )  ,  ]
ribo = raw_ribo[ match(keepGene, rownames(raw_ribo) )  ,  ]

rna[rna==0] = 1
ribo[ribo==0] = 1



print("Estimating differential expression of RNA:")

rna_label = data.frame(condition=rna_label)
rownames(rna_label) = colnames(rna)

dds_rna = DESeqDataSetFromMatrix(countData = rna, colData = rna_label, design = ~ condition)
dds_rna$condition = relevel(dds_rna$condition, ref=baseLevel)
dds_rna = DESeq(dds_rna)
rna_deseq2 = results(dds_rna)


print("Estimating differential expression of Ribo:")

ribo_label = data.frame(condition=ribo_label)
rownames(ribo_label) = colnames(ribo)

dds_ribo = DESeqDataSetFromMatrix(countData = ribo, colData = ribo_label, design = ~ condition)
dds_ribo$condition = relevel(dds_ribo$condition, ref=baseLevel)
dds_ribo = DESeq(dds_ribo)
ribo_deseq2 = results(dds_ribo)

print("Estimating differential translation efficiency:")

result_riborex = riborex(rna, ribo, rna_label, ribo_label)

res_rna = data.frame(rna_deseq2)
res_ribo = data.frame(ribo_deseq2)
res_TE = data.frame(result_riborex)

result_TEgene = cbind( res_rna[match(keepGene,rownames(res_rna)),], 
                       res_ribo[match(keepGene,rownames(res_ribo)),],
                       res_TE[match(keepGene,rownames(res_TE)),]  )

colnames(result_TEgene) = c( paste(colnames(res_rna),"_RNA",sep=""),
                             paste(colnames(res_ribo),"_Ribo",sep=""),
                             paste(colnames(res_TE),"_TE",sep="")  )


return(result_TEgene)

}
