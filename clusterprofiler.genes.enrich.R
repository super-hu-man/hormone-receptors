library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)

print('please paste your genes :')
read_input = function(l=readline(),ins=''){
                      if(l!=''&ins==''){ins=l;read_input(l=readline(),ins)}else if(l!=''&ins!=''){
                      ins=paste(ins,l,sep=';');read_input(l=readline(),ins)
                      }else{print(ins)}}
input = strsplit(read_input(),';',fixed=T)[[1]]
print('please enter the index of the gene start to enrich analysis : ( int type if needed !)')
genes = unique(input[readline():length(input)])

#
print('enrich GO BP')
bp = enrichGO(gene = genes , OrgDb = org.Hs.eg.db , keyType = 'SYMBOL' , readable = TRUE ,
              ont = 'BP' , pAdjustMethod = 'BH' , pvalueCutoff = 0.05 , qvalueCutoff  = 0.05 )
#
print('enrich GO CC')
cc = enrichGO(gene = genes , OrgDb = org.Hs.eg.db , keyType = 'SYMBOL' , readable = TRUE ,
              ont = 'CC' , pAdjustMethod = 'BH' , pvalueCutoff = 0.05 , qvalueCutoff  = 0.05 )
#
print('enrich GO MF')
mf = enrichGO(gene = genes , OrgDb = org.Hs.eg.db , keyType = 'SYMBOL' , readable = TRUE ,
              ont = 'MF' , pAdjustMethod = 'BH' , pvalueCutoff = 0.05 , qvalueCutoff  = 0.05 )
#
print('enrich KEGG')
kegg = enrichKEGG(gene = genes, organism = 'hsa' , pvalueCutoff = 0.05 , qvalueCutoff = 0.05 , 
                  use_internal_data=F )
