library(pheatmap)
setwd('D:/GitHub/Hormone-Receptors/')

infor = read.csv('D:/GitHub/Files.useful/GPL13534-11288.txt',header=T,sep='\t',row.names=1,stringsAsFactors=F,skip=37)
rq = 'n' 
files = list.files('D:/GitHub/Hormone-Receptors/result/Meth.correlation','cor.mRNA_Meth.txt',full.names = T)
tmp = data.frame(setNames(replicate(1,numeric(0), simplify = F), c("Meth_IDs")))
cor = data.frame(setNames(replicate(1,numeric(0), simplify = F), c("Meth_IDs")))
for(f in files){
  name = strsplit(strsplit(strsplit(f,'/',fixed=T)[[1]][6],'.',fixed=T)[[1]][1],'_',fixed = T)[[1]][1]
  print(name) ;
  df = read.csv(f,header = T,sep='\t',stringsAsFactors = F)
   cancer_cor = df[,c('Meth_IDs','logE_M_pearson_cor')] ; colnames(cancer_cor)[2] = name
   cor = merge(cor,cancer_cor,by=1,all=T)
  #print(date())
  df1 = df[df$Gene_name!='',]
  if (rq == 'y'){
    df1$logE_M_pearson_q = p.adjust( df1$logE_M_pearson_p ,'BH')
    df1$E_M_pearson_q    = p.adjust( df1$E_M_pearson_p ,'BH')
    df1$logE_M_spearman_q = p.adjust( df1$logE_M_spearman_p ,'BH')
    df1$E_M_spearman_q    = p.adjust( df1$E_M_spearman_p ,'BH')  }
  df2 = df1[df1$logE_M_spearman_q<=0.05,] ; #df2 = df2[abs(df2$logE_M_pearson_cor)>=0.4,]
  print(nrow(df2))
  x = df2[,c("Meth_IDs","Gene_name","Group","logE_M_spearman_cor")]
  colnames(x)[4] = name
  if(nrow(df2)<1){print(name);print('no pathway !');tmp[,name]=NA;next}
  tmp = merge(tmp,x,all=T)
}
tmp = data.frame(tmp[,1:3],Number=apply(tmp,1,function(x){sum(!is.na(x[4:36]))}),tmp[,4:36])
tmp = tmp[order(tmp$Number,decreasing=T),] ; rownames(tmp) = tmp$Meth_IDs
Probe_number = as.character(c(rep(NA,4),colSums(!is.na(tmp[,5:37]))))
tmp.out = rbind(Probe_number=Probe_number, tmp)
write.table(tmp.out,'./result/Meth.correlation/cor.cancers.summary.txt',col.names=T,row.names=F,quote=F,sep='\t')

tmp=read.csv('./result/Meth.correlation/cor.cancers.summary.txt',header=T,sep='\t',stringsAsFactors=F)
tmp = tmp[-1,] ; rownames(tmp) = tmp$Meth_IDs
tmp[is.na(tmp)]=0 ; tmp[,5:37][tmp[,5:37]>0]=1 ; tmp[,5:37][tmp[,5:37]<0]=-1
for(i in 1){plot(hclust(dist(t(tmp[tmp$Number>i,5:37]))),xlab=paste('as least in ',i+1,' cancer',sep='')) }
#pheatmap(tmp[tmp$Number>5,5:37],show_rownames = F,cluster_rows=F)

h=tmp[tmp$Number>1,5:37] 
# h=tmp[1:5000,5:37]
h[is.na(h)]=0 ;
res1 = pheatmap(h,cutree_rows = 2,cutree_cols = 2,show_rownames = F,cluster_rows=F)
tiff('./result/Meth.correlation/cor.sig.tiff',width=16,height=24,units='cm',res=600,compression='lzw')
res1
dev.off()
h[h>0] =1 ; h[h<0] = -1
res2 = pheatmap(h,cutree_rows = 2,cutree_cols = 2,show_rownames = F)
tiff('./result/Meth.correlation/cor.binomal.tiff',width=16,height=24,units='cm',res=600,compression='lzw')
res2
dev.off()

extract = function(x,infor){
  hc = as.data.frame(cutree(x$tree_row,k=2)) ; colnames(hc) = 'Cluster'
  cor.pos = rownames(hc)[hc$Cluster==2]
  cor.neg = rownames(hc)[hc$Cluster==1]
  #
  prb.pos = infor[cor.pos,c('Name','UCSC_RefGene_Name','UCSC_RefGene_Accession','UCSC_RefGene_Group','UCSC_CpG_Islands_Name')]
  prb.pos1 = prb.pos[rep(1:nrow(prb.pos),sapply(prb.pos$UCSC_RefGene_Name,function(x){a=length(strsplit(x,';',fixed=T)[[1]]);return(a)})),]
  prb.pos1$UCSC_RefGene_Name = unlist(sapply(prb.pos$UCSC_RefGene_Name,function(x){a=unlist(strsplit(x,';',fixed=T));return(a)}))
  prb.pos1$UCSC_RefGene_Group = unlist(sapply(prb.pos$UCSC_RefGene_Group,function(x){a=unlist(strsplit(x,';',fixed=T));return(a)}))
  prb.pos2 = unique(prb.pos1) ; rm('prb.pos1')
  #
  prb.neg = infor[cor.neg,c('Name','UCSC_RefGene_Name','UCSC_RefGene_Accession','UCSC_RefGene_Group','UCSC_CpG_Islands_Name')]
  prb.neg1 = prb.neg[rep(1:nrow(prb.neg),sapply(prb.neg$UCSC_RefGene_Name,function(x){a=length(strsplit(x,';',fixed=T)[[1]]);return(a)})),]
  prb.neg1$UCSC_RefGene_Name = unlist(sapply(prb.neg$UCSC_RefGene_Name,function(x){a=unlist(strsplit(x,';',fixed=T));return(a)}))
  prb.neg1$UCSC_RefGene_Group = unlist(sapply(prb.neg$UCSC_RefGene_Group,function(x){a=unlist(strsplit(x,';',fixed=T));return(a)}))
  prb.neg2 = unique(prb.neg1) ; rm('prb.neg1')
  res = list(prb.pos2, prb.neg2)
  #
  return(res)
}

x1 = extract(res1,infor)
write.table(x1[[1]],'./result/Meth.correlation/files1.pos.genes.txt',col.names=T,sep='\t',quote=F,row.names=F)
write.table(x1[[2]],'./result/Meth.correlation/files1.neg.genes.txt',col.names=T,sep='\t',quote=F,row.names=F)
x2 = extract(res2,infor)
write.table(x2[[1]],'./result/Meth.correlation/files2.pos.genes.txt',col.names=T,sep='\t',quote=F,row.names=F)
write.table(x2[[2]],'./result/Meth.correlation/files2.neg.genes.txt',col.names=T,sep='\t',quote=F,row.names=F)

onecluster = h[,names(cutree(res2$tree_col,k=2)[cutree(res2$tree_col,k=2)==2])]
onecluster = onecluster[rowSums(abs(onecluster)==1)>0,]
pheatmap(onecluster,cutree_rows = 2)