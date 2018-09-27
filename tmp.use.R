setwd('D:/GitHub/download.pan-cancer-altas/gene')

files=list.files('./','cor.mRNA_Meth.gene.txt',full.names = T)
dat = read.csv(files[1],header = T,sep='\t',row.names = 1,stringsAsFactors = F)
x =data.frame()
for(f in files){
  name = strsplit(strsplit(f,'/',fixed=T)[[1]][2],'.',fixed=T)[[1]][1]
  dat = read.csv(f,header = T,sep='\t',row.names = 1,stringsAsFactors = F)
  line = dat['AR',] ; rownames(line)=name
  x= rbind(x,line)
}
dat$logE_M_pearson_cor=as.numeric(dat$logE_M_pearson_cor)
dat$logE_M_pearson_q=as.numeric(dat$logE_M_pearson_q)

input = 'D:/GitHub/Hormone-Receptors/result/Meth.correlation/BLCA.cor.mRNA_Meth.txt'
probe = read.csv(input,header = T,sep='\t',row.names = 1,stringsAsFactors = F)
plot(probe$logE_M_pearson_cor,log(probe$logE_M_pearson_q,10))
plot(probe$logE_M_pearson_cor,-log(probe$logE_M_pearson_q,10),pch=20,cex=0.1)
a=density(-log(as.numeric(probe$logE_M_pearson_q,10))) ; plot(a)

setwd('D:/GitHub/Hormone-Receptors/')
path='D:/GitHub/Hormone-Receptors/result/clusterprofiler/f3'
tag = 'BP.simplify'
files = list.files(path,paste(tag,'.txt',sep=''),full.names = T)
tmp = data.frame()
for(f in files){
  name = strsplit(strsplit(strsplit(f,'/',fixed=T)[[1]][7],'.',fixed=T)[[1]][1],'_',fixed = T)[[1]][1]
  df = read.csv(f,header = T,sep='\t',stringsAsFactors = F)
  if(nrow(df)<1){print(name);print('no pathway !');tmp[,name]=NA;next}
  for(l in 1:nrow(df)){
    tmp[df[l,]$Description,name] = df[l,]$qvalue
    }
}
tmp = data.frame(Name=rownames(tmp),Number=apply(tmp,1,function(x){sum(!is.na(x))}),tmp)
tmp = tmp[order(tmp$Number,decreasing=T),] ; rownames(tmp) = paste('( ',tmp$Number,' cancers )  ',rownames(tmp),sep='')
use = tmp[tmp$Number>1,3:35] ; use[!is.na(use)] =1 ; use[is.na(use)] =0
#p = pheatmap(use,cluster_rows = T,fontsize_row =7.5)
tiff(paste('./result/clusterprofiler/files.summary/',tag,'.summary.tiff',sep=''),width=36,height=18,units='cm',res=600,compression='lzw')
p = pheatmap(use,cluster_rows = T,fontsize_row =7,cutree_rows = 2) ;dev.off()
use.out = data.frame(Cluster=as.data.frame(cutree(p$tree_row,k=2))[rownames(use),],tmp[rownames(use),])
use.out = rbind(Paths_number=as.character(c(rep(NA,3),colSums(!is.na(use.out[,4:35])))), use.out)
write.table(use.out,paste('./result/clusterprofiler/files.summary/',tag,'.summary.enrich.result.txt',sep=''),col.names = T,row.names = F,sep='\t',quote = F)
