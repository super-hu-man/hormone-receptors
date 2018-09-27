## setwd('D:/GitHub/Hormone-Receptors')
library(doSNOW)
library(foreach)
library(parallel)
args=commandArgs( trailingOnly = T)
n = as.numeric(args)
gene = 'AR'
files = list.dirs('../download.cbioprotal',full.names = T)[2:34]
# define test function
test = function(one,metd,datd1,datd2,gene){
  if ( sum(is.na(one)) > ncol(metd)*0.2 ) {x= rep(NA,13)}else{
  pl = cor.test(datd2, as.numeric(one) , method='pearson')
  sl = cor.test(datd2, as.numeric(one) , method='spearman')
  x= c(Gene = gene, logE_M_pearson_cor=pl$estimate , logE_M_pearson_p=pl$p.value , logE_M_pearson_q=NA ,
                    logE_M_spearman_cor=sl$estimate , logE_M_spearman_p=sl$p.value , E_M_spearman_q=NA ) }
  return(x)
}

for (f in files){
  name = toupper(strsplit(strsplit(f,'/',fixed=T)[[1]][3],'_',fixed=T)[[1]][1])
  print(name)
  dat = read.csv(paste(f,'/data_expression_merged.txt',sep=''),header = T,sep='\t',stringsAsFactors = F)
  rownames(dat) = paste(dat$Hugo_Symbol,dat$Entrez_Gene_Id,sep='|')
  colnames(dat) = gsub('[.]','-',colnames(dat))
  AR = dat[which(dat$Hugo_Symbol==gene),3:ncol(dat)]
  # test cor
  dat=log(dat[,3:ncol(dat)]+1,2) ; dat = dat[which(apply(dat,1,function(x){sum(x==0)<length(x)*0.3})),] 
  datd1 = as.numeric(AR) ; datd2 = as.numeric(log(AR+1,2))
  res = apply(dat, 1 , test , metd=dat , datd1=datd1 , datd2=datd2 , gene=gene)
  red = as.data.frame(t(res),stringsAsFactors=F)
  colnames(red) = c('Gene','logE_M_pearson_cor','logE_M_pearson_p','logE_M_pearson_q',
                           'logE_M_spearman_cor','logE_M_spearman_p','logE_M_spearman_q')
  # adjust p value
  red = red[is.na(red$logE_M_pearson_cor)=='FALSE',]
  red$logE_M_pearson_q = p.adjust( red$logE_M_pearson_p ,'BH')
  red$logE_M_spearman_q = p.adjust( red$logE_M_spearman_p ,'BH')
  
  res.out = data.frame(Meth_ID=rownames(red),red)
  write.table(res.out , paste('./result/mRNA.correlation/',name,'.cor.mRNA_Meth.txt',sep=''),col.names=T,row.names=F,sep='\t',quote=F)
}
print('done')

corf = list.files('./result/mRNA.correlation','.cor.mRNA_Meth.txt',full.names = T)
plot(c(-1,1.5)~c(0,0),type='l',xlim=c(-1,1),ylim=c(0,5), lty=2)
tmp = data.frame()
for (c in corf ){
  names = strsplit(strsplit(c,'/')[[1]][length(strsplit(c,'/')[[1]])],'[.]')[[1]][1] ; print(names)
  cors = read.csv(c,header = T,sep='\t',stringsAsFactors = F,row.names=1)
  corfilter = cors[ cors$logE_M_pearson_q<=0.05 ,]
  for(l in 1:nrow(corfilter)){
    line = corfilter[l,]
    tmp[rownames(corfilter)[l],names] = corfilter$logE_M_pearson_cor[l]
  }
  print(nrow(corfilter))
}
num = apply(tmp,1,function(x){sum(!is.na(x))})
tmps = data.frame(number=num,tmp) 
tmps = tmps[order(tmps$number,decreasing = T),]; tmps =tmps[-1,]
h = tmps[tmps$number>1,] ; h[is.na(h)]=0
pheatmap(h[,-1])


