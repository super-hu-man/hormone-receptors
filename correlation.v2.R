##
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
  pc = cor.test(datd1, as.numeric(one) , method='pearson')
  sl = cor.test(datd2, as.numeric(one) , method='spearman')
  sc = cor.test(datd1, as.numeric(one) , method='spearman')
  x= c(Gene = gene, logE_M_pearson_cor=pl$estimate , logE_M_pearson_p=pl$p.value , logE_M_pearson_q=NA ,
                    E_M_pearson_cor=pc$estimate    , E_M_pearson_p=pc$p.value    , E_M_pearson_q=NA ,
                    logE_M_spearman_cor=sl$estimate , logE_M_spearman_p=sl$p.value , E_M_spearman_q=NA ,
                    E_M_spearman_cor=pc$estimate    , E_M_spearman_p=pc$p.value    , E_M_spearman_q=NA ) }
  return(x)
}

for (f in files[n]){
  name = toupper(strsplit(strsplit(f,'/',fixed=T)[[1]][3],'_',fixed=T)[[1]][1])
  print(name)
  dat = read.csv(paste(f,'/data_expression_merged.txt',sep=''),header = T,sep='\t',stringsAsFactors = F)
  datu = dat[which(dat$Hugo_Symbol==gene),order(dat[which(dat$Hugo_Symbol==gene),3:ncol(dat)])+2]
  colnames(datu) = gsub('[.]','-',colnames(datu))
  # meth 
  met = read.csv(paste('../download.pan-cancer-atlas/Meth.Cancer/',name,'.Meth.txt',sep=''),header=T,row.names=1,stringsAsFactor=F,sep='\t')
  colnames(met) = gsub('[.]','-',colnames(met))
  share = intersect(colnames(datu),colnames(met))
  metd = met[,share]
  datd = datu[,share]
  # test cor
  datd1 = as.numeric(datd) ; datd2 = as.numeric(log(datd+1,2))
  res = apply(metd, 1 , test , metd=metd , datd1=datd1 , datd2=datd2 , gene=gene)
  red = as.data.frame(t(res),stringsAsFactors=F)
  colnames(red) = c('Gene','logE_M_pearson_cor','logE_M_pearson_p','logE_M_pearson_q',
                           'E_M_pearson_cor','E_M_pearson_p','E_M_pearson_q',
                           'logE_M_spearman_cor','logE_M_spearman_p','logE_M_spearman_q',
                           'E_M_spearman_cor','E_M_spearman_p','E_M_spearman_q')
  # adjust p value
  red = red[is.na(red$Gene)=='FALSE',]
  red$logE_M_pearson_q = p.adjust( red$logE_M_pearson_p ,'BH')
  red$E_M_pearson_q    = p.adjust( red$E_M_pearson_p ,'BH')
  red$logE_M_spearman_q = p.adjust( red$logE_M_spearman_p ,'BH')
  red$E_M_spearman_q    = p.adjust( red$E_M_spearman_p ,'BH')
  
  res.out = data.frame(Meth_ID=rownames(red),red)
  
  write.table(res.out , paste('./result/Meth.correlation/',name,'.cor.mRNA_Meth.txt',sep=''),col.names=T,row.names=F,sep='\t',quote=F)
}

print('done')
