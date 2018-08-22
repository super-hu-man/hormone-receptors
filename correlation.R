##
args=commandArgs( trailingOnly = T)
n = as.numeric(args)
gene = 'AR'
files = list.dirs('../download.cbioprotal',full.names = T)[2:34]
for(f in files[n]){
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
  test.res=data.frame(Gene=numeric(0), Meth_ID=numeric(0), 
                      logE_M_pearson_cor=numeric(0) , logE_M_pearson_p=numeric(0) , logE_M_pearson_q=numeric(0) ,
                      E_M_pearson_cor=numeric(0) , E_M_pearson_p=numeric(0) , E_M_pearson_q=numeric(0) ,
                      logE_M_spearman_cor=numeric(0) , logE_M_spearman_p=numeric(0) , logE_M_spearman_q=numeric(0) ,
                      E_M_spearman_cor=numeric(0) , E_M_spearman_p=numeric(0) , E_M_spearman_q=numeric(0) )
  for ( one in rownames(metd) ){
    if ( sum(is.na(metd[one,])) > ncol(metd)*0.2 ) {next}
    pl = cor.test(as.numeric(log(datd+1,2)), as.numeric(metd[one,]) , method='pearson')
    pc = cor.test(as.numeric(datd), as.numeric(metd[one,]) , method='pearson')
    sl = cor.test(as.numeric(log(datd+1,2)), as.numeric(metd[one,]) , method='spearman')
    sc = cor.test(as.numeric(datd), as.numeric(metd[one,]) , method='spearman')
    test.res[one,] = c(gene ,one , pl$estimate , pl$p.value , NA ,
				   pc$estimate , pc$p.value , NA ,
				   sl$estimate , sl$p.value , NA ,
				   pc$estimate , pc$p.value , NA ) }
  test.res$logE_M_pearson_q = p.adjust( test.res$logE_M_pearson_p ,'BH')
  test.res$E_M_pearson_q = p.adjust( test.res$E_M_pearson_p ,'BH')
  test.res$logE_M_spearman_q = p.adjust( test.res$logE_M_spearman_p ,'BH')
  test.res$E_M_spearman_q = p.adjust( test.res$E_M_spearman_p ,'BH')
  
  write.table(test.res,paste('./result/Meth.correlation/',name,'.cor.RNAvsMeth.txt',sep=''),col.names=T,row.names=F,sep='\t',quote=F)
}
print('done')
