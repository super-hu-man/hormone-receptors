files = list.files('./result/Meth.correlation','cor.mRNA_Meth.txt',full.names = T)
for(f in files){
  name = strsplit(strsplit(f,'/',fixed=T)[[1]][length(strsplit(f,'/',fixed=T)[[1]])],'.',fixed=T)[[1]][1]
  print(name)
  dat = read.csv(f,header=T,sep='\t',stringsAsFactors=F,row.names=1)
  png(paste('./result/Meth.correlation/',name,'.cor-q.png',sep=''))
  plot(dat$loglogE_M_spearman_cor,-log(dat$logE_M_spearman_q,10))
  #plot(dat$logloglogE_M_spearman_cor,dat$loglogE_M_spearman_q)
  abline(h=1.30103) ;n=sum(dat$E_M_pearson_q<0.05);x=dat[dat$logE_M_spearman_q<0.05,]
  r1=range(x[x$loglogE_M_spearman_cor<0,'loglogE_M_spearman_cor']);r2=range(x[x$loglogE_M_spearman_cor>0,'loglogE_M_spearman_cor'])
  text(0,1.30103,paste(paste(r1,collapse=' - '),paste(r2,collapse=' - '),n,sep=' '))
  print(r1);print(r2);print(n)
  abline(h=2)
  dev.off()
  
}

file = list.files('./','AR.txt',full.names = T)
res = data.frame()
for(a in file){
  print(a)
  name = strsplit(strsplit(a,'/',fixed=T)[[1]][length(strsplit(a,'/',fixed=T)[[1]])],'.',fixed=T)[[1]][1]
  dat = read.csv(a,header=T,sep='\t',stringsAsFactors=F,row.names=1)
  #colnames(dat) = gsub(colnames(dat),'.','-')
  rna = read.csv(paste('../../download.cbioprotal/',tolower(name),'_tcga_pan_can_atlas_2018/data_expression.txt',sep=''),header=T,sep='\t',stringsAsFactors=F)
  share = intersect(colnames(dat),colnames(rna))
  AR = sort(rna[which(rna$Hugo_Symbol=='AR'),share])
  dat = dat[,colnames(AR)]
  #plot(as.numeric(AR[1,]))
  #for(i in 1:nrow(dat)){c=cor(as.numeric(AR),as.numeric(dat[i,]))}
  c = as.data.frame(apply(dat, 1, cor ,y=as.numeric(AR),use='pairwise.complete.obs'))
  colnames(c)=name
  if(ncol(res)==0){res=c}else{res=cbind(res,c)}
}
