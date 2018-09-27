gene = 'AR'
files = list.dirs('../download.cbioprotal',full.names = T)[2:34]
for (f in files ){
  name = toupper(strsplit(strsplit(f,'/',fixed=T)[[1]][3],'_',fixed=T)[[1]][1])
  "
  dat = read.csv(paste(f,'/data_expression_merged.txt',sep=''),sep='\t',header=T,stringsAsFactors = F)
  dat1 = dat[which(dat$Hugo_Symbol==gene),3:ncol(dat)]
  par(mfrow=c(1,2))
  x1 = hist(as.numeric(dat1),breaks = 50)
  x2 = hist(log(as.numeric(dat1)+1,2),breaks = 20)
  Sys.sleep(5)
  print(f)
  "
  met = read.csv(paste('../download.pan-cancer-altas/Meth.Cancer/',name,'.AR.txt',sep=''),header=T,row.names=1,stringsAsFactor=F,sep='\t')
  colnames(met) = gsub('[.]','-',colnames(met))
  met2 = colMeans(met)
  y=density(met2)
  plot(y)
  for(i in 1:nrow(met)){
    x=density(as.numeric(met[i,]),na.rm = T)
    plot(x)
    scan()
  }
}


