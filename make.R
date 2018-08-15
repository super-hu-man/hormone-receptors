library(survival)

setwd('D:/GitHub/hormone-receptors/cli')
files = list.files('./','na')
for(f in files){
  #f=files[1]
  dat =read.csv(f,row.names = 1,header = T,sep='\t',stringsAsFactors = F)
  dat=dat[-c(1,2),]
  name = strsplit(f,'_',fixed=T)[[1]][4]
  tags = colnames(dat)
  h=data.frame(c('Tag','number_NA','other'))
  write.table(t(h),name,sep='\t',quote=F,row.names=F,col.names=F)
  del = read.csv('../del.tags',header = F,stringsAsFactors = F)[,1]
  for(t in tags){
    if (t %in% del){next}
    na.n = length(which(is.na(dat[,t])))
    if(na.n>=nrow(dat)*0.9 || max(table(dat[,t])==nrow(dat))){next}
    na.f = round(100*na.n/nrow(dat),2)
    mat = t(data.frame(table(dat[,t])))
    o1 = c(t,paste(na.n,'[',na.f,'] | N/A',sep=''))
    if(ncol(mat)<=1){
      fre = round(100*as.numeric(mat[2])/nrow(dat),2)
      o2=paste(mat[2],'[',fre,'%] | ',mat[1],sep='') }else{
        mat = mat[,order(mat[2,],decreasing=T)]
        fre = round(100*as.numeric(mat[2,])/nrow(dat),2)
        o2=paste(mat[2,],'[',fre,'%] | ',mat[1,],sep='') }
    out = t(data.frame(c(o1,o2)))
    write.table(out,name,sep='\t',quote=F,append=T,col.names=F,row.names=F)
  }
}