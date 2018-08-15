a=read.csv('y.txt',header=F,sep='\t',stringsAsFactors=F)
b=read.csv('keep.txt',header=F,sep='\t',stringsAsFactors=F)

tag = unique(b[,1:2])
x=NULL
for(i in c(1:nrow(tag))){
  l=a[intersect(which(a[,1]==tag[i,1]),which(a[,2]==tag[i,2])),]
  x = rbind(x,l)
}

write.table(x,'un.cox.ad.txt',row.names = F,col.names = F,sep='\t',quote = F)
