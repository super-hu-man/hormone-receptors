# 30% cut off diff analysis
library(limma)
library(edgeR)
library(DESeq2)
argv = commandArgs(trailingOnly = T)
fiho = sort(list.files('../download.firehose','uncv2.mRNAseq_raw_counts.txt'))
ucsc = sort(list.dirs('../download.cbioprotal')[2:34])

gene = 'AR|367'    ##
cutoff=30      ##
stap = argv ; print(stap)
file.create(paste('./result/diff/sample.compare.start-',stap,'.txt',sep=''))
file.create('./result/diff/sample.share.txt')
for(n in c(stap:stap)){
#for(n in c(stap:33)){
pro = strsplit(fiho[n],'.',fixed=T)[[1]][1]  ##

## input files
cli = read.csv('./ucsc/clinical.literature.final_use.txt',header = T,sep='\t',stringsAsFactors = F)
dat = read.csv(paste(ucsc[n],'/data_expression_merged.txt',sep=''),header = T,sep='\t',stringsAsFactors = F)
rownames(dat) = paste(dat$Hugo_Symbol,dat$Entrez_Gene_Id,sep='|')
datu = dat[,order(dat[which(rownames(dat)==gene),3:ncol(dat)])+2]
rsem = read.csv(paste('../download.firehose/',pro,'.uncv2.mRNAseq_raw_counts.txt',sep=''),header=T,sep='\t',row.names=1,stringsAsFactors=F)
share = intersect(colnames(datu),colnames(rsem))
b = sort(unique(substr(share,1,12)))
u =sort(unique(substr(colnames(datu),1,12)))[sort(unique(substr(colnames(datu),1,12)))%in%sort(unique(substr(colnames(rsem),1,12)))=='FALSE']
r =sort(unique(substr(colnames(rsem),1,12)))[sort(unique(substr(colnames(rsem),1,12)))%in%sort(unique(substr(colnames(datu),1,12)))=='FALSE']
write.table(paste(c('Both',b),collapse = '\t'),paste('./result/diff/sample.share.start-',stap,'.txt',sep=''),col.names = F,row.names = F,quote = F,sep='\t',append = T)
write.table(paste(c('UCSC',u),collapse = '\t'),paste('./result/diff/sample.share.start-',stap,'.txt',sep=''),col.names = F,row.names = F,quote = F,sep='\t',append = T)
write.table(paste(c('RSEM',r),collapse = '\t'),paste('./result/diff/sample.share.start-',stap,'.txt',sep=''),col.names = F,row.names = F,quote = F,sep='\t',append = T)

## DESeq2 diff raw count
print(paste(pro,n,'1',sep=' * '))
group.a = data.frame(sample=colnames(rsem),gp='T',stringsAsFactors = F)
group.a$gp[which(substr(colnames(rsem),14,14)!=0)] = 'N'
rsem = rsem[rowSums(rsem)>2*ncol(rsem),]
rsem = rsem[apply(rsem,1,function(x) sum(x=='0')<length(x)/2),]
all = DESeqDataSetFromMatrix(round(rsem),group.a,design=~1)
allnct = counts(estimateSizeFactors(all), normalized=T)
allnct = allnct[,order(allnct[which(rownames(allnct)==gene),])]
allnct = allnct[,substr(colnames(allnct),14,14)==0]
num = ncol(allnct) * cutoff/100
matrix = rsem[,c(colnames(allnct)[1:ceiling(num)],colnames(allnct)[(ncol(allnct)+1-ceiling(num)):ncol(allnct)])]
group = factor(c(rep('Low',ceiling(num)),rep('High',ceiling(num))),levels=c('Low','High'))
group = data.frame(names=colnames(matrix),group=group)
des = DESeqDataSetFromMatrix(countData = round(matrix),colData = group,design = ~group)
des2 = DESeq(des)
res = results(des2)
resOrdered = res[order(res$padj),]
resOrdered=as.data.frame(resOrdered)
out = data.frame(Gene=rownames(resOrdered),resOrdered)
lgn = colnames(allnct)[1:ceiling(num)]
hgn = colnames(allnct)[(ncol(allnct)+1-ceiling(num)):ncol(allnct)]
write.table(paste(c('Lgroup',lgn),collapse='\t'),paste('./result/diff/',pro,'.deseq.f1.txt',sep=''),
            col.names=F,row.names=F,sep='\t',quote = F)
write.table(paste(c('Hgroup',hgn),collapse='\t'),paste('./result/diff/',pro,'.deseq.f1.txt',sep=''),
            col.names=F,row.names=F,sep='\t',quote = F,append=T)
write.table(paste(colnames(out),collapse='\t'),paste('./result/diff/',pro,'.deseq.f1.txt',sep=''),
            col.names=F,row.names=F,sep='\t',quote=F,append = T)
write.table(out,paste('./result/diff/',pro,'.deseq.f1.txt',sep=''),col.names=F,row.names=F,sep='\t',quote = F,append=T)
x = c(colnames(allnct)[1:ceiling(num)],'.',colnames(allnct)[(ncol(allnct)-floor(num)):ncol(allnct)])

## t.test use ucsc normalized log2 RSEM+1
datu = na.omit(datu)
datu = datu[,substr(colnames(datu),14,14)==0]
datu = datu[rowSums(datu)>2*ncol(datu),]
datu = datu[apply(datu,1,function(x) sum(x=='0')<length(x)/2),]
num = ncol(datu) * cutoff/100
Lname = colnames(datu)[1:ceiling(num)]
Hname = colnames(datu)[(ncol(datu)+1-ceiling(num)):ncol(datu)]
print(paste(pro,n,'2',sep=' * '))
Lg = datu[,Lname]
Hg = datu[,Hname]
write.table(paste(c('Lgroup',Lname),collapse='\t'),paste('./result/diff/',pro,'.wilc.t.f2.txt',sep=''),
            col.names=F,row.names=F,sep='\t',quote = F)
write.table(paste(c('Hgroup',Hname),collapse='\t'),paste('./result/diff/',pro,'.wilc.t.f2.txt',sep=''),
            col.names=F,row.names=F,sep='\t',quote = F,append=T)
p.value = NULL
for(i in 1:nrow(datu)){
  Lr = paste(round(range(log(as.numeric(Lg[i,])+1,2)),5),collapse = '-')
  Lm = mean(log(as.numeric(Lg[i,])+1,2))
  Hr = paste(round(range(log(as.numeric(Hg[i,])+1,2)),5),collapse = '-')
  Hm = mean(log(as.numeric(Hg[i,])+1,2))
  fc = 2**(Hm-Lm)
  pv = t.test(log(as.numeric(Lg[i,])+1,2),log(as.numeric(Hg[i,])+1,2))$p.value
  pv2 = wilcox.test(log(as.numeric(Lg[i,])+1,2),log(as.numeric(Hg[i,])+1,2),exact=T)$p.value
  p.value = rbind(p.value,data.frame(
    Gene=rownames(datu)[i],L.range=Lr,L.mean=Lm,H.range=Hr,H.mean=Hm,fc=fc,t.test=pv,q.value.t=NA,wilcox=pv2,q.value.w=NA))
}
p.value$q.value.t = p.adjust(p.value$t.test,method = 'BH',nrow(p.value))
p.value$q.value.w = p.adjust(p.value$wilcox,method = 'BH',nrow(p.value))
p.value = p.value[order(p.value$q.value.w),]
write.table(paste(colnames(p.value),collapse='\t'),paste('./result/diff/',pro,'.wilc.t.f2.txt',sep=''),
            col.names=F,row.names=F,sep='\t',quote=F,append=T)
write.table(p.value,paste('./result/diff/',pro,'.wilc.t.f2.txt',sep=''),col.names=F,row.names=F,sep='\t',quote=F,append=T)

## DESeq with order
print(paste(pro,n,'3',sep=' * '))
rsem = rsem[,intersect(c(Lname,Hname),share)]
rsem = round(as.matrix(rsem))
group = factor(c(rep('Low',sum(Lname%in%colnames(rsem))),rep('High',sum(Hname%in%colnames(rsem)))),levels=c('Low','High'))
group = data.frame(names=colnames(rsem),group=group)
dds = DESeqDataSetFromMatrix(countData = rsem,colData = group,design = ~group)
dds2 <- DESeq(dds)
res <-  results(dds2)
resOrdered = res[order(res$padj),]
resOrdered=as.data.frame(resOrdered)
out = data.frame(Gene=rownames(resOrdered),resOrdered)
write.table(paste(c(paste('Lgroup',length(intersect(Lname,share)),sep=':'),intersect(Lname,share)),collapse='\t'),
            paste('./result/diff/',pro,'.deseq.f3.txt',sep=''),col.names=F,row.names=F,sep='\t',quote = F)
write.table(paste(c(paste('Lgroup',length(intersect(Hname,share)),sep=':'),intersect(Hname,share)),collapse='\t'),
            paste('./result/diff/',pro,'.deseq.f3.txt',sep=''),col.names=F,row.names=F,sep='\t',quote = F,append=T)
write.table(paste(colnames(out),collapse='\t'),paste('./result/diff/',pro,'.deseq.f3.txt',sep=''),
            col.names=F,row.names=F,sep='\t',quote=F,append = T)
write.table(out,paste('./result/diff/',pro,'.deseq.f3.txt',sep=''),col.names = T,row.names = F,sep='\t',quote = F)
y = c(Lname,'.',Hname)

xy =rbind(x,y)
write.table(pro,paste('./result/diff/sample.compare.start-',stap,'.txt',sep=''),col.names = F,row.names = F,sep='\t',quote = F,append = T)
write.table(xy,paste('./result/diff/sample.compare.start-',stap,'.txt',sep=''),col.names = F,row.names = F,sep='\t',quote = F,append = T)
}

print('work done !')
