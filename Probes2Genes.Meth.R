#setwd('C:/Users/root/Desktop/Meth_ID_match/')
#setwd('D:/GitHub/download.pan-cancer-altas/Meth.Cancer')
dat = read.csv('../Files.useful/GPL13534-11288.txt',header = T,row.names = 1,sep='\t',stringsAsFactors = F,skip = 37)

dat2 = dat[dat$UCSC_RefGene_Name!='',]
x=sapply(dat2$UCSC_RefGene_Name,function(i){x=unlist(strsplit(i,';',fixed=T));y=length(x);return(y)})
N=unlist(sapply(dat2$UCSC_RefGene_Name,function(i){x=unlist(strsplit(i,';',fixed=T));return(x)}))
A=unlist(sapply(dat2$UCSC_RefGene_Accession,function(i){x=unlist(strsplit(i,';',fixed=T));return(x)}))
G=unlist(sapply(dat2$UCSC_RefGene_Group,function(i){x=unlist(strsplit(i,';',fixed=T));return(x)}))
dat3 = dat2[rep(1:nrow(dat2),x),]
dat3$UCSC_RefGene_Name = N ; dat3$UCSC_RefGene_Group = G

use.regin = c('TSS1500' , 'TSS200' , '5\'UTR' , '1stExon' )
promoter = dat3[ which(dat3$UCSC_RefGene_Group %in% use.regin) , ]
promoter2 = promoter[!duplicated(promoter[,c('Name','UCSC_RefGene_Name','UCSC_RefGene_Accession','UCSC_RefGene_Group')]),]
tmp = promoter2[,c('Name','UCSC_RefGene_Name')]
dated = promoter2[!duplicated(tmp),]

rm(ls()[-which(ls()=='dated')]) ; gc()
region = unique(unlist(sapply(dated$UCSC_RefGene_Group,function(x){x=strsplit(x,';',fixed=T)[[1]];return(x)})))
genes  = unique(unlist(sapply(dated$UCSC_RefGene_Name,function(x){x=strsplit(x,';',fixed=T)[[1]];return(x)})))
probes = unique(dated$Name)

## argvs is the start and end number of files order , like ( 1:33 )
argvs = commandArgs(trailingOnly = T)
n = as.numeric(argvs)
files = list.files('../download.pan-cancer-altas/Meth.Cancer','Meth.txt.gz',full.names=T)
out.path = './result/Meth_gene.correlation/'
if(!dir.exists(out.path)){dir.create(out.path)}
avage = function(x,dated,txt){
  probe = dated[dated$UCSC_RefGene_Name == x,'Name']
  print(length(probe))
  genedat = txt[rownames(txt) %in% probe,]
  m = colMeans(genedat,na.rm = T)
  index = which(apply(genedat,2,function(x){sum(!is.na(x))<6 & sum(!is.na(x))<0.75*length(probe)}))
  m[index] = NA
  n = sum(is.na(m))
  return(m)
}
for (f in files[n]){
  name = strsplit(strsplit(f,'/',fixed=T)[[1]][length(strsplit(f,'/',fixed=T)[[1]])],'.',fixed=T)[[1]][1]
  txt = read.csv(gzfile(f),header = T,row.names = 1,sep='\t',stringsAsFactors = F)
  colnames(txt) = gsub('[.]','-',colnames(txt))
  res = lapply(genes[1:100],avage,dated=dated,txt=txt)
  res = t(data.frame(res)) ; rownames(res) = genes[1:100]
  res = res[rowSums(is.na(res))<ncol(txt)*0.5,]
  sums= table(unlist(apply(res,1,function(x){sum(is.na(x))})))
  print(f);print(sums)
  res = data.frame(Gene=rownames(res),res)
  write.table(res,paste(out.path,name,'.Meth.gene.txt',sep=''),col.names=T,row.names=F,sep='\t',quote=F)
}

