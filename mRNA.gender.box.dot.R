##
args='mRNA.RSEM'  # commandArgs( trailingOnly = T)  ## 'RSEM' or 'Protein.RPPA' 
library('ggplot2')
library(ggbeeswarm)
all.pro = c('ACC','BLCA','BRCA','CESC','CHOL','COAD','DLBC','ESCA','GBM','HNSC',
            'KICH','KIRC','KIRP','LAML','LGG','LIHC','LUAD','LUSC','MESO','OV',
            'PAAD','PCPG','PRAD','READ','SARC','SKCM','STAD','TGCT','THCA','THYM',
            'UCEC','UCS','UVM')
gender.pro = c('ACC','BLCA','BRCA','CHOL','COAD','DLBC','ESCA','GBM','HNSC',
               'KICH','KIRC','KIRP','LAML','LGG','LIHC','LUAD','LUSC','MESO',
               'PAAD','PCPG','READ','SARC','SKCM','STAD','THCA','THYM','UVM')
pro.list = c('BLCA','BRCA','COAD','HNSC','KICH','KIRC','KIRP','LIHC','LUAD','LUSC','PRAD','STAD','THCA','UCEC')
cli = read.table('./ucsc/clinical.literature.final_use.txt',header = T,sep='\t',row.names = 1,stringsAsFactors = F)
table = read.table(paste('./ucsc/',args,'.combined.txt',sep=''),header = T,row.names=1,sep='\t',stringsAsFactors = F)
table$Cancer = cli[table$bcr_patient_barcode,'type'] 
table$Gender = cli[table$bcr_patient_barcode,'gender']
table = table[is.na(table$Cancer)=='FALSE',]
t.e = table[substr(rownames(table),14,14)==0,]
n.e = table[substr(rownames(table),14,14)!=0,]

#dat = read.table('mRNA.RSEM.Tumors.matchinfor.txt',header = T,row.names=1,sep='\t',stringsAsFactors = F)
dat = t.e
d=NULL
for(gene in c('AR','ESR1','ESR2','PGR')) {
    df=dat[,c('bcr_patient_barcode','Cancer',gene,'gender')];
    colnames(df)[3]='exp' ; df$Gene=gene
    df$index=which(gene==c('AR','ESR1','ESR2','PGR'))
    d=rbind(d,df)
}
d[d$gender=='FEMALE','index']=d[d$gender=='FEMALE','index']-0.2
d[d$gender=='MALE','index']=d[d$gender=='MALE','index']+0.2
p = ggplot(data=d ,aes(x=Gene , y=exp, fill = gender))
p = p + geom_boxplot(alpha = 0,outlier.alpha = F,width=0.8,notch=F,notchwidth=0.1)
p = p + geom_quasirandom(aes(x=index),colour='grey85',size=1.8,shape=21,width=0.2)
p = p + geom_boxplot(colour='grey30',alpha = 0.3,outlier.alpha = F,width=0.8,notch=F,notchwidth=0.1,size=0.7)
tiff(filename = paste('All',"mRNA.gender.tiff",sep='.'),width =24,height = 16,units ="cm",compression="lzw",bg="white",res=1000)
plot(p)
dev.off()


df = NULL ; p.values = NULL
for(pro in gender.pro){
    t.p = t.e[t.e[,'Cancer']==pro,]
    t.u = t.p
    index = which(gender.pro==pro)
    t.u$type = 'Tumor'
    e.m = t.u[t.u$Gender=='MALE',gop]
    e.f = t.u[t.u$Gender=='FEMALE',gop]
    if(length(e.m)<=1 || length(e.f)<=1 ){p.v='.'} else {
      p.v = round(t.test(e.m,e.f)$p.value,3)
      p.w = round(wilcox.test(e.m,e.f)$p.value,3)}
    p.values = rbind(p.values , c(pro, p.v, p.w))
    data = t.u[,c(3,5,9,10)] ; data[data[,3]=='male',5]=index+0.15 ;data[data[,3]=='female',5]=index-0.15
    colnames(data) = c('exp','project','gender','Samples','point')
    df = rbind(df , data)
}
colnames(p.values) = c('pro','p.value') ; rownames(p.values) = p.values[,1]
colour = c("darkred","#3300FFFF")

p = ggplot(data=df ,aes(x=project , y=exp, fill = gender))
p = p + geom_boxplot(alpha = 0,outlier.alpha = F,width=0.6,notch=F,notchwidth=0.1,size=0.7)
p = p + geom_quasirandom(aes(x=point),colour='grey80',size=1.6,shape=21,width=0.12)
p = p + geom_boxplot(colour='grey33',alpha = 0.3,outlier.alpha = F,width=0.6,notch=F,notchwidth=0.1,size=0.7)
if (length(qv[qv<0.05])){
    p = p + annotate('text',x=gender.pro[which(qv<0.05)], y=rep(max(df[,1])+1,length(which(qv<0.05))), label=paste('',format(qv)[which(qv<0.05)],sep=''),size=3,colour='red') }
p = p + labs(title=paste(args[1],"differential expression between male and female",sep=' '),x="Cancer Type",y ="Expression of mRNA(RSEM log2 )")
p = p + theme(plot.title = element_text(hjust = 0.5))

tiff(filename = paste(args[1],"mRNA.gender.tiff",sep='.'),width =35,height = 18,units ="cm",compression="lzw",bg="white",res=1000)
plot(p)
dev.off()


## all
args=commandArgs( trailingOnly = T)
library(ggplot2)
library(ggbeeswarm)

dat = read.table('mRNA.RSEM.Tumors.matchinfor.txt',header = T,row.names=1,sep='\t',stringsAsFactors = F)
d=NULL
for(gene in c('AR','ESR1','ESR2','PGR')) {
    df=dat[,c('bcr_patient_barcode','type',gene,'gender')];
    colnames(df)[3]='exp' ; df$Gene=gene
    df$index=which(gene==c('AR','ESR1','ESR2','PGR'))
    d=rbind(d,df)
}
d[d$gender=='FEMALE','index']=d[d$gender=='FEMALE','index']-0.2
d[d$gender=='MALE','index']=d[d$gender=='MALE','index']+0.2
d=d[d$exp!='None',]
d$exp = as.numeric(d$exp)

p = ggplot(data=d ,aes(x=Gene , y=exp, fill = gender))
p = p + geom_boxplot(alpha = 0,outlier.alpha = F,width=0.8,notch=F,notchwidth=0.1)
p = p + geom_quasirandom(aes(x=index),colour='grey85',size=1.8,shape=21,width=0.2)
p = p + geom_boxplot(colour='grey30',alpha = 0.3,outlier.alpha = F,width=0.8,notch=F,notchwidth=0.1,size=0.7)
p = p + labs(title=paste('All cancer',"differential expression between male and female",sep=' '),x="Gene",y ="Expression of mRNA(RSEM log2 )")
p = p + theme(plot.title = element_text(hjust = 0.5))

tiff(filename = paste('All',"mRNA.gender.tiff",sep='.'),width =24,height = 16,units ="cm",compression="lzw",bg="white",res=300)
plot(p)
dev.off()
