library('ggplot2')
library(ggbeeswarm)
library(survival)
library(gridExtra)
library(survminer)

dat = read.csv('../download.firehose/gdac.broadinstitute.org_LIHC.mRNAseq_Preprocess.Level_3.2016012800.0.0/LIHC.uncv2.mRNAseq_RSEM_normalized_log2.txt',header=T,row.names=1,sep='\t',stringsAsFactors=F)
colnames(dat) = gsub('[.]','-',colnames(dat))
gene= dat['TRIP13|9319',] ; gene = as.data.frame(t(gene)) ; colnames(gene) = 'Gene'
gene[which(substr(rownames(gene),14,14)==0),'Type'] = 'Tumor (n = 373)' 
gene[which(substr(rownames(gene),14,14)!=0),'Type'] = 'Normal (n = 50)'
gene$Type = factor(gene$Type,levels = c('Tumor (n = 373)','Normal (n = 50)'))
wp = wilcox.test(as.numeric(gene[gene$Type=='Tumor (n = 373)','Gene']),as.numeric(gene[gene$Type=='Normal (n = 50)','Gene']))$p.value
tp = t.test(as.numeric(gene[gene$Type=='Tumor (n = 373)','Gene']),as.numeric(gene[gene$Type=='Normal (n = 50)','Gene']))$p.value

p = ggplot(data=gene ,aes(x=Type , y=Gene ,fill = Type))
p = p + geom_quasirandom(aes(x=Type),colour='grey90',size=1.8,shape=21,width=0.07)
p = p + stat_boxplot(colour='grey33',alpha = 0.5,outlier.alpha = F,width=0.24)
p = p + annotate('text',x=1.5, y=12, label=paste('wilcox-test p-value : ',round(wp,30),sep=''))
p = p + annotate('text',x=1.5, y=13, label=paste('t-test p-value : ',round(tp,40),sep=''))
p = p + labs(title="LIHC expression of tumor and normal",x="Cancer Type",y ="Expression of TRIP13 mRNA log2(RSEM + 1)",size=10)
p = p + theme(plot.title = element_text(hjust = 0.5), axis.text=element_text(size=12),axis.title=element_text(size=13),
              legend.text = element_text(size = 12),
              panel.border=element_rect(fill=NA),panel.background = element_rect(fill='white'))
tiff('LIHC.expression.tiff',width =15,height = 15,bg="white",units='cm',compression='lzw' ,res=500)
p ; dev.off()

cli = read.csv('./ucsc/clinical.literature.final_use.txt',header = 1,row.names=1,sep='\t',stringsAsFactors = F)
tsample = gene[which(substr(rownames(gene),14,14)==0),] ; tsample = tsample[order(rownames(tsample)),]
tsample2 = tsample[!duplicated(substr(rownames(tsample),1,12)),] ; rownames(tsample2) = substr(rownames(tsample2),1,12)
lihc = cbind(tsample2,cli[rownames(tsample2),])
lihc$Mean = 'L' ; lihc[lihc$Gene>mean(lihc$Gene),'Mean'] = 'H'
lihc$Median = 'L' ; lihc[lihc$Gene>median(lihc$Gene),'Median'] = 'H'
lihc$OS.time = round(lihc$OS.time /30,2) ; lihc$PFI.time = round(lihc$PFI.time /30,2)
out = data.frame(Sample=row.names(lihc),lihc)
write.table(out,'LIHC.CLI.txt',row.names=F,col.names=T,sep='\t',quote=F)
# 1
os = survfit(formula= Surv(lihc$OS.time, lihc$OS=='1')~Mean ,data=lihc)
pv1 = survdiff(formula= Surv(lihc$OS.time, lihc$OS=='1')~Mean ,data=lihc)
do = ggsurvplot(os,legend.title="Group: ",legend.labs=c('High','Low'),size=0.6,xlab='Survival Months')
d1 = do$plot + labs(title=paste('TRIP13','  OS',sep=''))
d1 = d1 + annotate('text',x=107, y=0.52, label='TRIP13  low  (n = 185)',size=4)
d1 = d1 + annotate('text',x=106, y=0.16, label='TRIP13  high  (n = 185)',size=4)
d1 = d1 + annotate('text',x=30, y=0.05, label=paste('Log-rank test, p = ',round(1-pchisq(pv1$chisq, length(pv1$n)-1),8),sep=''),size=5)
d1 = d1 + theme(axis.text=element_text(size=15),axis.title=element_text(size=50),plot.title=element_text(hjust = 0.5),
                legend.position='top',legend.justification=1)
tiff('survival.mean.OS.tiff',width =16,height = 16,bg="white",units='cm',compression='lzw' ,res=500)
d1 ; dev.off()
# 2
os = survfit(formula= Surv(lihc$PFI.time, lihc$PFI=='1')~Mean ,data=lihc)
pv1 = survdiff(formula= Surv(lihc$PFI.time, lihc$PFI=='1')~Mean ,data=lihc)
do = ggsurvplot(os,pval=F,legend.title="Group: ",legend.labs=c('High','Low'),pval.size=3,pval.coord=c(0,0.1),size=0.6,xlab='Survival Months')
d1 = do$plot+labs(title=paste('TRIP13','  PFI',sep=''))+theme(plot.title=element_text(hjust = 0.5),legend.position='top',legend.justification=1)
d1 = d1 + annotate('text',x=107, y=0.29, label='TRIP13  low  (n = 185)',size=4)
d1 = d1 + annotate('text',x=106, y=0.15, label='TRIP13  high  (n = 185)',size=4)
d1 = d1 + annotate('text',x=30, y=0.05, label=paste('Log-rank test, p = ',round(1-pchisq(pv1$chisq, length(pv1$n)-1),8),sep=''),size=5)
tiff('survival.mean.PFI.tiff',width =16,height = 16,bg="white",units='cm',compression='lzw' ,res=500)
d1 ; dev.off()
# 3
os = survfit(formula= Surv(lihc$OS.time, lihc$OS=='1')~Median ,data=lihc)
pv1 = survdiff(formula= Surv(lihc$OS.time, lihc$OS=='1')~Median ,data=lihc)
do = ggsurvplot(os,pval=F,legend.title="Group: ",legend.labs=c('High','Low'),pval.size=3,pval.coord=c(0,0.1),size=0.6,xlab='Survival Months')
d1 = do$plot+labs(title=paste('TRIP13','  OS',sep=''))+theme(plot.title=element_text(hjust = 0.5),legend.position='top',legend.justification=1)
d1 = d1 + annotate('text',x=107, y=0.52, label='TRIP13  low  (n = 184)',size=4)
d1 = d1 + annotate('text',x=106, y=0.16, label='TRIP13  high  (n = 186)',size=4)
d1 = d1 + annotate('text',x=30, y=0.05, label=paste('Log-rank test, p = ',round(1-pchisq(pv1$chisq, length(pv1$n)-1),8),sep=''),size=5)
tiff('survival.median.OS.tiff',width =16,height = 16,bg="white",units='cm',compression='lzw' ,res=500)
d1 ; dev.off()
# 4
os = survfit(formula= Surv(lihc$PFI.time, lihc$PFI=='1')~Median ,data=lihc)
pv1 = survdiff(formula= Surv(lihc$PFI.time, lihc$PFI=='1')~Median ,data=lihc)
do = ggsurvplot(os,pval=F,legend.title="Group: ",legend.labs=c('High','Low'),pval.size=3,pval.coord=c(0,0.1),size=0.6,xlab='Survival Months')
d1 = do$plot+labs(title=paste('TRIP13','  PFI',sep=''))+theme(plot.title=element_text(hjust = 0.5),legend.position='top',legend.justification=1)
d1 = d1 + annotate('text',x=107, y=0.29, label='TRIP13  low  (n = 184)',size=4)
d1 = d1 + annotate('text',x=106, y=0.15, label='TRIP13  high  (n = 186)',size=4)
d1 = d1 + annotate('text',x=30, y=0.05, label=paste('Log-rank test, p = ',round(1-pchisq(pv1$chisq, length(pv1$n)-1),8),sep=''),size=5)
tiff('survival.median.PFI.tiff',width =16,height = 16,bg="white",units='cm',compression='lzw' ,res=500)
d1 ; dev.off()

