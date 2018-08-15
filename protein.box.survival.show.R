args = commandArgs( trailingOnly = T)
input = './ucsc/Protein.RPPA.combined.txt'
gop = args
if (file.exists('./result')){print('dir result exist !')}else{dir.create('./result')}
if (file.exists(paste('./result/Protein.',gop,sep=''))){
  print('dir gene exist !')}else{dir.create(paste('./result/Protein.',gop,sep=''))}

library('ggplot2')
library(ggbeeswarm)
library(survival)
library(gridExtra)
library(survminer)

cli = read.table('./ucsc/clinical.literature.final_use.txt',header = T,sep='\t',row.names = 1,stringsAsFactors = F)
table = read.table(input,header = T,row.names=1,sep='\t',stringsAsFactors = F)
table$Cancer = cli[table$bcr_patient_barcode,'type']
table$Gender = cli[table$bcr_patient_barcode,'gender']
table$OS = cli[table$bcr_patient_barcode,'OS']
table$OS.time = cli[table$bcr_patient_barcode,'OS.time']
table$PFI = cli[table$bcr_patient_barcode,'PFI']
table$PFI.time = cli[table$bcr_patient_barcode,'PFI.time']
table = table[is.na(table$Cancer)=='FALSE',]
pro.list = sort(unique(table$Cancer))
t.e = table[substr(rownames(table),14,14)==0,]
n.e = table[substr(rownames(table),14,14)!=0,]

## step 1 : draw tumot of each cancer samples rsem ##
#####################################################
df = t.e 
df$exp = t.e[,gop]
p = ggplot(data=df ,aes(x=Cancer , y=exp,fill = Cancer))
p = p + geom_boxplot(alpha = 0,outlier.alpha = F,width=0.48,notch=F,notchwidth=0.1,size=0.3)
p = p + geom_quasirandom(aes(x=Cancer),colour='grey93',size=1.1,shape=21,width=0.1)
p = p + geom_boxplot(colour='grey40',alpha = 0.3,outlier.alpha = F,width=0.48,notch=F,notchwidth=0.1,size=0.3)
p = p + labs(title=paste(gop,"protein differential expression of tumor",sep=' '),x="Cancer Type",y ="Expression of protein")
p = p + theme(plot.title = element_text(hjust = 0.5), axis.text.x=element_text(size=7),legend.position ='none')
ggsave(filename = paste('./result/Protein.',gop,'/',gop,".protein.RPPA.Tumor.pdf",sep=''),p,width =30,height = 9,units ="cm",dpi=1000)

## step 2 : tumor between genders ##
####################################
pro.gender = pro.list
df = NULL ; p.values = NULL
for(pro in pro.gender){
  t.p = t.e[t.e[,'Cancer']==pro,] ; t.u = t.p
  index = which(pro.gender==pro)
  if (length(unique(t.p$Gender))==1){print(paste(pro,' only one gender !',sep='')) ;pro.gender = pro.gender[-index]}else{
    t.u$Type = 'Tumor'
    e.m = t.u[t.u$Gender=='MALE',gop]
    e.f = t.u[t.u$Gender=='FEMALE',gop]
    tvg = t.test(as.numeric(e.m),as.numeric(e.f))$p.value
    wvg = wilcox.test(as.numeric(e.m),as.numeric(e.f))$p.value
    p.values = rbind(p.values , c(pro, tvg, wvg))
    data = t.u[,c(gop,'Cancer','Type','Gender')]
    data[data[,'Gender']=='MALE','point']=index+0.185 ;data[data[,'Gender']=='FEMALE','point']=index-0.185
    df = rbind(df , data)
  }}
pva = as.data.frame(p.values,stringsAsFactors=F)
pva[,2] = as.numeric(pva[,2]) ; pva[,3] = as.numeric(pva[,3])
colnames(df)[1] = 'exp'
p = ggplot(data=df ,aes(x=Cancer , y=exp, fill = Gender))
p = p + geom_boxplot(alpha = 0,outlier.alpha = F,width=0.74,notch=F,notchwidth=0.1,size=0.3)
p = p + geom_quasirandom(aes(x=point),colour='grey93',size=1.1,shape=21,width=0.1)
p = p + geom_boxplot(colour='grey40',alpha = 0.3,outlier.alpha = F,width=0.74,notch=F,notchwidth=0.1,size=0.3)
p = p + annotate('text',x=pro.gender[which(pva[,2]<0.05)], y=rep(max(as.numeric(df[,1]))*1.1,length(which(pva[,2]<0.05))), label=paste('t=',format(round(pva[,2],3))[which(pva[,2]<0.05)],sep=''),size=2,colour='red')
p = p + annotate('text',x=pro.gender[which(pva[,3]<0.05)], y=rep(max(as.numeric(df[,1]))*1.2,length(which(pva[,3]<0.05))), label=paste('w=',format(round(pva[,3],3))[which(pva[,3]<0.05)],sep=''),size=2,colour='red')
p = p + labs(title=paste(gop,"differential expression between male and female",sep=' '),x="Cancer Type",y ="Expression of protein")
p = p + theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = paste('./result/Protein.',gop,'/',gop,".Protein.RPPA.gender.pdf",sep=''),p,width =30,height = 9,units ="cm",dpi=1000)

## step 4 : survival of high and low with median ##
###################################################
all.sur = t.e
col = c("#E7B800", "#2E9FDF")
all.sur[as.numeric(all.sur[,gop])<=median(as.numeric(all.sur[,gop])),'Group']='L'
all.sur[as.numeric(all.sur[,gop])>median(as.numeric(all.sur[,gop])),'Group']='H'
os = survfit(formula= Surv(all.sur$OS.time, all.sur$OS=='1')~Group ,data=all.sur)
do1 = ggsurvplot(os,pval=TRUE,conf.int=F,palette=col,legend.title="Group: ",legend.labs=c('High','Low'),pval.size=3, pval.coord=c(0,0.1),size=0.6,xlab='Survival Days')
do1 = do1$plot+labs(title=paste(gop,'OS',sep='   '))+theme_bw()+theme(plot.title=element_text(hjust = 0.5),
                                                                      legend.position='top',legend.justification=1)
pfi = survfit(formula= Surv(all.sur$PFI.time, all.sur$PFI=='1')~Group ,data=all.sur)
do2 = ggsurvplot(pfi,pval=TRUE,conf.int=F,palette=col,legend.title="Group: ",legend.labs=c('High','Low'),
                 pval.size=3,pval.coord=c(0,0.1),size=0.6,xlab='Survival Days')
do2 = do2$plot+labs(title=paste(gop,'PFI',sep='   '))+theme_bw()+theme(plot.title=element_text(hjust = 0.5),
                                                                       legend.position='top',legend.justification=1)
draw = arrangeGrob(do1,do2,nrow = 1,ncol = 2)
ggsave(paste('./result/Protein.',gop,'/',gop,'.all.HvsL.Protein.pdf',sep=''),draw,width =32,height = 12,units ='cm',dpi=1200)
x = list()
for (pro in pro.list){
  pro.data = t.e[t.e$Cancer==pro,]
  pro.data[as.numeric(pro.data[,gop])<=median(as.numeric(pro.data[,gop])),'Group']='L'
  pro.data[as.numeric(pro.data[,gop])>median(as.numeric(pro.data[,gop])),'Group']='H'
  # OS
  os = survfit(formula= Surv(pro.data$OS.time, pro.data$OS=='1')~Group ,data=pro.data)
  do1 = ggsurvplot(os,pval=TRUE,conf.int=F,palette=col,legend.title="Group: ",legend.labs=c('High','Low'),
                   pval.size=3,pval.coord=c(0,0.1),size=0.6,xlab='Survival Days')
  do1 = do1$plot+labs(title=paste(gop,'OS',sep='   '))+theme_bw()+theme(plot.title=element_text(hjust = 0.5),
                                                                        legend.position='top',legend.justification=1)
  x[[which(pro.list%in%pro)]] = do1
  # PFI
  if ('FALSE'%in%is.na(pro.data$PFI)){
    pfi = survfit(formula= Surv(pro.data$PFI.time, pro.data$PFI=='1')~Group ,data=pro.data)
    do2 = ggsurvplot(pfi,pval=TRUE,conf.int=F,palette=col,legend.title="Group: ",legend.labs=c('High','Low'),
                     pval.size=3,pval.coord=c(0,0.1),size=0.6,xlab='Survival Days')
    do2 = do2$plot+labs(title=paste(gop,'PFI',sep='   '))+theme_bw()+theme(plot.title=element_text(hjust = 0.5),
                                                                           legend.position='top',legend.justification=1)
  }else{do2 = ggplot()+theme(panel.background = element_blank())}
  x[[which(pro.list%in%pro)+33]] = do2
  draw = arrangeGrob(do1,do2,nrow = 1,ncol = 2)
  ggsave(paste('./result/Protein.',gop,'/',gop,'.',pro,'.HvsL.protein.pdf',sep=''),draw,width =32,height = 12,units ='cm',dpi=1200)
}


