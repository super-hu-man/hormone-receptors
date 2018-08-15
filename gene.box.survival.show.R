args = commandArgs( trailingOnly = T)
input = './ucsc/mRNA.RSEM.combined.txt'
gop = args
if (file.exists('./result')){print('dir result exist !')}else{dir.create('./result')}
if (file.exists(paste('./result/mRNA.',gop,sep=''))){
  print('dir gene exist !')}else{dir.create(paste('./result/mRNA.',gop,sep=''))}

library('ggplot2')
library(ggbeeswarm)
library(survival)
library(gridExtra)
library(survminer)

cli = read.table('./ucsc/clinical.literature.final_use.txt',header = T,sep='\t',row.names = 1,stringsAsFactors = F)
pro.list = unique(cli$type)
table = read.table(input,header = T,row.names=1,sep='\t',stringsAsFactors = F)
table$Cancer = cli[table$bcr_patient_barcode,'type']
table$Gender = cli[table$bcr_patient_barcode,'gender']
table$OS = cli[table$bcr_patient_barcode,'OS']
table$OS.time = cli[table$bcr_patient_barcode,'OS.time']
table$PFI = cli[table$bcr_patient_barcode,'PFI']
table$PFI.time = cli[table$bcr_patient_barcode,'PFI.time']
table = table[is.na(table$Cancer)=='FALSE',]
t.e = table[substr(rownames(table),14,14)==0,]
n.e = table[substr(rownames(table),14,14)!=0,]

## step 1 : draw T vs N all samples rsem ##
###########################################
df = NULL ; p.values = NULL
df.all = NULL ;p.vall = NULL
for(pro in pro.list){
  t.p = t.e[t.e[,'Cancer']==pro,]
  n.p = n.e[n.e[,'Cancer']==pro,]
  index = which(pro.list==pro)
  ## tumor vs normal
  if(nrow(n.p)<= 10){
    t.p$type = 'Tumor' ; t.p$point = index+0.185
    p.vall = rbind(p.vall,c(pro,1,1))
    data = rbind(t.p[,c(gop,'Cancer','type','point')],c(0,pro,'Normal',index-0.185))
    df.all = rbind(df.all , data)
  }else{
    t.p$type = 'Tumor' ; t.p$point = index+0.185
    n.p$type = 'Normal'; n.p$point = index-0.185
    wil = wilcox.test(as.numeric(t.p[,gop]) , as.numeric(n.p[,gop]))$p.value
    tup = t.test(as.numeric(t.p[,gop]) , as.numeric(n.p[,gop]))$p.value
    p.vall = rbind(p.vall,c(pro,tup,wil))
    data = rbind(t.p[,c(gop,'Cancer','type','point')],n.p[,c(gop,'Cancer','type','point')])
    df.all = rbind(df.all , data)
  }}
colnames(df.all)[1] = 'exp' ; colnames(df.all)[3] = 'Type'
df.all[,4]=as.numeric(df.all[,4]) ; df.all[,1]=as.numeric(df.all[,1])
pva = as.data.frame(p.vall,stringsAsFactors=F)
pva[,2] = as.numeric(pva[,2]) ; pva[,3] = as.numeric(pva[,3])
p = ggplot(data=df.all ,aes(x=Cancer , y=exp,fill = Type))
p = p + geom_boxplot(alpha = 0,outlier.alpha = F,width=0.75,notch=F,notchwidth=0.1,size=0.3)
p = p + geom_quasirandom(aes(x=point),colour='grey93',size=1.1,shape=21,width=0.1)
p = p + geom_boxplot(colour='grey40',alpha = 0.3,outlier.alpha = F,width=0.75,notch=F,notchwidth=0.1,size=0.3)
p = p + annotate('text',x=pro.list[which(pva[,2]<0.05)], y=rep(max(as.numeric(df.all[,1]))*1.05,length(which(pva[,2]<0.05))), label=paste('t=',format(round(pva[,2],3))[which(pva[,2]<0.05)],sep=''),size=2,colour='red')
p = p + annotate('text',x=pro.list[which(pva[,3]<0.05)], y=rep(max(as.numeric(df.all[,1]))*1.1,length(which(pva[,3]<0.05))), label=paste('w=',format(round(pva[,3],3))[which(pva[,3]<0.05)],sep=''),size=2,colour='red')
p = p + labs(title=paste(gop,"differential expression of tumor and normal",sep=' '),x="Cancer Type",y ="Expression of mRNA (log2 RSEM+1 )")
p = p + theme(plot.title = element_text(hjust = 0.5), axis.text.x=element_text(size=7))
#tiff(filename = paste(gop,"mRNA.matched.q.tiff",sep='.'),width =30,height = 12,units ="cm",compression="lzw",bg="white",res=1000)
ggsave(filename = paste('./result/mRNA.',gop,'/',gop,".mRNA.RSEM.TvsN.pdf",sep=''),p,width =30,height = 9,units ="cm",dpi=1000)

## step 2 : matched tumoe vs normal ##
######################################
pro.match = pro.list
df = NULL ; p.values = NULL
df.all = NULL ;p.vall = NULL
for(pro in pro.match){
  t.p = t.e[t.e[,'Cancer']==pro,]
  n.p = n.e[n.e[,'Cancer']==pro,]
  index = which(pro.match==pro)
  if (nrow(n.p)<=10){print(paste(pro,' sample is not enough !',sep='')) ; pro.match = pro.match[-index]}else{
    t.p=t.p[sort(rownames(t.p)),] ; t.p = t.p[duplicated(t.p[,1])=='FALSE',]
    n.p=n.p[sort(rownames(n.p)),] ; n.p = n.p[duplicated(n.p[,1])=='FALSE',]
    share = intersect(t.p[,1],n.p[,1])
    t.u = t.p[t.p[,1]%in%share,] ; t.u$type = 'Tumor' ; t.u$point = index+0.125
    n.u = n.p[n.p[,1]%in%share,] ; n.u$type = 'Normal'; n.u$point = index-0.125
    wpa = wilcox.test(as.numeric(t.u[,gop]) , as.numeric(n.u[,gop]),paired = T)$p.value
    tpa = t.test(as.numeric(t.u[,gop]) , as.numeric(n.u[,gop]), paired=T)$p.value
    p.values = rbind(p.values , c(pro, tpa, wpa))
    data = rbind(t.u[,c(gop,'Cancer','type','point')],n.u[,c(gop,'Cancer','type','point')])
    df = rbind(df , data)
  } }
colnames(df)[1] = 'exp' ; colnames(df)[3] = 'Type'
df[,4]=as.numeric(df[,4]) ; df[,1]=as.numeric(df[,1])
pva = as.data.frame(p.values,stringsAsFactors=F)
pva[,2] = as.numeric(pva[,2]) ; pva[,3] = as.numeric(pva[,3])
p = ggplot(data=df ,aes(x=Cancer , y=exp,fill = Type))
p = p + geom_boxplot(alpha = 0,outlier.alpha = F,width=0.5,notch=F,notchwidth=0.1,size=0.3)
p = p + geom_quasirandom(aes(x=point),colour='grey93',size=1.1,shape=21,width=0.1)
p = p + geom_boxplot(colour='grey40',alpha = 0.3,outlier.alpha = F,width=0.5,notch=F,notchwidth=0.1,size=0.3)
p = p + annotate('text',x=pro.match[which(pva[,2]<0.05)], y=rep(max(as.numeric(df[,1]))*1.05,length(which(pva[,2]<0.05))), label=paste('t=',format(round(pva[,2],3))[which(pva[,2]<0.05)],sep=''),size=2,colour='red')
p = p + annotate('text',x=pro.match[which(pva[,3]<0.05)], y=rep(max(as.numeric(df[,1]))*1.1,length(which(pva[,3]<0.05))), label=paste('w=',format(round(pva[,3],3))[which(pva[,3]<0.05)],sep=''),size=2,colour='red')
p = p + labs(title=paste(gop,"differential expression of matched samples",sep=' '),x="Cancer Type",y ="Expression of mRNA (log2 RSEM+1 )")
p = p + theme(plot.title = element_text(hjust = 0.5), axis.text.x=element_text(size=7))
#tiff(filename = paste(gop,"mRNA.matched.q.tiff",sep='.'),width =30,height = 12,units ="cm",compression="lzw",bg="white",res=1000)
ggsave(filename = paste('./result/mRNA.',gop,'/',gop,".mRNA.RSEM.matched.pdf",sep=''),p,width =30,height = 9,units ="cm",dpi=1000)

## step 3 : tumor between genders ##
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
p = p + annotate('text',x=pro.gender[which(pva[,2]<0.05)], y=rep(max(as.numeric(df[,1]))*1.05,length(which(pva[,2]<0.05))), label=paste('t=',format(round(pva[,2],3))[which(pva[,2]<0.05)],sep=''),size=2,colour='red')
p = p + annotate('text',x=pro.gender[which(pva[,3]<0.05)], y=rep(max(as.numeric(df[,1]))*1.1,length(which(pva[,3]<0.05))), label=paste('w=',format(round(pva[,3],3))[which(pva[,3]<0.05)],sep=''),size=2,colour='red')
p = p + labs(title=paste(args[1],"differential expression between male and female",sep=' '),x="Cancer Type",y ="Expression of mRNA(RSEM log2 )")
p = p + theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = paste('./result/mRNA.',gop,'/',gop,".mRNA.RSEM.gender.pdf",sep=''),p,width =30,height = 9,units ="cm",dpi=1000)

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
ggsave(paste('./result/mRNA.',gop,'/',gop,'.all.HvsL.mRNA.pdf',sep=''),draw,width =32,height = 12,units ='cm',dpi=1200)
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
  ggsave(paste('./result/mRNA.',gop,'/',gop,'.',pro,'.HvsL.mRNA.pdf',sep=''),draw,width =32,height = 12,units ='cm',dpi=1200)
}
#draw = marrangeGrob(x,nrow=33,ncol=2,top='All Cancers Projects')
#ggsave(paste('./result/',gop,'.project.HvsL.mRNA.tiff',sep=''),draw,width =32,height = 400,limitsize = F,units ='cm',dpi=1200)
