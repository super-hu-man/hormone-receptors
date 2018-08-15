##
#setwd('C:/Users/root/Desktop/')
library(survival)
library(gridExtra)
library(survminer)
args=commandArgs( trailingOnly = T)
filelist = c('mRNA.RSEM.Tumors.matchinfor.txt', 'Protein.RPPA.Tumors.matchinfor.txt')

input = filelist[2]
list.p =c('AR','ERALPHA','ERALPHA_pS118','PR')
list.m = c('AR','ESR1','ESR2','PGR')
type.tag = strsplit(input,'.',fixed=TRUE)[[1]][1]

data=read.table(paste('./ucsc/',input,sep=''),header = T,sep='\t',row.names = 1,stringsAsFactors = F)
if (type.tag != 'mRNA'){ list = list.p }else{ list = list.m }

all.pro = unique(data$type)
data = data[sort(rownames(data)),]
data = data[duplicated(data$bcr_patient_barcode)=='FALSE',]

for(p in all.pro ){
    ## exp diff
    pro.data = data[data$type==p,]
    col = c("#E7B800", "#2E9FDF")
    x = list()
    for(i in list){
        index = which(i==list)
        df = pro.data ; df[,i]=gsub('None','-100000',df[,i])
        df[as.numeric(df[,i]) <= median(as.numeric(df[,i])),'group']='Low' ; df[as.numeric(df[,i])>median(as.numeric(df[,i])),'group']='High'
        print(c( p, i, table(df$group) ))
        ## OS survival  ## c("#E7B800", "#2E9FDF")
        name = paste(i,'OS',sep='     ')
        os = survfit(formula= Surv(df$OS.time, df$OS=='1')~group ,data=df)
        do = ggsurvplot(os,pval=TRUE,conf.int=T,palette=col,legend.title="Group: ",legend.labs=c('High','Low'),pval.size=3,pval.coord=c(0,0.1),size=0.6,xlab='Survival Days')
        do = do$plot+labs(title=name)+theme_bw()+theme(plot.title=element_text(hjust = 0.5),legend.position='top',legend.justification=1)
        x[[index]] = do
        ## PFS survival
        if (nrow(df[is.na(df$PFI)=='FALSE'&&is.na(df$PFI.time)=='FALSE',]) < 10 ){ x[[index+4]] = ggplot()+theme(panel.background = element_blank()) }else{
        name = paste(i,'PFS',sep='     ')
        pfs = survfit(formula= Surv(df$PFI.time, df$PFI=='1')~group ,data=df)
        dp = ggsurvplot(pfs,pval=TRUE,conf.int=T,palette=col,legend.title="Group: ",legend.labs=c('High','Low'),pval.size=3,pval.coord=c(0,0.1),size=0.6,xlab='Survival Days')
        dp = dp$plot+labs(title=name)+theme_bw()+theme(plot.title=element_text(hjust = 0.5),legend.position='top',legend.justification=1)
        x[[index+4]] = dp}
    }
    #tiff(filename = paste(p,paste(type.tag,".exp.tiff",sep=''),sep='.'),width =32,height = 48,units ="cm",compression="lzw",bg="white",res=1000)
    draw = marrangeGrob(x,nrow=4,ncol=2,top=p)
    ggsave(paste(paste(p,paste(type.tag,".exp.pdf",sep=''),sep='.')),draw,width =32,height = 48,units ='cm',dpi=1200)
    #dev.off()
}
x = list()
for(i in list){
    index = which(i==list)
    df = data ; df[,i]=gsub('None','-100000',df[,i])
    df[as.numeric(df[,i])<=median(as.numeric(df[,i])),'group']='Low' ; df[as.numeric(df[,i])>median(as.numeric(df[,i])),'group']='High'
    name = paste(i,'OS',sep='     ')
    if(nrow(df)<5) {df=data.frame(PFI=c('0','0'),PFI.time=c(0,0),group=c('High','Low'))}
    os = survfit(formula= Surv(df$OS.time, df$OS=='1')~group ,data=df)
    do = ggsurvplot(os,pval=TRUE,conf.int=T,palette=col,legend.title="Group: ",legend.labs=c('High','Low'),pval.size=3,pval.coord=c(0,0.1),size=0.6,xlab='Survival Days')
    do = do$plot+labs(title=name)+theme_bw()+theme(plot.title=element_text(hjust = 0.5),legend.position='top',legend.justification=1)
    x[[index]] = do
    name = paste(i,'PFS',sep='     ')
    pfs = survfit(formula= Surv(df$PFI.time, df$PFI=='1')~group ,data=df)
    dp = ggsurvplot(pfs,pval=TRUE,conf.int=T,palette=col,legend.title="Group: ",legend.labs=c('High','Low'),pval.size=3,pval.coord=c(0,0.1),size=0.6,xlab='Survival Days')
    dp = dp$plot+labs(title=name)+theme_bw()+theme(plot.title=element_text(hjust = 0.5),legend.position='top',legend.justification=1)
    x[[index+4]] = dp
}
#tiff(filename = paste(p,paste(type.tag,".exp.tiff",sep=''),sep='.'),width =32,height = 48,units ="cm",compression="lzw",bg="white",res=1000)
draw = marrangeGrob(x,nrow=4,ncol=2,top='All Cancers')
ggsave(paste(paste('All',paste(type.tag,".exp.pdf",sep=''),sep='.')),draw,width =32,height = 48,units ='cm',dpi=1200)
#dev.off()

#### gender survival diff ####
data=read.table('./ucsc/clinical.literature.final_use.txt',header = T,sep='\t',row.names = 1,stringsAsFactors = F)
pros = unique(data$type)
for(i in pros ){
    df = data[data$type==i,]
    m = Surv(df$OS.time, df$OS=='1')
    name = paste(i,'OS',sep='     ')
    d1 = ggsurvplot(survfit(formula = m ~ gender ,data=df),pval=TRUE,conf.int=T,palette=col,legend.title="Group: ",pval.size=3,pval.coord=c(0,0.1),size=0.6,xlab='Survival Days')
    d1 = d1$plot+labs(title=name)+theme_bw()+theme(plot.title=element_text(hjust = 0.5),legend.position='top',legend.justification=1)
    n = Surv(df$PFI.time, df$PFI=='1')
    name = paste(i,'PFS',sep='     ')
    if (nrow(df[is.na(df$PFI)=='FALSE'&&is.na(df$PFI.time)=='FALSE',]) < 10 ){ x[[index+4]] = ggplot()+theme(panel.background = element_blank()) }else{
    d2 = ggsurvplot(survfit(formula = n ~ gender ,data=df),pval=TRUE,conf.int=T,palette=col,legend.title="Group: ",pval.size=3,pval.coord=c(0,0.1),size=0.6,xlab='Survival Days')
    d2 = d2$plot+labs(title=name)+theme_bw()+theme(plot.title=element_text(hjust = 0.5),legend.position='top',legend.justification=1) }
    tiff(filename = paste(i,'gender.tiff',sep='.'),width =48,height = 24,units ="cm",compression="lzw",bg="white",res=1000)
    grid.arrange(d1,d2,nrow=1,ncol=2)
    dev.off()
}
