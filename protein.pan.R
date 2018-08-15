##
args=commandArgs( trailingOnly = T)
library('ggplot2')
library(ggbeeswarm)

protein = c('AR','ERALPHA','ERALPHA_pS118','PR')
cli = read.table('combined.cli.txt',header = T,sep='\t',row.names = 1,stringsAsFactors = F)
input = read.table('TCPA.PANCAN19.L4.txt',header = T,stringsAsFactors = F)               # pan19
#input = read.table('TCPA.32Cancer.Tumor.txt',header = T,stringsAsFactors = F)           # 32 cancers
pan = NULL
for(p in protein){
    x=cbind(input[,1:4],input[,p])
    colnames(x)[5]='Exp' ; x$Protein=p
    pan = rbind(pan, x)
}
pan$gender = cli[substr(pan[,1],1,12),4]
pan = pan[is.na(pan$gender)=='FALSE',]
pan.pros = sort(unique(pan$Cancer_Type))

## all pan 19
pos.inf = data.frame(protein,c(-0.3,-0.1,0.1,0.3))
rownames(pos.inf) = pos.inf[,1]
df =NULL
for(p in pan.pros){
    tp = pan[pan$Cancer_Type==p,]
    tp$index = which(pan.pros==p) + pos.inf[tp$Protein,2]
    df=rbind(df,tp)
}
p = ggplot(data=df ,aes(x=Cancer_Type , y=Exp, fill = Protein))
p = p + geom_boxplot(alpha = 0,outlier.alpha = F,width=0.8,notch=F,notchwidth=0.1)
p = p + geom_quasirandom(aes(x=index),colour='grey85',size=1.3,shape=21,width=0.12)
p = p + geom_boxplot(alpha = 0.4,outlier.alpha = F,width=0.8,notch=F,notchwidth=0.1)
p = p + labs(title=paste('AR/ESR1/PGR',"protein expression",sep=' '),x="Cancer Type",y ="Expression of protein")
p = p + theme(plot.title = element_text(hjust = 0.5))
tiff(filename = paste('Genes',"protein.gender.tiff",sep='.'),width =15+length(pan.pros),height = 15,units ="cm",compression="lzw",bg="white",res=1000)
plot(p)
dev.off()

## gender
for (gene in protein){
gender.pro = pan.pros
df =NULL ; p.vas = NULL
df.all = NULL ; p.all = NULL
pan.g =pan[pan$Protein==gene,]
for(p in gender.pro){
    tp = pan.g[pan.g$Cancer_Type==p,] ; ind = which(gender.pro==p) ; tp$index = which(gender.pro==p)
    tp[tp$gender=='female','index'] = tp[tp$gender=='female','index'] - 0.15
    tp[tp$gender=='male','index'] = tp[tp$gender=='male','index'] + 0.15
    if(length(tp[tp$gender=='male','gender']) <2 || length(tp[tp$gender=='female','gender']) <2 ){
        print(p);gender.pro=gender.pro[-ind] } else {
        df=rbind(df,tp)
        p.va = t.test(tp[tp$gender=='male','Exp'],tp[tp$gender=='female','Exp'])$p.value
        p.vas = rbind(p.vas,c(p,p.va))}
    df.all = rbind(df.all, tp)
    p.vall = t.test(tp$Exp, pan.g$Exp)$p.value
    p.all = rbind(p.all,c(p, p.vall))
}
q.all = data.frame(p.all[,1],round(p.adjust(p.all[,2]),3))
q.vall = q.all[q.all[,2]<=0.05,]
p = ggplot(data=df.all ,aes(x=Cancer_Type , y=Exp,fill=Cancer_Type))
p = p + geom_quasirandom(colour='grey90',size=1.8,shape=21,width=0.25)
p = p + geom_boxplot(colour='black',alpha = 0.35,outlier.alpha = F,width=0.36,size=0.5)
p = p + labs(title=paste(gene,"expression of tumor samples",sep=' '),x="Cancer Type",y ="Expression of protein ")
p = p + theme(plot.title = element_text(hjust = 0.5),axis.text.x=element_text(size=8),legend.position='none')
p = p + geom_hline(yintercept = as.numeric(mean(df.all$Exp)),colour='red',linetype='dashed',size=0.5)
p = p + annotate('text',x=q.vall[,1], y=rep(max(df.all$Exp)+1,nrow(q.vall)), label=paste('',format(q.all[,2])[q.all[,2]<0.05],sep=''),size=2,colour='red')
tiff(filename = paste(gene,"protein.exp.tiff",sep='.'),width =15+length(unique(df.all$Cancer_Type)),height = 15,units ="cm",compression="lzw",bg="white",res=1000)
plot(p)
dev.off()
q.vs = data.frame(p.vas[,1],round(p.adjust(p.vas[,2]),3))
qv = q.vs[q.vs[,2]<=0.05,]
p = ggplot(data=df ,aes(x=Cancer_Type , y=Exp, fill = gender))
p = p + geom_boxplot(alpha = 0,outlier.alpha = F,width=0.6,notch=F,notchwidth=0.1)
p = p + geom_quasirandom(aes(x=index),colour='grey85',size=1.8,shape=21,width=0.12)
p = p + geom_boxplot(colour='grey30',alpha = 0.3,outlier.alpha = F,width=0.6,notch=F,notchwidth=0.1)
p = p + annotate('text',x=as.character(qv[,1]), y=rep(max(df$Exp)+1,nrow(qv)), label=paste('',format(qv[,2]),sep=''),size=3,colour='red')
p = p + labs(title=paste(gene,"protein expression between male and female",sep=' '),x="Cancer Type",y ="Expression of protein")
p = p + theme(plot.title = element_text(hjust = 0.5))
tiff(filename = paste(gene,"protein.gender.tiff",sep='.'),width =15+length(unique(df$Cancer_Type)),height = 15,units ="cm",compression="lzw",bg="white",res=1000)
plot(p)
dev.off()
}


## all
for(a in protein){pan[pan$Protein==a,'index']=which(a==protein)-0.2}
pan[pan$gender=='male','index']=pan[pan$gender=='male','index']+0.4
p = ggplot(data=pan ,aes(x=Protein , y=Exp, fill = gender))
p = p + geom_boxplot(alpha = 0,outlier.alpha = F,width=0.8,notch=F,notchwidth=0.1)
p = p + geom_quasirandom(aes(x=index),colour='grey85',size=1.8,shape=21,width=0.2)
p = p + geom_boxplot(colour='grey30',alpha = 0.3,outlier.alpha = F,width=0.8,notch=F,notchwidth=0.1,size=0.7)
tiff(filename = paste('All',"protein.gender.tiff",sep='.'),width =28,height = 18,units ="cm",compression="lzw",bg="white",res=1000)
plot(p)
dev.off()
