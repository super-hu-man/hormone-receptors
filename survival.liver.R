library(survival)

cli = read.csv('./cli.from.TCGA-CDR/LIHC.filter-NA.txt',header = T,sep='\t',row.names = 1,stringsAsFactors = F)
rna = read.csv('../download.cbioprotal/lihc_tcga_pan_can_atlas_2018/data_expression.txt',header=T,sep='\t',stringsAsFactors=F)
ar = rna[which(rna$Hugo_Symbol=='AR'),-c(1,2)]
colnames(ar) = substr(gsub('[.]','-',colnames(ar)),1,12)
ar = t(ar) ; ar = data.frame(sn=rownames(ar),ar)
ar$group = 'L' ; ar[ar$X1074>median(as.numeric(ar$X1074)),'group'] ='H'
cli2 = data.frame(ar,cli,stringsAsFactors = F)
cli3=cli2[apply(cli2,1,function(x){sum(is.na(x))==0}),]

y = Surv(cli3$OS.time,cli3$OS==1)
y = Surv(cli2$OS.time,cli2$OS==1)
p = survfit(y~group , data=cli2)
d = survdiff(y~group , data=cli2)
