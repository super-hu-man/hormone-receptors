library(ComplexHeatmap)

#mc = read.csv('./mutation/mc3.v0.2.8.PUBLIC.xena.txt',header = T,sep='\t',stringsAsFactors = F)
#cli = read.csv('./files/clinical.literature.final_use.txt',header = T,sep='\t',stringsAsFactors = F)
#pros = sort(unique(cli$type))

gene = 'AR'
cutoff = 30
list = c("Missense_Mutation","In_Frame_Del","Nonsense_Mutation","Splice_Site","In_Frame_Ins",
         "Frame_Shift_Del","Nonstop_Mutation","Frame_Shift_Ins","Translation_Start_Site")
files = list.dirs('../download.cbioprotal',full.names = T)[2:34] 

for(fd in files){
  # sort samples 
  name = toupper(strsplit(strsplit(fd,'/',fixed=T)[[1]][3],'_',fixed=T)[[1]][1])
  print(name)
  dat = read.csv(paste(fd,'/data_expression_merged.txt',sep=''),header = T,sep='\t',stringsAsFactors = F)
  datu = dat[,order(dat[which(dat$Hugo_Symbol==gene),3:ncol(dat)])+2]
  colnames(datu) = gsub('[.]','-',colnames(datu))
  num = ncol(datu) * cutoff/100
  Lname = colnames(datu)[1:ceiling(num)]
  Hname = colnames(datu)[(ncol(datu)-floor(num)):ncol(datu)]
### mutations ###
#################
  mut = read.csv(paste(fd,'data_mutations_extended.txt',sep='/'),header=T,sep='\t',stringsAsFactors=F)
  mut = mut[mut$Tumor_Sample_Barcode %in% colnames(datu),]
  # filter 1000 exon mutation
  exon = mut[mut$EXON!='.',]
  info = as.data.frame(table(exon$Tumor_Sample_Barcode),stringsAsFactors=F)
  mut1 = mut[mut$Tumor_Sample_Barcode %in% info[info$Freq<=1000,'Var1'],]
  # prepare matrix  ^as.data.frame(matrix(numeric(0),ncol=4))^
  outp = datu[which(dat$Hugo_Symbol==gene),colnames(datu) %in% mut1$Tumor_Sample_Barcode]
  rownames(outp) = 'Exp' ; outp['Group',] = 'M'
  outp['Group',colnames(outp) %in% Lname] = 'Low'
  outp['Group',colnames(outp) %in% Hname] = 'High'
  # filter less than 5% non-silents mutation samples
  mut2 = mut1[mut1$Variant_Classification %in% list,]
  tmp = data.frame(mut2$Hugo_Symbol,mut2$Tumor_Sample_Barcode)
  tmp = unique(tmp)
  info = as.data.frame(table(tmp$mut2.Hugo_Symbol),stringsAsFactors=F)
  numbersample = length(unique(mut2$Tumor_Sample_Barcode))
  mut2 = mut2[mut2$Hugo_Symbol %in% info[info$Freq>=0.05*numbersample,'Var1'],]
  print(sum(info$Freq>=0.05*numbersample))
  # write mutation matrix
  out = outp
  for(l in 1:nrow(mut2)){
    line = mut2[l,]
    if(is.null(out[line$Hugo_Symbol,line$Tumor_Sample_Barcode])||is.na(out[line$Hugo_Symbol,line$Tumor_Sample_Barcode])){
      out[line$Hugo_Symbol,line$Tumor_Sample_Barcode]=line$Variant_Classification  } else {
      out[line$Hugo_Symbol,line$Tumor_Sample_Barcode]=paste(out[line$Hugo_Symbol,line$Tumor_Sample_Barcode],line$Variant_Classification,sep=';')
    }
  }
  # statistics information
  out['Group','LvsH_number'] = paste(sum(Lname%in%colnames(out)),sum(Hname%in%colnames(out)),sep=' | ')
  for(g in rownames(out)[3:nrow(out)]){
    Lg = out[g, Lname[Lname%in%colnames(out)]]
    Hg = out[g, Hname[Hname%in%colnames(out)]]
    X2 = data.frame(c(sum(is.na(Lg)),length(Lname[Lname%in%colnames(out)])-sum(is.na(Lg))) ,
                    c(sum(is.na(Hg)),length(Hname[Hname%in%colnames(out)])-sum(is.na(Hg))))
    f.test =  fisher.test(X2)
    out[g,'LvsH_number'] = paste(sum(is.na(Lg)=='FALSE'),sum(is.na(Hg)=='FALSE'),sep=' | ')
    out[g,'p_value'] = f.test$p.value
  }
  out$q_value[3:nrow(out)] = p.adjust(out$p_value[3:nrow(out)],method = 'BH',n=length(out$p_value[3:nrow(out)]))
  #out = out[,c('LvsH_number','p_value','q_value',Lname[Lname%in%colnames(out)],Hname[Hname%in%colnames(out)])]
  out = out[,c('LvsH_number', 'p_value','q_value', colnames(outp))]
  out = data.frame(Gene=rownames(out),out,stringsAsFactors = F)
  write.table(out,paste('./result/Mutation/',name,'.matrix.L-H.txt',sep=''),row.names=F,col.names=T,sep='\t',quote=F)
}

source('complexheatmap.R')
