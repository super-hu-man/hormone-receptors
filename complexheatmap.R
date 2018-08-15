library(ComplexHeatmap)
library(circlize)

files = list.files('./result/Mutation','matrix',full.names = T)
cli = read.csv('./ucsc/clinical.literature.final_use.txt',header = T,sep='\t',stringsAsFactors = F,row.names = 1)

for(f in files){
  name = strsplit(strsplit(f,'/',fixed=T)[[1]][4],'.',fixed=T)[[1]][1]
  print(name)
  mat = read.csv(f,header = T,row.names = 1,sep='\t',stringsAsFactors = F)
  colnames(mat)  = gsub('[.]','-',colnames(mat))
  mat = mat[,c(1:3,which(mat['Group',]=='Low'),which(mat['Group',]=='High'))]
  # group make
  group = mat[1:2,4:ncol(mat)]
  group['Gender',] = cli[substr(colnames(group),1,12),'gender']
  group['Age',] = cli[substr(colnames(group),1,12),'age_at_initial_pathologic_diagnosis']
  group = as.data.frame(t(group),stringsAsFactors = F) ; 
  group$Age = as.numeric(group$Age) ; group$Exp = as.numeric(log(as.numeric(group$Exp)+1,2))
  # 
  mat = mat[-c(1:2),]
  mat = mat[mat$q_value<=0.05,]
  if(nrow(mat)==0){print('no statistical significance !');next}
  mat = mat[order(mat$q_value),]
  # right bar make
  barpt = mat$LvsH_number ; bargp=as.data.frame(sapply(barpt,strsplit,split=' | ',fixed=T),stringsAsFactors=F) ; 
  rownames(bargp)=c('Low','High') ; colnames(bargp) = rownames(mat)
  mat = mat[,4:ncol(mat)]

#prepare
alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#EEE9E9", col = NA))},
  Missense_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#0496FF", col = NA))},
  In_Frame_Del = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#FFB4A2", col = NA))},
  Nonsense_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#8AC926", col = NA))},
  Splice_Site = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#BA55D3", col = NA))},
  In_Frame_Ins = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#008B00", col = NA))},
  Frame_Shift_Del = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#FFBC42", col = NA))},
  Nonstop_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#00FF7F", col = NA))},
  Frame_Shift_Ins = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#F6511D", col = NA))},
  Translation_Start_Site = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#1C1C1C", col = NA))}
)
col = c('Missense_Mutation'='#0496FF' , 'In_Frame_Del'='#FFB4A2'    , 'Nonsense_Mutation'='#8AC926' ,
        'Splice_Site'='#BA55D3'       , 'In_Frame_Ins'='#008B00'    , 'Frame_Shift_Del'='#FFBC42',
        'Nonstop_Mutation'='#00FF7F'  , 'Frame_Shift_Ins'='#F6511D' , 'Translation_Start_Site'='#1C1C1C')

# top groups annotations
tgp = HeatmapAnnotation(Exp = anno_barplot(group$Exp,gp=gpar(fill='#FFD400',col='#FFD400' ),border=F,axis=T,axis_side='right'),
                        Group = group$Group , Gender = group$Gender , Age = group$Age ,
                        gap=unit(0.8,'mm'),annotation_name_offset = unit(3,'mm'), 
                        name=colnames(group), show_annotation_name=TRUE,annotation_name_side='left',
                        annotation_legend_param = list(legend_height=unit(0.1,'cm'),legend_direction='horizontal'),
                        col=list( Group=c("Low"="blue", "High"="red"),
                                  Age=colorRamp2(c(min(group$Age,na.rm=T),max(group$Age,na.rm=T)), c("white","#FFA500")),
                                  Gender=c('FEMALE'='#FF6161','MALE'='#1592A0')) ,
                        annotation_height = c(1.8, 1, 1, 1))
# right annotation : number of each gene
rgp = as.matrix(t(apply(bargp,2,as.numeric)))
rgp[,1] = as.numeric(rgp[,1])/sum(group$Group=='Low')*100 ; rgp[,2] = as.numeric(rgp[,2])/sum(group$Group=='High')*100
anno_bar_beside = function(index) {
  n = length(index)
  data_scale = range(rowSums(rgp, na.rm = TRUE), na.rm = TRUE)
  data_scale = data_scale + c(-0.05, 0.05)*(data_scale[2] - data_scale[1])
  pushViewport(viewport(xscale = c(0,max(rgp)*1.5), yscale = c(0.5, n + 0.5)))
  width = rgp[index, 1] ; x_coor = width/2
  grid.rect(x = x_coor, y = n - seq_along(index) + 1.15, width = abs(width), height = 0.3, default.units = "native",gp=gpar(fill='blue',col=NA))
  width = rgp[index, 2] ; x_coor = width/2
  grid.rect(x = x_coor, y = n - seq_along(index) + 0.85, width = abs(width), height = 0.3, default.units = "native",gp=gpar(fill='red',col=NA))
  at = c(0,grid.pretty(data_scale)) ; label = at
  grid.xaxis(main = FALSE,gp=gpar(fontsize=7), at=at, label=label/100)
  grid.lines(x=unit(0,'npc')) ; grid.lines(x=unit(0,'npc'),y=unit(c(0,0.5),'mm'),gp=gpar(col='white'))
  upViewport()
}
gnr = rowAnnotation(rgp=anno_bar_beside,width = unit(3, "cm"))
# left gene name
gname =rownames(mat)
gnl = rowAnnotation(gn=row_anno_text(gname,just='right',offset=unit(20,'mm')),width = unit(2, "cm"))

# center landscape
pic = oncoPrint(mat, alter_fun=alter_fun, col=col,show_row_barplot = F,row_title_side='left',show_pct=F,show_row_names = F,
          top_annotation = tgp,top_annotation_height = unit(3.5,'cm'),
          get_type = function(x){strsplit(x, ";")[[1]]}, 
          column_order = NULL, row_order = NULL,heatmap_legend_param = list(title='Alterations',nrow=2,height=unit(3,'cm')))

tiff(paste('./result/Mutation/',name,'.tiff'),width =6+0.12*ncol(mat),height = 8+0.6*nrow(mat),units = 'cm',compression = 'lzw',res=800)
draw(gnl+pic+gnr, heatmap_legend_side='top',annotation_legend_side='top')
dev.off()
}
