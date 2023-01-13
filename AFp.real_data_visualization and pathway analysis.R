#-------------------------------------------------visualization

rm(list=ls())
library(ggplot2)
data.dir="data/"
out.dir=paste0(data.dir,"/out/")
setwd(data.dir)
load("LungData.all.Rdata")
Data.Gene<-LungData.all$gene
Data.Clinical<-LungData.all$clinic
library(dplyr)
library(gplots)
#"fev6prd1"
selected.variable<-c("fev1pd1a","fvcprd1","ratiopre","WBCDIFF1","WBCDIFF4")
#confounder.variable<-c("pkyrs","Cigarette_smoking","BMI","AGE","GENDER")
confounder.variable<-c("BMI","AGE","GENDER")
Data.confounder<-dplyr::select(Data.Clinical,confounder.variable)
Y<-dplyr::select(Data.Clinical,selected.variable)

###############filter genes

dim(Data.Gene)
index1<-which(apply(Data.confounder,1,function(x){sum(is.na(x))})!=0)
index2<-which(apply(Y,1,function(x){sum(is.na(x))})!=0)
index_delete<-unique(c(index1,index2))
Data.Gene<-Data.Gene[-index_delete,]
Data.Clinical<-Data.Clinical[-index_delete,]
Data.confounder<-Data.confounder[-index_delete,]
Y<-Y[-index_delete,]

X<-Data.Gene
Z<-Data.confounder
name1<-colnames(X)

#--------------------------------------read summarized result

load("lung_disease_result_all.Rdata")
load("Varibility_index.Rdata")
load("comembership_Caleb.AFp.Rdata")


index<-which((final.res$AFp$AFp.pvalue*15966)<0.05)
Coef<-final.res$input.stat$coef[index,]
Pvalue<-final.res$input.stat$pvalue[index,]
variability<-variability.index$AFp[index,]

weight.matrix<-final.res$AFp$W[index,]
X1<-X[,index]
Distance.matrix<-Distance.matrix.AFp[index,index]
name1<-name1[index]
colnames(Distance.matrix)<-rownames(Distance.matrix)<-rep("",dim(Distance.matrix)[1])
Distance.matrix<-1-Distance.matrix
rm(Distance.matrix.AFp)
#
library(tightClust)
##res<-tight.clust(Distance.matrix, 7, 30, random.seed=12315)
##save(res,file="tight_clustering_K=7.Rdata")
load("tight_clustering_K=7.Rdata")
index1<-which(res$cluster==1)
index2<-which(res$cluster==2)
index3<-which(res$cluster==3)
index4<-which(res$cluster==4)
index5<-which(res$cluster==6)
index6<-which(res$cluster==5)
index7<-which(res$cluster==7)
index<-c(index1,index2,index3,index4,index5,index6,index7)
colors<-c(rep("1",length(index1)),rep("2",length(index2)),rep("3",length(index3)),rep("4",length(index4)),
          rep("5",length(index5)),rep("6",length(index6)),rep("7",length(index7)))


pdf(paste0(out.dir,"heatmap_tightclust_K=7.pdf"),height = 7,width = 7)
heatmap.2(1-Distance.matrix[index,index],trace = "none",col="greenred",Rowv = NULL,
          Colv=NULL,RowSideColors = colors,ColSideColors = colors,
          density.info = "none",key.title = NA,key.xlab = "Co-membership",key.ylab = NA,
          keysize=1, key.par = list(cex=1))



legend("top",legend = paste("C",1:7,sep=""),
       col = as.character(1:7), 
       cex=1,lty=1,lwd=3,horiz = T
)
dev.off()

data.heatmap<-(-log10(Pvalue)*sign(Coef))[index,]
data.heatmap[data.heatmap>10]<-10
data.heatmap[data.heatmap<(-10)]<-(-10)
rownames(data.heatmap)<-rep("",nrow(data.heatmap))
colnames(data.heatmap)<-c("FEV1", "FVC", "ratiopre", "WBCDIFF1","WBCDIFF4")

pdf(paste0(out.dir,"heatmap_Pvalue_tightclust_K=7.pdf"),width = 14,height = 10)
heatmap.2(data.heatmap,trace = "none",col="greenred",Rowv = NULL,
          Colv=NULL,RowSideColors = colors,
          density.info = "none",key.title = NA,key.xlab = "-log10(pvalue)",key.ylab = NA,
          keysize=1.5, key.par = list(cex=1.5),margins=c(10,10),cexCol=2)
legend("top",legend = paste("C",1:7,sep=""),
       col = as.character(1:7), 
       cex=1.3,lty=1,lwd=3,horiz = T
)
dev.off()

#--------------------------

bt_list=list()
index1<-1:length(which(res$cluster==1))
library(ComplexHeatmap)

data1<-X1[,index1]
Coef1<-Coef[index1,]
variability1<-variability[index1,]
weight.matrix1<-weight.matrix[index1,]
#---------------------
data1<-t(data1)
order1<-order(apply(variability1,1,sum),decreasing = F)
data1<-data1[order1,]
variability1<-variability1[order1,]

Coef1<-Coef1[order1,]
weight.matrix1<-weight.matrix1[order1,]
Sign_weight1<-weight.matrix1*sign(Coef1)

data1<-t(scale(t(data1)))

Y$fvcprd1<-(Y$fvcprd1-mean(Y$fvcprd1))/sd(Y$fvcprd1)
Y$WBCDIFF1<-(Y$WBCDIFF1-mean(Y$WBCDIFF1))/sd(Y$WBCDIFF1)


Y$WBCDIFF4<-(Y$WBCDIFF4-mean(Y$WBCDIFF4))/sd(Y$WBCDIFF4)
col_fun1 = circlize::colorRamp2(c(0, 0.5, 1), c("green", "black", "red"),space="RGB")
column_ha = HeatmapAnnotation(fev1pd1a=Y$fev1pd1a,fvcprd1=Y$fvcprd1,
                              ratiopre=Y$ratiopre,WBCDIFF1=Y$WBCDIFF1,
                              WBCDIFF4=Y$WBCDIFF4,
                              annotation_name_side = "left",
                              show_legend = F,
                              height=unit(1, "cm"),
                              simple_anno_size = unit(0.2, "cm"),
                              annotation_name_gp= gpar(fontsize = 8),
                              col = list(fev1pd1a = col_fun1,fvcprd1=circlize::colorRamp2(c(min(Y$fvcprd1), 0, max(Y$WBCDIFF1)), c("green", "black", "red"),space="RGB"),
                                         ratiopre=col_fun1,WBCDIFF1=circlize::colorRamp2(c(min(Y$WBCDIFF1), 0, max(Y$WBCDIFF1)), c("green", "black", "red"),space="RGB"),
                                         WBCDIFF4=circlize::colorRamp2(c(min(Y$WBCDIFF4), 0, max(Y$WBCDIFF4)), c("green", "black", "red"),space="RGB")))

library(circlize)
library(tidyr)
row_ha = as.dendrogram(hclust(dist(data1)))
col_fun = circlize::colorRamp2(c(min(data1), 0, max(data1)), c("green", "black", "red"),space="RGB")
p1<-Heatmap(data1, top_annotation = column_ha,
            cluster_rows=row_ha,cluster_columns=T,show_row_names=F,show_column_names=F,
            col=col_fun,name="Expression",show_row_dend = F,show_column_dend = F,show_heatmap_legend = F,width=unit(12, "cm"),
            row_title_side="left")
FEV1=sapply(1:dim(data1)[1],function(x) cor(data1[x,],Y$fev1pd1a,method = "spearman"))
FVC=sapply(1:dim(data1)[1],function(x) cor(data1[x,],Y$fvcprd1,method = "spearman"))
ratiopre=sapply(1:dim(data1)[1],function(x) cor(data1[x,],Y$ratiopre,method = "spearman"))
WBCDIFF1=sapply(1:dim(data1)[1],function(x) cor(data1[x,],Y$WBCDIFF1,method = "spearman"))
WBCDIFF4=sapply(1:dim(data1)[1],function(x) cor(data1[x,],Y$WBCDIFF4,method = "spearman"))

box_df=data.frame(FEV1=FEV1,FVC=FVC,ratiopre=ratiopre,WBCDIFF1=WBCDIFF1,WBCDIFF4=WBCDIFF4)
box_df2<-pivot_longer(box_df,everything(),names_to="Phenotypes",values_to="Correlations")
p4<-ggplot(box_df2, aes(x=Phenotypes,y=Correlations)) + xlab("Phenotypes")+ylab("Correlations")+theme_bw()+
  geom_boxplot(fill="red",lwd=1)+theme(legend.position = "right",
                                         plot.title = element_text(size = rel(2), hjust = 0.5,face="bold"),
                                         axis.title = element_text(size = rel(2),face="bold"),axis.text = element_text(size=20,face = "bold"),
                                         legend.text = element_text(size=15),panel.border = element_rect(colour = "black",size=2))+
  ggtitle("Cluster 1")+scale_y_continuous(limits = c(-.6, .6))+geom_hline(yintercept=0, linetype="dashed", color = "red")

bt_list[[1]]=p4


col_fun1 = circlize::colorRamp2(c(0, 0.5, 1), c("green", "black", "red"),space="RGB")
colnames(variability1)<-c("FEV1","FVC","Ratiopre","WBCDIFF1","WBCDIFF4")
p2<-Heatmap(variability1,cluster_rows=row_ha,cluster_columns=F,show_row_names=F,show_column_names=T,show_row_dend = F,
            col=col_fun1,name="Variability index",show_heatmap_legend = T,width=unit(4, "cm"),column_names_side="top",
            column_names_gp = grid::gpar(fontsize = 12)
)

col_fun2 = structure(c("yellow","black","cyan2"), names = c("-1", "0", "1"))

colnames(Sign_weight1)<-c("FEV1","FVC","Ratiopre","WBCDIFF1","WBCDIFF4")

p3<-Heatmap(Sign_weight1,cluster_rows=row_ha,cluster_columns=F,show_row_names=F,show_column_names=T,show_row_dend = F,
            col=col_fun2,name="Signed AW-Fisher weights",show_heatmap_legend = F,width=unit(4, "cm"),column_names_side="top",
            column_names_gp = grid::gpar(fontsize = 12))

ht_list1<-p3
vt_list1<-p2
module1<-draw(ht_list1, ht_gap = unit(0.5, "cm"))
#----------------------------------------module2
#----------------------------------------
a1<-which(res$cluster==1)
a2<-which(res$cluster==2)
a3<-which(res$cluster==3)
a4<-which(res$cluster==4)
a5<-which(res$cluster==5)

a6<-which(res$cluster==6)
a7<-which(res$cluster==7)

index1<-(length(a1)+1):(length(a1)+length(a2))
library(ComplexHeatmap)

data1<-X1[,index1]
Coef1<-Coef[index1,]
variability1<-variability[index1,]
weight.matrix1<-weight.matrix[index1,]
#---------------------
data1<-t(data1)
order1<-order(apply(variability1,1,sum),decreasing = F)
data1<-data1[order1,]
variability1<-variability1[order1,]

Coef1<-Coef1[order1,]
weight.matrix1<-weight.matrix1[order1,]
Sign_weight1<-weight.matrix1*sign(Coef1)

data1<-t(scale(t(data1)))

Y$fvcprd1<-(Y$fvcprd1-mean(Y$fvcprd1))/sd(Y$fvcprd1)
Y$WBCDIFF1<-(Y$WBCDIFF1-mean(Y$WBCDIFF1))/sd(Y$WBCDIFF1)


Y$WBCDIFF4<-(Y$WBCDIFF4-mean(Y$WBCDIFF4))/sd(Y$WBCDIFF4)
col_fun1 = circlize::colorRamp2(c(0, 0.5, 1), c("green", "black", "red"),space="RGB")
column_ha = HeatmapAnnotation(fev1pd1a=Y$fev1pd1a,fvcprd1=Y$fvcprd1,
                              ratiopre=Y$ratiopre,WBCDIFF1=Y$WBCDIFF1,
                              WBCDIFF4=Y$WBCDIFF4,
                              annotation_name_side = "left",
                              show_legend = F,
                              height=unit(1, "cm"),
                              simple_anno_size = unit(0.2, "cm"),
                              annotation_name_gp= gpar(fontsize = 8),
                              col = list(fev1pd1a = col_fun1,fvcprd1=circlize::colorRamp2(c(min(Y$fvcprd1), 0, max(Y$WBCDIFF1)), c("green", "black", "red"),space="RGB"),
                                         ratiopre=col_fun1,WBCDIFF1=circlize::colorRamp2(c(min(Y$WBCDIFF1), 0, max(Y$WBCDIFF1)), c("green", "black", "red"),space="RGB"),
                                         WBCDIFF4=circlize::colorRamp2(c(min(Y$WBCDIFF4), 0, max(Y$WBCDIFF4)), c("green", "black", "red"),space="RGB")))

library(circlize)
row_ha = as.dendrogram(hclust(dist(data1)))
col_fun = circlize::colorRamp2(c(min(data1), 0, max(data1)), c("green", "black", "red"),space="RGB")
p1<-Heatmap(data1, top_annotation = column_ha,
            cluster_rows=row_ha,cluster_columns=T,show_row_names=F,show_column_names=F,
            col=col_fun,name="Expression",show_row_dend = F,show_column_dend = F,show_heatmap_legend = F,width=unit(12, "cm"),
            row_title_side="left")
FEV1=sapply(1:dim(data1)[1],function(x) cor(data1[x,],Y$fev1pd1a,method = "spearman"))
FVC=sapply(1:dim(data1)[1],function(x) cor(data1[x,],Y$fvcprd1,method = "spearman"))
ratiopre=sapply(1:dim(data1)[1],function(x) cor(data1[x,],Y$ratiopre,method = "spearman"))
WBCDIFF1=sapply(1:dim(data1)[1],function(x) cor(data1[x,],Y$WBCDIFF1,method = "spearman"))
WBCDIFF4=sapply(1:dim(data1)[1],function(x) cor(data1[x,],Y$WBCDIFF4,method = "spearman"))

box_df=data.frame(FEV1=FEV1,FVC=FVC,ratiopre=ratiopre,WBCDIFF1=WBCDIFF1,WBCDIFF4=WBCDIFF4)
box_df2<-pivot_longer(box_df,everything(),names_to="Phenotypes",values_to="Correlations")
p4<-ggplot(box_df2, aes(x=Phenotypes,y=Correlations)) + xlab("Phenotypes")+ylab("Correlations")+theme_bw()+
  geom_boxplot(fill="red",lwd=1)+theme(legend.position = "right",
                                       plot.title = element_text(size = rel(2), hjust = 0.5,face="bold"),
                                       axis.title = element_text(size = rel(2),face="bold"),axis.text = element_text(size=20,face = "bold"),
                                       legend.text = element_text(size=15),panel.border = element_rect(colour = "black",size=2))+
  ggtitle("Cluster 2")+scale_y_continuous(limits = c(-.6, .6))+geom_hline(yintercept=0, linetype="dashed", color = "red")

bt_list[[2]]=p4


col_fun1 = circlize::colorRamp2(c(0, 0.5, 1), c("green", "black", "red"),space="RGB")
colnames(variability1)<-c("FEV1","FVC","Ratiopre","WBCDIFF1","WBCDIFF4")
p2<-Heatmap(variability1,cluster_rows=row_ha,cluster_columns=F,show_row_names=F,show_column_names=T,show_row_dend = F,
            col=col_fun1,name="Variability index",show_heatmap_legend = T,width=unit(4, "cm"),column_names_side="top",
            column_names_gp = grid::gpar(fontsize = 8))

col_fun2 = structure(c("yellow","black","cyan2"), names = c("-1", "0", "1"))
colnames(Sign_weight1)<-c("Fev1pd1a","Fvcprd1","Ratiopre","WBCDIFF1","WBCDIFF4")
p3<-Heatmap(Sign_weight1,cluster_rows=row_ha,cluster_columns=F,show_row_names=F,show_column_names=T,show_row_dend = F,
            col=col_fun2,name="Signed AW-Fisher weights",show_heatmap_legend = F,width=unit(4, "cm"),column_names_side="top",
            column_names_gp = grid::gpar(fontsize = 8))

ht_list2<-p3
vt_list2<-p2
module2<-draw(ht_list2, ht_gap = unit(0.5, "cm"))
#----------------------------------------module3
#----------------------------------------
index1<-(length(a1)+length(a2)+1):(length(a1)+length(a2)+length(a3))
library(ComplexHeatmap)

data1<-X1[,index1]
Coef1<-Coef[index1,]
variability1<-variability[index1,]
weight.matrix1<-weight.matrix[index1,]
#---------------------
data1<-t(data1)
order1<-order(apply(variability1,1,sum),decreasing = F)
data1<-data1[order1,]
variability1<-variability1[order1,]

Coef1<-Coef1[order1,]
weight.matrix1<-weight.matrix1[order1,]
Sign_weight1<-weight.matrix1*sign(Coef1)

data1<-t(scale(t(data1)))

Y$fvcprd1<-(Y$fvcprd1-mean(Y$fvcprd1))/sd(Y$fvcprd1)
Y$WBCDIFF1<-(Y$WBCDIFF1-mean(Y$WBCDIFF1))/sd(Y$WBCDIFF1)


Y$WBCDIFF4<-(Y$WBCDIFF4-mean(Y$WBCDIFF4))/sd(Y$WBCDIFF4)
col_fun1 = circlize::colorRamp2(c(0, 0.5, 1), c("green", "black", "red"),space="RGB")
column_ha = HeatmapAnnotation(fev1pd1a=Y$fev1pd1a,fvcprd1=Y$fvcprd1,
                              ratiopre=Y$ratiopre,WBCDIFF1=Y$WBCDIFF1,
                              WBCDIFF4=Y$WBCDIFF4,
                              annotation_name_side = "left",
                              show_legend = F,
                              height=unit(1, "cm"),
                              simple_anno_size = unit(0.2, "cm"),
                              annotation_name_gp= gpar(fontsize = 8),
                              col = list(fev1pd1a = col_fun1,fvcprd1=circlize::colorRamp2(c(min(Y$fvcprd1), 0, max(Y$WBCDIFF1)), c("green", "black", "red"),space="RGB"),
                                         ratiopre=col_fun1,WBCDIFF1=circlize::colorRamp2(c(min(Y$WBCDIFF1), 0, max(Y$WBCDIFF1)), c("green", "black", "red"),space="RGB"),
                                         WBCDIFF4=circlize::colorRamp2(c(min(Y$WBCDIFF4), 0, max(Y$WBCDIFF4)), c("green", "black", "red"),space="RGB")))

library(circlize)
row_ha = as.dendrogram(hclust(dist(data1)))
col_fun = circlize::colorRamp2(c(min(data1), 0, max(data1)), c("green", "black", "red"),space="RGB")
p1<-Heatmap(data1, top_annotation = column_ha,
            cluster_rows=row_ha,cluster_columns=T,show_row_names=F,show_column_names=F,
            col=col_fun,name="Expression",show_row_dend = F,show_column_dend = F,show_heatmap_legend = F,width=unit(12, "cm"),
            row_title_side="left")
FEV1=sapply(1:dim(data1)[1],function(x) cor(data1[x,],Y$fev1pd1a,method = "spearman"))
FVC=sapply(1:dim(data1)[1],function(x) cor(data1[x,],Y$fvcprd1,method = "spearman"))
ratiopre=sapply(1:dim(data1)[1],function(x) cor(data1[x,],Y$ratiopre,method = "spearman"))
WBCDIFF1=sapply(1:dim(data1)[1],function(x) cor(data1[x,],Y$WBCDIFF1,method = "spearman"))
WBCDIFF4=sapply(1:dim(data1)[1],function(x) cor(data1[x,],Y$WBCDIFF4,method = "spearman"))

box_df=data.frame(FEV1=FEV1,FVC=FVC,ratiopre=ratiopre,WBCDIFF1=WBCDIFF1,WBCDIFF4=WBCDIFF4)
box_df2<-pivot_longer(box_df,everything(),names_to="Phenotypes",values_to="Correlations")
p4<-ggplot(box_df2, aes(x=Phenotypes,y=Correlations)) + xlab("Phenotypes")+ylab("Correlations")+theme_bw()+
  geom_boxplot(fill="red",lwd=1)+theme(legend.position = "right",
                                       plot.title = element_text(size = rel(2), hjust = 0.5,face="bold"),
                                       axis.title = element_text(size = rel(2),face="bold"),axis.text = element_text(size=20,face = "bold"),
                                       legend.text = element_text(size=15),panel.border = element_rect(colour = "black",size=2))+
  ggtitle("Cluster 3")+scale_y_continuous(limits = c(-.6, .6))+geom_hline(yintercept=0, linetype="dashed", color = "red")

bt_list[[3]]=p4


col_fun1 = circlize::colorRamp2(c(0, 0.5, 1), c("green", "black", "red"),space="RGB")
colnames(variability1)<-c("FEV1","FVC","Ratiopre","WBCDIFF1","WBCDIFF4")
p2<-Heatmap(variability1,cluster_rows=row_ha,cluster_columns=F,show_row_names=F,show_column_names=T,show_row_dend = F,
            col=col_fun1,name="Variability index",show_heatmap_legend = F,width=unit(4, "cm"),column_names_side="top",
            column_names_gp = grid::gpar(fontsize = 8))

col_fun2 = structure(c("yellow","black","cyan2"), names = c("-1", "0", "1"))
colnames(Sign_weight1)<-c("FEV1","FVC","Ratiopre","WBCDIFF1","WBCDIFF4")
p3<-Heatmap(Sign_weight1,cluster_rows=row_ha,cluster_columns=F,show_row_names=F,show_column_names=T,show_row_dend = F,
            col=col_fun2,name="Signed AW-Fisher weights",show_heatmap_legend = F,width=unit(4, "cm"),column_names_side="top",
            column_names_gp = grid::gpar(fontsize = 8))

ht_list3<-p3
vt_list3<-p2
module3<-draw(ht_list3, ht_gap = unit(0.5, "cm"))




#----------------------------------------module4
#----------------------------------------
index1<-(length(a1)+length(a2)+length(a3)+1):(length(a1)+length(a2)+length(a3)+length(a4))
library(ComplexHeatmap)

data1<-X1[,index1]
Coef1<-Coef[index1,]
variability1<-variability[index1,]
weight.matrix1<-weight.matrix[index1,]
#---------------------
data1<-t(data1)
order1<-order(apply(variability1,1,sum),decreasing = F)
data1<-data1[order1,]
variability1<-variability1[order1,]

Coef1<-Coef1[order1,]
weight.matrix1<-weight.matrix1[order1,]
Sign_weight1<-weight.matrix1*sign(Coef1)

data1<-t(scale(t(data1)))

Y$fvcprd1<-(Y$fvcprd1-mean(Y$fvcprd1))/sd(Y$fvcprd1)
Y$WBCDIFF1<-(Y$WBCDIFF1-mean(Y$WBCDIFF1))/sd(Y$WBCDIFF1)


Y$WBCDIFF4<-(Y$WBCDIFF4-mean(Y$WBCDIFF4))/sd(Y$WBCDIFF4)
col_fun1 = circlize::colorRamp2(c(0, 0.5, 1), c("green", "black", "red"),space="RGB")
column_ha = HeatmapAnnotation(fev1pd1a=Y$fev1pd1a,fvcprd1=Y$fvcprd1,
                              ratiopre=Y$ratiopre,WBCDIFF1=Y$WBCDIFF1,
                              WBCDIFF4=Y$WBCDIFF4,
                              annotation_name_side = "left",
                              show_legend = F,
                              height=unit(1, "cm"),
                              simple_anno_size = unit(0.2, "cm"),
                              annotation_name_gp= gpar(fontsize = 8),
                              col = list(fev1pd1a = col_fun1,fvcprd1=circlize::colorRamp2(c(min(Y$fvcprd1), 0, max(Y$WBCDIFF1)), c("green", "black", "red"),space="RGB"),
                                         ratiopre=col_fun1,WBCDIFF1=circlize::colorRamp2(c(min(Y$WBCDIFF1), 0, max(Y$WBCDIFF1)), c("green", "black", "red"),space="RGB"),
                                         WBCDIFF4=circlize::colorRamp2(c(min(Y$WBCDIFF4), 0, max(Y$WBCDIFF4)), c("green", "black", "red"),space="RGB")))

library(circlize)
row_ha = as.dendrogram(hclust(dist(data1)))
col_fun = circlize::colorRamp2(c(min(data1), 0, max(data1)), c("green", "black", "red"),space="RGB")
p1<-Heatmap(data1, top_annotation = column_ha,
            cluster_rows=row_ha,cluster_columns=T,show_row_names=F,show_column_names=F,
            col=col_fun,name="Expression",show_row_dend = F,show_column_dend = F,show_heatmap_legend = F,width=unit(12, "cm"),
            row_title_side="left")
FEV1=sapply(1:dim(data1)[1],function(x) cor(data1[x,],Y$fev1pd1a,method = "spearman"))
FVC=sapply(1:dim(data1)[1],function(x) cor(data1[x,],Y$fvcprd1,method = "spearman"))
ratiopre=sapply(1:dim(data1)[1],function(x) cor(data1[x,],Y$ratiopre,method = "spearman"))
WBCDIFF1=sapply(1:dim(data1)[1],function(x) cor(data1[x,],Y$WBCDIFF1,method = "spearman"))
WBCDIFF4=sapply(1:dim(data1)[1],function(x) cor(data1[x,],Y$WBCDIFF4,method = "spearman"))

box_df=data.frame(FEV1=FEV1,FVC=FVC,ratiopre=ratiopre,WBCDIFF1=WBCDIFF1,WBCDIFF4=WBCDIFF4)
box_df2<-pivot_longer(box_df,everything(),names_to="Phenotypes",values_to="Correlations")
p4<-ggplot(box_df2, aes(x=Phenotypes,y=Correlations)) + xlab("Phenotypes")+ylab("Correlations")+theme_bw()+
  geom_boxplot(fill="red",lwd=1)+theme(legend.position = "right",
                                       plot.title = element_text(size = rel(2), hjust = 0.5,face="bold"),
                                       axis.title = element_text(size = rel(2),face="bold"),axis.text = element_text(size=20,face = "bold"),
                                       legend.text = element_text(size=15),panel.border = element_rect(colour = "black",size=2))+
  ggtitle("Cluster 4")+scale_y_continuous(limits = c(-.6, .6))+geom_hline(yintercept=0, linetype="dashed", color = "red")

bt_list[[4]]=p4


col_fun1 = circlize::colorRamp2(c(0, 0.5, 1), c("green", "black", "red"),space="RGB")
colnames(variability1)<-c("FEV1","FVC","Ratiopre","WBCDIFF1","WBCDIFF4")
p2<-Heatmap(variability1,cluster_rows=row_ha,cluster_columns=F,show_row_names=F,show_column_names=T,show_row_dend = F,
            col=col_fun1,name="Variability index",show_heatmap_legend = F,width=unit(4, "cm"),column_names_side="top",
            column_names_gp = grid::gpar(fontsize = 8))

col_fun2 = structure(c("yellow","black","cyan2"), names = c("-1", "0", "1"))
colnames(Sign_weight1)<-c("FEV1","FVC","Ratiopre","WBCDIFF1","WBCDIFF4")
p3<-Heatmap(Sign_weight1,cluster_rows=row_ha,cluster_columns=F,show_row_names=F,show_column_names=T,show_row_dend = F,
            col=col_fun2,name="Signed AW-Fisher weights",show_heatmap_legend = F,width=unit(4, "cm"),column_names_side="top",
            column_names_gp = grid::gpar(fontsize = 8))

ht_list4<-p3
vt_list4<-p2
module4<-draw(ht_list4, ht_gap = unit(0.5, "cm"))

#----------------------------------------module5
#----------------------------------------
index1<-(length(a1)+length(a2)+length(a3)+length(a4)+1):(length(a1)+length(a2)+length(a3)+length(a4)+length(a5))
library(ComplexHeatmap)

data1<-X1[,index1]
Coef1<-Coef[index1,]
variability1<-variability[index1,]
weight.matrix1<-weight.matrix[index1,]
#---------------------
data1<-t(data1)
order1<-order(apply(variability1,1,sum),decreasing = F)
data1<-data1[order1,]
variability1<-variability1[order1,]

Coef1<-Coef1[order1,]
weight.matrix1<-weight.matrix1[order1,]
Sign_weight1<-weight.matrix1*sign(Coef1)

data1<-t(scale(t(data1)))

Y$fvcprd1<-(Y$fvcprd1-mean(Y$fvcprd1))/sd(Y$fvcprd1)
Y$WBCDIFF1<-(Y$WBCDIFF1-mean(Y$WBCDIFF1))/sd(Y$WBCDIFF1)


Y$WBCDIFF4<-(Y$WBCDIFF4-mean(Y$WBCDIFF4))/sd(Y$WBCDIFF4)
col_fun1 = circlize::colorRamp2(c(0, 0.5, 1), c("green", "black", "red"),space="RGB")
column_ha = HeatmapAnnotation(fev1pd1a=Y$fev1pd1a,fvcprd1=Y$fvcprd1,
                              ratiopre=Y$ratiopre,WBCDIFF1=Y$WBCDIFF1,
                              WBCDIFF4=Y$WBCDIFF4,
                              annotation_name_side = "left",
                              show_legend = F,
                              height=unit(1, "cm"),
                              simple_anno_size = unit(0.2, "cm"),
                              annotation_name_gp= gpar(fontsize = 8),
                              col = list(fev1pd1a = col_fun1,fvcprd1=circlize::colorRamp2(c(min(Y$fvcprd1), 0, max(Y$WBCDIFF1)), c("green", "black", "red"),space="RGB"),
                                         ratiopre=col_fun1,WBCDIFF1=circlize::colorRamp2(c(min(Y$WBCDIFF1), 0, max(Y$WBCDIFF1)), c("green", "black", "red"),space="RGB"),
                                         WBCDIFF4=circlize::colorRamp2(c(min(Y$WBCDIFF4), 0, max(Y$WBCDIFF4)), c("green", "black", "red"),space="RGB")))

library(circlize)
row_ha = as.dendrogram(hclust(dist(data1)))
col_fun = circlize::colorRamp2(c(min(data1), 0, max(data1)), c("green", "black", "red"),space="RGB")
p1<-Heatmap(data1, top_annotation = column_ha,
            cluster_rows=row_ha,cluster_columns=T,show_row_names=F,show_column_names=F,
            col=col_fun,name="Expression",show_row_dend = F,show_column_dend = F,show_heatmap_legend = F,width=unit(12, "cm"),
            row_title_side="left")
FEV1=sapply(1:dim(data1)[1],function(x) cor(data1[x,],Y$fev1pd1a,method = "spearman"))
FVC=sapply(1:dim(data1)[1],function(x) cor(data1[x,],Y$fvcprd1,method = "spearman"))
ratiopre=sapply(1:dim(data1)[1],function(x) cor(data1[x,],Y$ratiopre,method = "spearman"))
WBCDIFF1=sapply(1:dim(data1)[1],function(x) cor(data1[x,],Y$WBCDIFF1,method = "spearman"))
WBCDIFF4=sapply(1:dim(data1)[1],function(x) cor(data1[x,],Y$WBCDIFF4,method = "spearman"))

box_df=data.frame(FEV1=FEV1,FVC=FVC,ratiopre=ratiopre,WBCDIFF1=WBCDIFF1,WBCDIFF4=WBCDIFF4)
box_df2<-pivot_longer(box_df,everything(),names_to="Phenotypes",values_to="Correlations")
p4<-ggplot(box_df2, aes(x=Phenotypes,y=Correlations)) + xlab("Phenotypes")+ylab("Correlations")+theme_bw()+
  geom_boxplot(fill="red",lwd=1)+theme(legend.position = "right",
                                       plot.title = element_text(size = rel(2), hjust = 0.5,face="bold"),
                                       axis.title = element_text(size = rel(2),face="bold"),axis.text = element_text(size=20,face = "bold"),
                                       legend.text = element_text(size=15),panel.border = element_rect(colour = "black",size=2))+
  ggtitle("Cluster 6")+scale_y_continuous(limits = c(-.6, .6))+geom_hline(yintercept=0, linetype="dashed", color = "red")

bt_list[[6]]=p4


col_fun1 = circlize::colorRamp2(c(0, 0.5, 1), c("green", "black", "red"),space="RGB")
colnames(variability1)<-c("FEV1","FVC","Ratiopre","WBCDIFF1","WBCDIFF4")
p2<-Heatmap(variability1,cluster_rows=row_ha,cluster_columns=F,show_row_names=F,show_column_names=T,show_row_dend = F,
            col=col_fun1,name="Variability index",show_heatmap_legend = F,width=unit(4, "cm"),column_names_side="top",
            column_names_gp = grid::gpar(fontsize = 8))

col_fun2 = structure(c("yellow","black","cyan2"), names = c("-1", "0", "1"))
colnames(Sign_weight1)<-c("FEV1","FVC","Ratiopre","WBCDIFF1","WBCDIFF4")
p3<-Heatmap(Sign_weight1,cluster_rows=row_ha,cluster_columns=F,show_row_names=F,show_column_names=T,show_row_dend = F,
            col=col_fun2,name="Signed AW-Fisher weights",show_heatmap_legend = F,width=unit(4, "cm"),column_names_side="top",
            column_names_gp = grid::gpar(fontsize = 8))

ht_list5<-p3
vt_list5<-p2
module5<-draw(ht_list5, ht_gap = unit(0.5, "cm"))


#----------------------------------------module6
#----------------------------------------

index1<-(length(a1)+length(a2)+length(a3)+length(a4)+length(a5)+1):(length(a1)+length(a2)+length(a3)+length(a4)+length(a5)+length(a6))
library(ComplexHeatmap)

data1<-X1[,index1]
Coef1<-Coef[index1,]
variability1<-variability[index1,]
weight.matrix1<-weight.matrix[index1,]
#---------------------
data1<-t(data1)
order1<-order(apply(variability1,1,sum),decreasing = F)
data1<-data1[order1,]
variability1<-variability1[order1,]

Coef1<-Coef1[order1,]
weight.matrix1<-weight.matrix1[order1,]
Sign_weight1<-weight.matrix1*sign(Coef1)

data1<-t(scale(t(data1)))

Y$fvcprd1<-(Y$fvcprd1-mean(Y$fvcprd1))/sd(Y$fvcprd1)
Y$WBCDIFF1<-(Y$WBCDIFF1-mean(Y$WBCDIFF1))/sd(Y$WBCDIFF1)


Y$WBCDIFF4<-(Y$WBCDIFF4-mean(Y$WBCDIFF4))/sd(Y$WBCDIFF4)
col_fun1 = circlize::colorRamp2(c(0, 0.5, 1), c("green", "black", "red"),space="RGB")
column_ha = HeatmapAnnotation(fev1pd1a=Y$fev1pd1a,fvcprd1=Y$fvcprd1,
                              ratiopre=Y$ratiopre,WBCDIFF1=Y$WBCDIFF1,
                              WBCDIFF4=Y$WBCDIFF4,
                              annotation_name_side = "left",
                              show_legend = F,
                              height=unit(1, "cm"),
                              simple_anno_size = unit(0.2, "cm"),
                              annotation_name_gp= gpar(fontsize = 8),
                              col = list(fev1pd1a = col_fun1,fvcprd1=circlize::colorRamp2(c(min(Y$fvcprd1), 0, max(Y$WBCDIFF1)), c("green", "black", "red"),space="RGB"),
                                         ratiopre=col_fun1,WBCDIFF1=circlize::colorRamp2(c(min(Y$WBCDIFF1), 0, max(Y$WBCDIFF1)), c("green", "black", "red"),space="RGB"),
                                         WBCDIFF4=circlize::colorRamp2(c(min(Y$WBCDIFF4), 0, max(Y$WBCDIFF4)), c("green", "black", "red"),space="RGB")))

library(circlize)
row_ha = as.dendrogram(hclust(dist(data1)))
col_fun = circlize::colorRamp2(c(min(data1), 0, max(data1)), c("green", "black", "red"),space="RGB")
p1<-Heatmap(data1, top_annotation = column_ha,
            cluster_rows=row_ha,cluster_columns=T,show_row_names=F,show_column_names=F,
            col=col_fun,name="Expression",show_row_dend = F,show_column_dend = F,show_heatmap_legend = F,width=unit(12, "cm"),
            row_title_side="left")
FEV1=sapply(1:dim(data1)[1],function(x) cor(data1[x,],Y$fev1pd1a,method = "spearman"))
FVC=sapply(1:dim(data1)[1],function(x) cor(data1[x,],Y$fvcprd1,method = "spearman"))
ratiopre=sapply(1:dim(data1)[1],function(x) cor(data1[x,],Y$ratiopre,method = "spearman"))
WBCDIFF1=sapply(1:dim(data1)[1],function(x) cor(data1[x,],Y$WBCDIFF1,method = "spearman"))
WBCDIFF4=sapply(1:dim(data1)[1],function(x) cor(data1[x,],Y$WBCDIFF4,method = "spearman"))

box_df=data.frame(FEV1=FEV1,FVC=FVC,ratiopre=ratiopre,WBCDIFF1=WBCDIFF1,WBCDIFF4=WBCDIFF4)
box_df2<-pivot_longer(box_df,everything(),names_to="Phenotypes",values_to="Correlations")
p4<-ggplot(box_df2, aes(x=Phenotypes,y=Correlations)) + xlab("Phenotypes")+ylab("Correlations")+theme_bw()+
  geom_boxplot(fill="red",lwd=1)+theme(legend.position = "right",
                                       plot.title = element_text(size = rel(2), hjust = 0.5,face="bold"),
                                       axis.title = element_text(size = rel(2),face="bold"),axis.text = element_text(size=20,face = "bold"),
                                       legend.text = element_text(size=15),panel.border = element_rect(colour = "black",size=2))+
  ggtitle("Cluster 5")+scale_y_continuous(limits = c(-.6, .6))+geom_hline(yintercept=0, linetype="dashed", color = "red")

bt_list[[5]]=p4


col_fun1 = circlize::colorRamp2(c(0, 0.5, 1), c("green", "black", "red"),space="RGB")
colnames(variability1)<-c("FEV1","FVC","Ratiopre","WBCDIFF1","WBCDIFF4")
p2<-Heatmap(variability1,cluster_rows=row_ha,cluster_columns=F,show_row_names=F,show_column_names=T,show_row_dend = F,
            col=col_fun1,name="Variability index",show_heatmap_legend = F,width=unit(4, "cm"),column_names_side="top",
            column_names_gp = grid::gpar(fontsize = 8))

col_fun2 = structure(c("yellow","black","cyan2"), names = c("-1", "0", "1"))
colnames(Sign_weight1)<-c("FEV1","FVC","Ratiopre","WBCDIFF1","WBCDIFF4")
p3<-Heatmap(Sign_weight1,cluster_rows=row_ha,cluster_columns=F,show_row_names=F,show_column_names=T,show_row_dend = F,
            col=col_fun2,name="Signed AW-Fisher weights",show_heatmap_legend = F,width=unit(4, "cm"),column_names_side="top",
            column_names_gp = grid::gpar(fontsize = 8))
#ht_list6<-p1+p2+p3
ht_list6<-p3
vt_list6<-p2
module6<-draw(ht_list6, ht_gap = unit(0.5, "cm"))

#----------------------------------------module7
#----------------------------------------

index1<-(length(a1)+length(a2)+length(a3)+length(a4)+length(a5)+length(a6)+1):(length(a1)+length(a2)+length(a3)+length(a4)+length(a5)+length(a6)+length(a7))
library(ComplexHeatmap)

data1<-X1[,index1]
Coef1<-Coef[index1,]
variability1<-variability[index1,]
weight.matrix1<-weight.matrix[index1,]
#---------------------
data1<-t(data1)
order1<-order(apply(variability1,1,sum),decreasing = F)
data1<-data1[order1,]
variability1<-variability1[order1,]

Coef1<-Coef1[order1,]
weight.matrix1<-weight.matrix1[order1,]
Sign_weight1<-weight.matrix1*sign(Coef1)

data1<-t(scale(t(data1)))
Y$fvcprd1<-(Y$fvcprd1-mean(Y$fvcprd1))/sd(Y$fvcprd1)
Y$WBCDIFF1<-(Y$WBCDIFF1-mean(Y$WBCDIFF1))/sd(Y$WBCDIFF1)


Y$WBCDIFF4<-(Y$WBCDIFF4-mean(Y$WBCDIFF4))/sd(Y$WBCDIFF4)
col_fun1 = circlize::colorRamp2(c(0, 0.5, 1), c("green", "black", "red"),space="RGB")
column_ha = HeatmapAnnotation(fev1pd1a=Y$fev1pd1a,fvcprd1=Y$fvcprd1,
                              ratiopre=Y$ratiopre,WBCDIFF1=Y$WBCDIFF1,
                              WBCDIFF4=Y$WBCDIFF4,
                              annotation_name_side = "left",
                              show_legend = F,
                              height=unit(1, "cm"),
                              simple_anno_size = unit(0.2, "cm"),
                              annotation_name_gp= gpar(fontsize = 8),
                              col = list(fev1pd1a = col_fun1,fvcprd1=circlize::colorRamp2(c(min(Y$fvcprd1), 0, max(Y$WBCDIFF1)), c("green", "black", "red"),space="RGB"),
                                         ratiopre=col_fun1,WBCDIFF1=circlize::colorRamp2(c(min(Y$WBCDIFF1), 0, max(Y$WBCDIFF1)), c("green", "black", "red"),space="RGB"),
                                         WBCDIFF4=circlize::colorRamp2(c(min(Y$WBCDIFF4), 0, max(Y$WBCDIFF4)), c("green", "black", "red"),space="RGB")))

library(circlize)
row_ha = as.dendrogram(hclust(dist(data1)))
col_fun = circlize::colorRamp2(c(min(data1), 0, max(data1)), c("green", "black", "red"),space="RGB")
p1<-Heatmap(data1, top_annotation = column_ha,
            cluster_rows=row_ha,cluster_columns=T,show_row_names=F,show_column_names=F,
            col=col_fun,name="Expression",show_row_dend = F,show_column_dend = F,show_heatmap_legend = F,width=unit(12, "cm"),
            row_title_side="left")
FEV1=sapply(1:dim(data1)[1],function(x) cor(data1[x,],Y$fev1pd1a,method = "spearman"))
FVC=sapply(1:dim(data1)[1],function(x) cor(data1[x,],Y$fvcprd1,method = "spearman"))
ratiopre=sapply(1:dim(data1)[1],function(x) cor(data1[x,],Y$ratiopre,method = "spearman"))
WBCDIFF1=sapply(1:dim(data1)[1],function(x) cor(data1[x,],Y$WBCDIFF1,method = "spearman"))
WBCDIFF4=sapply(1:dim(data1)[1],function(x) cor(data1[x,],Y$WBCDIFF4,method = "spearman"))

box_df=data.frame(FEV1=FEV1,FVC=FVC,ratiopre=ratiopre,WBCDIFF1=WBCDIFF1,WBCDIFF4=WBCDIFF4)
box_df2<-pivot_longer(box_df,everything(),names_to="Phenotypes",values_to="Correlations")
p4<-ggplot(box_df2, aes(x=Phenotypes,y=Correlations)) + xlab("Phenotypes")+ylab("Correlations")+theme_bw()+
  geom_boxplot(fill="red",lwd=1)+theme(legend.position = "right",
                                       plot.title = element_text(size = rel(2), hjust = 0.5,face="bold"),
                                       axis.title = element_text(size = rel(2),face="bold"),axis.text = element_text(size=20,face = "bold"),
                                       legend.text = element_text(size=15),panel.border = element_rect(colour = "black",size=2))+
  ggtitle("Cluster 7")+scale_y_continuous(limits = c(-.6, .6))+geom_hline(yintercept=0, linetype="dashed", color = "red")

bt_list[[7]]=p4


col_fun1 = circlize::colorRamp2(c(0, 0.5, 1), c("green", "black", "red"),space="RGB")
colnames(variability1)<-c("FEV1","FVC","Ratiopre","WBCDIFF1","WBCDIFF4")
p2<-Heatmap(variability1,cluster_rows=row_ha,cluster_columns=F,show_row_names=F,show_column_names=T,show_row_dend = F,
            col=col_fun1,name="Variability index",show_heatmap_legend = F,width=unit(4, "cm"),column_names_side="top",
            column_names_gp = grid::gpar(fontsize = 8))

col_fun2 = structure(c("yellow","black","cyan2"), names = c("-1", "0", "1"))
colnames(Sign_weight1)<-c("FEV1","FVC","Ratiopre","WBCDIFF1","WBCDIFF4")
p3<-Heatmap(Sign_weight1,cluster_rows=row_ha,cluster_columns=F,show_row_names=F,show_column_names=T,show_row_dend = F,
            col=col_fun2,name="Signed weight estimation",show_heatmap_legend = T,width=unit(4, "cm"),column_names_side="top",
            column_names_gp = grid::gpar(fontsize = 8))

ht_list7<-p3
vt_list7<-p2
module7<-draw(ht_list7, ht_gap = unit(0.5, "cm"))







pdf(paste0(out.dir,"weights.pdf"))

ht_list=ht_list1 %v% ht_list2 %v% ht_list3 %v% ht_list4 %v% ht_list6 %v% ht_list5 %v% ht_list7 
draw(ht_list,heatmap_legend_side = "bottom")
dev.off()
pdf(paste0(out.dir,"vari.pdf"))

vt_list=vt_list1 %v% vt_list2 %v% vt_list3 %v% vt_list4 %v% vt_list6 %v% vt_list5 %v% vt_list7 
draw(vt_list,heatmap_legend_side = "bottom")
dev.off()
library(ggpubr)
pdf(paste0(out.dir,"cor_boxplot.pdf"),width = 17.5,height = 17.5)
p_boxplot=ggarrange(plotlist = bt_list,ncol=2,nrow=4)
p_boxplot
dev.off()

#------------------------------pathway analysis
Gene_set1<-name1[index1]
Gene_set2<-name1[index2]
Gene_set3<-name1[index3]
Gene_set4<-name1[index4]
Gene_set5<-name1[index5]
Gene_set6<-name1[index6]
Gene_set7<-name1[index7]

Gene_set<-list(module1=Gene_set1,module2=Gene_set2,module3=Gene_set3,
               module4=Gene_set4,module5=Gene_set5,module6=Gene_set6,
               module7=Gene_set7)


load("pathway.Yujia.Rdata")
background<-colnames(Data.Gene)
pathway_size_limit<-800
new_pathway_list<-pathway.list[which(unlist(lapply(pathway.list,length))>15 & unlist(lapply(pathway.list,length))<pathway_size_limit)]

#
pvalue.set1<-ifelse(background%in%Gene_set1,0.01,1)
pvalue.set2<-ifelse(background%in%Gene_set2,0.01,1)
pvalue.set3<-ifelse(background%in%Gene_set3,0.01,1)
pvalue.set4<-ifelse(background%in%Gene_set4,0.01,1)
pvalue.set5<-ifelse(background%in%Gene_set5,0.01,1)
pvalue.set6<-ifelse(background%in%Gene_set6,0.01,1)
pvalue.set7<-ifelse(background%in%Gene_set7,0.01,1)

data1<-data.frame(gene=background,set1=pvalue.set1,set2=pvalue.set2,
                  set3=pvalue.set3,set4=pvalue.set4,
                  set5=pvalue.set5,set6=pvalue.set6,set7=pvalue.set7)
write.csv(data1,file="data_pathway.csv",row.names = F)

#set 1, 2 and 3
output<-gsa.fisher.new(x=c(Gene_set1,Gene_set2,Gene_set3),background=background,pathway = new_pathway_list)
output$pvalue<-as.numeric(output$pvalue)
output$qvalue<-as.numeric(output$qvalue)
index<-order(as.numeric(output$pvalue),decreasing = F)
output<-output[index,]

output1<-output[1:10,]
write.csv(output,file="module1_2_3_significant_pathway.csv")

#set4
output<-gsa.fisher.new(x=Gene_set4,background=background,pathway = new_pathway_list)
output$pvalue<-as.numeric(output$pvalue)
output$qvalue<-as.numeric(output$qvalue)
index<-order(as.numeric(output$pvalue),decreasing = F)
output<-output[index,]

write.csv(output,file="module4_significant_pathway.csv")
output2<-output[1:10,]
#module set5
output<-gsa.fisher.new(x=Gene_set5,background=background,pathway = new_pathway_list)
output$pvalue<-as.numeric(output$pvalue)
output$qvalue<-as.numeric(output$qvalue)
index<-order(as.numeric(output$pvalue),decreasing = F)
output<-output[index,]

output3<-output[1:10,]
write.csv(output,file="module5_significant_pathway.csv")


#set 6,7
output<-gsa.fisher.new(x=c(Gene_set6,Gene_set7),background=background,pathway = new_pathway_list)
output$pvalue<-as.numeric(output$pvalue)
output$qvalue<-as.numeric(output$qvalue)
index<-order(as.numeric(output$pvalue),decreasing = F)
output<-output[index,]

output4<-output[1:10,]
write.csv(output,file="module6_7_significant_pathway.csv")


data1<-data.frame(module1=output1$pathway,module2=output2$pathway,module3=output3$pathway,
                  module4=output4$pathway)
data1<-pivot_longer(data1,everything(),names_to = "module",values_to = "pathway")
data2<-data.frame(module1=output1$pvalue,module2=output2$pvalue,module3=output3$pvalue,
                  module4=output4$pvalue)
data2<-pivot_longer(data2,everything(),names_to = "module",values_to = "pvalue")
data3<-data.frame(module1=output1$qvalue,module2=output2$qvalue,module3=output3$qvalue,
                  module4=output4$qvalue)
data3<-pivot_longer(data3,everything(),names_to = "module",values_to = "qvalue")
data1$pvalue<-data2$pvalue
data1$qvalue<-data3$qvalue
library(ggplot2)
data1$pvalue<-as.numeric(data1$pvalue)
data1$pathway<-as.character(data1$pathway)
data1$pathway[10]<-" KEGG MAPK signaling pathway"

data1$label<-ifelse(data1$qvalue<0.05,"*","")


library(forcats)


data1$pathway<-fct_reorder(data1$pathway,data1$pvalue)

ggplot(data1)+aes(x=pathway,y=pvalue,fill=module)+
  geom_bar(stat = "identity")+facet_wrap(~module,scales="free",nrow=4,ncol=1)+
  scale_fill_discrete(name="module",labels = c("module1 (C1+C2+C3)","module2 (C4)","module3 (C5)","module4 (C6+C7)"))+
  coord_flip()+ylab("-log10(pvalue)")+theme_bw()+theme(legend.position = "bottom",
                                                       plot.title = element_text(size = rel(2), hjust = 0.5),
                                                       axis.title = element_text(size = rel(2)),
                                                       legend.text = element_text(size=15),
                                                       axis.text.x = element_text(face="bold", 
                                                                                  size=15),
                                                       axis.text.y = element_text(face="bold", 
                                                                                  size=15),
                                                       strip.text.x = element_text(
                                                         size = 15, face = "bold"
                                                       ))


ggplot(data1)+aes(x=pathway,y=pvalue,fill=module)+
  geom_bar(stat = "identity")+facet_wrap(~module,scales="free",nrow=4,ncol=1)+
  scale_fill_discrete(name="module",labels = c("module1 (C1+C2+C3)","module2 (C4)","module3 (C5)","module4 (C6+C7)"))+
  coord_flip()+ylab("-log10(pvalue)")+theme_bw()+theme(legend.position = "bottom",
                                                       plot.title = element_text(size = rel(2), hjust = 0.5),
                                                       axis.title = element_blank(),
                                                       legend.text = element_text(size=15),
                                                       axis.text.x = element_text(face="bold", 
                                                                                  size=15),
                                                       axis.text.y = element_blank(),
                                                       strip.text.x = element_text(
                                                         size = 15, face = "bold"
                                                       ))

data.pathway<-rbind(data1[data1$module=="module1",],data1[data1$module=="module2",],data1[data1$module=="module3",],data1[data1$module=="module4",])
write.csv(data.pathway,file="data.pathway_top10.csv")
#------------------------------------------
gsa.fisher.new <- function(x, background, pathway) {
  ####x is the list of query genes
  ####backgroud is a list of background genes that query genes from 
  ####pathway is a list of different pathway genes
  names.pathway<-names(pathway)
  count_table<-matrix(0,2,2)
  x<-toupper(x)
  background<-toupper(background)
  index<-which(toupper(background) %in% toupper(x)==FALSE)
  background_non_gene_list<-background[index]
  x<-toupper(x)
  pathway<-lapply(pathway,function(x) intersect(toupper(background),toupper(x)))
  get.fisher <- function(path) {
    
    res <- NA
    #        count.table <- table(background %in% x, background %in% toupper(path))
    ####in the gene list and in the pathway
    count_table[1,1]<-sum(x %in% path)
    #count_table[1,1]<-sum(is.na(charmatch(x,path))==0)
    ####in the gene list but not in the pathway
    count_table[1,2]<-length(x)-count_table[1,1]
    ####not in the gene list but in the pathway
    count_table[2,1]<-sum(background_non_gene_list%in% path)
    #count_table[2,1]<-sum(is.na(charmatch(background_non_gene_list,path))==0)
    ####not in the gene list and not in the pathway
    count_table[2,2]<-length(background_non_gene_list)-count_table[2,1]       
    matched_gene<-x[x %in% path]
    
    match_num<-length(matched_gene)
    
    overlap_info<-array(0,dim=4)
    names(overlap_info)<-c("DE in Geneset","DE not in Genese","NonDE in Geneset","NonDE out of Geneset")
    overlap_info[1]=count_table[1,1]
    overlap_info[2]=count_table[1,2]
    overlap_info[3]=count_table[2,1]
    overlap_info[4]=count_table[2,2]
    if(length(count_table)==4){
      res <- fisher.test(count_table, alternative="greater")$p}
    return(list(p_value=res,match_gene=matched_gene,match_num=match_num,
                fisher_table=overlap_info))
  }
  p_val<-array(0,dim=length(pathway))
  
  match_gene_list<-list(length(pathway))
  
  num1<-array(0,dim=length(pathway))
  
  num2<-matrix(0,nrow=length(pathway),ncol=4)
  colnames(num2)<-c("DE in Geneset","DE not in Genese","NonDE in Geneset","NonDE out of Geneset")
  for(i in 1:length(pathway)){
    result<-get.fisher(pathway[[i]])
    p_val[i]<-result$p_value
    match_gene_list[[i]]<-result$match_gene
   
    num1[i]<-result$match_num
    num2[i,]<-result$fisher_table
  }
  
  names(p_val) <-names(pathway)
  q_val <- p.adjust(p_val, "BH")
  
  
  
  summary<-data.frame(pvalue=p_val,
                      qvalue=q_val,
                      pathway=names(pathway)
                      
  )
  
  
  a<-format(summary,digits=3)    
  
  
  

  return(a)
}




