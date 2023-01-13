#setwd("~/Box/yujia_pitt/pitt/research/AWFisher/lung disease")
rm(list=ls())
#out.dir="C:/Users/acdsee8602/iCloudDrive/AWFISHER_THEORY/GENES/code/real_data/code_result"
out.dir="/home/yuf31/Genes"
setwd(out.dir)
source("/home/yuf31/Genes/Sim_function_new.R")
load("LungData.all.Rdata")
Data.Gene<-LungData.all$gene
Data.Clinical<-LungData.all$clinic
library(dplyr)
#"fev6prd1"
selected.variable<-c("fev1pd1a","fvcprd1","ratiopre","WBCDIFF1","WBCDIFF4")
#confounder.variable<-c("pkyrs","Cigarette_smoking","BMI","AGE","GENDER")
confounder.variable<-c("BMI","AGE","GENDER")
Data.confounder<-dplyr::select(Data.Clinical,confounder.variable)
Y<-dplyr::select(Data.Clinical,selected.variable)

###############filter genes
# Index.Mean<-which(apply(Data.Gene,2,mean)>5)
# Data.Gene<-Data.Gene[,Index.Mean]
# Index.Sd<-which(apply(Data.Gene,2,function(x){sqrt(var(x))})>0.5)
# Data.Gene<-Data.Gene[,Index.Sd]
dim(Data.Gene)
index1<-which(apply(Data.confounder,1,function(x){sum(is.na(x))})!=0)
index2<-which(apply(Y,1,function(x){sum(is.na(x))})!=0)
index_delete<-unique(c(index1,index2))
Data.Gene<-Data.Gene[-index_delete,]
Data.Clinical<-Data.Clinical[-index_delete,]
Data.confounder<-Data.confounder[-index_delete,]
Y<-Y[-index_delete,]
#X<-cbind(Data.Gene,Data.confounder)
X<-Data.Gene
Z<-Data.confounder





# GetSortMatrix_exact<-function(K){
#   W<-do.call(expand.grid, rep(list(c(0, 1)), K))[-1,]
#   return(W)
# }
# Permutation<-function(X){
#   n<-dim(X)[1]
#   for(i in 1:ncol(X)){
#     index.permute<-sample(1:n,replace = F)
#     X[,i]<-X[index.permute,i]
#   }
#   return(X)
# }


num.bootstrap<-50
number.perm<-100
mc.cores=1
set.seed(12315)
library(parallel)
wcs<-mclapply(1:num.bootstrap,function(ii){
  print(ii)
  set.seed(ii)
  index.bootstrap<-sample(1:nrow(X),replace = T)
  X1<-X[index.bootstrap,]
  Y1<-Y[index.bootstrap,]
  Z1<-Data.confounder[index.bootstrap,]
 
  #res.origin<-Getpvalue(X = as.matrix(X1),Y = Y1,Z = as.matrix(Z1))
  Y.residual<-Y1
  for(i in 1:ncol(Y1)){
    data.temp<-cbind(data.frame(Y=Y1[,i]),Z1)
    model<-lm(Y~.,data=data.temp)
    Y.residual[,i]<-model$residuals
  }
  r1<-X1
  for(i in 1:ncol(X1)){
    data.temp<-cbind(data.frame(X=X1[,i]),Z1)
    model<-lm(X~.,data=data.temp)
    r1[,i]<-model$residuals
  }
  
  pvalue.bootstrap<-Getpvalue(X = as.matrix(X1),Y = Y1,Z = as.matrix(Z1))
  #save(pvalue.bootstrap,file=paste("bootstrap_",i,".Rdata",sep=""))
  Permutation.list<-GetPermutationPvalue(as.matrix(X1),Y1,number.perm,ncores=1)
  #save(Permutation.list,file=paste("bootstrap_permutation_",i,".Rdata",sep=""))
  Pvalue.permutation<-NULL
  for(i in 1:100){
    Pvalue.permutation<-rbind(Pvalue.permutation,Permutation.list$Pvalue[[i]])
  }
  
  #AFp
  mod.AFp<-AFp(pvalue = pvalue.bootstrap$pvalue,pvalue.perm = Pvalue.permutation)
  #AFp.linear
  mod.AFp.linear<-AFp.linear(pvalue.bootstrap$pvalue,pvalue.perm = Pvalue.permutation)
  #AFz
  mod.AFz<-AFz(pvalue = pvalue.bootstrap$pvalue,pvalue.perm = Pvalue.permutation)
  #AFz.linear
  mod.AFz.linear<-AFz.linear(pvalue =pvalue.bootstrap$pvalue,pvalue.perm = Pvalue.permutation)
  
  #minp
  minp<-apply(pvalue.bootstrap$pvalue,1,min)
  minp.perm<-apply(Pvalue.permutation,1,min)
  p.minp<-unlist(lapply(minp,function(x){sum(minp.perm<x)/length(minp.perm)}))
  #p.minp<-sum(minp.perm<minp)/100000
  #EW Fisher
  EW<-apply(pvalue.bootstrap$pvalue,1,function(x){sum(-log(x))})
  EW.perm<-apply(Pvalue.permutation,1,function(x){sum(-log(x))})
  p.EW<-unlist(lapply(EW,function(x){sum(EW.perm>x)/length(EW.perm)}))
  #p.EW<-sum(EW.perm>EW)/100000
  
  pvalue.manova<-c()
  for(i in 1:ncol(X)){
    mod<-manova(as.matrix(Y)~X[,i])
    mod1<-summary(mod, test="Wilks",intercept = T)
    pvalue.manova[i]<-mod1$stats[2,6]
  }
  
  final.res<-list(manova=pvalue.manova,AFp=mod.AFp,AFp.linear=mod.AFp.linear,AFz=mod.AFz,
                  AFz.linear=mod.AFz.linear,minp=p.minp,EW=p.EW,
                  input.stat=pvalue.bootstrap)
  save(final.res,file=paste("bootstrap_",ii,".Rdata",sep=""))
  return(final.res)
},mc.cores = mc.cores)


save(wcs,file="lung_disease_bootstrap.Rdata")
#Permutation.list<-GetPermutationPvalue(Data.Gene,Y,number.perm,Data.confounder,ncores=mc.cores)
#save(Permutation.list,file="Permutation.list_1000.Rdata")


