setwd("~/OneDrive - University of Pittsburgh/_BoxMigration/yujia_pitt/pitt/research/AWFisher/lung disease")
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
#Index.Mean<-which(apply(Data.Gene,2,mean)>5)
#Data.Gene<-Data.Gene[,Index.Mean]
#Index.Sd<-which(apply(Data.Gene,2,function(x){sqrt(var(x))})>0.5)
#Data.Gene<-Data.Gene[,Index.Sd]
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

res.origin<-Getpvalue(X = as.matrix(X),Y = Y,Z = as.matrix(Z))
Y.residual<-Y
for(i in 1:ncol(Y)){
  data.temp<-cbind(data.frame(Y=Y[,i]),Z)
  model<-lm(Y~.,data=data.temp)
  Y.residual[,i]<-model$residuals
}
r1<-X
for(i in 1:ncol(X)){
  data.temp<-cbind(data.frame(X=X[,i]),Z)
  model<-lm(X~.,data=data.temp)
  r1[,i]<-model$residuals
}

library(parallel)
res.perm<-GetPermutationPvalue(as.matrix(r1),Y.residual,number.perm=100,ncores = 1)
Pvalue.permutation<-NULL
for(i in 1:100){
  Pvalue.permutation<-rbind(Pvalue.permutation,res.perm$Pvalue[[i]])
}

#AFp
t1<-Sys.time()
mod.AFp<-AFp(pvalue = res.origin$pvalue,pvalue.perm = Pvalue.permutation)
t2<-Sys.time()
t2-t1
#AFp.linear
t1<-Sys.time()
mod.AFp.linear<-AFp.linear(pvalue = res.origin$pvalue,pvalue.perm = Pvalue.permutation)
t2<-Sys.time()
t2-t1
#AFz
mod.AFz<-AFz(pvalue = res.origin$pvalue,pvalue.perm = Pvalue.permutation)
#AFz.linear
mod.AFz.linear<-AFz.linear(pvalue = res.origin$pvalue,pvalue.perm = Pvalue.permutation)

#minp
minp<-apply(res.origin$pvalue,1,min)
minp.perm<-apply(Pvalue.permutation,1,min)
p.minp<-unlist(lapply(minp,function(x){sum(minp.perm<x)/length(minp.perm)}))
#p.minp<-sum(minp.perm<minp)/100000
#EW Fisher
EW<-apply(res.origin$pvalue,1,function(x){sum(-log(x))})
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
                input.stat=res.origin)
save(final.res,file="lung_disease_result_all.Rdata")

#GWAS Methods

# mod.aSPU.ind<-apply(X,2,function(x1){
#   GEEaSPU(traits = Y,geno = as.matrix(x1),Z=Z,model = "gaussian",corstr = "independence")[10]
# })
# mod.aSPU.ex<-apply(X,2,function(x1){
#   GEEaSPU(traits = Y,geno = as.matrix(x1),Z=Z,model = "gaussian",corstr = "exchangeable")[10]
# })

#GWAS Methods
X1<-apply(X,2,list)
mod.aSPU.ind<-mclapply(X1,function(x1){
  GEEaSPU(traits = Y,geno = as.matrix(x1[[1]]),Z=Z,model = "gaussian",corstr = "independence",n.sim = 1596600)[10]
},mc.cores = 200)
# mod.aSPU.ex<-mclapply(X1,function(x1){
#   GEEaSPU(traits = Y,geno = as.matrix(x1[[1]]),Z=Z,model = "gaussian",corstr = "exchangeable",n.sim = 10000000)[10]
# },mc.cores = 50)#speed 3
mod.aSPU.ex<-mclapply(X1,function(x1){
  GEEaSPU(traits = Y,geno = as.matrix(x1[[1]]),Z=Z,model = "gaussian",corstr = "exchangeable",n.sim = 1596600)[10]
},mc.cores = 200)#speed 5
#aSPU<-GEEaSPU(traits = y,geno = as.matrix(x[,1]),Z=NULL,model = "gaussian",corstr = "independence")
### TATES ###
mod.TATEs <-apply(res.origin$pvalue,1,function(x1){
  TATES(Y=Y,pval_matrix=matrix(x1,nrow=1),trait=5,n.snp = 1)[6]
})
# final.res<-list(mod.aSPU.ind=mod.aSPU.ind,mod.aSPU.ex=mod.aSPU.ex,mod.TATEs=mod.TATEs,
#                 input.stat=res.origin)
final.res<-list(mod.aSPU.ex=mod.aSPU.ex)
save(final.res,file="Real_data_aSPU_ex_Method_new.Rdata")

#invisible(capture.output(suppressMessages(pan<-GEEaSPU(traits = Y,geno = as.matrix(X1[[15669]][[1]]),Z=Z,model = "gaussian",corstr = "exchangeable",n.sim = 1000000))))

#-----------------------------------supplementary functions
# library(Rfast)
# GetPermutationPvalue<-function(r1,Y.residual,number.perm,ncores=1){
#   G=ncol(r1)
#   K<-ncol(Y)
#   Pvalue.permutation<-NULL
#   Tstat.permutation<-NULL
#   #Coef.permutation<-NULL
#   wcs<-mclapply(1:number.perm,function(i){
#     set.seed(i)
#     print(i)
#     Data.Gene.perm<-Permutation(r1)
#     
#     pvalue<-matrix(0,nrow=ncol(Data.Gene.perm),ncol=ncol(Y.residual))
#     t.stat<-matrix(0,nrow=ncol(Data.Gene.perm),ncol=ncol(Y.residual))
#     for(i in 1:ncol(Y.residual)){
#       res<-regression(x=Data.Gene.perm,y=Y.residual[,i])
#       pvalue[,i]<-res[,2]
#       t.stat[,i]<-res[,1]
#     }
#     res<-list(pvalue=pvalue,t.stat=t.stat)
#     return(res)
#   },mc.cores =ncores )
#   
#   for(i in 1:number.perm){
#     Pvalue.permutation[[i]]<-wcs[[i]]$pvalue
#     Tstat.permutation[[i]]<-wcs[[i]]$t.stat
#     #Coef.permutation[[i]]<-wcs[[i]]$coef
#   }
#   res<-list(Pvalue=Pvalue.permutation,Tstat=Tstat.permutation)
#   return(res)
# }
# 
# 
# AFp.linear<-function(pvalue,pvalue.perm){
#   W<-GetSortMatrix_linear(ncol(pvalue))
#   W<-t(W)
#   index.pvalue<-order(pvalue,decreasing = F)
#   pvalue.sort<-apply(pvalue,1,function(x){sort(x,decreasing = F)})
#   pvalue.perm.sort<-apply(pvalue.perm,1,function(x){sort(x,decreasing = F)})
#   pvalue.perm.sort<-t(pvalue.perm.sort)
#   TEMP1<-(-2*log(pvalue.perm.sort))%*%W
#   num.nonzero<-apply(W,2,function(x){sum(x>0)})
#   names(num.nonzero)<-1:ncol(pvalue)
#   num<-ifelse(t(pvalue.sort)<0.1,1,0)%*%W
#   colnames(num)<-1:ncol(pvalue)
#   
#   TStat<-(-2*log(t(pvalue.sort))%*%W)
#   AFp.linear<-matrix(NA,nrow=nrow(TStat),ncol=ncol(TStat))
#   for(j in 1:ncol(AFp.linear)){
#     print(j)
#     score.perm<-TEMP1[,j]
#     score<-TStat[,j]
#     index<-order(c(score,score.perm),decreasing = F)
#     index1<-match(1:nrow(TStat),index)
#     for(i in 1:nrow(AFp.linear)){
#       large.all<-(nrow(TEMP1)+nrow(TStat)-index1[i])
#       large.sample<-sum(index1>index1[i])
#       AFp.linear[i,j]<-(large.all-large.sample)/nrow(TEMP1)
#       if(AFp.linear[i,j]==0){{AFp.linear[i,j]<-1/nrow(TEMP1)}}
#     }
#   }
#   AFp.linear.Stat<-apply(AFp.linear,1,min)#get the test statistic
#   final.AFp.linear<-matrix(0,nrow=nrow(AFp.linear),ncol=length(pvalue))
#   
#   index<-apply(AFp.linear,1,function(x){which(x==min(x))})
#   num1<-list()
#   index1<-c()
#   for(i in 1:length(index)){
#     num1[[i]]<-num[i,][index[[i]]]
#     index1[i]<-index[[i]][which.max(num1[[i]])]
#   }
#   #num1<-lapply(index,function(x){num[x]})
#   #index1<-index[which.max(num1)]
#   final.AFp.linear<-W[,index1]
#   #final.AFp.linear[which(pvalue.sort>0.1)]<-0
#   sig.feature<-index.pvalue[final.AFp.linear==1]
#   #----------------For permuted data
#   AFp.linear.perm<-matrix(NA,nrow=nrow(TEMP1),ncol=ncol(TEMP1))
#   for(j in 1:ncol(AFp.linear.perm)){
#     print(j)
#     score.perm<-TEMP1[,j]
#     #score<-TStat[,j]
#     index<-match(score.perm,sort(score.perm,decreasing = F))
#     AFp.linear.perm[,j]<-(nrow(TEMP1)-index+1)/nrow(TEMP1)
#     #index1<-match(1:nrow(TStat),index)
#     
#   }
#   AFp.linear.Stat.perm<-apply(AFp.linear.perm,1,min)
#   Pvalue<-rep(NA,length(AFp.linear.Stat))
#   index<-order(c(AFp.linear.Stat,AFp.linear.Stat.perm),decreasing = F)
#   index1<-match(1:length(AFp.linear.Stat),index)
#   for(i in 1:length(AFp.linear.Stat)){
#     large.all<-(length(AFp.linear.Stat.perm)+length(AFp.linear.Stat)-index1[i])
#     large.sample<-sum(index1>index1[i])
#     Pvalue[i]<-1-(large.all-large.sample)/length(AFp.linear.Stat.perm)
#     
#   }
#   # pvalue.res<-c()
#   # for(i in 1:length(AFp.linear.Stat)){
#   #   pvalue.res[i]<-sum(AFp.linear.Stat.perm<AFp.linear.Stat[i])/nrow(pvalue.perm)
#   #   if(pvalue.res[i]==0){pvalue.res[i]<-1/nrow(pvalue.perm)}
#   # }
#   
#   #pvalue.res<-sum(AFp.linear.Stat.perm<AFp.linear.Stat)/nrow(pvalue.perm)
#   #if(pvalue.res==0){pvalue.res<-1/nrow(pvalue.perm)}
#   res<-list(index=sig.feature,AFp.linear.pvalue=pvalue.res)
#   return(res)
# }
