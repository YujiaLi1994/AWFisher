library(mvtnorm)
library(MASS)
library(parallel)
library(Rfast)
library(mvtnorm)
library(MASS)
library(parallel)
library(Rfast)
#library(aSPU)
library(invgamma)
library(statmod)
library(MSKAT)
library(CompQuadForm)
#library(MultiPhen)




#---------------simulation generating function
Sim1D<-function(N,sigma_x,sigma_mu,sigma_y){
  Mu<-matrix(rnorm(4*N,0,sigma_mu),ncol=N)
  y1<-y2<-y3<-y4<-y5<-y6<-y7<-y8<-y9<-y10<-c()
  x1<-x2<-x3<-x4<-x5<-matrix(nrow=N,ncol=50)
  for(i in 1:N){
    y1[i]<-rnorm(1,Mu[1,i],sigma_y[1])
    y2[i]<-rnorm(1,Mu[1,i],sigma_y[2])
    y3[i]<-rnorm(1,Mu[1,i],sigma_y[3])
    y4[i]<-rnorm(1,Mu[1,i],sigma_y[4])
    y5[i]<-rnorm(1,Mu[1,i]+Mu[2,i],sigma_y[5])
    y6[i]<-rnorm(1,Mu[1,i]+Mu[2,i],sigma_y[6])
    y7[i]<-rnorm(1,Mu[1,i]+Mu[2,i],sigma_y[7])
    y8[i]<-rnorm(1,Mu[1,i]+Mu[2,i],sigma_y[8])
    y9[i]<-rnorm(1,Mu[1,i]+Mu[2,i],sigma_y[9])
    y10[i]<-rnorm(1,Mu[3,i],sigma_y[10])
    x1[i,]<-rnorm(50,Mu[1,i],sigma_x[1])
    x2[i,]<-rnorm(50,Mu[2,i],sigma_x[2])
    x3[i,]<-rnorm(50,Mu[3,i],sigma_x[3])
    #x4[i,]<-rnorm(50,Mu[4,i],sigma_x[4])
    
  }
  
  
  x<-cbind(x1,x2,x3)
  y<-cbind(y1,y2,y3,y4,y5,y6,y7,y8,y9,y10)
  data<-list(x=x,y=y)
  return(data)
}

#---------------with covariate adjustment
Sim2D<-function(N,sigma_x,sigma_mu,sigma_y,sigma_c){
  z<-rnorm(N,0,sigma_c)#covariate add on Mu1
  Mu<-matrix(rnorm(4*N,0,sigma_mu),ncol=N)
  y1<-y2<-y3<-y4<-y5<-y6<-y7<-y8<-y9<-y10<-c()
  x1<-x2<-x3<-x4<-x5<-matrix(nrow=N,ncol=50)
  for(i in 1:N){
    y1[i]<-rnorm(1,Mu[1,i]+z[i],sigma_y[1])
    y2[i]<-rnorm(1,Mu[1,i]+z[i],sigma_y[2])
    y3[i]<-rnorm(1,Mu[1,i]+z[i],sigma_y[3])
    y4[i]<-rnorm(1,Mu[1,i]+z[i],sigma_y[4])
    y5[i]<-rnorm(1,Mu[1,i]+z[i]+Mu[2,i],sigma_y[5])
    y6[i]<-rnorm(1,Mu[1,i]+z[i]+Mu[2,i],sigma_y[6])
    y7[i]<-rnorm(1,Mu[1,i]+z[i]+Mu[2,i],sigma_y[7])
    y8[i]<-rnorm(1,Mu[1,i]+z[i]+Mu[2,i],sigma_y[8])
    y9[i]<-rnorm(1,Mu[1,i]+z[i]+Mu[2,i],sigma_y[9])
    y10[i]<-rnorm(1,Mu[3,i],sigma_y[10])
    x1[i,]<-rnorm(50,Mu[1,i]+z[i],sigma_x[1])
    x2[i,]<-rnorm(50,Mu[2,i],sigma_x[2])
    x3[i,]<-rnorm(50,Mu[3,i],sigma_x[3])
    #x4[i,]<-rnorm(50,Mu[4,i],sigma_x[4])
    
  }
  
  
  x<-cbind(x1,x2,x3)
  y<-cbind(y1,y2,y3,y4,y5,y6,y7,y8,y9,y10)
  data<-list(x=x,y=y,z=z)
  return(data)
}

# Sim_new_with_covariate<-function(N,sigma_x,sigma_mu,sigma_y,sigma_c){
#   z<-rnorm(N,0,sigma_c)#covariate add on Mu1
#   Mu<-matrix(rnorm(4*N,0,sigma_mu),ncol=N)
#   y1<-y2<-y3<-y4<-y5<-c()
#   
#   x1<-x2<-x3<-x4<-x5<-matrix(nrow=N,ncol=50)
#   for(i in 1:N){
#     y1[i]<-rnorm(1,Mu[1,i]+z[i],sigma_y)
#     y2[i]<-rnorm(1,Mu[1,i]+Mu[2,i]+z[i],sigma_y)
#     y3[i]<-rnorm(1,Mu[2,i]+Mu[3,i],sigma_y)
#     y4[i]<-rnorm(1,Mu[3,i]+Mu[4,i],sigma_y)
#     y5[i]<-rnorm(1,Mu[4,i],sigma_y)
#     x1[i,]<-rnorm(50,Mu[1,i]+z[i],sigma_x)
#     x2[i,]<-rnorm(50,Mu[2,i],sigma_x)
#     x3[i,]<-rnorm(50,Mu[3,i],sigma_x)
#     x4[i,]<-rnorm(50,Mu[4,i],sigma_x)
#     
#   }
#   
#   
#   x<-cbind(x1,x2,x3,x4)
#   y<-cbind(y1,y2,y3,y4,y5)
#   data<-list(x=x,y=y,z=z)
#   return(data)
# }

# Sim_new<-function(N,sigma_x,sigma_mu,sigma_y){
#   Mu<-matrix(rnorm(4*N,0,sigma_mu),ncol=N)
#   y1<-y2<-y3<-y4<-y5<-c()
#   x1<-x2<-x3<-x4<-x5<-matrix(nrow=N,ncol=50)
#   for(i in 1:N){
#     y1[i]<-rnorm(1,Mu[1,i],sigma_y)
#     y2[i]<-rnorm(1,Mu[1,i]+Mu[2,i],sigma_y)
#     y3[i]<-rnorm(1,Mu[2,i]+Mu[3,i],sigma_y)
#     y4[i]<-rnorm(1,Mu[3,i]+Mu[4,i],sigma_y)
#     y5[i]<-rnorm(1,Mu[4,i],sigma_y)
#     x1[i,]<-rnorm(50,Mu[1,i],sigma_x)
#     x2[i,]<-rnorm(50,Mu[2,i],sigma_x)
#     x3[i,]<-rnorm(50,Mu[3,i],sigma_x)
#     x4[i,]<-rnorm(50,Mu[4,i],sigma_x)
# 
#   }
#   
#   
#   x<-cbind(x1,x2,x3,x4)
#   y<-cbind(y1,y2,y3,y4,y5)
#   data<-list(x=x,y=y)
#   return(data)
# }


Getpvalue<-function(X,Y,Z){
  #X,Y and Z should be matrix
  #if(!is.matrix(X)){stop("X should be a matrix")}
  #if(!is.matrix(Z)){stop("Z should be a matrix")}
  pvalue<-matrix(0,nrow=ncol(X),ncol=ncol(Y))
  Coef<-matrix(0,nrow=ncol(X),ncol=ncol(Y))
  t.stat<-matrix(0,nrow=ncol(X),ncol=ncol(Y))
  Y<-as.matrix(Y)
  for(i in 1:ncol(Y)){
    for(j in 1:ncol(X)){
      #print(j)
      if(is.null(Z)){
        data.temp<-cbind(data.frame(Y=Y[,i],X=X[,j]))
      }else{
        data.temp<-cbind(data.frame(Y=Y[,i],X=X[,j]),Z)
      }
      model<-lm(Y~.,data=data.temp)
      Coef[j,i]<-model$coefficients[2]
      model1<-summary(model)
      t.stat[j,i]<-model1$coefficients[2,3]
      pvalue[j,i]<-model1$coefficients[2,4]
    }
  }
  res<-list(pvalue=pvalue,t.stat=t.stat,coef=Coef)
  return(res)
}
#Getpvalue(X=as.matrix(x1),Y=y,Z=as.matrix(x2))


GetSortMatrix_linear<-function(n){
  Sort.matrix<-matrix(0,nrow=n,ncol=n)
  for(i in 1:n){
    Sort.matrix[i,1:i]<-1
  }
  return(Sort.matrix)
}
#GetSortMatrix_linear(3)


GetSortMatrix_exact<-function(n){
  W<-do.call(expand.grid, rep(list(c(0, 1)), n))[-1,]
  return(W)
}

#GetSortMatrix_exact(3)

Permutation<-function(X){
  n<-dim(X)[1]
  index.permute<-sample(1:n,replace = F)
  X<-X[index.permute,]
  
  return(X)
}
#prove the procedure: multiple regression y~x1+x2 is the same as y~r1+x2, where r1 is the residual of x1 regressed from x2.
# Getpvalue(X=as.matrix(x1),Y=y,Z=as.matrix(x2))
# 
# model<-lm(x1~x2,data=x)
# Getpvalue(X=as.matrix(model$residuals),Y=y,Z=as.matrix(x2))
# 
# Y.residual<-Y
# for(i in 1:ncol(Y)){
#   data.temp<-cbind(data.frame(Y=Y[,i]),cbind(x$x2))
#   model<-lm(Y~.,data=data.temp)
#   Y.residual[,i]<-model$residuals
# }
# 
# Getpvalue(X=as.matrix(model$residuals),Y=Y.residual,Z=as.matrix(rep(1,100)))

GetPermutationPvalue<-function(Data.Gene,Y.residual,number.perm,ncores=1){
  G=ncol(Data.Gene)
  K<-ncol(Y.residual)
  Pvalue.permutation<-NULL
  Tstat.permutation<-NULL
  #Coef.permutation<-NULL
  wcs<-mclapply(1:number.perm,function(i){
    set.seed(i)
    print(i)
    Data.Gene.perm<-Permutation(Data.Gene)
    
    pvalue<-matrix(0,nrow=ncol(Data.Gene.perm),ncol=ncol(Y.residual))
    t.stat<-matrix(0,nrow=ncol(Data.Gene.perm),ncol=ncol(Y.residual))
    for(i in 1:ncol(Y.residual)){
      res<-regression(x=Data.Gene.perm,y=Y.residual[,i])
      pvalue[,i]<-res[,2]
      t.stat[,i]<-res[,1]
    }
    res<-list(pvalue=pvalue,t.stat=t.stat)
    return(res)
  },mc.cores =ncores )
  
  for(i in 1:number.perm){
    Pvalue.permutation[[i]]<-wcs[[i]]$pvalue
    Tstat.permutation[[i]]<-wcs[[i]]$t.stat
    #Coef.permutation[[i]]<-wcs[[i]]$coef
  }
  res<-list(Pvalue=Pvalue.permutation,Tstat=Tstat.permutation)
  return(res)
}





AFp<-function(pvalue,pvalue.perm){
  W<-GetSortMatrix_exact(ncol(pvalue))
  W<-t(W)
  TEMP1<-(-2*log(pvalue.perm))%*%W
  num.nonzero<-apply(W,2,function(x){sum(x>0)})
  #names(num.nonzero)<-1:10
  #num<-ifelse(pvalue<0.1,1,0)%*%W
  #colnames(num)<-1:10
  
  TStat<-(-2*log(pvalue)%*%W)
  
  
  
  AFp<-matrix(NA,nrow=nrow(TStat),ncol=ncol(TStat))
  for(j in 1:ncol(AFp)){
    print(j)
    score.perm<-TEMP1[,j]
    score<-TStat[,j]
    #index1<-match(score,sort(c(score,score.perm),decreasing = F))
    index<-order(c(score,score.perm),decreasing = F)
    index1<-match(1:length(score),index)
    for(i in 1:nrow(AFp)){
      large.all<-(nrow(TEMP1)+nrow(TStat)-index1[i])
      large.sample<-sum(index1>index1[i])
      AFp[i,j]<-(large.all-large.sample)/nrow(TEMP1)
      if(AFp[i,j]==0){{AFp[i,j]<-1/nrow(TEMP1)}}
    }
    #}
  }
  AFp.Stat<-apply(AFp,1,min)#get the test statistic
  final.AFp<-matrix(0,nrow=nrow(AFp),ncol=ncol(pvalue))
  for(i in 1:nrow(final.AFp)){
    index<-which(AFp[i,]==AFp.Stat[i])
    
    num1<-num.nonzero[index]
    Stat<-TStat[i,index]
    index.num<-which.max(num1)
    if(length(index.num)>1){index.num<-index.num[which.max(Stat[index.num])]}
    
    index1<-index[index.num]
    final.AFp[i,]<-W[,index1]
    index<-which(final.AFp[i,]==1)
    final.AFp[i,index]<-ifelse(pvalue[i,index]>0.1,0,1)
  }
  
  #final.AFp[which(pvalue>0.1)]<-0
  #----------------For permuted data
  AFp.perm<-matrix(NA,nrow=nrow(TEMP1),ncol=ncol(TEMP1))
  for(j in 1:ncol(AFp.perm)){
    print(j)
    score.perm<-TEMP1[,j]
    #index<-order(c(score.perm),decreasing = F)
    index<-match(score.perm,sort(score.perm,decreasing = F))
    AFp.perm[,j]<-(nrow(TEMP1)-index+1)/nrow(TEMP1)
    #index1<-match(1:nrow(TStat),index)
    
  }
  AFp.Stat.perm<-apply(AFp.perm,1,min)
  pvalue.res<-unlist(lapply(AFp.Stat,function(x){sum(AFp.Stat.perm<x)/nrow(pvalue.perm)}))
  #pvalue.res<-sum(AFp.Stat.perm<AFp.Stat)/nrow(pvalue.perm)
  pvalue.res<-unlist(lapply(pvalue.res,function(x){
    if(x==0){
      x<-1/nrow(pvalue.perm)
    }
    return(x)
    }))
  
  #if(pvalue.res==0){pvalue.res<-1/nrow(pvalue.perm)}
  res<-list(W=final.AFp,AFp.pvalue=pvalue.res)
  return(res)
}


AFz<-function(pvalue,pvalue.perm){
  W<-GetSortMatrix_exact(ncol(pvalue))
  W<-t(W)
  TEMP1<-(-2*log(pvalue.perm))%*%W
  TStat<-(-2*log(pvalue)%*%W)
  AFz<-matrix(NA,nrow=nrow(TStat),ncol=ncol(TStat))
  for(j in 1:ncol(AFz)){
    AFz[,j]<-abs((TStat[,j]-mean(TEMP1[,j]))/sd(TEMP1[,j]))
  }
  AFz.Stat<-apply(AFz,1,max)#get the test statistic
  final.AFz<-matrix(NA,nrow=nrow(TStat),ncol=ncol(pvalue))
  for(i in 1:nrow(AFz)){
    final.AFz[i,]<-W[,which.max(AFz[i,])]
    index<-which(final.AFz[i,]==1)
    final.AFz[i,index]<-ifelse(pvalue[i,index]>0.1,0,1)
  }
  
  #----------------For permuted data
  AFz.perm<-matrix(NA,nrow=nrow(TEMP1),ncol=ncol(TEMP1))
  for(j in 1:ncol(AFz.perm)){
    print(j)
    AFz.perm[,j]<-abs((TEMP1[,j]-mean(TEMP1[,j]))/sd(TEMP1[,j]))
    #index1<-match(1:nrow(TStat),index)
    
  }
  
  AFz.Stat.perm<-apply(AFz.perm,1,max)
  pvalue.res<-unlist(lapply(AFz.Stat,function(x){sum(AFz.Stat.perm>x)/nrow(pvalue.perm)}))
  #pvalue.res<-sum(AFp.Stat.perm<AFp.Stat)/nrow(pvalue.perm)
  pvalue.res<-unlist(lapply(pvalue.res,function(x){
    if(x==0){
      x<-1/nrow(pvalue.perm)
    }
    return(x)
  }))
  
  #pvalue.res<-sum(AFz.Stat.perm>AFz.Stat)/nrow(pvalue.perm)
  #if(pvalue.res==0){pvalue.res<-1/nrow(pvalue.perm)}
  res<-list(W=final.AFz,AFz.pvalue=pvalue.res)
  return(res)
}

AFz.linear<-function(pvalue,pvalue.perm){
  W<-GetSortMatrix_linear(ncol(pvalue))
  W<-t(W)
  index.pvalue<-matrix(NA,nrow=nrow(pvalue),ncol=ncol(pvalue))
  pvalue.sort<-matrix(NA,nrow=nrow(pvalue),ncol=ncol(pvalue))
  for(i in 1:nrow(index.pvalue)){
    index.pvalue[i,]<-order(pvalue[i,],decreasing = F)
    pvalue.sort[i,]<-sort(pvalue[i,],decreasing = F)
  }
  
  pvalue.perm.sort<-apply(pvalue.perm,1,function(x){sort(x,decreasing = F)})
  pvalue.perm.sort<-t(pvalue.perm.sort)
  
  TEMP1<-(-2*log(pvalue.perm.sort))%*%W
  TStat<-(-2*log(pvalue.sort)%*%W)
  AFz.linear<-matrix(NA,nrow=nrow(TStat),ncol=ncol(TStat))
  for(j in 1:ncol(AFz.linear)){
    AFz.linear[,j]<-abs((TStat[,j]-mean(TEMP1[,j]))/sd(TEMP1[,j]))
  }
  AFz.Stat.linear<-apply(AFz.linear,1,max)#get the test statistic
  
  final.AFz.linear<-matrix(NA,nrow=nrow(TStat),ncol=ncol(pvalue))
  for(i in 1:nrow(AFz.linear)){
    final.AFz.linear[i,]<-W[,which.max(AFz.linear[i,])]
    sig.feature<-index.pvalue[i,][final.AFz.linear[i,]==1]
    final.AFz.linear[i,]<-as.numeric((1:ncol(pvalue))%in%sig.feature)
    index<-which(final.AFz.linear[i,]==1)
    final.AFz.linear[i,index]<-ifelse(pvalue[i,index]>0.1,0,1)
  }
  
  
  
  #final.AFz.linear<-W[,which.max(AFz.linear)]
  #sig.feature<-index.pvalue[final.AFz.linear==1]
  #----------------For permuted data
  AFz.perm.linear<-matrix(NA,nrow=nrow(TEMP1),ncol=ncol(TEMP1))
  for(j in 1:ncol(AFz.perm.linear)){
    print(j)
    AFz.perm.linear[,j]<-abs((TEMP1[,j]-mean(TEMP1[,j]))/sd(TEMP1[,j]))
    #index1<-match(1:nrow(TStat),index)
    
  }
  
  AFz.Stat.perm.linear<-apply(AFz.perm.linear,1,max)
  
  pvalue.res<-unlist(lapply(AFz.Stat.linear,function(x){sum(AFz.Stat.perm.linear>x)/nrow(pvalue.perm)}))
  #pvalue.res<-sum(AFp.Stat.perm<AFp.Stat)/nrow(pvalue.perm)
  pvalue.res<-unlist(lapply(pvalue.res,function(x){
    if(x==0){
      x<-1/nrow(pvalue.perm)
    }
    return(x)
  }))
  res<-list(W=final.AFz.linear,AFz.linear.pvalue=pvalue.res)
  return(res)
}



AFp.linear<-function(pvalue,pvalue.perm){
  W<-GetSortMatrix_linear(ncol(pvalue))
  W<-t(W)
  index.pvalue<-matrix(NA,nrow=nrow(pvalue),ncol=ncol(pvalue))
  pvalue.sort<-matrix(NA,nrow=nrow(pvalue),ncol=ncol(pvalue))
  for(i in 1:nrow(index.pvalue)){
    index.pvalue[i,]<-order(pvalue[i,],decreasing = F)
    pvalue.sort[i,]<-sort(pvalue[i,],decreasing = F)
  }
  
  pvalue.perm.sort<-apply(pvalue.perm,1,function(x){sort(x,decreasing = F)})
  pvalue.perm.sort<-t(pvalue.perm.sort)
  
  TEMP1<-(-2*log(pvalue.perm.sort))%*%W
  num.nonzero<-apply(W,2,function(x){sum(x>0)})
  #names(num.nonzero)<-1:10
  #num<-ifelse(pvalue.sort<0.1,1,0)%*%W
  #colnames(num)<-1:10
  
  TStat<-(-2*log(pvalue.sort)%*%W)
  AFp.linear<-matrix(NA,nrow=nrow(TStat),ncol=ncol(TStat))
  for(j in 1:ncol(AFp.linear)){
    print(j)
    score.perm<-TEMP1[,j]
    score<-TStat[,j]
    index<-order(c(score,score.perm),decreasing = F)
    index1<-match(1:nrow(TStat),index)
    for(i in 1:nrow(AFp.linear)){
      large.all<-(nrow(TEMP1)+nrow(TStat)-index1[i])
      large.sample<-sum(index1>index1[i])
      AFp.linear[i,j]<-(large.all-large.sample)/nrow(TEMP1)
      if(AFp.linear[i,j]==0){{AFp.linear[i,j]<-1/nrow(TEMP1)}}
    }
  }
  AFp.linear.Stat<-apply(AFp.linear,1,min)#get the test statistic
  final.AFp.linear<-matrix(0,nrow=nrow(AFp.linear),ncol=ncol(pvalue))
  
  for(i in 1:nrow(final.AFp.linear)){
    index<-which(AFp.linear[i,]==AFp.linear.Stat[i])
    
    num1<-num.nonzero[index]
    Stat<-TStat[i,index]
    index.num<-which.max(num1)
    if(length(index.num)>1){index.num<-index.num[which.max(Stat[index.num])]}
    
    index1<-index[index.num]
    final.AFp.linear[i,]<-W[,index1]
    index<-which(final.AFp.linear[i,]==1)
    final.AFp.linear[i,index]<-ifelse(pvalue.sort[i,index]>0.1,0,1)
    
    sig.feature<-index.pvalue[i,final.AFp.linear[i,]==1]
    final.AFp.linear[i,]<-as.numeric((1:ncol(pvalue))%in%sig.feature)
  }
  
  
  
  # index<-which(AFp.linear==min(AFp.linear))
  # num1<-num[index]
  # index1<-index[which.max(num1)]
  # final.AFp.linear<-W[,index1]
  #final.AFp.linear[which(pvalue.sort>0.1)]<-0
  
  #----------------For permuted data
  AFp.linear.perm<-matrix(NA,nrow=nrow(TEMP1),ncol=ncol(TEMP1))
  for(j in 1:ncol(AFp.linear.perm)){
    print(j)
    score.perm<-TEMP1[,j]
    #score<-TStat[,j]
    index<-match(score.perm,sort(score.perm,decreasing = F))
    AFp.linear.perm[,j]<-(nrow(TEMP1)-index+1)/nrow(TEMP1)
    #index1<-match(1:nrow(TStat),index)
    
  }
  AFp.linear.Stat.perm<-apply(AFp.linear.perm,1,min)
  pvalue.res<-unlist(lapply(AFp.linear.Stat,function(x){sum(AFp.linear.Stat.perm<x)/nrow(pvalue.perm)}))
  #pvalue.res<-sum(AFp.Stat.perm<AFp.Stat)/nrow(pvalue.perm)
  pvalue.res<-unlist(lapply(pvalue.res,function(x){
    if(x==0){
      x<-1/nrow(pvalue.perm)
    }
    return(x)
  }))
  
  res<-list(W=final.AFp.linear,AFp.linear.pvalue=pvalue.res)
  return(res)
}
