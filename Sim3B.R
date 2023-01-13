#speed4 279890
library(truncnorm)

Getpvalue<-function(X,Y,Z,type="continuous"){
  #X,Y and Z should be matrix
  #if(!is.matrix(X)){stop("X should be a matrix")}
  #if(!is.matrix(Z)){stop("Z should be a matrix")}
  pvalue<-matrix(0,nrow=ncol(X),ncol=ncol(Y))
  Coef<-matrix(0,nrow=ncol(X),ncol=ncol(Y))
  t.stat<-matrix(0,nrow=ncol(X),ncol=ncol(Y))
  Y<-as.matrix(Y)
  if(type=="continuous"){
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
  }else if(type=="count"){
    for(i in 1:ncol(Y)){
      for(j in 1:ncol(X)){
        #print(j)
        if(is.null(Z)){
          data.temp<-cbind(data.frame(Y=Y[,i],X=X[,j]))
        }else{
          data.temp<-cbind(data.frame(Y=Y[,i],X=X[,j]),Z)
        }
        model<-glm(Y~.,data=data.temp,family ="poisson")
        Coef[j,i]<-model$coefficients[2]
        model1<-summary(model)
        t.stat[j,i]<-model1$coefficients[2,3]
        pvalue[j,i]<-model1$coefficients[2,4]
      }
    }
  }else if(type=="binary"){
    for(i in 1:ncol(Y)){
      for(j in 1:ncol(X)){
        #print(j)
        if(is.null(Z)){
          data.temp<-cbind(data.frame(Y=Y[,i],X=X[,j]))
        }else{
          data.temp<-cbind(data.frame(Y=Y[,i],X=X[,j]),Z)
        }
        model<-glm(Y~.,data=data.temp,family ="binomial")
        Coef[j,i]<-model$coefficients[2]
        model1<-summary(model)
        t.stat[j,i]<-model1$coefficients[2,3]
        pvalue[j,i]<-model1$coefficients[2,4]
      }
    }
  }else{
    stop("no distributional family found")
  }
  
  res<-list(pvalue=pvalue,t.stat=t.stat,coef=Coef)
  return(res)
}
#Getpvalue(X=as.matrix(x1),Y=y,Z=as.matrix(x2))




GetPermutationPvalue<-function(Data.Gene,Y.residual,number.perm,ncores=1,type="continuous"){
  G=ncol(Data.Gene)
  K<-ncol(Y.residual)
  Pvalue.permutation<-NULL
  Tstat.permutation<-NULL
  #Coef.permutation<-NULL
  if(type=="continuous"){
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
  }else if(type=="count"){
    wcs<-mclapply(1:number.perm,function(i){
      set.seed(i)
      print(i)
      Data.Gene.perm<-Permutation(Data.Gene)
      
      pvalue<-matrix(0,nrow=ncol(Data.Gene.perm),ncol=ncol(Y.residual))
      t.stat<-matrix(0,nrow=ncol(Data.Gene.perm),ncol=ncol(Y.residual))
      for(i in 1:ncol(Y.residual)){
        res<-score.glms(x=Data.Gene.perm,y=Y.residual[,i],oiko = "poisson")
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
  }else if(type=="binary"){
    wcs<-mclapply(1:number.perm,function(i){
      set.seed(i)
      print(i)
      Data.Gene.perm<-Permutation(Data.Gene)
      
      pvalue<-matrix(0,nrow=ncol(Data.Gene.perm),ncol=ncol(Y.residual))
      t.stat<-matrix(0,nrow=ncol(Data.Gene.perm),ncol=ncol(Y.residual))
      for(i in 1:ncol(Y.residual)){
        res<-score.glms(x=Data.Gene.perm,y=Y.residual[,i],oiko = "binomial")
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
  }else{
    stop("no distributional family found")
  }
  res<-list(Pvalue=Pvalue.permutation,Tstat=Tstat.permutation)
  return(res)
}

Sim3D<-function(N,sigma_x,sigma_mu,sigma_y){
  #z<-rnorm(N,0,sigma_c)#covariate add on Mu1
  Mu<-matrix(rnorm(3*N,0,sigma_mu),ncol=N)
  Mu<-exp(Mu)
  #z<-rtruncnorm(N,a=0,mean=1,sigma_c)#covariate add on Mu1
  #Mu<-matrix(rtruncnorm(4*N,a=1,mean=2,sd=sigma_mu),ncol=N)
  #rtruncnorm(ngenes,a=eff_a/2,mean=eff_a,sd=1)
  #Mu<-matrix(rpois(4*N,mean_mu),ncol=N)
  y1<-y2<-y3<-y4<-y5<-y6<-y7<-y8<-y9<-y10<-c()
  
  x1<-x2<-x3<-x4<-x5<-matrix(nrow=N,ncol=50)
  for(i in 1:N){
    y1[i]<-rpois(1,Mu[1,i])
    y2[i]<-rpois(1,Mu[1,i])
    y3[i]<-rpois(1,Mu[1,i])
    y4[i]<-rpois(1,Mu[1,i])
    # y5[i]<-rpois(1,exp(Mu[1,i]+Mu[2,i]))
    # y6[i]<-rpois(1,exp(Mu[1,i]+Mu[2,i]))
    # y7[i]<-rpois(1,exp(Mu[1,i]+Mu[2,i]))
    # y8[i]<-rpois(1,exp(Mu[1,i]+Mu[2,i]))
    # y9[i]<-rpois(1,exp(Mu[1,i]+Mu[2,i]))
    y5[i]<-rnorm(1,Mu[1,i]+Mu[2,i],sigma_y[1])
    y6[i]<-rnorm(1,Mu[1,i]+Mu[2,i],sigma_y[2])
    y7[i]<-rnorm(1,Mu[1,i]+Mu[2,i],sigma_y[3])
    y8[i]<-rnorm(1,Mu[1,i]+Mu[2,i],sigma_y[4])
    y9[i]<-rnorm(1,Mu[1,i]+Mu[2,i],sigma_y[5])
    # y5[i]<-rbinom(1,size=1,p=(exp(Mu[1,i]+Mu[2,i]))/(1+exp(Mu[1,i]+Mu[2,i])))
    # y6[i]<-rbinom(1,size=1,p=(exp(Mu[1,i]+Mu[2,i]))/(1+exp(Mu[1,i]+Mu[2,i])))
    # y7[i]<-rbinom(1,size=1,p=(exp(Mu[1,i]+Mu[2,i]))/(1+exp(Mu[1,i]+Mu[2,i])))
    # y8[i]<-rbinom(1,size=1,p=(exp(Mu[1,i]+Mu[2,i]))/(1+exp(Mu[1,i]+Mu[2,i])))
    # y9[i]<-rbinom(1,size=1,p=(exp(Mu[1,i]+Mu[2,i]))/(1+exp(Mu[1,i]+Mu[2,i])))
    
    y10[i]<-rnorm(1,Mu[3,i],sigma_y[6])
    x1[i,]<-rnorm(50,Mu[1,i],sigma_x)
    x2[i,]<-rnorm(50,Mu[2,i],sigma_x)
    x3[i,]<-rnorm(50,Mu[3,i],sigma_x)
    #x4[i,]<-rnorm(50,Mu[4,i],sigma_x)
    
  }
  
  
  x<-cbind(x1,x2,x3)
  y<-cbind(y1,y2,y3,y4,y5,y6,y7,y8,y9,y10)
  data<-list(x=x,y=y)
  return(data)
}


N<-100
sigma_x<-rep(0.5,3)
sigma_y<-c(0.01,2,2,2,2,1)

library(mvtnorm)
library(MASS)
library(parallel)
library(Rfast)
library(mvtnorm)
library(MASS)
library(parallel)
library(Rfast)
library(aSPU)
library(invgamma)
library(statmod)
library(MSKAT)
library(CompQuadForm)
library(MultiPhen)
source("tates.R")
source("Sim_function_new.R")
sigma_mu_vector<-c(0,0.4,0.6)
#num.sim<-500
num.cores<-30

for(ind.sigma.mu in 1:length(sigma_mu_vector)){
  if(ind.sigma.mu==1){
    num.sim<-500
  }else{
    num.sim<-500
  }
  sigma_mu<-sigma_mu_vector[ind.sigma.mu]
  res<-mclapply(1:num.sim,function(ind.sim){
    set.seed(ind.sim)
    data<-Sim3D(N,sigma_x = sigma_x,sigma_mu = sigma_mu,sigma_y=sigma_y)
    x<-data$x
    y<-data$y
    x<-as.data.frame(x)
    #z<-data$z
    colnames(x)<-paste("x",1:ncol(x),sep="")
    
    #quantile(cor(y[,1],x[,1:50]))
    #MANOVA
    # pvalue.manova<-rep(NA,ncol(x))
    # for(i in 1:ncol(x)){
    #   mod<-manova(y~x[,i])
    #   mod1<-summary(mod, test="Wilks",intercept = T)
    #   pvalue.manova[i]<-mod1$stats[2,6]
    # }
    #quantile(pvalue.manova)
    res.origin.count<-Getpvalue(X = as.matrix(x),Y = as.matrix(y[,1:4]),Z = NULL,type = "count")
    #res.origin.binary<-Getpvalue(X = as.matrix(x),Y = as.matrix(y[,5:9]),Z = NULL,type="binary")
    res.origin.continuous<-Getpvalue(X = as.matrix(x),Y = as.matrix(y[,5:10]),Z = NULL,type="continuous")
    
    #pvalue<-cbind(res.origin.count$pvalue,res.origin.binary$pvalue,res.origin.continuous$pvalue)
    pvalue<-cbind(res.origin.count$pvalue,res.origin.continuous$pvalue)
    
    res.perm.count<-GetPermutationPvalue(as.matrix(x),as.matrix(y[,1:4]),number.perm=500,ncores=1,type = "count")
    Pvalue.permutation.count<-NULL
    for(i in 1:500){
      Pvalue.permutation.count<-rbind(Pvalue.permutation.count,res.perm.count$Pvalue[[i]])
    }
    
    # res.perm.binary<-GetPermutationPvalue(as.matrix(x),as.matrix(y[,5:9]),number.perm=500,ncores=1,type = "binary")
    # Pvalue.permutation.binary<-NULL
    # for(i in 1:500){
    #   Pvalue.permutation.binary<-rbind(Pvalue.permutation.binary,res.perm.binary$Pvalue[[i]])
    # }
    
    res.perm.continuous<-GetPermutationPvalue(as.matrix(x),as.matrix(y[,5:10]),number.perm=500,ncores=1,type = "continuous")
    Pvalue.permutation.continuous<-NULL
    for(i in 1:500){
      Pvalue.permutation.continuous<-rbind(Pvalue.permutation.continuous,res.perm.continuous$Pvalue[[i]])
    }
    #Pvalue.permutation<-cbind(Pvalue.permutation.count,Pvalue.permutation.binary,Pvalue.permutation.continuous)
    Pvalue.permutation<-cbind(Pvalue.permutation.count,Pvalue.permutation.continuous)
    
    #GWAS
    # mod.aSPU.ind<-apply(x,2,function(x1){
    #   GEEaSPU(traits = y,geno = as.matrix(x1),Z=z,model = "gaussian",corstr = "independence",n.sim = 1000)[10]
    # })
    # mod.aSPU.ex<-apply(x,2,function(x1){
    #   GEEaSPU(traits = y,geno = as.matrix(x1),Z=z,model = "gaussian",corstr = "exchangeable",n.sim = 1000)[10]
    # })
    
    
    #AFp
    mod.AFp<-AFp(pvalue = pvalue,pvalue.perm = Pvalue.permutation)
    #AFp.linear
    #mod.AFp.linear<-AFp.linear(pvalue = res.origin$pvalue,pvalue.perm = Pvalue.permutation)
    #AFz
    mod.AFz<-AFz(pvalue = pvalue,pvalue.perm = Pvalue.permutation)
    #AFz.linear
    #mod.AFz.linear<-AFz.linear(pvalue = res.origin$pvalue,pvalue.perm = Pvalue.permutation)
    
    #minp
    minp<-apply(pvalue,1,min)
    minp.perm<-apply(Pvalue.permutation,1,min)
    p.minp<-unlist(lapply(minp,function(x){sum(minp.perm<x)/length(minp.perm)}))
    #p.minp<-sum(minp.perm<minp)/100000
    #sum(p.minp<0.05)/length(p.minp)
    #EW Fisher
    EW<-apply(pvalue,1,function(x){sum(-log(x))})
    EW.perm<-apply(Pvalue.permutation,1,function(x){sum(-log(x))})
    p.EW<-unlist(lapply(EW,function(x){sum(EW.perm>x)/length(EW.perm)}))
    #sum(p.EW<0.05)/length(p.EW)
    
    # index1<-which(p.minp<0.05)
    # index2<-which(p.EW<0.05)
    # index3<-which(mod.AFp$AFp.pvalue<0.05)
    # index4<-which(pvalue.manova<0.05)
    #p.EW<-sum(EW.perm>EW)/100000
    ### TATES ###
    mod.TATEs <-apply(pvalue,1,function(x1){
      TATES(Y=y,pval_matrix=matrix(x1,nrow=1),trait=5,n.snp = 1)[6]
    })
    
    #-----------------------------------truncated at 10^(-8)
    # pvalue<-res.origin$pvalue
    # pvalue<-ifelse(pvalue<10^(-6),10^(-6),pvalue)
    # #AFp
    # mod.AFp.truc<-AFp(pvalue = pvalue,pvalue.perm = Pvalue.permutation)
    # #AFp.linear
    # mod.AFp.linear.truc<-AFp.linear(pvalue = pvalue,pvalue.perm = Pvalue.permutation)
    # #AFz
    # mod.AFz.truc<-AFz(pvalue = pvalue,pvalue.perm = Pvalue.permutation)
    # #AFz.linear
    # mod.AFz.linear.truc<-AFz.linear(pvalue = pvalue,pvalue.perm = Pvalue.permutation)
    # 
    # #minp
    # minp<-apply(pvalue,1,min)
    # minp.perm<-apply(Pvalue.permutation,1,min)
    # p.minp.truc<-unlist(lapply(minp,function(x){sum(minp.perm<x)/length(minp.perm)}))
    # #p.minp<-sum(minp.perm<minp)/100000
    # #EW Fisher
    # EW<-apply(pvalue,1,function(x){sum(-log(x))})
    # EW.perm<-apply(Pvalue.permutation,1,function(x){sum(-log(x))})
    # p.EW.truc<-unlist(lapply(EW,function(x){sum(EW.perm>x)/length(EW.perm)}))
    # #p.EW<-sum(EW.perm>EW)/100000
    # ### TATES ###
    # mod.TATEs.truc <-apply(pvalue,1,function(x1){
    #   TATES(Y=y,pval_matrix=matrix(x1,nrow=1),trait=5,n.snp = 1)[6]
    # })
    # 
    # 
    final.res<-list(
                    AFp=mod.AFp,AFz=mod.AFz,minp=p.minp,EW=p.EW,mod.TATEs=mod.TATEs,
                    input.stat=pvalue)
    
    # final.res<-list(manova=pvalue.manova,mod.aSPU.ind=mod.aSPU.ind,mod.aSPU.ex=mod.aSPU.ex,
    #                 AFp=mod.AFp,AFp.linear=mod.AFp.linear,AFz=mod.AFz,
    #                 AFz.linear=mod.AFz.linear,minp=p.minp,EW=p.EW,mod.TATEs=mod.TATEs,
    #                 AFp.truc=mod.AFp.truc,AFp.linear.truc=mod.AFp.linear.truc,AFz.truc=mod.AFz.truc,
    #                 AFz.linear.truc=mod.AFz.linear.truc,minp.truc=p.minp.truc,EW.truc=p.EW.truc,
    #                 mod.TATEs.truc=mod.TATEs.truc,
    #                 input.stat=res.origin)
    return(final.res)
  },mc.cores = num.cores)
  file_name<-paste("Sim3E_hierarchical_sigma_mu=",sigma_mu,".Rdata",sep="")
  save(res,file=file_name)
}
  

