#speed3 279890


N<-100
sigma_x<-rep(0.5,3)
sigma_y<-c(rep(2,9),1)
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
num.cores<-50

for(ind.sigma.mu in 1:length(sigma_mu_vector)){
  if(ind.sigma.mu==1){
    num.sim<-500
  }else{
    num.sim<-500
  }
  sigma_mu<-sigma_mu_vector[ind.sigma.mu]
  res<-mclapply(1:num.sim,function(ind.sim){
    set.seed(ind.sim)
    data<-Sim1D(N,sigma_x = sigma_x,sigma_mu = sigma_mu,sigma_y = sigma_y)
    x<-data$x
    y<-data$y
    x<-as.data.frame(x)
    colnames(x)<-paste("x",1:ncol(x),sep="")
    
    #quantile(cor(y[,1],x[,1:50]))
    #MANOVA
    pvalue.manova<-rep(NA,ncol(x))
    for(i in 1:ncol(x)){
      mod<-manova(y~x[,i])
      mod1<-summary(mod, test="Wilks",intercept = T)
      pvalue.manova[i]<-mod1$stats[2,6]
    }
    #quantile(pvalue.manova)
    res.origin<-Getpvalue(X = as.matrix(x),Y = y,Z = NULL)
    res.perm<-GetPermutationPvalue(as.matrix(x),y,number.perm=500,ncores=1)
    Pvalue.permutation<-NULL
    for(i in 1:500){
      Pvalue.permutation<-rbind(Pvalue.permutation,res.perm$Pvalue[[i]])
    }
    
    #GWAS
    mod.aSPU.ind<-apply(x,2,function(x1){
      GEEaSPU(traits = y,geno = as.matrix(x1),Z=NULL,model = "gaussian",corstr = "independence",n.sim = 1000)[10]
    })
    mod.aSPU.ex<-apply(x,2,function(x1){
      GEEaSPU(traits = y,geno = as.matrix(x1),Z=NULL,model = "gaussian",corstr = "exchangeable",n.sim = 1000)[10]
    })
    
    
    #AFp
    mod.AFp<-AFp(pvalue = res.origin$pvalue,pvalue.perm = Pvalue.permutation)
    #AFp.linear
    #mod.AFp.linear<-AFp.linear(pvalue = res.origin$pvalue,pvalue.perm = Pvalue.permutation)
    #AFz
    mod.AFz<-AFz(pvalue = res.origin$pvalue,pvalue.perm = Pvalue.permutation)
    #AFz.linear
    #mod.AFz.linear<-AFz.linear(pvalue = res.origin$pvalue,pvalue.perm = Pvalue.permutation)
    
    #minp
    minp<-apply(res.origin$pvalue,1,min)
    minp.perm<-apply(Pvalue.permutation,1,min)
    p.minp<-unlist(lapply(minp,function(x){sum(minp.perm<x)/length(minp.perm)}))
    #p.minp<-sum(minp.perm<minp)/100000
    #sum(p.minp<0.05)/length(p.minp)
    #EW Fisher
    EW<-apply(res.origin$pvalue,1,function(x){sum(-log(x))})
    EW.perm<-apply(Pvalue.permutation,1,function(x){sum(-log(x))})
    p.EW<-unlist(lapply(EW,function(x){sum(EW.perm>x)/length(EW.perm)}))
    #sum(p.EW<0.05)/length(p.EW)
    
    # index1<-which(p.minp<0.05)
    # index2<-which(p.EW<0.05)
    # index3<-which(mod.AFp$AFp.pvalue<0.05)
    # index4<-which(pvalue.manova<0.05)
    #p.EW<-sum(EW.perm>EW)/100000
    ### TATES ###
    mod.TATEs <-apply(res.origin$pvalue,1,function(x1){
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
    final.res<-list(manova=pvalue.manova,mod.aSPU.ind=mod.aSPU.ind,mod.aSPU.ex=mod.aSPU.ex,
                    AFp=mod.AFp,AFz=mod.AFz,minp=p.minp,EW=p.EW,mod.TATEs=mod.TATEs,
                    input.stat=res.origin)
    
    # final.res<-list(manova=pvalue.manova,mod.aSPU.ind=mod.aSPU.ind,mod.aSPU.ex=mod.aSPU.ex,
    #                 AFp=mod.AFp,AFp.linear=mod.AFp.linear,AFz=mod.AFz,
    #                 AFz.linear=mod.AFz.linear,minp=p.minp,EW=p.EW,mod.TATEs=mod.TATEs,
    #                 AFp.truc=mod.AFp.truc,AFp.linear.truc=mod.AFp.linear.truc,AFz.truc=mod.AFz.truc,
    #                 AFz.linear.truc=mod.AFz.linear.truc,minp.truc=p.minp.truc,EW.truc=p.EW.truc,
    #                 mod.TATEs.truc=mod.TATEs.truc,
    #                 input.stat=res.origin)
    return(final.res)
  },mc.cores = num.cores)
  file_name<-paste("Sim1D_hierarchical_sigma_mu=",sigma_mu,".Rdata",sep="")
  save(res,file=file_name)
}
  

