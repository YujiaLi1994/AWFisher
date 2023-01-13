#setwd("~/Desktop/AWFisher_02_12_2021/real data")
rm(list=ls())
library(stringr)
#out.dir="C:/Users/acdsee8602/iCloudDrive/AWFISHER_THEORY/GENES/code/real_data/code_result"
out.dir="/home/yuf31/Genes"
setwd(out.dir)

load("lung_disease_bootstrap.Rdata")
load("lung_disease_result_all.Rdata")
variability.index.AFp<-matrix(NA,nrow(wcs[[1]]$AFp$W),5)
temp<-c()
for(i in 1:nrow(variability.index.AFp)){
  #print(i)
  for(j in 1:ncol(variability.index.AFp)){
    for(k in (1:50)){
      temp[k]<-wcs[[k]]$AFp$W[i,j]
    }
    variability.index.AFp[i,j]<-4*(var(temp)*49/50)
  }
  
  
}

variability.index.AFp.linear<-matrix(NA,nrow(wcs[[1]]$AFp.linear$W),5)
temp<-c()
for(i in 1:nrow(variability.index.AFp.linear)){
  #print(i)
  for(j in 1:ncol(variability.index.AFp.linear)){
    for(k in (1:50)){
      temp[k]<-wcs[[k]]$AFp.linear$W[i,j]
    }
    variability.index.AFp.linear[i,j]<-4*(var(temp)*49/50)
  }
  
  
}

variability.index.AFz<-matrix(NA,nrow(wcs[[1]]$AFz$W),5)
temp<-c()
for(i in 1:nrow(variability.index.AFz)){
 # print(i)
  for(j in 1:ncol(variability.index.AFz)){
    for(k in (1:50)){
      temp[k]<-wcs[[k]]$AFz$W[i,j]
    }
    variability.index.AFz[i,j]<-4*(var(temp)*49/50)
  }
  
  
}

variability.index.AFz.linear<-matrix(NA,nrow(wcs[[1]]$AFz.linear$W),5)
temp<-c()
for(i in 1:nrow(variability.index.AFz.linear)){
  print(i)
  for(j in 1:ncol(variability.index.AFz.linear)){
    for(k in (1:50)){
      temp[k]<-wcs[[k]]$AFz.linear$W[i,j]
    }
    variability.index.AFz.linear[i,j]<-4*(var(temp)*49/50)
  }
  
  
}

variability.index<-list(AFp=variability.index.AFp,AFp.linear=variability.index.AFp.linear,
                        AFz=variability.index.AFz,AFz.linear=variability.index.AFz.linear)
save(variability.index,file="Varibility_index.Rdata")
#--------------------------------------
#Caleb's way to categorization
#--------------------------------------
Distance.matrix.AFp<-matrix(0,nrow(wcs[[1]]$AFp$W),nrow(wcs[[1]]$AFp$W))
for(i in 1:50){
  print(i)
  temp<-wcs[[i]]
  matrix1<-temp$AFp$W*sign(temp$input.stat$coef)
  res.string<-apply(matrix1,1,function(x){
    str_c(x,collapse = "")
  })
  res.comember<-sapply(res.string,function(x){
    as.numeric(x==res.string)
  })
  Distance.matrix.AFp<-Distance.matrix.AFp+res.comember
}
Distance.matrix.AFp<-Distance.matrix.AFp/50
#Distance.matrix[1:5,1:5]  

Distance.matrix.AFp.linear<-matrix(0,nrow(wcs[[1]]$AFp.linear$W),nrow(wcs[[1]]$AFp.linear$W))
for(i in 1:50){
  print(i)
  temp<-wcs[[i]]
  matrix1<-temp$AFp.linear$W*sign(temp$input.stat$coef)
  res.string<-apply(matrix1,1,function(x){
    str_c(x,collapse = "")
  })
  res.comember<-sapply(res.string,function(x){
    as.numeric(x==res.string)
  })
  Distance.matrix.AFp.linear<-Distance.matrix.AFp.linear+res.comember
}
Distance.matrix.AFp.linear<-Distance.matrix.AFp.linear/50


Distance.matrix.AFz<-matrix(0,nrow(wcs[[1]]$AFz$W),nrow(wcs[[1]]$AFz$W))
for(i in 1:50){
  print(i)
  temp<-wcs[[i]]
  matrix1<-temp$AFz$W*sign(temp$input.stat$coef)
  res.string<-apply(matrix1,1,function(x){
    str_c(x,collapse = "")
  })
  res.comember<-sapply(res.string,function(x){
    as.numeric(x==res.string)
  })
  Distance.matrix.AFz<-Distance.matrix.AFz+res.comember
}
Distance.matrix.AFz<-Distance.matrix.AFz/50

#rm(Distance.matrix.AFp)
#rm(Distance.matrix.AFp.linear)
#rm(Distance.matrix.AFz)

Distance.matrix.AFz.linear<-matrix(0,nrow(wcs[[1]]$AFz.linear$W),nrow(wcs[[1]]$AFz.linear$W))
for(i in 1:50){
  print(i)
  temp<-wcs[[i]]
  matrix1<-temp$AFz.linear$W*sign(temp$input.stat$coef)
  res.string<-apply(matrix1,1,function(x){
    str_c(x,collapse = "")
  })
  res.comember<-sapply(res.string,function(x){
    as.numeric(x==res.string)
  })
  Distance.matrix.AFz.linear<-Distance.matrix.AFz.linear+res.comember
}
Distance.matrix.AFz.linear<-Distance.matrix.AFz.linear/50

save(Distance.matrix.AFp,file="comembership_Caleb.AFp.Rdata") 
save(Distance.matrix.AFp.linear,file="comembership_Caleb.AFp.linear.Rdata") 
save(Distance.matrix.AFz,file="comembership_Caleb.AFz.Rdata") 
save(Distance.matrix.AFz.linear,file="comembership_Caleb.AFz.linear.Rdata") 




#-------------------------------------------------------
#Way 1 to clustering.
#--------------------------------------
Distance.all<-matrix(NA,nrow(wcs[[1]]$weight.matrix.bootstrap),10)
Distance<-lapply(1:5,function(j){
  Distance.temp<-matrix(NA,nrow(wcs[[1]]$weight.matrix.bootstrap),2)
  temp<-c()
  for(i in 1:nrow(Distance.temp)){
    print(i)
    
    
    for(k in (1:50)){
      temp[k]<-(wcs[[k]]$weight.matrix.bootstrap[i,j])*(sign(wcs[[k]]$coef.matrix[i,j]))
      Distance.temp[i,1]<-sum(temp==1)
      Distance.temp[i,2]<-sum(temp==-1)
      #Distance.temp[i,3]<-sum(temp==0)
    }
    #Distance.all<-cbind(Distance.all,Distance.temp)
  }
  return(Distance.temp)
})


Distance.all[,1:2]<-Distance[[1]]
Distance.all[,3:4]<-Distance[[2]]
Distance.all[,5:6]<-Distance[[3]]
Distance.all[,7:8]<-Distance[[4]]
Distance.all[,9:10]<-Distance[[5]]
Distance.all[1:5,]
save(Distance.all,file="Distance_1_Yujia.Rdata")