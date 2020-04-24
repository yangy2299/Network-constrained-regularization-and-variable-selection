library(MASS)
library(glmnet)
library(ncvreg)


# define values of beta for 4 models
p=110
n=100
beta1=c(5,rep(5/sqrt(10),10),-5,rep(-5/sqrt(10),10),3,rep(3/sqrt(10),10),
        -3,rep(-3/sqrt(10),10),rep(0,p-44))

beta2=c(5,rep(-5/sqrt(10),3),rep(5/sqrt(10),7),-5,rep(5/sqrt(10),3),
        rep(-5/sqrt(10),7),3,rep(-3/sqrt(10),3),rep(3/sqrt(10),7),-3,
        rep(3/sqrt(10),3),rep(-3/sqrt(10),7),rep(0,p-44))

beta3=c(5,rep(0.5,10),-5,rep(-0.5,10),3,rep(0.3,10),-3,rep(-0.3,10),rep(0,p-44))

beta4=c(5,rep(-0.5,3),rep(0.5,7),-5,rep(0.5,3),rep(-0.5,7),3,rep(-0.3,3),
        rep(0.3,7),-3,rep(0.3,3),rep(-0.3,7),rep(0,p-44))


#define Laplacian matrice of graph
lcell=diag(1,11)
lcell[1,2:11]=-1/sqrt(10)
lcell[2:11,1]=-1/sqrt(10)
l=matrix(0,p,p)
for (i in 1:(p/11)){
  l[(11*(i-1)+1):(i*11),(11*(i-1)+1):(i*11)]=lcell
}


#spectral decomposition of L matrix
s=eigen(l)$vectors%*%diag(sqrt(eigen(l)$values))


# sigma for noise
sig=1


###########################
# function to generate data of 4 models
mod.fcn=function(n,p){
  x=matrix(0,n,p)
  for (i in 1:(p/11)){
    xtf=rnorm(n,0,1)
    x[,11*(i-1)+1]=xtf
    x[,(11*(i-1)+2):(11*i)]=rnorm(10*n,0.7*xtf,0.51)
  }
  er1=rnorm(n,0,sig) 
  er2=rnorm(n,0,sig)
  er3=rnorm(n,0,sig)
  er4=rnorm(n,0,sig)
  y1=x%*%beta1+er1  
  y2=x%*%beta2+er2
  y3=x%*%beta3+er3
  y4=x%*%beta4+er4
  y0=cbind(y1,y2,y3,y4)
  return(list(x=x,y=y0))
}


#  define functions
# True positive
tp.fcn=function(x){
  sum(x[1:44]!=0)/44
}

# True negative
tn.fcn=function(x){
  n=length(x)
  sum(x[45:n]==0)/(n-44)
}


#  CD for LASSO(mode1), elastic net(mode2),network-constrained(mode3)  
s.fcn=function(z,l){
  sign(z)*(abs(z)-l)*((abs(z)-l)>0)
} 

cd.fcn=function(xtr,ytr,xt,yt,lmd1,lmd2,mode){
  if (mode==1){
    p=ncol(xtr)
    n1=nrow(xtr)
    n2=nrow(xt)
    beta = c(rep(0,p))
    beta.next = c(rep(0.5,p))
    z=0
      while (max(abs(beta.next-beta))> 1e-2){
        beta = beta.next  
        for (k in 1:p){
          z = t(ytr-xtr[,-k]%*%beta.next[-k])%*%xtr[,k]/n1
          beta.next[k]=s.fcn(z,lmd1)/(t(xtr[,k])%*%xtr[,k]/n1)
        }
      }
      pred.t = sum((yt-xt%*%beta.next)^2)/n2
      return(list(beta=beta.next,mse=pred.t))
  }else{      
      if (mode==2){
        sm=diag(1,ncol(xtr))
      }else{
        sm=l
      }
      xtr=rbind(xtr,sqrt(lmd2)*sm)/sqrt(1+lmd2)
      ytr=c(ytr,rep(0,ncol(xt)))
      lmd=lmd1/sqrt(1+lmd2)
      
      p=ncol(xtr)
      n1=nrow(xtr)
      n2=nrow(xt)
      beta = c(rep(0,p))
      beta.next = c(rep(0.5,p))
      z=0
      while (max(abs(beta.next-beta))> 1e-2){
        beta = beta.next  
        for (k in 1:p){
          z = t(ytr-xtr[,-k]%*%beta.next[-k])%*%xtr[,k]/n1
          beta.next[k]=s.fcn(z,lmd)/(t(xtr[,k])%*%xtr[,k]/n1)
        }
      }
      betahat=beta.next/sqrt(1+lmd2)
      pred.t = sum((yt-xt%*%betahat)^2)/n2
      return(list(beta=betahat,mse=pred.t))
  }
}  


###########################
lmdv1=seq(0.05,1.05,length=5)
lmdv2=seq(0.1,10.1,length=5)
msel1=rep(0,5)
msel2=matrix(rep(0,25),5)
msel3=matrix(rep(0,25),5)
mat=NULL

for (k in 1:50){ # 50 replicates of data simulation
  set.seed(k)
  data_train=mod.fcn(n,p)
  xtrain=data_train$x
  ytrain0=data_train$y
  xtr=xtrain[1:50,]
  xt=xtrain[51:100,]
  ytr0=ytrain0[1:50,]
  yt0=ytrain0[51:100,]
  
  data_test=mod.fcn(n,p)
  xtest=data_test$x
  ytest0=data_test$y

  for (m in 1:4){      # iterations for 4 models
    ytr=ytr0[,m]
    yt=yt0[,m]
    ytest=ytest0[,m]
    ytrain=ytrain0[,m]
    for (i in 1:5){   # iterations for lambda1
      lmd1=lmdv1[i]
      fit01=cd.fcn(xtr,ytr,xt,yt,lmd1,mode=1)
      msel1[i]=fit01$mse
      for (j in 1:5){ # iterations for 4 lambda2
        lmd2=lmdv2[j]
        fit02=cd.fcn(xtr,ytr,xt,yt,lmd1,lmd2,mode=2)
        fit03=cd.fcn(xtr,ytr,xt,yt,lmd1,lmd2,mode=3)
        msel2[i,j]=fit02$mse
        msel3[i,j]=fit03$mse
      }
    }
    lmd11=lmdv1[which(msel1 == min(msel1))] #choose lambda
    lmd12=lmdv1[which(msel2 == min(msel2), arr.ind=T)[1]]
    lmd22=lmdv2[which(msel2 == min(msel2), arr.ind=T)[2]]
    lmd13=lmdv1[which(msel3 == min(msel3), arr.ind=T)[1]]
    lmd23=lmdv2[which(msel3 == min(msel3), arr.ind=T)[2]]
    fit1=cd.fcn(xtrain,ytrain,xtest,ytest,lmd11,mode=1)# refit data
    fit2=cd.fcn(xtrain,ytrain,xtest,ytest,lmd12,lmd22,mode=2)
    fit3=cd.fcn(xtrain,ytrain,xtest,ytest,lmd13,lmd23,mode=3)
    mate=c(tp.fcn(fit1$beta),tp.fcn(fit2$beta),tp.fcn(fit3$beta),
           tn.fcn(fit1$beta),tn.fcn(fit2$beta),tn.fcn(fit3$beta),
           fit1$mse,fit2$mse,fit3$mse)
    mat=rbind(mat,mate)
  }
}


# rearrange data
mod1=matrix(rep(0,450),50)
mod2=mod1
mod3=mod1
mod4=mod1
table=NULL
for (i in 1:50){
  mod1[i,]=mat[((i-1)*4+1),]
  mod2[i,]=mat[((i-1)*4+2),]
  mod3[i,]=mat[((i-1)*4+3),]
  mod4[i,]=mat[((i-1)*4+4),]
}
row1=colMeans(mod1)
row2=apply(mod1,2,sd)/sqrt(50)
row3=colMeans(mod2)
row4=apply(mod2,2,sd)/sqrt(50)
row5=colMeans(mod3)
row6=apply(mod3,2,sd)/sqrt(50)
row7=colMeans(mod4)
row8=apply(mod4,2,sd)/sqrt(50)
table=rbind(row1,row2,row3,row4,row5,row6,row7,row8)
colnames(table)=c("Ltp","Etp","Ntp","Ltn","Etn","Ntn","Lmse","Emse","Nmse")
rownames(table)=c("m1avg","m1se","m2avg","m2se","m3avg","m3se","m4avg","m4se")
write.csv(table,file = "tab.csv")
table

  