#################################################
# Covariance-on-covariance regression
#################################################

library("mvtnorm")

rm(list=ls())

source("COCReg.R")

#################################################
# simulation parameter setting
p<-10
q<-5

# X
set.seed(100)
gamma.mat0<-matrix(runif(p),nrow=p,ncol=p)
gamma.mat<-qr.Q(qr(gamma.mat0))
for(j in 1:p)
{
  if(gamma.mat[which.max(abs(gamma.mat[,j])),j]<0)
  {
    gamma.mat[,j]<-(-gamma.mat[,j])
  }
}
Gamma1<-gamma.mat
# t(Gamma1)%*%Gamma1
# Gamma1%*%t(Gamma1)

x.eigen.m<-exp(seq(1,-2,length.out=p))
x.eigen.sd<-0.5

# Y
set.seed(500)
gamma.mat0<-matrix(runif(q),nrow=q,ncol=q)
gamma.mat<-qr.Q(qr(gamma.mat0))
for(j in 1:q)
{
  if(gamma.mat[which.max(abs(gamma.mat[,j])),j]<0)
  {
    gamma.mat[,j]<-(-gamma.mat[,j])
  }
}
Gamma2<-gamma.mat
# t(Gamma2)%*%Gamma2
# Gamma2%*%t(Gamma2)

y.eigen.m<-exp(seq(1,-2,length.out=q))
y.eigen.sd<-0.5

# index of correlated eigenvalues
x.idx<-c(1,3)
y.idx<-c(2,4)

# regression model parameters
alpha<-c(3,2)
beta<-cbind(c(1,-1),c(-1,1))
colnames(beta)<-paste0("C",1:ncol(beta))
#################################################

#################################################
# generate data
n<-100
nTx<-100
nTy<-100

uvec<-rep(nTx,n)      # # of time points of X
vvec<-rep(nTy,n)      # # of time points of Y

set.seed(200)

# generate data
# covariates
W<-cbind(rep(1,n),rbinom(n,size=1,prob=0.5))

# covariance matrices
Delta<-array(NA,c(p,p,n))
Sigma<-array(NA,c(q,q,n))
Lambda1<-matrix(NA,n,p)             # eigenvalues of X
for(j in 1:p)
{
  Lambda1[,j]<-exp(rnorm(n=n,mean=log(x.eigen.m[j]),sd=x.eigen.sd))
}
Lambda2<-matrix(NA,n,q)             # eigenvalues of Y
for(k in 1:q)
{
  ftmp<-which(y.idx==k)
  if(length(ftmp)>0)
  {
    Lambda2[,k]<-exp(alpha[ftmp]*log(Lambda1[,x.idx[ftmp]])+W%*%beta[,ftmp])
  }else
  {
    Lambda2[,k]<-exp(rnorm(n=n,mean=log(y.eigen.m[k]),sd=y.eigen.sd))
  }
}
for(i in 1:n)
{
  # X
  Delta[,,i]<-Gamma1%*%diag(Lambda1[i,])%*%t(Gamma1)
  
  # Y
  Sigma[,,i]<-Gamma2%*%diag(Lambda2[i,])%*%t(Gamma2)
}

# X
X<-vector("list",length=n)
for(i in 1:n)
{
  X[[i]]<-rmvnorm(n=uvec[i],mean=rep(0,p),sigma=Delta[,,i])
}

# Y
Y<-vector("list",length=n)
for(i in 1:n)
{
  Y[[i]]<-rmvnorm(n=vvec[i],mean=rep(0,q),sigma=Sigma[,,i])
}
#################################################

#################################################
# method parameters
max.itr<-1000
tol<-1e-4

trace<-TRUE
score.return<-TRUE

Hy<-diag(rep(1,q))
Hx<-diag(rep(1,p))

DfD.thred<-2
nD<-5

verbose<-TRUE

boot<-TRUE
sims<-500
boot.ci.type<-"se"
conf.level<-0.95
outlier.remove<-FALSE
#################################################

#################################################
# run Covariance-on-Covariance regression

# use DfD criterion
re1<-NULL
try(re1<-COCReg(Y,X,W,stop.crt="DfD",DfD.thred=DfD.thred,Hy=Hy,Hx=Hx,max.itr=max.itr,tol=tol,trace=FALSE,score.return=score.return,verbose=verbose))

# specify # of components
re2<-NULL
try(re2<-COCReg(Y,X,W,stop.crt="nD",nD=nD,Hy=Hy,Hx=Hx,max.itr=max.itr,tol=tol,trace=FALSE,score.return=score.return,verbose=verbose))

# inference
re.boot<-vector("list",length=ncol(re1$gamma))
for(jj in 1:ncol(re1$gamma))
{
  try(re.boot[[jj]]<-COCReg.coef.boot(Y,X,W,gamma=re1$gamma[,jj],theta=re1$theta[,jj],boot=boot,sims=sims,boot.ci.type=boot.ci.type,conf.level=conf.level,
                                      max.itr=max.itr,tol=tol,trace=FALSE,score.return=score.return,verbose=verbose))
  
  print(jj)
}
#################################################


