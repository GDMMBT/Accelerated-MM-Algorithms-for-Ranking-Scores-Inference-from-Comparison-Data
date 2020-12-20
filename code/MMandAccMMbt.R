##################################################################################
###MM algorithm for Bradley-Terry. w is the winning matrix.

Bayesian.MM=function(init,tiny,alpha,beta,w){
  n=dim(w)[1]
  iterationmatrix=matrix(0,100000,n)
  i=0
  x1=init
  x0=x1+rep(10*tiny,length(x1))
  u=c()
  l=c()
  while(sum(abs(x0-x1)>tiny)>0) {
    x0=x1
    for(j in 1:n){
      for(k in 1:n){
        if(k != j){
          u[k]=w[j,k]
          l[k]=(w[j,k]+w[k,j])/(exp(x0[j])+exp(x0[k]))
        }
        else{u[k]=0
        l[k]=0
        }
      }
      x1[j]=log((alpha-1+sum(u))/(beta+sum(l)))
    }
    i=i+1
    iterationmatrix[i,]=x1
  }
  iterationmatrix=iterationmatrix[1:i,1:n]
  return(list(Iterationmatrix=iterationmatrix))
}

##################################################################################
###MM algorithm for Bradley-Terry with normalisation at each step such that \sum_{i=1}^n exp(w_i)=1.

Bayesian.MM.Norm=function(init,tiny,alpha,beta,w){
  n=dim(w)[1]
  iterationmatrix=matrix(0,100000,n)
  i=0
  x1=init
  x0=x1+rep(10*tiny,length(x1))
  u=c()
  l=c()
  while(sum(abs(x0-x1)>tiny)>0 && i<3000) {
    x0=x1
    for(j in 1:n){
      for(k in 1:n){
        if(k != j){
          u[k]=w[j,k]
          l[k]=(w[j,k]+w[k,j])/(exp(x0[j])+exp(x0[k]))
        }
        else{u[k]=0
        l[k]=0
        }
      }
      x1[j]=log((alpha-1+sum(u))/(beta+sum(l)))
    }
    i=i+1
    x1=log(exp(x1)/sum(exp(x1)))  
    iterationmatrix[i,]=x1
  }
  iterationmatrix=iterationmatrix[1:i,1:n]
  return(list(Iterationmatrix=iterationmatrix))
}

##################################################################################
###AccMM algorithm for Bradley-Terry. w is the winning matrix.

Acc.Bayesian.MM=function(init,tiny,alpha,beta,w){
  n=dim(w)[1]
  iterationmatrix=matrix(0,100000,n)
  i=0
  x1=init
  x0=x1+rep(10*tiny,length(x1))
  u=c()
  l=c()
  while(sum(abs(x0-x1)>tiny)>0) {
    x0=x1
    for(j in 1:n){
      for(k in 1:n){
        if(k != j){
          u[k]=w[j,k]
          l[k]=(w[j,k]+w[k,j])/(exp(x0[j])+exp(x0[k]))
        }
        else{u[k]=0
        l[k]=0
        }
      }
      x1[j]=log((alpha-1+sum(u))/(beta+sum(l)))
    }
    c=log(n*(alpha-1)/beta)-log(sum(exp(x1)))
    x1=x1+c
    i=i+1
    iterationmatrix[i,]=x1
  }
  iterationmatrix=iterationmatrix[1:i,1:n]
  return(list(Iterationmatrix=iterationmatrix))
}

##################################################################################

#value settings for experiments MM and AccMM:

init=rep(0,n)
tiny=10^-4

#only for sample case
#w=w0
#n=n0 
#m.max=mm.max

ImatrixMM=Bayesian.MM(init,tiny,11,10,w)$Iterationmatrix
ImatrixMM.21=Bayesian.MM(init,tiny,2,1,w)$Iterationmatrix
ImatrixMM.1=Bayesian.MM(init,tiny,1.1,0.1,w)$Iterationmatrix
ImatrixMM.01=Bayesian.MM(init,tiny,1.01,0.01,w)$Iterationmatrix
ImatrixMM0=Bayesian.MM(init,tiny,1,0,w)$Iterationmatrix

dim(ImatrixMM)
dim(ImatrixMM.21)
dim(ImatrixMM.1)
dim(ImatrixMM.01)
dim(ImatrixMM0)

##################################################################################
####For GIFGIF Sample datasets
ImatrixMM.Norm=Bayesian.MM.Norm(init,tiny,11,10,w)$Iterationmatrix
ImatrixMM.Norm.21=Bayesian.MM.Norm(init,tiny,2,1,w)$Iterationmatrix
ImatrixMM.Norm.1=Bayesian.MM.Norm(init,tiny,1.1,0.1,w)$Iterationmatrix
ImatrixMM.Norm.01=Bayesian.MM.Norm(init,tiny,1.01,0.01,w)$Iterationmatrix
ImatrixMM.Norm0=Bayesian.MM.Norm(init,tiny,1,0,w)$Iterationmatrix

dim(ImatrixMM.Norm)
dim(ImatrixMM.Norm.21)
dim(ImatrixMM.Norm.1)
dim(ImatrixMM.Norm.01)
dim(ImatrixMM.Norm0)

##################################################################################
AccImatrixMM=Acc.Bayesian.MM(init,tiny,11,10,w)$Iterationmatrix
AccImatrixMM.21=Acc.Bayesian.MM(init,tiny,2,1,w)$Iterationmatrix
AccImatrixMM.1=Acc.Bayesian.MM(init,tiny,1.1,0.1,w)$Iterationmatrix
AccImatrixMM.01=Acc.Bayesian.MM(init,tiny,1.01,0.01,w)$Iterationmatrix

dim(AccImatrixMM)
dim(AccImatrixMM.21)
dim(AccImatrixMM.1)
dim(AccImatrixMM.01)
