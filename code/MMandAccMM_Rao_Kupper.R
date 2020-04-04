###MM algorithm for Rao_Kupper. c=theta
###s is comparison matrix with s_ij as the number of observed paired comparisons of items i and j such that either i wins j or there is tie outcome.
Bayesian.MM=function(init,tiny,alpha,beta,s,c){
  n=dim(s)[1]
  iterationmatrix=matrix(0,100000,n)
  t=0
  x1=init
  x0=x1+rep(10*tiny,n)
  while(sum(abs(x0-x1)>tiny)>0) {
    x0=x1
    for(i in 1:n){
      u=c()
      l=c()
      for(j in 1:n){
        if(j != i){
          u[j]=s[i,j]
          l[j]=((s[i,j])/(exp(x0[i])+c*exp(x0[j])))+((c*s[j,i])/(exp(x0[j])+c*exp(x0[i])))
        }
        else{u[j]=0
        l[j]=0
        }
      }
      x1[i]=log((alpha-1+sum(u))/(beta+sum(l)))
    }
    t=t+1
    iterationmatrix[t,]=x1
  }
  iterationmatrix=iterationmatrix[1:t,1:n]
  return(list(Iterationmatrix=iterationmatrix))
}

###AccMM algorithm for Rao_Kupper. c=theta
###s is comparison matrix with s_ij as the number of observed paired comparisons of items i and j such that either i wins j or there is tie outcome.
Acc.Bayesian.MM=function(init,tiny,alpha,beta,s,c){
  n=dim(s)[1]
  iterationmatrix=matrix(0,100000,n)
  t=0
  x1=init
  x0=x1+rep(10*tiny,n)
  while(sum(abs(x0-x1)>tiny)>0) {
    x0=x1
    for(i in 1:n){
      u=c()
      l=c()
      for(j in 1:n){
        if(j != i){
          u[j]=s[i,j]
          l[j]=((s[i,j])/(exp(x0[i])+c*exp(x0[j])))+((c*s[j,i])/(exp(x0[j])+c*exp(x0[i])))
        }
        else{u[j]=0
        l[j]=0
        }
      }
      x1[i]=log((alpha-1+sum(u))/(beta+sum(l)))
    }
    c1=log(n*(alpha-1)/beta)-log(sum(exp(x1)))
    x1=x1+c1
    t=t+1
    iterationmatrix[t,]=x1
  }
  iterationmatrix=iterationmatrix[1:t,1:n]
  return(list(Iterationmatrix=iterationmatrix))
}

#####################################################################################################################
#value settings for experiments MM and AccMM:

init=rep(0,n)
tiny=10^(-4)
c=sqrt(2)

#only for sample case
#s=s0
#n=n0

ImatrixMM=Bayesian.MM(init,tiny,11,10,s,c)$Iterationmatrix
ImatrixMM.21=Bayesian.MM(init,tiny,2,1,s,c)$Iterationmatrix
ImatrixMM.1=Bayesian.MM(init,tiny,1.1,0.1,s,c)$Iterationmatrix
ImatrixMM.01=Bayesian.MM(init,tiny,1.01,0.01,s,c)$Iterationmatrix
ImatrixMM0=Bayesian.MM(init,tiny,1,0,s,c)$Iterationmatrix

dim(ImatrixMM)
dim(ImatrixMM.21)
dim(ImatrixMM.1)
dim(ImatrixMM.01)

dim(ImatrixMM0)

AccImatrixMM=Acc.Bayesian.MM(init,tiny,11,10,s,c)$Iterationmatrix
AccImatrixMM.21=Acc.Bayesian.MM(init,tiny,2,1,s,c)$Iterationmatrix
AccImatrixMM.1=Acc.Bayesian.MM(init,tiny,1.1,0.1,s,c)$Iterationmatrix
AccImatrixMM.01=Acc.Bayesian.MM(init,tiny,1.01,0.01,s,c)$Iterationmatrix

dim(AccImatrixMM)
dim(AccImatrixMM.21)
dim(AccImatrixMM.1)
dim(AccImatrixMM.01)

