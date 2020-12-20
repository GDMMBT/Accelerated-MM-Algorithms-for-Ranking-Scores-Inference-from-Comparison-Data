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