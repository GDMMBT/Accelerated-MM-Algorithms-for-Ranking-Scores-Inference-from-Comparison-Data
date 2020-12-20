###We use following function to compute the log aposteriori probablity. Input x is parameter value w.

log.aposteriori.prob=function(x){
  l=c()
  l1=c()
  for(i in 1:n){
    for(j in 1:n){
      if (j!=i) {l1[j]=w1[i,j]*log(exp(x[i])/(exp(x[i])+exp(x[j])))}
      else{l1[j]=0}
    }
    l[i]=sum(l1)+(alpha-1)*x[i]-beta*exp(x[i])
  }
  return(sum(l))
}