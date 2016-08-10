#####--- Theoretical prediction based on complete randomness

## Functions to deal with overflows calculating factorials of numbers > 100


#functions factorials
stirling_log <- function(n) n*log(n)-n+(0.5*log(2*pi*n)) # approximation of Log(n!), only for large numbers
log_comb_stirling <- function(n,r) stirling_log(n) - ( (stirling_log(n-r)) + stirling_log(r))
L=51;N=L^2;
#  Pn11[i,j,k] = (log_comb_stirling(sumA[j],n11[k])*log_comb_stirling(N-sumA[j],sumP[i]-n11[k]))/(log_comb_stirling(N,sumP[i])*log_comb_stirling(N,sumA[j]));
Pn11= vector("list",n)
n11 = vector("list",n)
pn11 = vector("list",n)
pn10 = vector("list",n)
pn01 = vector("list",n)
pn00 = vector("list",n)
nmin=matrix(0,n,m);
nmax=matrix(0,n,m);
Mutual_theo = vector("list",n)

sumP <- sapply(P, function(x) sum(x>0))
sumA <- sapply(A, function(x) sum(x>0))

for(i in 1:n){
 Pn11[[i]] = vector("list",m);
 n11[[i]] = vector("list",m);
 pn11[[i]] = vector("list",m);
 pn10[[i]] = vector("list",m);
 pn01[[i]] = vector("list",m);
 pn00[[i]] = vector("list",m);
 Mutual_theo[[i]] = vector("list",m)
 for(j in 1:m){ 
  nmin[i,j] = max(0,sumP[i]+sumA[j]-L*L);
  nmax[i,j] = min(sumP[i],sumA[j]);
  n11[[i]][[j]] =  seq(nmin[i,j],nmax[i,j]);
  for(k in 1:length(n11)){
   if(sumP[i] > sumA[j]){
      Pn11[[i]][[j]] = sapply(n11[[i]][[j]], function(x)  (log_comb_stirling(sumP[i],x)*log_comb_stirling(N-sumP[i],sumA[j]-x))/(log_comb_stirling(N,sumP[i])*log_comb_stirling(N,sumA[j])))
   }else{
      Pn11[[i]][[j]] = sapply(n11[[i]][[j]], function(x)  (log_comb_stirling(sumA[j],x)*log_comb_stirling(N-sumA[j],sumP[i]-x))/(log_comb_stirling(N,sumP[i])*log_comb_stirling(N,sumA[j])))
   }
   # Mutual information
   pn11[[i]][[j]]=sapply(n11[[i]][[j]], function(x) x/N);
   pn10[[i]][[j]]=sapply(n11[[i]][[j]], function(x) (sumP[i] - x)/N);
   pn01[[i]][[j]]=sapply(n11[[i]][[j]], function(x) (sumA[j] - x)/N);
   pn00[[i]][[j]]= 1 - pn01[[i]][[j]] - pn10[[i]][[j]] - pn11[[i]][[j]];
   Mutual_theo[[i]][[j]]= pn11[[i]][[j]]*log(N*N*pn11[[i]][[j]]/(sumP[i]*sumA[j])) + pn10[[i]][[j]]*log(N*N*pn10[[i]][[j]]/(sumP[i]*(N-sumA[j]))) + pn01[[i]][[j]]*log(N*N*pn01[[i]][[j]]/((N-sumP[i])*sumA[j])) + pn00[[i]][[j]]*log(N*N*pn00[[i]][[j]]/((N-sumP[i])*(N-sumA[j])));
 }
}
}

save(Pn11,file="prob_n11.RData")
save(pn11,pn10,pn01,pn11,file="prob_mutual.RData")
save(Mutual_theo,file="mutual_theo.RData")



mean_pn11=matrix(0,n,m)
var_pn11=matrix(0,n,m)
##  sum over all possibles n11 to obtain average value from expected random overlap
Pn11_2=vector("list",n)
Mutual_theo2=vector("list",n)

for(i in 1:n){
 for(j in 1:m){
    Pn11[[i]][[j]] <- sapply(Pn11[[i]][[j]], function(x) if(is.nan(x)==1){x=0}else{x=x} );
    Mutual_theo[[i]][[j]] <- sapply(Mutual_theo[[i]][[j]], function(x) if(is.nan(x)==1){x=0}else{x=x});
    mean_pn11[i,j] <-sum(Pn11[[i]][[j]]*Mutual_theo[[i]][[j]])
    var_pn11[i,j]<- sum( (Pn11[[i]][[j]]*Mutual_theo[[i]][[j]] -  mean_pn11[i,j])^2 )
 }
}

save(mean_pn11,file="mean_mutual_info_theo.RData")
