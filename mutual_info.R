#######------------------ Mutual information ----------------########
## Data of Plant-AMF spatially explicit community (J. Klironomos)
##------% Monte-Carlo and theoretical estimates %-------##


## Observed mutual information
# Calculate the number of (0,0), (0,1), (1,0) and (1,1)
n=18;m=15;L=51
p11=matrix(0,n,m)
p10=matrix(0,n,m)
p01=matrix(0,n,m)
p00=matrix(0,n,m)

for(i in 1:n){
 for(j in 1:m){
   #p11[i,j]=0;p01[i,j]=0;p10[i.j]=0;p00[i,j]=0;
   for(k in 1:L){
     for(kk in 1:L){
       if(P[[i]][k,kk]==1 && A[[j]][k,kk]==1) p11[i,j]=p11[i,j]+1
       if(P[[i]][k,kk]==1 && A[[j]][k,kk]==0) p10[i,j]=p10[i,j]+1
       if(P[[i]][k,kk]==0 && A[[j]][k,kk]==1) p01[i,j]=p01[i,j]+1
       if(P[[i]][k,kk]==0 && A[[j]][k,kk]==0) p00[i,j]=p00[i,j]+1
     }
   }
 }
}

sumP=rep(0,n); sumA=rep(0,m);
for(i in 1:n) sumP[i]=sum(P[[i]]>0)
for(j in 1:m) sumA[j]=sum(A[[j]]>0)
for(i in 1:n) sumP[i]=sum(P[[i]]>0)/(L*L)
for(j in 1:m) sumA[j]=sum(A[[j]]>0)/(L*L)
# Mutual information

p10=p10/(L*L);
p01=p01/(L*L);
p00=p00/(L*L);
p11=p11/(L*L);
pp1=(p10+p11)/(L*L);
pp0=(p00+p01)/(L*L);
pa1=(p01+p11)/(L*L);
pa0=(p00+p10)/(L*L);
N=L*L
Mutual=matrix(0,n,m);
Mutual2=matrix(0,n,m);
HX = matrix(0,n,m);
HY = matrix(0,n,m);
HXY = matrix(0,n,m);


for(i in 1:n){
 for(j in 1:m){
  Mutual[i,j]= p11[i,j]*log2(N*N*p11[i,j]/(sumP[i]*sumA[j])) + p10[i,j]*log2(N*N*p10[i,j]/(sumP[i]*(N-sumA[j]))) + p01[i,j]*log2(N*N*p01[i,j]/((N-sumP[i])*sumA[j])) + p00[i,j]*log2((N*N*p00[i,j]/((N-sumP[i])*(N-sumA[j]))));
  # Entropy X - plants
  HX[i,j] = - (pp1[i,j]*log2(pp1[i,j]) + pp0[i,j]*log2(pp0[i,j])) 
  # Enropy Y AMF
  HY[i,j] = - (pa1[i,j]*log2(pa1[i,j]) + pa0[i,j]*log2(pa0[i,j])) 
  # Entropy of both
  HXY[i,j] = -(p11[i,j]*log2(p11[i,j]) + p10[i,j]*log2(p10[i,j]) + p01[i,j]*log2(p01[i,j]) + p00[i,j]*log2(p00[i,j]))
  
  Mutual2[i,j] = HX[i,1]+HY[1,j] - HXY[i,j];
 }
}

save(HX,HY,HXY,file="entropy_obs.RData")
save(Mutual,file="mutual_obs.RData")
save(Mutual2,file="mutual_entropy.RData")

bubble<- Mutual2/HXY;
bubble1 <- matrix(bubble,270,1)
bubble3 <- matrix(Mutual2,270,1);

bubble2<-apply(bubble1,1, function(x) if(is.na(x)==1){x=0}else{x=x})

symbols(coords,circles=bubble2,bg=rgb(red=0, green=0, blue=255, alpha=75, max=255), inches=0.1 )

## Symmetric dependence (Dx,y)
n = 18; m = 15;
Dxy = matrix(0,n,m);
Dx = matrix(0,n,m);
Dy = matrix(0,n,m);

for(i in 1:n){
   for(j in 1:m){
       Dxy[i,j] = Mutual[i,j]/(sqrt(HX[i,j]*HY[i,j]));
       Dx[i,j] = Mutual[i,j]/HY[i,j]; # Knowledge of Y tells me about X; i.e. information of AMF species tells me about Plant species
       Dy[i,j] = Mutual[i,j]/HX[i,j]; # Knowledge of X tells me about Y; i.e. information of plant species tells me about AMF species
       if(is.na(Dxy[i,j])==1){
         Dxy[i,j]=0;
       }
       if(is.na(Dx[i,j])==1){
         Dx[i,j]=0;
       }
       if(is.na(Dy[i,j])==1){
         Dy[i,j]=0;
       }
    }
}

asym = abs(Dx - Dy); # Asymmetry of species dependence, new way of estimating asymmetry
save(asym,Dx,Dy,Dxy,file="symmetric_dependence.RData")

## Plot of the interaction matrix between AMF-plants for species mutual dependence.
# it would be interesting to plot Dx|y and Dy|x
library(fields);
library(lattice);
whiteb= colorRampPalette(c("white", "gray", "black"), space = "Lab")
whitered= colorRampPalette(c("white", "orange", "red"), space = "Lab")
blueored= colorRampPalette(c("white", "yellow", "red"), space = "Lab")
brewer.div = colorRampPalette(c("red", "yellow", "blue"), space = "Lab")
## Plot of 4 matrices for mutual information Dx, Dy, Dxy, and asymmetry

levelplot(Dxy,col.regions = whiteb(100), cuts = 99, xlab="", ylab="AMF")
tp2 <- levelplot(Dx,col.regions = whiteb(100), cuts = 99, xlab="Plants", ylab="AMF")
tp3 <- levelplot(Dy,col.regions = whiteb(100), cuts = 99, xlab="Plants", ylab="AMF")
tp4 <- levelplot(asym,col.regions = whiteb(100), cuts = 99, xlab="Plants", ylab="AMF")

tp1 <- levelplot(Dxy,col.regions = whiteb(100), cuts = 99, xlab="", ylab="AMF", scales = list(alternating = 0, at=0))
tp2 <- levelplot(Dx,col.regions = whiteb(100), cuts = 99, xlab="Plants", ylab="AMF", scales = list(alternating = 0, at=0))
tp3 <- levelplot(Dy,col.regions = whiteb(100), cuts = 99, xlab="Plants", ylab="AMF", scales = list(alternating = 0, at=0))
tp4 <- levelplot(asym,col.regions = whiteb(100), cuts = 99, xlab="Plants", ylab="AMF", scales = list(alternating = 0, at=0))


plot(tp1, split = c(1, 1, 2, 2)); 
plot(tp2, split = c(1, 2, 2, 2), newpage = FALSE);
plot(tp3, split = c(2, 1, 2, 2), newpage = FALSE);
plot(tp4, split = c(2, 2, 2, 2), newpage = FALSE);

postscript(file="asym_sym.eps")
jpeg(file="asym_dep.jpg") 
plot(tp1, split = c(1, 1, 1, 2)); 
plot(tp4, split = c(1, 2, 1, 2), newpage = FALSE);
dev.off()

coords=c(0,0);

for(i in 1:n){
  for(j in 1:m){
    if(is.na(Dxy[i,j])==0){
        coords=rbind(coords,c(i,j));
    }
  }
}
coords=coords[-1,]



bubble = matrix(Dxy,270,1);
coor <- expand.grid(Plants,AMF);
symbols(coords,circles=bubble*100,bg=rgb(red=0, green=0, blue=255, alpha=75, max=255), inches=FALSE )


##### Monte Carlo mutual information ----------------------------------------
load("Pnullmatrices.RData");load("anullmatrices.RData");

## Observed mutual information
# Calculate the number of (0,0), (0,1), (1,0) and (1,1)
n=18;m=15;L=51; sims=1000;

p11mc = vector("list", sims);
p10mc = vector("list", sims);
p01mc = vector("list", sims);
p00mc = vector("list", sims);

for(l in 1:sims){
p11mc[[l]]=matrix(0,n,m);
p10mc[[l]]=matrix(0,n,m);
p01mc[[l]]=matrix(0,n,m);
p00mc[[l]]=matrix(0,n,m);
for(i in 1:n){
 for(j in 1:m){
   #p11[i,j]=0;p01[i,j]=0;p10[i.j]=0;p00[i,j]=0;
   for(k in 1:L){
     for(kk in 1:L){
       if(Peq[[i]][[l]][k,kk]==1 && Aeq[[j]][[l]][k,kk]==1) p11mc[[l]][i,j]=p11mc[[l]][i,j]+1
       if(Peq[[i]][[l]][k,kk]==1 && Aeq[[j]][[l]][k,kk]==0) p10mc[[l]][i,j]=p10mc[[l]][i,j]+1
       if(Peq[[i]][[l]][k,kk]==0 && Aeq[[j]][[l]][k,kk]==1) p01mc[[l]][i,j]=p01mc[[l]][i,j]+1
       if(Peq[[i]][[l]][k,kk]==0 && Aeq[[j]][[l]][k,kk]==0) p00mc[[l]][i,j]=p00mc[[l]][i,j]+1
     }
   }
 }
}
}

NP <- vector("list",n);
Na <- vector("list",m);
for(i in 1:n) NP[i] <- sapply(Peq[[i]], function(x) sum(x>0)/(L*L));
for(j in 1:m) Na[j] <- sapply(Aeq[[j]], function(x) sum(x>0)/(L*L));

# Mutual information
p11mc2 <- sapply(p11mc,function(x) x/(51*51))
p10mc2 <- sapply(p10mc,function(x) x/(51*51))
p01mc2 <- sapply(p01mc,function(x) x/(51*51))
p00mc2 <- sapply(p00mc,function(x) x/(51*51))

Mutual=matrix(0,n,m);

for(i in 1:n){
 for(j in 1:m){
  Mutual[i,j]= p11[i,j]*log(p11[i,j]/(sumP[i]*sumA[j])) + p10[i,j]*log(p10[i,j]/(sumP[i]*sumA[j])) + p01[i,j]*log(p01[i,j]/(sumP[i]*sumA[j])) + p00[i,j]*log(p00[i,j]/(sumP[i]*sumA[j]));
 }
}

#####--- Theoretical prediction based on complete randomness

## Functions to deal with overflows calculating factorials of numbers > 100

stirling <- function(n) sqrt(2*pi*n)*(n/exp(1))^n 
stirling_log <- function(n) n*log(n)-n+(0.5*log(2*pi*n)) # approximation of Log(n!), only for large numbers
# factorial() only works until factorial(170)
#log_comb <- function(n,r) log(factorial(n))- ( (log(factorial(n-r)) + log(factorial(r))))
log_comb_stirling <- function(n,r) stirling_log(n) - ( (stirling_log(n-r)) + stirling_log(r))
# if factorial > 100
#comb2 <- function(n,r) stirling(n)/(stirling(n-r)*stirling(r))
# if factorial < 100
#comb <- function(n,r) factorial(n)/(factorial(n-r)*factorial(r));

#functions factorials
stirling_log <- function(n) n*log(n)-n+(0.5*log(2*pi*n)) 
log_comb_stirling <- function(n,r) stirling_log(n) - ( (stirling_log(n-r)) + stirling_log(r))
L=51;N=L^2; n=18; m=15;
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
Mutual_theo_2 =vector("list",n)
HXTtheo = vector("list",n)

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
 Mutual_theo_2[[i]] =vector("list",m)
 HXTtheo[[i]] = vector("list",m)
 
 for(j in 1:m){ 
  nmin[i,j] = max(0,sumP[i]+sumA[j]-L*L);
  nmax[i,j] = min(sumP[i],sumA[j]);
  n11[[i]][[j]] =  seq(nmin[i,j],nmax[i,j]);
  for(k in 1:length(n11)){
   if(sumP[i] > sumA[j]){
      Pn11[[i]][[j]] = sapply(n11[[i]][[j]], function(x)  exp(log_comb_stirling(sumP[i],x) + log_comb_stirling(N-sumP[i],sumA[j]-x) - (log_comb_stirling(N,sumA[j]))))
   }else{
      Pn11[[i]][[j]] = sapply(n11[[i]][[j]], function(x)  exp(log_comb_stirling(sumA[j],x) + log_comb_stirling(N-sumA[j],sumP[i]-x) - (log_comb_stirling(N,sumP[i]))))
   }
   # Mutual information
   pn11[[i]][[j]]=sapply(n11[[i]][[j]], function(x) x/N);
   pn10[[i]][[j]]=sapply(n11[[i]][[j]], function(x) (sumP[i] - x)/N);
   pn01[[i]][[j]]=sapply(n11[[i]][[j]], function(x) (sumA[j] - x)/N);
   pn00[[i]][[j]]= 1 - pn01[[i]][[j]] - pn10[[i]][[j]] - pn11[[i]][[j]];
   
   Mutual_theo[[i]][[j]]= pn11[[i]][[j]]*log2(N*N*pn11[[i]][[j]]/(sumP[i]*sumA[j])) + pn10[[i]][[j]]*log2(N*N*pn10[[i]][[j]]/(sumP[i]*(N-sumA[j]))) + pn01[[i]][[j]]*log2(N*N*pn01[[i]][[j]]/((N-sumP[i])*sumA[j])) + pn00[[i]][[j]]*log2(N*N*pn00[[i]][[j]]/((N-sumP[i])*(N-sumA[j])));
   
   HXTtheo[[i]][[j]] = - ( pn11[[i]][[j]]*log2(pn11[[i]][[j]]) + pn10[[i]][[j]]*log2(pn10[[i]][[j]]) + pn01[[i]][[j]]*log2(pn01[[i]][[j]]) + pn00[[i]][[j]]*log2(pn00[[i]][[j]]));
   
   #Mutual_theo_2[[i]][[j]] = HX[i] + HY[j] - HXTtheo[[i]][[j]];
   
 }
}
}

save(Pn11,file="prob_n11.RData")
save(pn11,pn10,pn01,pn11,file="prob_mutual.RData")
save(Mutual_theo,file="mutual_theo.RData")
save(Mutual_theo_2,HXTtheo,file="mutual1_entropy.RData")

mean_pn11_2=matrix(0,n,m)
mean_pn11=matrix(0,n,m)
var_pn11=matrix(0,n,m)
##  sum over all possibles n11 to obtain average value from expected random overlap
Pn11=vector("list",n)
Mutual_theo2=vector("list",n)

for(i in 1:n){
 for(j in 1:m){
    Pn11[[i]][[j]] <- sapply(Pn11[[i]][[j]], function(x) if(is.nan(x)==1){x=0}else{x=x} );
    Mutual_theo[[i]][[j]] <- sapply(Mutual_theo[[i]][[j]], function(x) if(is.nan(x)==1){x=0}else{x=x});
    #Mutual_theo_2[[i]][[j]] <- sapply(Mutual_theo_2[[i]][[j]], function(x) if(is.nan(x)==1){x=0}else{x=x});
     #mean_pn11_2[i,j] <-sum(Pn11[[i]][[j]]*Mutual_theo_2[[i]][[j]]);
    mean_pn11[i,j] <-sum(Pn11[[i]][[j]]*Mutual_theo[[i]][[j]]);
    var_pn11[i,j]<- sum( (Pn11[[i]][[j]]*Mutual_theo[[i]][[j]] -  mean_pn11[i,j])^2 )
 }
}

save(mean_pn11,file="mean_mutual_info_theo.RData")

# Expected Variance: sum((Pn11[[i]][[j]]*Mutual_theo[[i]][[j]] -  mean_pn11[i,j])^2)

pmatrix <- function(N1,N2,N,n11){
  
  n01=N2-n11;
  n10=N1-n11;
  n00=N-n10-n01-n11;
  pP1=(n11+n10)/N;
  pP0=(n01+n00)/N;
  pA1=(n11+n01)/N;
  pA0=(n10+n00)/N;
  
  return(c(pP1,pP0,pA1,pA0));
}

test <-seq(1,2601);
N=2601;

testx <- seq(1,100);
test2 <- sapply(test, function(x) log_comb_stirling(N,x))
test3 <- sapply(testx, function(x) factorial(100)/(factorial(100-x)*factorial(x)))
