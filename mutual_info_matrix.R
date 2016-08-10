###### Calculating confidence intervals for each interaction ##
## This script use the theoretical probabilities calculated by mutual_info.R of spatial overlap between AMF-Plants

## alpha = 0.05, lower_tail = 0.025; upper_tail = 0.975 of the cumulative prob distribution

n=18;m=15; N=51*51;

load("prob_n11.RData") #Pn11
load("prob_mutual.RData") #pn11, etc
load("n11_theo.RData")
for(i in 1:n){
 for(j in 1:m){
    Pn11[[i]][[j]] <- sapply(Pn11[[i]][[j]], function(x) if(is.nan(x)==1){x=0}else{x=x} );
  }
}

# First obtain cumsum of prob functions
cum_pn11 = vector("list",n)

for(i in 1:n){
  cum_pn11[[i]]=vector("list",m)
  for(j in 1:m){
     cum_pn11[[i]][[j]] = cumsum(Pn11[[i]][[j]]);
  }
}



# Second, find overlap values for lower_tail and upper_tail of each function (x| N1, N2, N) and sum prob

probs_low=vector("list",n)
probs_upper=vector("list",n)
for(i in 1:n){
  probs_low[[i]]=vector("list",m)
  probs_upper[[i]]=vector("list",m)
  for(j in 1:m){
   for(k in 1:length(cum_pn11[[i]][[j]])){
    if(cum_pn11[[i]][[j]][k] < 0.025){
       probs_low[[i]][[j]] = c(probs_low[[i]][[j]],k);
    }
    if(cum_pn11[[i]][[j]][k] > 0.975){
      probs_upper[[i]][[j]] = c(probs_upper[[i]][[j]],k);
    }
   }
  }
}
min_overlap=matrix(0,n,m)
max_overlap=matrix(0,n,m)
for(i in 1:n){
  for(j in 1:m){
     if(is.null(probs_low[[i]][[j]])==1){ probs_low[[i]][[j]]==0;}
     if(is.null(probs_upper[[i]][[j]])==1){ probs_upper[[i]][[j]]==length(probs_upper[[i]][[j]]);}
     min_overlap[i,j]=max(probs_low[[i]][[j]])
     max_overlap[i,j]=min(probs_upper[[i]][[j]])
  }
}
save(max_overlap,min_overlap,file="maxmin_overlaps_theo.RData")
# Third Evaluate whether obs < lower_tail or obs > upper_tail
# Observed overlaps
n11_obs=matrix(0,n,m)
 for(i in 1:n){
    for(j in 1:m){
        n11_obs[i,j] = sum(S[,,i,j]);
    }
 }

save(n11_obs,file="n11_observed.RData")

## Compare observed versus theoretical quantiles
for(i in 1:n){
  for(j in 1:m){
    if(is.finite(max_overlap[i,j])==0){ 
       max_overlap[i,j]=length(cum_pn11[[i]][[j]])
      }
  }
}
 
inter_matrix=matrix(0,n,m)
mutualism=matrix(0,n,m)
seggregation=matrix(0,n,m)

for(i in 1:n){
  for(j in 1:m){
    if(n11_obs[i,j]<min_overlap[i,j]){
        inter_matrix[i,j]= 1;
        seggregation[i,j]=1;
    }else if(n11_obs[i,j]>max_overlap[i,j]){
        inter_matrix[i,j]= 2; 
        mutualism[i,j]=1; 
    }else{
        inter_matrix[i,j]= 0;
    }
  }
}
 
library(bipartite)
nested(mutualism,method="NODF") # 8.3
nested(seggregation,method="NODF") # 10.06
  save(inter_matrix,mutualism,seggregation,file="matrices_sig_MI_theo.RData")
  
sortweb(mutualism,sort.order="dec")

