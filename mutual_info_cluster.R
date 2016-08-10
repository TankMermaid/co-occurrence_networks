
##### Monte Carlo mutual information ----------------------------------------
load("Pnullmatrices.RData");load("anullmatrices.RData");

## Observed mutual information
# Calculate the number of (0,0), (0,1), (1,0) and (1,1)
n=18;m=15;L=51; sims=1000;
N=L*L;

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

## more efficient way !
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
            p11mc[[l]][i,j] = Peq[[i]][[l]][k,kk]*Aeq[[j]][[l]][k,kk] + p11mc[[l]][i,j];
     }
    }
    p10mc[[l]][i,j] = NP[[i]][[l]]- p11mc[[l]][i,j];
    p01mc[[l]][i,j] = Na[[j]][[l]] - p11mc[[l]][i,j];
    p00mc[[l]][i,j] = N - NP[[i]][[l]] - Na[[j]][[l]];
  }
}
}

save(p11mc,p10mc,p01mc,p00mc,file="montecarlo_overlap.RData")

NP <- vector("list",n);
Na <- vector("list",m);
for(i in 1:n) NP[[i]] <- sapply(Peq[[i]], function(x) sum(x>0));
for(j in 1:m) Na[[j]] <- sapply(Aeq[[j]], function(x) sum(x>0));
for(i in 1:n) NP[[i]] <- sapply(Peq[[i]], function(x) sum(x>0)/(L*L));
for(j in 1:m) Na[[j]] <- sapply(Aeq[[j]], function(x) sum(x>0)/(L*L));

save(NP,file="totalsumP.RData")
save(Na,file="totalsumA.RData")
for(i in 1:n) NP[[i]] <- sapply(NP[[i]], function(x) x*N);
for(j in 1:m) Na[[j]] <- sapply(Na[[j]], function(x) x*N);
# Mutual information
p11mc2 <- vector("list",sims);
p10mc2 <- vector("list",sims);
p01mc2 <- vector("list",sims);
p00mc2 <- vector("list",sims);
for(l in 1:sims){
    p11mc2[[l]] = sapply(p11mc[[l]],function(x) x/(51*51))
    p10mc2[[l]] = sapply(p10mc[[l]],function(x) x/(51*51))
    p01mc2[[l]] = sapply(p01mc[[l]],function(x) x/(51*51))
    p00mc2[[l]] = sapply(p00mc[[l]],function(x) x/(51*51))
}
save(p11mc2,p10mc2,p01mc2,p00mc2,file="pmc2.RData")

for(l in 1:sims){
     p11mc2[[l]] = matrix(p11mc2[[l]],18,15)
     p10mc2[[l]] = matrix(p10mc2[[l]],18,15)
     p01mc2[[l]] = matrix(p01mc2[[l]],18,15)
     p00mc2[[l]] = matrix(p00mc2[[l]],18,15)
}     

Mutual_mc = vector("list",sims)
HXY_mc = vector("list",sims)
HXmc = vector("list",sims)
HYmc = vector("list",sims)
Mmc = vector("list",sims)

for(l in 1:sims){
  Mutual_mc[[l]]=matrix(0,n,m);
  HXY_mc[[l]]=matrix(0,n,m);
  HXmc[[l]] = matrix(0,n,m);
  HYmc[[l]] = matrix(0,n,m);
  Mmc[[l]]=matrix(0,n,m);
  for(i in 1:n){
    for(j in 1:m){
      Mutual_mc[[l]][i,j] =  p11mc2[[l]][i,j]*log2(p11mc2[[l]][i,j]/(NP[[i]][[l]]*Na[[j]][[l]])) + p10mc2[[l]][i,j]*log2(p10mc2[[l]][i,j]/(NP[[i]][[l]]*(1-Na[[j]][[l]]))) + p01mc2[[l]][i,j]*log2(p01mc2[[l]][i,j]/((1-NP[[i]][[l]])*Na[[j]][[l]])) + p00mc2[[l]][i,j]*log2(p00mc2[[l]][i,j]/((1-NP[[i]][[l]])*(1-Na[[j]][[l]])));
      #if(is.na(Mutual_mc[[l]][i,j])==1) Mutual_mc[[l]][i,j]=0;
      # Entropy of both 
      HXY_mc[[l]][i,j] = - (p11mc2[[l]][i,j]*log2(p11mc2[[l]][i,j]) + p10mc2[[l]][i,j]*log2(p10mc2[[l]][i,j]) + p01mc2[[l]][i,j]*log2(p01mc2[[l]][i,j]) + p00mc2[[l]][i,j]*log2(p00mc2[[l]][i,j]) )
      
      #Entropy plants
      HXmc[[l]][i,j] = - (NP[[i]][[l]]*log2(NP[[i]][[l]]) + (1-NP[[i]][[l]])*log2(1-NP[[i]][[l]]) )
      HYmc[[l]][i,j] = - (Na[[j]][[l]]*log2(Na[[j]][[l]]) + (1-Na[[j]][[l]])*log2(1-Na[[j]][[l]]) )
      
      Mmc[[l]][i,j] = HXmc[[l]][i,j] + HYmc[[l]][i,j] - HXY_mc[[l]][i,j];
 }
}
}

save(Mutual_mc,file="mutual_info_mc.RData")
load("mutual_info_mc.RData")

save(HXY_mc,HYmc,HXmc,file="entropy_mc")

meanMutual_mc=array(0,dim=c(18,15,1000))
mutual_plot <- array(meanMutual_mc,dim=c(270,1,1000))
entropyXY = array(0,dim=c(18,15,1000))
mut_entro= array(0,dim=c(18,15,1000))

Yplot =  meanMutual_mc[1,1,]
Xplot = rep(1,1000);
   plot.new()
   plot.window(xlim=c(1,270), xaxs="i", ylim=c(0,0.007), yaxs="i")
   for(k in 1:270){
    points(rep(k,sims),mut_entro[k,1,]);
   }     
 axis(1)
 axis(2)
 box()

  for(i in 1:n){
    for(j in 1:m){
      for(l in 1:sims){
        meanMutual_mc[i,j,l]=Mutual_mc[[l]][i,j];
	entropyXY[i,j,l]=HXY_mc[[l]][i,j];
	mut_entro[i,j,l]= meanMutual_mc[i,j,l]/entropyXY[i,j,l];
      }
    }
   }
  
  save(mut_entro,file="mut_entro.RData")
  
  
  meanmcmut=array(0,dim=c(18,15))
  sdmcmut=array(0,dim=c(18,15))
  medianmcmut=array(0,dim=c(18,15))
  for(i in 1:n){
    for(j in 1:m){
       meanmcmut[i,j]=mean(meanMutual_mc[i,j,],na.rm=TRUE);
       sdmcmut[i,j]=sd(meanMutual_mc[i,j,],na.rm=TRUE);
       medianmcmut[i,j]=median(meanMutual_mc[i,j,],na.rm=TRUE);
    }
  }

save(sdmcmut,meanmcmut,file="sd_meanmcmut.RData")

## Matrix plot bubble (mean)
coords=c(0,0);
for(i in 1:n){
  for(j in 1:m){
      coords = rbind(coords,c(i,j));
   }
}

matrix_plot <- expand.grid(1:n,1:m);
coords <- cbind(matrix_plot[,1],matrix_plot[,2])

bubble <- matrix(meanmcmut,270,1);
# bubble plot with transparency (alpha)
symbols(coords,circles=bubble,bg=rgb(red=0, green=0, blue=255, alpha=25, max=255) )

obs=Mutual;


z_score=matrix(0,n,m);
binary_mc=matrix(0,n,m);
binary_mc_sig=matrix(0,n,m);
pvalue_mc=matrix(0,n,m);
for(i in 1:n){
  for(j in 1:m){
      z_score[i,j]=(obs[i,j]-meanmcmut[i,j])/sdmcmut[i,j];
      pvalue_mc[i,j]=pnorm(-abs(z_score[i,j]))
      if(is.nan(obs[i,j])==1) obs[i,j]=0; 
      if(obs[i,j]>meanmcmut[i,j]){
	  binary_mc[i,j]=1
      }
      if(pvalue_mc[i,j]<0.05/(15*18)){  #bonferroni correction
	  binary_mc_sig[i,j]=1
      }
  }
}

save(z_score,pvalue_mc,binary_mc_sig,file="results_mutualinfo.RData")

## Z-scores calculated with corrected Mutual Information: MI/HXY
obs=Mutual2;


for(i in 1:n){
  for(j in 1:m){
     if(is.na(obs[i,j])==1) obs[i,j]=0
  }
}

z_score_2=matrix(0,n,m);
binary_mc_2=matrix(0,n,m);
binary_mc_sig_2=matrix(0,n,m);
pvalue_mc_2=matrix(0,n,m);
for(i in 1:n){
  for(j in 1:m){
      z_score_2[i,j]=(obs[i,j]-meanmcmut[i,j])/sdmcmut[i,j];
      pvalue_mc_2[i,j]=pnorm(-abs(z_score_2[i,j]))
      if(is.nan(obs[i,j])==1) obs[i,j]=0; 
      if(obs[i,j]>meanmcmut[i,j]){
	  binary_mc_2[i,j]=1
      }
      if(pvalue_mc_2[i,j]<0.05/(15*18)){  #bonferroni correction
	  binary_mc_sig_2[i,j]=1
      }
  }
}

save(z_score_2,pvalue_mc_2,binary_mc_sig_2,file="results_mutualinfo_2.RData")

coord2=c(0,0);
for(i in 1:n){
  for(j in 1:m){
    if(binary_mc_sig_2[i,j]==1){
        coord2=rbind(coord2,c(i,j));
    }
  }
}
coord2=coord2[-1,]


matrix_plot <- expand.grid(1:n,1:m);
coords <- cbind(matrix_plot[,2],matrix_plot[,1])
bin <- sortweb(binary_mc_sig_2,sort.order="inc")
coord2=c(0,0);
for(i in 1:n){
  for(j in 1:m){
    if(binary_mc_sig_2[i,j]==1){
        coord2=rbind(coord2,c(i,j));
    }
  }
}

squares1 = matrix(bin,270,1)

symbols(coord2,squares=rep(1,92),bg=rgb(red=0, green=0, blue=255, alpha=125, max=255),inches=0.2)