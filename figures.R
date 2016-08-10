
## C - SCORES  PLOTS
## Significant interactions with Bonferroni correction for aggregated (1) and segregated (2) occurrences

# Data:
save(Z_cs_eq, Z_cs_fix, Z_cs_swap, Z_cs_env, pvalc_eq, pvalc_fix, pvalc_swap, pvalc_env, file="C_scores_Z_scores.RData")
load("C_scores_Z_scores.RData")

n=18;m=15;
bin_sig1 = matrix(0,n,m)
bin_sig2 = matrix(0,n,m)
bin_sig3 = matrix(0,n,m)
bin_sig4 = matrix(0,n,m)

for(i in 1:n){
  for(j in 1:m){
     # Poisson CSR null model
     if(pvalc_eq[i,j]< 0.05/270 && Z_cs_eq[i,j]< 0 ){
           bin_sig1[i,j] = 1
     }else if( pvalc_eq[i,j]< 0.05/270 && Z_cs_eq[i,j]> 0 ){
           bin_sig1[i,j] = 2;
     }
     # Fixed null model
     if(pvalc_fix[i,j]< 0.05/270 && Z_cs_fix[i,j]< 0 ){
           bin_sig2[i,j] = 1
     }else if( pvalc_fix[i,j]< 0.05/270 && Z_cs_fix[i,j]> 0 ){
           bin_sig2[i,j] = 2;
     }
     # swap null model
     if(pvalc_swap[i,j]< 0.05/270 && Z_cs_swap[i,j]< 0 ){
           bin_sig4[i,j] = 1
     }else if( pvalc_swap[i,j]< 0.05/270 && Z_cs_swap[i,j]> 0 ){
           bin_sig4[i,j] = 2;
     }
     # Environmentally constrained null model
     if(pvalc_env[i,j]< 0.05/270 && Z_cs_env[i,j]< 0 ){
           bin_sig3[i,j] = 1
     }else if( pvalc_env[i,j]< 0.05/270 && Z_cs_env[i,j]> 0 ){
           bin_sig3[i,j] = 2;
     }
   }
}

# Nestedness CSR, Fixed, swap and ENV.
bin1=matrix(0,n,m)
bin2=matrix(0,n,m)
bin3=matrix(0,n,m)
bin4=matrix(0,n,m)
library("bipartite")
for(i in 1:n){
  for(j in 1:m){
   if(bin_sig1[i,j]==1){
       bin1[i,j]=1;
   }
   if(bin_sig2[i,j]==1){
       bin2[i,j]=1;
   }
   if(bin_sig3[i,j]==1){
       bin3[i,j]=1;
   }
   if(bin_sig4[i,j]==1){
       bin4[i,j]=1;
   }
  }
}



save(bin_sig1,bin_sig2,bin_sig3,bin_sig4,file="binary_sig_cscores.RData")
load("binary_sig_cscores.RData")


postscript(file="matrix_cscore_2.eps",width=15,height=10,horizontal=TRUE)
#par(mfcol=c(1,3),pin=c(3,3))

nf <- layout(matrix(c(1,2,3),1,3,byrow=TRUE), c(4,4,4), c(4,4,4), TRUE)


par(mar=c(4,4,1,1))
#par(mfcol=c(2,2), mar=c(4,4,1,1), oma=c(1.5,4,1,1))
 # Null model equip Poisson CSR
coord1e=c(0,0);
for(i in 1:n){
  for(j in 1:m){
    if(bin_sig1[i,j]==1){
        coord1e=rbind(coord1e,c(i,j));
    }
  }
}
coord1e=coord1e[-1,]

# aggregation blue squares
symbols(coord1e,squares=rep(1,dim(coord1e)[1]),fg="white",bg=rgb(red=0, green=0, blue=255, alpha=255, max=255),inches=0.17, xlim=c(0,18),ylim=c(0,15),xlab="Plants",ylab="AMF", cex.lab=1.3,main="CSR",xaxt='n', yaxt='n')
# seggregation red squares
coord2e=c(0,0);
for(i in 1:n){
  for(j in 1:m){
    if(bin_sig1[i,j]==2){
        coord2e=rbind(coord2e,c(i,j));
    }
  }
}
coord2e=coord2e[-1,]

symbols(coord2e,squares=rep(1,dim(coord2e)[1]),fg="white",bg=rgb(red=255, green=0, blue=0, alpha=255, max=255),inches=0.17, xlim=c(0,18),ylim=c(0,15), add=TRUE,xlab="Plants",ylab="AMF", main="CSR",xaxt='n', yaxt='n')


 
par(mar=c(4,4,1,1)) 
# Null model fixed
coord1=c(0,0);
for(i in 1:n){
  for(j in 1:m){
    if(bin_sig2[i,j]==1){
        coord1=rbind(coord1,c(i,j));
    }
  }
}
coord1=coord1[-1,]

# aggregation blue squares
symbols(coord1,squares=rep(1,dim(coord1)[1]),fg="white",bg=rgb(red=0, green=0, blue=255, alpha=255, max=255),inches=0.17,xlab="Plants",ylab="AMF",cex.lab=1.3, xlim=c(0,18),ylim=c(0,15),main="Fixed",xaxt='n', yaxt='n')
# seggregation red squares
coord2=c(0,0);
for(i in 1:n){
  for(j in 1:m){
    if(bin_sig2[i,j]==2){
        coord2=rbind(coord2,c(i,j));
    }
  }
}
coord2=coord2[-1,]

symbols(coord2,squares=rep(1,dim(coord2)[1]),fg="white",bg=rgb(red=255, green=0, blue=0, alpha=255, max=255),inches=0.17, xlim=c(0,18),ylim=c(0,15),xlab="Plants",ylab="AMF",cex.lab=1.3, add=TRUE, main="Fixed",xaxt='n', yaxt='n')

 par(mar=c(4,4,1,1))
# Null model environm.
coord3=c(0,0);
for(i in 1:n){
  for(j in 1:m){
    if(bin_sig3[i,j]==1){
        coord3=rbind(coord3,c(i,j));
    }
  }
}
coord3=coord3[-1,]

# aggregation blue squares
symbols(coord3,squares=rep(1,dim(coord3)[1]),fg="white",bg=rgb(red=0, green=0, blue=255, alpha=255, max=255),inches=0.17, xlim=c(0,18),ylim=c(0,15),xlab="Plants",ylab="AMF",cex.lab=1.3,main="Env. constrained",xaxt='n', yaxt='n')
# seggregation red squares
coord4=c(0,0);
for(i in 1:n){
  for(j in 1:m){
    if(bin_sig3[i,j]==2){
        coord4=rbind(coord4,c(i,j));
    }
  }
}
coord4=coord4[-1,]

symbols(coord4,squares=rep(1,dim(coord4)[1]),fg="white",bg=rgb(red=255, green=0, blue=0, alpha=255, max=255),inches=0.17, xlim=c(0,18),ylim=c(0,15), add=TRUE,xlab="Plants",ylab="AMF",cex.lab=1.3, main="Env. constrained",xaxt='n', yaxt='n')

## Swap null model
par(mar=c(4,4,1,1))
# Null model environm.
coord5=c(0,0);
for(i in 1:n){
  for(j in 1:m){
    if(bin_sig4[i,j]==1){
        coord5=rbind(coord5,c(i,j));
    }
  }
}
coord5=coord5[-1,]

# aggregation blue squares
symbols(coord5,squares=rep(1,dim(coord5)[1]),fg="white",bg=rgb(red=0, green=0, blue=255, alpha=255, max=255),inches=0.17, xlim=c(0,18),ylim=c(0,15),xlab="Plants",ylab="AMF",cex.lab=1.3,main="Swap",xaxt='n', yaxt='n')
# seggregation red squares
coord6=c(0,0);
for(i in 1:n){
  for(j in 1:m){
    if(bin_sig4[i,j]==2){
        coord6=rbind(coord6,c(i,j));
    }
  }
}
coord6=coord6[-1,]

symbols(coord6,squares=rep(1,dim(coord6)[1]),fg="white",bg=rgb(red=255, green=0, blue=0, alpha=255, max=255),inches=0.17, xlim=c(0,18),ylim=c(0,15), add=TRUE,xlab="Plants",ylab="AMF",cex.lab=1.3, main="Swap",xaxt='n', yaxt='n')



dev.off()

postscript(file="test.eps",width=10,height=10)
nf <- layout(matrix(c(1,2,3),1,3,byrow=TRUE), c(5,5,5), c(5,5,5), TRUE)
par(mar=c(3,3,1,1))
hist(runif(100),ylab="mm",xlab="pp")
par(mar=c(3,3,1,1))
hist(runif(100),ylab="mm",xlab="pp")
par(mar=c(3,3,1,1))
hist(runif(100),ylab="mm",xlab="pp")
dev.off()

par(fig=c(0,0.4,0.2,0.6))
hist(runif(100),ylab="mm",xlab="pp")


#### Matrices plot for Mutual information observed vs. predicted


### Matrices plot for fij observed vs. predicted

save(binaryfij_eq,binaryfij_fix,binaryfij_env,binaryfij_swap,file="sig_matrices_fij.RData")
load("sig_matrices_fij.RData")


postscript(file="matrix_fij2.eps",width=10,height=10,horizontal=FALSE)
par(mfrow=c(3,1))
postscript(file="matrix_fij.eps",width=10,height=10,horizontal=FALSE)
par(mfrow=c(2,2))

postscript(file="csr_fij.eps",width=10,height=10,horizontal=FALSE)
# CSR
coord1fij=c(0,0);
for(i in 1:n){
  for(j in 1:m){
    if(binaryfij_eq[i,j]==1){
        coord1fij=rbind(coord1fij,c(i,j));
    }
  }
}
coord1fij=coord1fij[-1,]
symbols(coord1fij,squares=rep(1,dim(coord1fij)[1]),fg="white",bg=rgb(red=0, green=0, blue=255, alpha=255, max=255),inches=0.18, xlim=c(0,18),ylim=c(0,15),xlab="Plants",ylab="AMF",cex.lab=1.3,cex.main=1.5,main="CSR",xaxt='n', yaxt='n')
dev.off()
# fixed

coord2fij=c(0,0);
for(i in 1:n){
  for(j in 1:m){
    if(binaryfij_fix[i,j]==1){
        coord2fij=rbind(coord2fij,c(i,j));
    }
  }
}
coord2fij=coord2fij[-1,]
symbols(coord2fij,squares=rep(1,dim(coord2fij)[1]),fg="white",bg=rgb(red=0, green=0, blue=255, alpha=255, max=255),inches=0.18, xlim=c(0,18),ylim=c(0,15),xlab="Plants",ylab="AMF",cex.lab=1.3,cex.main=1.5,main="Fixed",xaxt='n', yaxt='n')

# Swap
postscript(file="swap_fij.eps",width=10,height=10,horizontal=FALSE)
coord3fij=c(0,0);
for(i in 1:n){
  for(j in 1:m){
    if(binaryfij_swap[i,j]==1){
        coord3fij=rbind(coord3fij,c(i,j));
    }
  }
}
coord3fij=coord3fij[-1,]
symbols(coord3fij,squares=rep(1,dim(coord3fij)[1]),fg="white",bg=rgb(red=0, green=0, blue=255, alpha=255, max=255),inches=0.18, xlim=c(0,18),ylim=c(0,15),xlab="Plants",ylab="AMF",cex.lab=1.3,cex.main=1.5,main="Swap",xaxt='n', yaxt='n')
dev.off()
# Env
postscript(file="env_fij.eps",width=10,height=10,horizontal=FALSE)
coord4fij=c(0,0);
for(i in 1:n){
  for(j in 1:m){
    if(binaryfij_env[i,j]==1){
        coord4fij=rbind(coord4fij,c(i,j));
    }
  }
}
coord4fij=coord4fij[-1,]
symbols(coord4fij,squares=rep(1,dim(coord4fij)[1]),fg="white",bg=rgb(red=0, green=0, blue=255, alpha=255, max=255),inches=0.18, xlim=c(0,18),ylim=c(0,15),xlab="Plants",ylab="AMF",cex.lab=1.3,cex.main=1.5,main="Env. const.",xaxt='n', yaxt='n')

dev.off()


## Correlograms  for null models
# Load Data from simulations
load("Pnullmatrices.RData")
valbin <- matrix(Peq[[1]][[1]],51*51,1, byrow=TRUE)
valbinswap <- matrix(Pswap[[1]][[1]],51*51,1, byrow=TRUE)
valbinfix <- matrix(Pfix[[1]][[1]],51*51,1, byrow=TRUE)
valbinenv <- matrix(Penv[[1]][[1]],51*51,1, byrow=TRUE)


coordsfij<-cell2nb(51, 51, type="rook", torus=FALSE)

nmcor <- sp.correlogram(coordsfij,valbin,order=10, method="I",style = "W")
plot.spcor(nmcor)

nmcorswap <- sp.correlogram(coordsfij,valbinswap,order=10, method="I",style = "W")
plot.spcor(nmcorswap)

nmcorenv<- sp.correlogram(coordsfij,valbinenv,order=10, method="I",style = "W")
plot.spcor(nmcorenv)

## CSR correlogram Plant 1 Bromus inermis
valbincsr=vector("list",sims)
nmcorcsr=vector("list",sims)
moran_csr=vector("list",sims)
for(l in 1:sims){
    valbincsr[[l]] = matrix(Peq[[1]][[l]],51*51,1, byrow=TRUE);
    nmcorcsr[[l]] <- sp.correlogram(coordsfij,valbincsr[[l]],order=10, method="I",style = "W", randomisation=FALSE)
    moran_csr[[l]]=nmcorcsr[[l]]$res[,1]
}

save(nmcorcsr,moran_csr,valbincsr,file="csr_corr.RData")
 
## Env correlogram Plant 1 Bromus inermis
valbinenv=vector("list",sims)
nmcorenv=vector("list",sims)
moran_env=vector("list",sims)
for(l in 1:sims){
    valbinenv[[l]] = matrix(Penv[[1]][[l]],51*51,1, byrow=TRUE);
    nmcorenv[[l]] <- sp.correlogram(coordsfij,valbinenv[[l]],order=10, method="I",style = "W", randomisation=FALSE)
    moran_env[[l]]=nmcorenv[[l]]$res[,1]
}

save(nmcorenv,moran_env,valbinenv,file="env_corr.RData")


## SWAP correlogram Plant 1 Bromus inermis

plot.spcor(nmcorswap$res[,1],type="l")
valbinswap=vector("list",sims)
nmcorswap=vector("list",sims)
moran=vector("list",sims)
for(l in 1:sims){
    valbinswap[[l]] = matrix(Pswap[[1]][[l]],51*51,1, byrow=TRUE);
    nmcorswap[[l]] <- sp.correlogram(coordsfij,valbinswap[[l]],order=10, method="I",style = "W", randomisation=FALSE)
    moran[[l]]=nmcorswap[[l]]$res[,1]
}

# 40 values of Moran I from simulated data of swap null model
save(nmcorswap,moran,valbinswap,file="swap_corr.RData")

#Observed Plant 1 Bromus inermis

valbinobs<- matrix(P[[1]],51*51,1, byrow=TRUE)

nmcorobs<- sp.correlogram(coordsfij,valbinobs,order=10, method="I",style = "W")
plot.spcor(nmcorobs)


# observed versus swap plant 1 correlogram
plot(moran[[1]],type="l",col="grey", xlab="Lag", ylab="Moran's I")
for(l in 1:35){
 lines(moran[[l]],type="l",col="grey")
}
lines(nmcorobs$res[,1],type="l",col="red",lwd=1.8)


## Correlograms
postscript(file="correlograms_nm.eps",width=10,height=10,horizontal=FALSE)
par(mfrow=c(2,2))

plot.spcor(nmcor, xlab="Lag", ylab="Moran's I", main="CSR", xlim=c(0,11), ylim=c(-0.05,0.15))

#nmcorfix <- sp.correlogram(coordsfij,valbinfix,order=10, method="I",style = "W")
plot.spcor(nmcorfix,main="Fixed",xlim=c(0,11), ylim=c(-0.05,0.15))

plot.spcor(nmcorenv, main="Env. const.",xlim=c(0,11), ylim=c(-0.05,0.15))

plot(moran[[1]],type="l",col="grey", xlab="lags", ylab="Moran's I",main="Swap",xlim=c(0,11), ylim=c(-0.05,0.15))
for(l in 1:35){
 lines(moran[[l]],type="l",col="grey")
}
lines(nmcorobs$res[,1],type="l",col="red",lwd=1.8)
lines(seq(0,11),rep(0,12))
dev.off()


#### Mutual information observed versus predicted

save(inter_matrix,mutualism,seggregation,file="matrices_sig_MI_theo.RData")


labels_amf <- c("A.denticulata", "A.morrowiae", "A.spinosa", "E.colombiana", "G.claroideum", "G.etunicatum", "G.intraradices", "G.macrocarpum", "G.mosseae", "G.gigantea", "G.margarita", "S.calospora", "S.dipurpurascens", "S.heterogama" ,"S.pellucia");

labels_plants <- c("B.inermis", "F.virginiana", "P.lanceolata", "S.graminifolia", "C.leucanthemum", "P.compresa", "H.pratense", "A.novae-angliae", "H.perforatum", "E.vulgare", "S.canadensis", "P.pratensis", "D.carota", "A.millefolium", "E.strigosus", "C.arvensis", "R.serotina", "P.vulgaris")


load("matrices_sig_MI_theo.RData")

coordmut=c(0,0);
coordseg=c(0,0);
for(i in 1:n){
  for(j in 1:m){
    if(mutualism[i,j]==1){
        coordmut=rbind(coordmut,c(i,j));
    }
    if(seggregation[i,j]==1){
        coordseg=rbind(coordseg,c(i,j));
    }
  }
}
coordmut=coordmut[-1,]
coordseg=coordseg[-1,]

postscript(file="mutual_info_matrix2.eps",width=15,height=15,horizontal=FALSE)
layout(1,15,15, TRUE)
par(mar=c(7, 7, 2, 2)) 
symbols(coordmut,squares=rep(1,dim(coordmut)[1]),fg="white",bg=rgb(red=0, green=0, blue=255, alpha=255, max=255),inches=0.18, xlim=c(0,18),ylim=c(0,15),xlab="",ylab="",cex.lab=1.3,cex.main=1.5,main="Mutual Information",xaxt='n', yaxt='n')
symbols(coordseg,squares=rep(1,dim(coordseg)[1]),fg="white",bg=rgb(red=255, green=0, blue=0, alpha=255, max=255),inches=0.18, xlim=c(0,18),ylim=c(0,15),xlab="Plants",ylab="AMF",cex.lab=1.3,cex.main=1.5,xaxt='n', yaxt='n', add=TRUE)

axis(1, at= NULL,labels = FALSE,xaxt='n', yaxt='n')
text(seq(1, 18, by=1), par("usr")[3]-1.5, labels = labels_plants, srt = 90, pos = 1, xpd = TRUE, cex=0.85)
axis(2, at= NULL,labels = FALSE,xaxt='n', yaxt='n')
text(rep(-3,15),seq(1.4, 15.4, by=1), par("usr")[3]-0.3, labels = labels_amf, srt = 0, pos = 1, xpd = TRUE,cex=0.85)

dev.off()

## Coincidences of signif. interactions between C_score and MI

# CSR mutualisms
bin_sig1
mut_c=matrix(0,n,m);

for(i in 1:n){
   for(j in 1:m){
       if(bin_sig1[i,j]==1){
          mut_c[i,j]=1;
       }
   }
}
coinc=rep(0,n)
for(i in 1:n){
         coinc[i] = sum(mutualism[i,]*mut_c[i,]);
}

sum(coinc)

seg_c=matrix(0,n,m);

for(i in 1:n){
   for(j in 1:m){
       if(bin_sig1[i,j]==2){
          seg_c[i,j]=1;
       }
   }
}
coinc_seg=rep(0,n)
for(i in 1:n){
         coinc_seg[i] = sum(seggregation[i,]*seg_c[i,]);
}

## Coincidences of signif. interactions between C_score and MI and fij

binaryfij_eq

coinc_mutfij=rep(0,n)
for(i in 1:n){
         coinc_mutfij[i] = sum(binaryfij_eq[i,]*mut_c[i,]);
}

coinc_mutmi=rep(0,n)
for(i in 1:n){
         coinc_mutmi[i] = sum(binaryfij_eq[i,]*mutualism[i,]);
}
