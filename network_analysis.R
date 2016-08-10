## ---------- NETWORK STRUCTURE OF PLANT-AMF COMMUNITIES ---------- ##

## Data: 4 null models , 1000 null matrices for each species

## Structure script:
# Nestedness: Brualdi, NODF, Temperature (using VEGAN package)
# Conectance
# nesteddisc(comm)
# nestedtemp(comm, ...)


## Calculate Nestedness (NODF) calling library bipartite 

# Load Data from simulations
load("binary_eq.RData");
load("binary_fix.RData");
load("binary_swap.RData");
load("binary_env.RData");
load("nodf_observed.RData");

th = seq(0,1,0.01);

n=18;m=15;sims=1000;

library(bipartite);
library(VEGAN);

nodf_eq = matrix(0,length(th),sims);
nodf_fix = matrix(0,length(th),sims);
nodf_swap = matrix(0,length(th),sims);
nodf_env = matrix(0,length(th),sims);

temp_eq = matrix(0,length(th),sims);
temp_fix = matrix(0,length(th),sims);
temp_swap = matrix(0,length(th),sims);
temp_env = matrix(0,length(th),sims);

# Nestedness for  interaction strength NODF
 for(h in 1:length(th)){
  for(l in 1:sims){
    nodf_eq[h,l] = nested(Bineq[[h]][[l]],method="NODF");
    nodf_fix[h,l] = nested(Binfix[[h]][[l]],method="NODF");
    nodf_swap[h,l] = nested(Binswap[[h]][[l]],method="NODF");
    nodf_env[h,l] = nested(Binenv[[h]][[l]],method="NODF");
    temp_eq[h,l] = nested(Bineq[[h]][[l]],method="binmatnest2");
    temp_fix[h,l] = nested(Binfix[[h]][[l]],method="binmatnest2");
    temp_swap[h,l] = nested(Binswap[[h]][[l]],method="binmatnest2");
    temp_env[h,l] = nested(Binenv[[h]][[l]],method="binmatnest2");
 }
}

save(nodf_eq,nodf_fix,nodf_swap,nodf_env,file="nodf_sim_1000.RData")
load("nodf_sim_1000.RData")
save(nodf_swap,file="nodf_swap_2000.RData")

# Mean + SD nestedness for simulations at each threshold value MP

  meannodf_eq <- apply(nodf_eq,1,mean)
  stdnodf_eq <- apply(nodf_eq,1,sd)
  meannodf_fix <- apply(nodf_fix,1,mean)
  stdnodf_fix <- apply(nodf_fix,1,sd)
  meannodf_swap <- apply(nodf_swap,1,mean)
  stdnodf_swap <- apply(nodf_swap,1,sd)
  meannodf_env <- apply(nodf_env,1,mean)
  stdnodf_env <- apply(nodf_env,1,sd)
  
 

 for(h in 1:length(thp)){
      error_P[h] = qnorm(0.975)*stdnodfP[h]/sqrt(sims);
      error_P2[h] = qnorm(0.975)*stdnodfP2[h]/sqrt(sims);
      upper_clP[h] = meannodfP[h] + error_P[h];
      upper_clP2[h] = meannodfP2[h] + error_P2[h];
      lower_clP[h] = meannodfP[h] - error_P[h];
      lower_clP2[h] = meannodfP2[h] - error_P2[h];
 }

## Confidence intervals

error_nodf_eq <- apply(as.array(stdnodf_eq),1, function(x) qnorm(0.975)*x);
error_nodf_fix <- apply(as.array(stdnodf_fix),1, function(x) qnorm(0.975)*x);
error_nodf_swap <- apply(as.array(stdnodf_swap),1, function(x) qnorm(0.975)*x);
error_nodf_env <- apply(as.array(stdnodf_env),1, function(x) qnorm(0.975)*x);
# 

 upper_eq <- rbind(meannodf_eq,error_nodf_eq);
 upper_nodf_eq <- apply(upper_eq,2,sum);
 upper_fix <- rbind(meannodf_fix,error_nodf_fix);
 upper_nodf_fix <- apply(upper_fix,2,sum);
 upper_swap <- rbind(meannodf_swap,error_nodf_swap);
 upper_nodf_swap <- apply(upper_swap,2,sum);
 upper_env <- rbind(meannodf_env,error_nodf_env);
 upper_nodf_env <- apply(upper_env,2,sum);
 lower_nodf_eq <- rep(0,1,length(th));
 lower_nodf_fix <- rep(0,1,length(th));
 lower_nodf_swap <- rep(0,1,length(th))
 lower_nodf_env <- rep(0,1,length(th))
 
 for(h in 1:length(th)){
    lower_nodf_eq[h] = meannodf_eq[h] - error_nodf_eq[h];
    lower_nodf_fix[h] = meannodf_fix[h] - error_nodf_fix[h];
    lower_nodf_swap[h] = meannodf_swap[h] - error_nodf_swap[h];
    lower_nodf_env[h] = meannodf_env[h] - error_nodf_env[h];
 }
 
 save(upper_nodf_swap,lower_nodf_swap,meannodf_swap,stdnodf_swap,file="stat_swap_nodf.RData")
 
 load("stat_swap_nodf.RData")
## Z-scores for NODF, Temp 
## Plotting 2*2 figures (4 null models)
par(mfrow=c(2,2))
plot(th,meannodf_eq,col="grey",xlab="Threshold", ylab="NODF",main="Null model equiprobable")
apply(nodf_eq,2,points,col="grey",x=th);
lines(th,meannodf_eq, lwd=3); 
lines(th,nodf_obs,lwd=2,type="l",col="red")
lines(th,upper_nodf_eq,col="blue"); lines(th,lower_nodf_eq,col="blue")
 
plot(th,meannodf_fix,col="grey",xlab="Threshold", ylab="NODF",main="Null model fixed");
apply(nodf_fix,2,points,col="grey",x=th); 
lines(th,meannodf_fix, lwd=3); 
lines(th,nodf_obs,lwd=2,type="l",col="red")
lines(th,upper_nodf_fix,col="blue"); lines(th,lower_nodf_fix,col="blue")

plot(th,meannodf_swap,col="grey",xlab="Threshold", ylab="NODF",main="Null model swap",);
apply(nodf_swap,2,points,col="grey",x=th); 
lines(th,meannodf_swap, lwd=3);
lines(th,nodf_obs,lwd=2,type="l",col="red")
lines(th,upper_nodf_swap,col="blue"); lines(th,lower_nodf_swap,col="blue")

plot(th,meannodf_env,col="grey",xlab="Threshold", ylab="NODF",main="Null model environmentally constrained");
apply(nodf_env,2,points,col="grey",x=th); 
lines(th,meannodf_env,lwd=3); 
lines(th,nodf_obs,lwd=2,type="l",col="red")
lines(th,upper_nodf_env,col="blue"); lines(th,lower_nodf_env,col="blue")

## Z-scores for NODF, Temp 
## Plotting 3*1 figures (3 null models)



#layout(matrix(c(1,2,3),1,3,byrow=TRUE),c(4,4,4),c(4,4,4),TRUE)
dev.new(width=9,height=4)
postscript(file="nested_threshold2.eps",width=9,height=4)
par(mfrow=c(1,3),mar=c(4,4,3,1),oma = c( 0, 0, 2, 0 ))
plot(th,meannodf_eq,col="grey",cex.lab=1.3,xlab="", ylab="NODF",main="CSR")
apply(nodf_eq,2,points,col="grey",x=th);
lines(th,meannodf_eq, lwd=3); 
lines(th,nodf_obs,lwd=2,type="l",col="red")
lines(th,upper_nodf_eq,col="blue"); lines(th,lower_nodf_eq,col="blue")

plot(th,meannodf_swap,col="grey",cex.lab=1.5,xlab="Threshold", ylab="",main="Swap");
apply(nodf_swap,2,points,col="grey",x=th); 
lines(th,meannodf_swap, lwd=3); 
lines(th,nodf_obs,lwd=2,type="l",col="red")
lines(th,upper_nodf_swap,col="blue"); lines(th,lower_nodf_swap,col="blue")

plot(th,meannodf_env,col="grey",cex.lab=1.3,xlab="", ylab="",main="Env. constrained");
apply(nodf_env,2,points,col="grey",x=th); 
lines(th,meannodf_env,lwd=3); 
lines(th,nodf_obs,lwd=2,type="l",col="red")
lines(th,upper_nodf_env,col="blue"); lines(th,lower_nodf_env,col="blue")

dev.off()

## Plot of (Obs - mean_nm)

dif1 <- nodf_obs-meannodf_env;
dif2 <- nodf_obs-meannodf_eq;
dif3 <- nodf_obs-meannodf_swap;
xx <- seq(0,0.6,by=0.01)

postscript("obs_exp_nodf.eps")

plot(xx,dif1[1:61], type='o', pch= 20, lty=40, col='blue', lwd=2, xlab="Threshold", ylab="Observed - Expected") 
lines(xx,dif2[1:61], type='o', pch= 20, lty=14, col='green', lwd=2) 
lines(xx,dif3[1:61], type='o', pch= 20, col='red', lwd=2) 

legend(0.5, max(dif1), c('ENV','CSR','SS'), cex=0.8, col=c('blue','green','red'), pch=c(20,20,20), lty=c(40,14,1), lwd=2,bty='n');
dev.off()

plot(xx,nodf_obs-meannodf_env, type='o', pch= 20, col= gray(0.65),lwd=2, xlim=c(0,0.6), ylim=c(-5,25)) 
lines(xx,nodf_obs-meannodf_eq, type='o', pch= 20,  col=gray(0.35), lwd=2,xlim=c(0,0.6), ylim=c(-5,25)) 
lines(xx,nodf_obs-meannodf_swap, type='o', pch= 20, col='black', lwd=2, xlim=c(0,0.6), ylim=c(-5,25)) 


# Z_scores= (obs-exp)/sd

Z_eq =  rep(0,1,length(th));
Z_fix =  rep(0,1,length(th));
Z_swap =  rep(0,1,length(th));
Z_env =  rep(0,1,length(th));
pval_eq = rep(0,1,length(th));
pval_fix = rep(0,1,length(th));
pval_swap = rep(0,1,length(th));
pval_env = rep(0,1,length(th));


for(h in 1:length(th)){
    Z_eq[h] = (nodf_obs[h]-meannodf_eq[h])/stdnodf_eq[h];
    Z_fix[h] = (nodf_obs[h]-meannodf_fix[h])/stdnodf_fix[h];
    Z_swap[h] = (nodf_obs[h]-meannodf_swap[h])/stdnodf_swap[h];
    Z_env[h] = (nodf_obs[h]-meannodf_env[h])/stdnodf_env[h];
    pval_eq[h] = 2*pnorm(-abs(Z_eq[h]))
    pval_fix[h] = 2*pnorm(-abs(Z_fix[h])) 
    pval_swap[h] = 2*pnorm(-abs(Z_swap[h]))
    pval_env[h] = 2*pnorm(-abs(Z_env[h]))
}
save(Z_eq,Z_fix,Z_swap,Z_env,pval_eq,pval_fix,pval_swap,pval_env,file="Z_scores.RData")


##  Co-occurrence fij values from null models - Calculate Z-scores of fij


save(MEANeq,Ceq,Asymeq,file="co_occurrence_eq.RData")
save(MEANfix,Cfix,Asymfix,file="co_occurrence_fix.RData")
save(MEANswap,Cswap,Asymswap,file="co_occurrence_swap.RData")
save(MEANenv,Cenv,Asymenv,file="co_occurrence_env.RData")
load("co_occurrence_eq.RData")
load("co_occurrence_fix.RData")
load("co_occurrence_swap.RData")
load("co_occurrence_env.RData")

save(MEAN,file="observed_mean_fij.RData")
load("observed_mean_fij.RData")

# Z_scores= (obs-exp)/sd of fij values from null models (no threshold) mean of mean_fij


fij_eq2=array(0,c(n,m,sims))
fij_fix2=array(0,c(n,m,sims))
fij_swap2=array(0,c(n,m,sims))
fij_env2=array(0,c(n,m,sims))

for(i in 1:n){
  for(j in 1:m){
    for(l in 1:sims){
        fij_eq2[i,j,l]= MEANeq[[l]][i,j];
	fij_fix2[i,j,l]=MEANfix[[l]][i,j];
	fij_swap2[i,j,l]=MEANswap[[l]][i,j];
	fij_env2[i,j,l]=MEANenv[[l]][i,j];
     }
   }
}


mean2fij_eq=matrix(0,n,m);
mean2fij_fix=matrix(0,n,m);
mean2fij_swap=matrix(0,n,m);
mean2fij_env=matrix(0,n,m);
sdfij_eq=matrix(0,n,m);
sdfij_fix=matrix(0,n,m);
sdfij_swap=matrix(0,n,m);
sdfij_env=matrix(0,n,m);

for(i in 1:n){
  for(j in 1:m){
      mean2fij_eq[i,j]=mean(fij_eq2[i,j,]);
      mean2fij_fix[i,j]=mean(fij_fix2[i,j,]);
      mean2fij_swap[i,j]=mean(fij_swap2[i,j,]);
      mean2fij_env[i,j]=mean(fij_env2[i,j,]);
      sdfij_eq[i,j]=sd(fij_eq2[i,j,]);
      sdfij_fix[i,j]=sd(fij_fix2[i,j,]);
      sdfij_swap[i,j]=sd(fij_swap2[i,j,]);
      sdfij_env[i,j]=sd(fij_env2[i,j,]);
   }
}


binaryfij_eq=matrix(0,n,m);
binaryfij_fix=matrix(0,n,m);
binaryfij_swap=matrix(0,n,m);
binaryfij_env=matrix(0,n,m);
Z_fij_eq=matrix(0,n,m);
Z_fij_fix=matrix(0,n,m);
Z_fij_swap=matrix(0,n,m);
Z_fij_env=matrix(0,n,m);
pvalfij_eq=matrix(0,n,m);
pvalfij_fix=matrix(0,n,m);
pvalfij_swap=matrix(0,n,m);
pvalfij_env=matrix(0,n,m);



for(i in 1:n){
  for(j in 1:m){
     Z_fij_eq[i,j]=(MEAN[i,j]-mean2fij_eq[i,j])/sdfij_eq[i,j];
     Z_fij_fix[i,j]=(MEAN[i,j]-mean2fij_fix[i,j])/sdfij_fix[i,j];
     Z_fij_swap[i,j]=(MEAN[i,j]-mean2fij_swap[i,j])/sdfij_swap[i,j];
     Z_fij_env[i,j]=(MEAN[i,j]-mean2fij_env[i,j])/sdfij_env[i,j];
     pvalfij_eq[i,j]=pnorm(-abs(Z_fij_eq[i,j]))
     pvalfij_fix[i,j]=pnorm(-abs(Z_fij_fix[i,j]))
     pvalfij_swap[i,j]=pnorm(-abs(Z_fij_swap[i,j]))
     pvalfij_env[i,j]=pnorm(-abs(Z_fij_env[i,j]))
 
     if(pvalfij_eq[i,j] < 0.05/(18*15)){
        binaryfij_eq[i,j]=1;
     }
     if(pvalfij_fix[i,j]< 0.05/(18*15)){
	binaryfij_fix[i,j]=1;
     }
     if(pvalfij_swap[i,j]< 0.05/(18*15)){
	binaryfij_swap[i,j]=1;
     }
     if(pvalfij_env[i,j]< 0.05/(18*15)){
	binaryfij_env[i,j]=1;
     }
   }
}

save(binaryfij_eq,binaryfij_fix,binaryfij_env,binaryfij_swap,file="sig_matrices_fij.RData")

## Z-scores for C_scores



Z_cs_eq=matrix(0,n,m);
Z_cs_fix=matrix(0,n,m);
Z_cs_swap=matrix(0,n,m);
Z_cs_env=matrix(0,n,m);
binaryc_eq=matrix(0,n,m);
binaryc_fix=matrix(0,n,m);
binaryc_swap=matrix(0,n,m);
binaryc_env=matrix(0,n,m);

Ceq2=array(0,c(n,m,sims))
Cfix2=array(0,c(n,m,sims))
Cswap2=array(0,c(n,m,sims))
Cenv2=array(0,c(n,m,sims))

for(i in 1:n){
  for(j in 1:m){
    for(l in 1:sims){
        Ceq2[i,j,l]= Ceq[[l]][i,j];
	Cfix2[i,j,l]=Cfix[[l]][i,j];
	Cswap2[i,j,l]=Cswap[[l]][i,j];
	Cenv2[i,j,l]=Cenv[[l]][i,j];
     }
   }
}
 
save(Ceq2,Cfix2,Cswap2,Cenv2,file="C_scores_mc2.RData");

meanCeq=matrix(0,n,m);
meanCfix=matrix(0,n,m);
meanCswap=matrix(0,n,m);
meanCenv=matrix(0,n,m);
sdCeq=matrix(0,n,m);
sdCfix=matrix(0,n,m);
sdCswap=matrix(0,n,m);
sdCenv=matrix(0,n,m);

for(i in 1:n){
  for(j in 1:m){
      meanCeq[i,j]=mean(Ceq2[i,j,]);
      meanCfix[i,j]=mean(Cfix2[i,j,]);
      meanCswap[i,j]=mean(Cswap2[i,j,]);
      meanCenv[i,j]=mean(Cenv2[i,j,]);
      sdCeq[i,j]=sd(Ceq2[i,j,]);
      sdCfix[i,j]=sd(Cfix2[i,j,]);
      sdCswap[i,j]=sd(Cswap2[i,j,]);
      sdCenv[i,j]=sd(Cenv2[i,j,]);
   }
}


binaryc_eq=matrix(0,n,m);
binaryc_fix=matrix(0,n,m);
binaryc_swap=matrix(0,n,m);
binaryc_env=matrix(0,n,m);
Z_cs_eq=matrix(0,n,m);
Z_cs_fix=matrix(0,n,m);
Z_cs_swap=matrix(0,n,m);
Z_cs_env=matrix(0,n,m);
pvalc_eq=matrix(0,n,m);
pvalc_fix=matrix(0,n,m);
pvalc_swap=matrix(0,n,m);
pvalc_env=matrix(0,n,m);

save(Cswap2, meanCswap, sdCswap, file="swap_cscore.RData")

for(i in 1:n){
  for(j in 1:m){
     Z_cs_eq[i,j]=(C_score[i,j]-meanCeq[i,j])/sdCeq[i,j];
     Z_cs_fix[i,j]=(C_score[i,j]-meanCfix[i,j])/sdCfix[i,j];
     Z_cs_swap[i,j]=(C_score[i,j]-meanCswap[i,j])/sdCswap[i,j];
     Z_cs_env[i,j]=(C_score[i,j]-meanCenv[i,j])/sdCenv[i,j];
     pvalc_eq[i,j]=pnorm(-abs(Z_cs_eq[i,j]))
     pvalc_fix[i,j]=pnorm(-abs(Z_cs_fix[i,j]))
     pvalc_swap[i,j]=pnorm(-abs(Z_cs_swap[i,j]))
     pvalc_env[i,j]=pnorm(-abs(Z_cs_env[i,j]))
 
     if(C_score[i,j]> meanCeq[i,j]){
        binaryc_eq[i,j]=1;
     }
     if(C_score[i,j]> meanCfix[i,j]){
	binaryc_fix[i,j]=1;
     }
     if(C_score[i,j]> meanCswap[i,j]){
	binaryc_swap[i,j]=1;
     }
     if(C_score[i,j]> meanCenv[i,j]){
	binaryc_env[i,j]=1;
     }
   }
}
save(Z_cs_swap,pvalc_swap,file="zscore_cscore_swap.RData")
critic_pval=0.05/(n*m)    # Bonferroni correction
save(Z_cs_eq, Z_cs_fix, Z_cs_swap, Z_cs_env, pvalc_eq, pvalc_fix, pvalc_swap, pvalc_env, file="C_scores_Z_scores.RData")

## Significant interactions with Bonferroni correction for aggregated (1) and segregated (2) occurrences
n=18;m=15;
bin_sig = matrix(0,n,m)

for(i in 1:n){
  for(j in 1:m){
     if(pvalc_eq[i,j]< 0.05/270 && Z_cs_eq[i,j]< 0 ){
           bin_sig[i,j] = 1
     }else if( pvalc_eq[i,j]< 0.05/270 && Z_cs_eq[i,j]> 0 ){
           bin_sig[i,j] = 2;
     }
   }
}


## Binary matrices null models with different stat significance 0.05, 0.01, critic_pval
binaryc05_eq=matrix(0,n,m);
binaryc05_fix=matrix(0,n,m);
binaryc05_swap=matrix(0,n,m);
binaryc05_env=matrix(0,n,m);
binaryc01_eq=matrix(0,n,m);
binaryc01_fix=matrix(0,n,m);
binaryc01_swap=matrix(0,n,m);
binaryc01_env=matrix(0,n,m);
for(i in 1:n){
  for(j in 1:m){
    if(pvalc_eq[i,j]<0.05){
        binaryc05_eq[i,j]=1;
     }
    if(pvalc_eq[i,j]<0.01){
        binaryc01_eq[i,j]=1;
     }
    if(pvalc_fix[i,j]<0.05){
        binaryc05_fix[i,j]=1;
     }
     if(pvalc_fix[i,j]<0.01){
        binaryc01_fix[i,j]=1;
     }
     if(pvalc_swap[i,j]<0.05){
        binaryc05_swap[i,j]=1;
     }
     if(pvalc_swap[i,j]<0.01){
        binaryc01_swap[i,j]=1;
     }
     if(pvalc_env[i,j]<0.05){
        binaryc05_env[i,j]=1;
     }
     if(pvalc_env[i,j]<0.01){
        binaryc01_env[i,j]=1;
     }
  }
}

## Connectance

con_eq = matrix(0,length(thp),sims);
con_fix = matrix(0,length(thp),sims);
con_swap = matrix(0,length(thp),sims);
con_env = matrix(0,length(thp),sims);

 
# Nestedness for mean interaction strength
 for(h in 1:length(thp)){
  for(l in 1:sims){
    con_eq[h,l] = sum(Bineq[[h]][[l]]>0)/(n*m)
    con_fix[h,l] = sum(Binfix[[h]][[l]]>0)/(n*m)
    con_swap[h,l] = sum(Binswap[[h]][[l]]>0)/(n*m)
    con_env[h,l] = sum(Binenv[[h]][[l]]>0)/(n*m)
 }
}

# Mean nestedness for simulations at each threshold value
 meannodf = rep(0,1,length(thp));
 meancon = rep(0,1,length(thp));
 meannodf2 = rep(0,1,length(thp));
  meancon2 = rep(0,1,length(thp));
 stdnodf = rep(0,1,length(thp));
 stdcon = rep(0,1,length(thp));
 stdnodf2 = rep(0,1,length(thp));
 stdcon2 = rep(0,1,length(thp));

 for(h in 1:length(thp)){
     meannodf[h] = mean(nestedness[h,]);
     meancon[h] = mean(connectance[h,]);
     stdnodf[h] = sd(nestedness[h,]);
  }
  
 for(h in 1:length(thp)){
     meannodf2[h] = mean(nestedness2[h,]);
     meancon2[h] = mean(connectance2[h,]);
     stdnodf2[h] = sd(nestedness2[h,]);
     stdcon2[h] = sd(connectance2[h,]);
  }
error = rep(0,1,length(thp));
error2 = rep(0,1,length(thp));
  upper_cl= rep(0,1,length(thp));
  upper_cl2= rep(0,1,length(thp));
  lower_cl = rep(0,1,length(thp));
   lower_cl2 = rep(0,1,length(thp));

 for(h in 1:length(thp)){
      error[h] = qnorm(0.975)*stdnodf[h]/sqrt(sims);
      upper_cl[h] = meannodf[h] + error[h];
      lower_cl[h] = meannodf[h] - error[h];
 }
for(h in 1:length(thp)){
      error2[h] = qnorm(0.975)*stdnodf2[h]/sqrt(sims);
      upper_cl2[h] = meannodf2[h] + error2[h];
      lower_cl2[h] = meannodf2[h] - error2[h];
 }
 library(plotrix);
 
plotCI(thp,meannodf,error,main="Basic")
plotCI(thp,meannodf,error,2*error,lwd=2,col="red",scol="blue",
main="Add colors to the points and error bars")
plotCI(thp,meannodf,error,main="Basic")
plotCI(thp,meannodf,error,2*error,lwd=2,col="red",scol="blue",main="Nestedness I")
plotCI(thp,meannodf2,error2,2*error2,lwd=2,col="red",scol="blue",main="Nestedness II")

 plot(thp,meannodf2,"b", xlabel="Threshold value", ylabel="Nestedness"); 
 lines(thp,upper_cl2); lines(thp,lower_cl2)
 
 plot(nestedness[1,1]);
 
library(fields);

jpeg('image_nullmodel1_amf.jpg')
image.plot(nestedness, main="Nestedness null model (I) AMF-Plant Network")
dev.off()
 
 jpeg('nullmodel1_nodf_amf.jpg')
 plotCI(thp,meannodf,error,2*error,lwd=2,col="red",scol="blue",main="Null model I")
 dev.off()
 
 jpeg('nullmodel2_nodf_amf.jpg')
 plotCI(thp,meannodf2,error2,2*error2,lwd=2,col="red",scol="blue",main="Null model II"); 
 dev.off()
 
 ## ----------------------------- Coords for modularity analysis ------------------------  ###

 ## Total number of matrices to run: 3(index: C, F, I)* 1 null model (CSR) + 2(index: C, F)* 2 null models ( ENV, SS) = 7 matrices
 ## Athough, we can skip calculation of mod for the SS, so it will be: 3 (CSR) + 2 (ENV) = 5 matrices

# load matrices:

load("sig_matrices_fij.RData")
load("binary_sig_cscores.RData")
n=18; m=15;
## --------- C-score matrices -------- ##

## CSR
coord1e=c(0,0);
for(i in 1:n){
  for(j in 1:m){
    if(bin_sig1[i,j]==1){
        coord1e=rbind(coord1e,c(i,j));
    }
  }
}
coord1e=coord1e[-1,]
write(coord1e,"coords_csr_c.dat",ncolumns=2) # file for SA program

## ENV
coord_env_c=c(0,0);
for(i in 1:n){
  for(j in 1:m){
    if(bin_sig3[i,j]==1){
        coord_env_c=rbind(coord_env_c,c(i,j));
    }
  }
}
coord_env_c=coord_env_c[-1,]
write(coord_env_c,"coords_env_c.dat",ncolumns=2) # file for SA program

## --------- Fij matrices -------- ##


## CSR
coord_csr_f=c(0,0);
for(i in 1:n){
  for(j in 1:m){
    if(binaryfij_eq[i,j]==1){
        coord_csr_f=rbind(coord_csr_f,c(i,j));
    }
  }
}
coord_csr_f=coord_csr_f[-1,]
write(coord_csr_f,"coords_csr_f.dat",ncolumns=2) # file for SA program

## ENV
coord_env_f=c(0,0);
for(i in 1:n){
  for(j in 1:m){
    if(binaryfij_env[i,j]==1){
        coord_env_f=rbind(coord_env_f,c(i,j));
    }
  }
}
coord_env_f=coord_env_f[-1,]
write(coord_env_f,"coords_env_f.dat",ncolumns=2) # file for SA program


## --------- Mutual information -------- ##

load("matrices_sig_MI_theo.RData")


## CSR
coord_csr_i=c(0,0);
for(i in 1:n){
  for(j in 1:m){
    if(inter_matrix[i,j]==1){
        coord_csr_i=rbind(coord_csr_i,c(i,j));
    }
  }
}
coord_csr_i=coord_csr_i[-1,]
write(coord_csr_i,"coords_csr_mutualinfo.dat",ncolumns=2) # file for SA program



## ---------------------------------
x=0;y=0;
for(i in 1:n){
   for(j in 1:m){
     if(binaryc_swap[i,j]==1){
         x=c(x,paste("p",i)); 
	 y=c(y,j);
     }
   }
 }
 coords <- cbind(x,y)
 write.table(coords,"coords.dat") # file for SA program
