##  R script for AMF-Plant null models development ##

## Data from John Klironomos

## Structure of the code:
# 1) Importing spatial matrix data for all species
# 2) Implementing functions of the null models
# 2.1) Null model I: Equiprobable
# 2.2) Null model II: Fixed rows
# 2.3) Null model III: sequential swap 
# 2.4) Null model IV: Environ. constrained

# For each null model 1000 null matrices will be produced for a total of 4000 null matrices
# Z-scores will calculated and upper-lower tail from null model distributions of C-score and fij (and asymmetry)

# 3) Calculate co-occurrence matrices (n*m) with step function and cont. function 

# 4) Threshold loop analysis, th->[0,1] and generate binary matrices (adjacency matrices)

# 5) Finally: estimate nestedness( Brualdi, Temp, NODF), connectance and modularity

## Outputs: C.scores, fij distributions and stats; nestedness, connectance and modularity profile with thresholds
## Plot networks, profiles, binary matrices, histograms
## 

# files needed: *.dat from plants and AMF, postlda.RData

## Importing data

# Plants data

 p1 <- matrix(scan("p1.dat", n = 51*51,sep="\t",fill=TRUE,blank.lines.skip=FALSE), 51, 51, byrow = TRUE)
 p2 <- matrix(scan("p2.dat", n = 51*51,sep="\t",fill=TRUE,blank.lines.skip=FALSE), 51, 51, byrow = TRUE)
 p3 <- matrix(scan("p3.dat", n = 51*51,sep="\t",fill=TRUE,blank.lines.skip=FALSE), 51, 51, byrow = TRUE)
 p4 <- matrix(scan("p4.dat", n = 51*51,sep="\t",fill=TRUE,blank.lines.skip=FALSE), 51, 51, byrow = TRUE)
 p5 <- matrix(scan("p5.dat", n = 51*51,sep="\t",fill=TRUE,blank.lines.skip=FALSE), 51, 51, byrow = TRUE)
 p6 <- matrix(scan("p6.dat", n = 51*51,sep="\t",fill=TRUE,blank.lines.skip=FALSE), 51, 51, byrow = TRUE)
 p7 <- matrix(scan("p7.dat", n = 51*51,sep="\t",fill=TRUE,blank.lines.skip=FALSE), 51, 51, byrow = TRUE)
 p8 <- matrix(scan("p8.dat", n = 51*51,sep="\t",fill=TRUE,blank.lines.skip=FALSE), 51, 51, byrow = TRUE)
 p9 <- matrix(scan("p9.dat", n = 51*51,sep="\t",fill=TRUE,blank.lines.skip=FALSE), 51, 51, byrow = TRUE)
 p10 <- matrix(scan("p10.dat", n = 51*51,sep="\t",fill=TRUE,blank.lines.skip=FALSE), 51, 51, byrow = TRUE)
 p11 <- matrix(scan("p11.dat", n = 51*51,sep="\t",fill=TRUE,blank.lines.skip=FALSE), 51, 51, byrow = TRUE)
 p12 <- matrix(scan("p12.dat", n = 51*51,sep="\t",fill=TRUE,blank.lines.skip=FALSE), 51, 51, byrow = TRUE)
 p13 <- matrix(scan("p13.dat", n = 51*51,sep="\t",fill=TRUE,blank.lines.skip=FALSE), 51, 51, byrow = TRUE)
 p14 <- matrix(scan("p14.dat", n = 51*51,sep="\t",fill=TRUE,blank.lines.skip=FALSE), 51, 51, byrow = TRUE)
 p15 <- matrix(scan("p15.dat", n = 51*51,sep="\t",fill=TRUE,blank.lines.skip=FALSE), 51, 51, byrow = TRUE)
 p16 <- matrix(scan("p16.dat", n = 51*51,sep="\t",fill=TRUE,blank.lines.skip=FALSE), 51, 51, byrow = TRUE)
 p17 <- matrix(scan("p17.dat", n = 51*51,sep="\t",fill=TRUE,blank.lines.skip=FALSE), 51, 51, byrow = TRUE)
 p18 <- matrix(scan("p18.dat", n = 51*51,sep="\t",fill=TRUE,blank.lines.skip=FALSE), 51, 51, byrow = TRUE)
 
  # Mycorrizae data
 a1 <- matrix(scan("a1.dat", n = 51*51,sep="\t",fill=TRUE,blank.lines.skip=FALSE), 51, 51, byrow = TRUE)
 a2 <- matrix(scan("a2.dat", n = 51*51,sep="\t",fill=TRUE,blank.lines.skip=FALSE), 51, 51, byrow = TRUE)
 a3 <- matrix(scan("a3.dat", n = 51*51,sep="\t",fill=TRUE,blank.lines.skip=FALSE), 51, 51, byrow = TRUE)
 a4 <- matrix(scan("a4.dat", n = 51*51,sep="\t",fill=TRUE,blank.lines.skip=FALSE), 51, 51, byrow = TRUE)
 a5 <- matrix(scan("a5.dat", n = 51*51,sep="\t",fill=TRUE,blank.lines.skip=FALSE), 51, 51, byrow = TRUE)
 a6 <- matrix(scan("a6.dat", n = 51*51,sep="\t",fill=TRUE,blank.lines.skip=FALSE), 51, 51, byrow = TRUE)
 a7 <- matrix(scan("a7.dat", n = 51*51,sep="\t",fill=TRUE,blank.lines.skip=FALSE), 51, 51, byrow = TRUE)
 a8 <- matrix(scan("a8.dat", n = 51*51,sep="\t",fill=TRUE,blank.lines.skip=FALSE), 51, 51, byrow = TRUE)
 a9 <- matrix(scan("a9.dat", n = 51*51,sep="\t",fill=TRUE,blank.lines.skip=FALSE), 51, 51, byrow = TRUE)
 a10 <- matrix(scan("a10.dat", n = 51*51,sep="\t",fill=TRUE,blank.lines.skip=FALSE), 51, 51, byrow = TRUE)
 a11 <- matrix(scan("a11.dat", n = 51*51,sep="\t",fill=TRUE,blank.lines.skip=FALSE), 51, 51, byrow = TRUE)
 a12 <- matrix(scan("a12.dat", n = 51*51,sep="\t",fill=TRUE,blank.lines.skip=FALSE), 51, 51, byrow = TRUE)
 a13 <- matrix(scan("a13.dat", n = 51*51,sep="\t",fill=TRUE,blank.lines.skip=FALSE), 51, 51, byrow = TRUE)
 a14 <- matrix(scan("a14.dat", n = 51*51,sep="\t",fill=TRUE,blank.lines.skip=FALSE), 51, 51, byrow = TRUE)
 a15 <- matrix(scan("a15.dat", n = 51*51,sep="\t",fill=TRUE,blank.lines.skip=FALSE), 51, 51, byrow = TRUE)
 
 P = list(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18);
 #P = append(p18, P);
 A = list(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15);
 
## Parameters:

L=51;
n=18; m=15;

sims=1000;
thp = seq(0,1,0.01);
tha = seq(0,1,0.01);
## Null models

# Null model I: Equiprobable

nullequip <- function(M){
   
   dimensions = dim(M);
   x = dimensions[1];
   y = dimensions[2];
  # set.seed(2008)
    sumM<-sum(M>0)
    relsum = sumM/(x*y);
    
    Mrand = matrix(0,x,y);
    
    for (k in 1:x){
     for (l in 1:y){
         if( relsum > runif(1) ){
	     Mrand[k,l] = 1;
	  }
	}
    }
    
    relsum2 = sum(Mrand>0)/(x*y)
    
     return(Mrand)
}

# Null model II: fixed rows

nullfixed <- function(M){
   
   dimensions = dim(M);
   x = dimensions[1];
   y = dimensions[2];
   #set.seed(2008)
   sum = 0;
   sumrow=rep(0,x);
   
   for (k in 1:x){
     for (l in 1:y){
         if( M[k,l] == 1){
	     sum = sum +1;
	  }
	}
     sumrow[k]=sum(M[k,]);
    }
    relsum = sum/(x*y);
    
    Mrand = matrix(0,x,y);
    
    
  for(i in 1:x){
      j=1;
      while( sumrow[i]!= sum(Mrand[i,]) && j<=y ){
        if( runif(1) < relsum){
            Mrand[i,j]=1;
         }
         if(j==y){
           j=1;
         }
         j=j+1;
     }
   }
   
   relsum2 = sum(Mrand>0)/(x*y)
    
     return(Mrand)
}


# Null model III: sequential swap

swapnull <- function(M,SWAPS){

   dimensions = dim(M);
   x = dimensions[1];
   y = dimensions[2];

   #set.seed(2008)

   relsum=(sum(M>0))/(x*y)
   
   krow = rep(0,SWAPS)

   lrow = rep(0,SWAPS)

   kcol = rep(0,SWAPS)

   lcol = rep(0,SWAPS)
   
   # Create vector of positions for swappings
   for(i in 1:SWAPS){
     krow[i] = floor(runif(1)*100);

     while(krow[i]<1 || krow[i]>x){
       krow[i] = floor(runif(1)*100);
     }

     lrow[i] = floor(runif(1)*100);

     while(lrow[i]<1 || lrow[i]>x){   
      lrow[i] = floor(runif(1)*100);
     }

     kcol[i] = floor(runif(1)*100);

     while(kcol[i]<1 || kcol[i]>y){  
      kcol[i] = floor(runif(1)*100);
     }

     lcol[i] = floor(runif(1)*100);
     while(lcol[i]<1 || lcol[i]>y){  
       lcol[i] = floor(runif(1)*100);
     }
   }

   Mr = M;
     # Loop for swappings 
     for( s in 1:SWAPS){
      
         kr =krow[s];
         lr = lrow[s];
         kc = kcol[s];
         lc = lcol[s];  

       if( (sum(Mr[kr,kc] + Mr[kr,lc])==1) && (sum(Mr[lr,lc]+Mr[lr,kc])==1) && (sum(Mr[lr,lc]+Mr[kr,lc])==1)  && (sum(Mr[lr,kc]+Mr[kr,kc])==1) ){
         # temp variables
          p= Mr[kr,kc]; q = Mr[kr,lc]; r=Mr[lr,lc]; t=Mr[lr,kc];
         # swapping
          Mr[kr,kc] = q; Mr[kr,lc]=p;
          Mr[lr,lc] = t; Mr[lr,kc]=r;
        } 
     }  
     
     return(Mr)

}


# Null model IV: environ. constrained

load("post_lda_plants.RData") # postzp,postza

envnull <- function(post,L){  # post vector of site-specific distribution

    null<-sapply(post, function(x) if(x>runif(1)){ x=1}else{x=0})
    null1<-null[1:L];   # transforming vector of presences into matrix L*L
    for(n in 2:L-1)  null1= rbind(null1,null[(n*L+1):(n*L+L)])  

    return(null1)   # return matrix
}



## Creating lists for null matrices

# Plant's list null matrices
Peq = vector("list", n);
Pfix = vector("list", n);
Pswap = vector("list", n);
Penv = vector("list", n);

# Animal's list null matrices
Aeq = vector("list", m);
Afix = vector("list", m);
Aswap = vector("list", m);
Aenv = vector("list", m);

# Create empty matrix of null model spatial matrices
for( i in 1:n){
     Peq[[i]] = vector("list", sims);
     Pfix[[i]]= vector("list", sims);
     Pswap[[i]]= vector("list", sims);
     Penv[[i]]= vector("list", sims);
     
}

for( j in 1:m){
     Aeq[[j]] = vector("list", sims);
     Afix[[j]] = vector("list", sims);
     Aswap[[j]] = vector("list", sims);
     Aenv[[j]] = vector("list", sims);
}

SWAPS=500
# Create matrices from null model plants
for( i in 1:n){ 
  for( l in 1:sims){
        Peq[[i]][[l]] = nullequip(P[[i]]);
        Pfix[[i]][[l]] = nullfixed(P[[i]])
	Pswap[[i]][[l]] = swapnull(P[[i]],SWAPS)
	Penv[[i]][[l]] = envnull(postzp[[i]],L)
  }
}
  # Create matrices from null model animals
for( j in 1:m){
  for( l in 1:sims){
       Aeq[[j]][[l]] = nullequip(A[[j]]);
       Afix[[j]][[l]] = nullfixed(A[[j]]);
       Aswap[[j]][[l]] = swapnull(A[[j]],SWAPS);
       Aenv[[j]][[l]] = envnull(postza[[j]],L);
  }
}

# Matrices created !! ---- time consumed ?

# Saving files of null model spatial matrices

save(Peq,Pfix,Pswap,Penv,file="Pnullmatrices.RData")
save(Aeq,Afix,Aswap,Aenv,file="anullmatrices.RData")


## Create matrices S's for each simulation, product and sum of presences

Seq = vector("list", sims);
Sfix = vector("list", sims);
Sswap = vector("list", sims);
Senv = vector("list", sims);

for (l in 1:sims){
    Seq[[l]] = array(0, c(1,ncol(p1),n,m)); 
    Sfix[[l]] = array(0, c(1,ncol(p1),n,m)); 
    Sswap[[l]] = array(0, c(1,ncol(p1),n,m)); 
    Senv[[l]] = array(0, c(1,ncol(p1),n,m)); 
}
 
for (l in 1:sims){
   for(i in 1:n){
      for(j in 1:m){
	 for(k in 1:ncol(p1)){
	      Seq[[l]][,k,i,j] = sum(Peq[[i]][[l]][,k]*Aeq[[j]][[l]][,k]);
              Sfix[[l]][,k,i,j]= sum(Pfix[[i]][[l]][,k]*Afix[[j]][[l]][,k]);
	      Sswap[[l]][,k,i,j]= sum(Pswap[[i]][[l]][,k]*Aswap[[j]][[l]][,k]);
	      Senv[[l]][,k,i,j]= sum(Penv[[i]][[l]][,k]*Aenv[[j]][[l]][,k]);
	 }
      }
   }
}

save(Seq,Sfix,Sswap,Senv,file="S_overlap.RData")

# Create MP's and MA's matrices , MA -> fji=Nij/Nj --- MP -> fij=Nij/Ni  

MPeq = vector("list",sims);
MPfix = vector("list",sims);
MPswap = vector("list",sims);
MPenv = vector("list",sims);
MAeq = vector("list",sims);
MAfix = vector("list",sims);
MAswap = vector("list",sims);
MAenv = vector("list",sims);


for(l in 1:sims){
   MPeq[[l]] = matrix(0,n,m);
   MPfix[[l]] =  matrix(0,n,m);
   MPswap[[l]] =  matrix(0,n,m);
   MPenv[[l]] =  matrix(0,n,m);
}

for(l in 1:sims){
   MAeq[[l]] = matrix(0,m,n);
   MAfix[[l]] = matrix(0,m,n);
   MAswap[[l]] = matrix(0,m,n);
   MAenv[[l]] = matrix(0,m,n);
}
## fij and fji -----> fij=Nij/Ni ---- fji=Nij/Nj ----- S[[l]][,,i,j] = Nij --> co-occurrence
for (l in 1:sims){
   for(i in 1:n){
      for(j in 1:m){
           MPeq[[l]][i,j] = sum(Seq[[l]][,,i,j])/sum(Peq[[i]][[l]]);
           MPfix[[l]][i,j] = sum(Sfix[[l]][,,i,j])/sum(Pfix[[i]][[l]]);
	   MPswap[[l]][i,j] = sum(Sswap[[l]][,,i,j])/sum(Pswap[[i]][[l]]);
	   MPenv[[l]][i,j] = sum(Sswap[[l]][,,i,j])/sum(Penv[[i]][[l]]);
       }
   }
}

for (l in 1:sims){
   for(j in 1:m){
      for(i in 1:n){
           MAeq[[l]][j,i] = sum(Seq[[l]][,,i,j])/sum(Aeq[[j]][[l]]);
           MAfix[[l]][j,i] = sum(Sfix[[l]][,,i,j])/sum(Afix[[j]][[l]]);
	   MAswap[[l]][j,i] = sum(Sswap[[l]][,,i,j])/sum(Aswap[[j]][[l]]);
	   MAenv[[l]][j,i] = sum(Senv[[l]][,,i,j])/sum(Aenv[[j]][[l]]);
       }
   }
}

save(MAeq,MAfix,MAswap,MAenv,file="overlap_p.RData");
save(MPeq,MPfix,MPswap,MPenv,file="overlap_a.RData");



## Threshold loop for all simulations

## Undirected interaction stregnth -- Mean between MP and MA
 
MEANeq = vector("list", sims);
MEANfix = vector("list", sims); 
MEANswap = vector("list", sims);
MEANenv = vector("list", sims);

# C-scores
Ceq = vector("list", sims);
Cfix = vector("list", sims);
Cswap = vector("list", sims);
Cenv = vector("list", sims);

# Asymmetry

Asymeq = vector("list", sims);
Asymfix = vector("list", sims);
Asymswap = vector("list", sims);
Asymenv = vector("list", sims);

for( l in 1:sims){
   MEANeq[[l]] = matrix(0,n,m);
   MEANfix[[l]] = matrix(0,n,m);
   MEANswap[[l]] = matrix(0,n,m);
   MEANenv[[l]] = matrix(0,n,m);
   Ceq[[l]] =  matrix(0,n,m);
   Cfix[[l]] =  matrix(0,n,m);
   Cswap[[l]] =  matrix(0,n,m);
   Cenv[[l]] =  matrix(0,n,m);
   Asymeq[[l]] = matrix(0,n,m);
   Asymfix[[l]] = matrix(0,n,m);
   Asymswap[[l]] = matrix(0,n,m);
   Asymenv[[l]] = matrix(0,n,m);
}

## mean  Fij's co-occurrence, C-scores and asymmetries (C-scores have to be normalized by plants and animals totals)
for(l in 1:sims){
 for(i in 1:n){
   for(j in 1:m){
       MEANeq[[l]][i,j] = (MPeq[[l]][i,j] + MAeq[[l]][j,i])/2;
       Ceq[[l]][i,j] = (sum(P[[i]]) - sum(Seq[[l]][,,i,j]))*(sum(A[[j]]) - sum(Seq[[l]][,,i,j]))/(n*m);
       Asymeq[[l]][i,j] = abs(MPeq[[l]][i,j] - MAeq[[l]][j,i])/(max(MPeq[[l]][i,j],MAeq[[l]][j,i]))
       MEANfix[[l]][i,j] = (MPfix[[l]][i,j] + MAfix[[l]][j,i])/2; 
       Cfix[[l]][i,j] = (sum(P[[i]]) - sum(Sfix[[l]][,,i,j]))*(sum(A[[j]]) - sum(Sfix[[l]][,,i,j]))/(n*m)
       Asymfix[[l]][i,j] = abs(MPfix[[l]][i,j] - MAfix[[l]][j,i])/(max(MPfix[[l]][i,j],MAfix[[l]][j,i]))
       MEANswap[[l]][i,j] = (MPswap[[l]][i,j] + MAswap[[l]][j,i])/2; 
       Cswap[[l]][i,j] = (sum(P[[i]]) - sum(Sswap[[l]][,,i,j]))*(sum(A[[j]]) - sum(Sswap[[l]][,,i,j]))/(n*m)
       Asymswap[[l]][i,j] = abs(MPswap[[l]][i,j] - MAswap[[l]][j,i])/(max(MPswap[[l]][i,j],MAswap[[l]][j,i]))
       MEANenv[[l]][i,j] = (MPenv[[l]][i,j] + MAenv[[l]][j,i])/2; 
       Cenv[[l]][i,j] = (sum(P[[i]]) - sum(Senv[[l]][,,i,j]))*(sum(A[[j]]) - sum(Senv[[l]][,,i,j]))/(n*m)
       Asymenv[[l]][i,j] = abs(MPenv[[l]][i,j] - MAenv[[l]][j,i])/(max(MPenv[[l]][i,j],MAenv[[l]][j,i]))
   }
 }
}

save(MEANeq,Ceq,Asymeq,file="co_occurrence_eq.RData")
save(MEANfix,Cfix,Asymfix,file="co_occurrence_fix.RData")
save(MEANswap,Cswap,Asymswap,file="co_occurrence_swap.RData")
save(MEANenv,Cenv,Asymenv,file="co_occurrence_env.RData")

## Create empty binary matrices n=100 binary matrices
Bineq = vector("list", length(thp));
Binfix = vector("list", length(thp));
Binswap = vector("list", length(thp));
Binenv = vector("list", length(thp));

for(h in 1:length(thp)){
    Bineq[[h]] = vector("list",sims);
    Binfix[[h]] = vector("list",sims);
    Binswap[[h]] = vector("list",sims);
    Binenv[[h]] = vector("list",sims);
}
# Mean binary matrices for all simulations and thresholds
for(h in 1:length(thp)){
   for(l in 1:sims){
      Bineq[[h]][[l]] = matrix(0,n,m);
      Binfix[[h]][[l]] = matrix(0,n,m);
      Binswap[[h]][[l]] = matrix(0,n,m);
      Binenv[[h]][[l]] = matrix(0,n,m);
      for(i in 1:n){
        for(j in 1:m){
	    if(MEANeq[[l]][i,j]>thp[h]){
	        Bineq[[h]][[l]][i,j]=1;
	    }
	    if(MEANfix[[l]][i,j]>thp[h]){
	        Binfix[[h]][[l]][i,j]=1;
	    }
	    if(MEANswap[[l]][i,j]>thp[h]){
	        Binswap[[h]][[l]][i,j]=1;
	    }
	    if(MEANenv[[l]][i,j]>thp[h]){
	        Binenv[[h]][[l]][i,j]=1;
	    }
        }
      }
    }
 }

save(Bineq,file="binary_eq.RData");
save(Binfix,file="binary_fix.RData");
save(Binswap,file="binary_swap.RData");
save(Binenv,file="binary_env.RData");
