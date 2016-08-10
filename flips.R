###### ---------------------- F L I P S -------------------- #####

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
 
 ### Create empty lists for the flips (only applied to AMF)
 A_flip = list(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15);
 
## Parameters:

  L=51;
  n=18; m=15;
  p <- rep(0,m);
  sets=10;
 
## Apply flips to AMF data

   for(i in 1:m){
      p[i]=0;
      for(l in 1:L){
         for(ll in 1:L){
            
          
            # Randomly change presence of AMF species    
              if(runif(1)>0.75){
                  if(A_flip[[i]][l,ll] == 1 ){
                     A_flip[[i]][l,ll] = 0;
                  }else{ 
                     A_flip[[i]][l,ll] = 1;
                  }
                  p[i]=p[i]+1;
              }
            
          } 
       }
       flip=A_flip[[i]]
       save(flip, file=sprintf("A_flip_%d.RData", i))
    }   
    


A_flip=vector("list",m);


## Save multiple files of A flip sets (100).
 for(i in 1:m){
   A_flip[[i]] = vector("list",sets);
   for(s in 1:sets){
      A_flip[[i]][[s]] = A_flip1[[i]];
      p[i]=0;
      for(l in 1:L){
         for(ll in 1:L){
            
            # Randomly change presence of AMF species    
              if(runif(1)>0.75){
                  if(A_flip[[i]][[s]][l,ll] == 1 ){
                     A_flip[[i]][[s]][l,ll] = 0;
                  }else{ 
                     A_flip[[i]][[s]][l,ll] = 1;
                  }
                  p[i]=p[i]+1;
              }
            
          } 
       }
    }
   aflip = A_flip[[i]]
   save(aflip, file=sprintf("A_flip_%d.RData", i))   
}


### This works!! FLIPS 

flips=250
p=0;
total=1000
posx = floor(runif(total,min=1,max=52))
posx2 = floor(runif(total,min=1,max=52))
posy = floor(runif(total,min=1,max=52))
posy2 = floor(runif(total,min=1,max=52))
  while(p < flips){   
       for(k in 1:total){   
            # Randomly change presence of AMF species    
              if(A_flip[[1]][posx[k],posy[k]] == 1 && A_flip[[1]][posx2[k],posy2[k]]==0){
                     A_flip[[1]][posx[k],posy[k]] = 0;
                     A_flip[[1]][posx2[k],posy2[k]] = 1;
                  
                  p=p+1;
               }
        }   
}
   
posx= vector("list",m)
posx2 = vector("list",m)
posy = vector("list",m)
posy2 = vector("list",m)

total=10000

ab = rep(0,m);
ab25 = rep(0,m);

for(i in 1:m){
  posx[[i]]= floor(runif(total,min=1,max=52));
  posx2[[i]]= floor(runif(total,min=1,max=52));
  posy[[i]]= floor(runif(total,min=1,max=52));
  posy2[[i]]= floor(runif(total,min=1,max=52));
  ab[i]=sum(A_flip[[i]]>0);
  ab25[i] = floor(ab[i]*0.15);
 p=0; flips=ab25[i];
 while(p < flips){   
       for(k in 1:total){   
            # Randomly change presence of AMF species    
              if(A_flip[[i]][posx[[i]][k],posy[[i]][k]] == 1 && A_flip[[i]][posx2[[i]][k],posy2[[i]][k]]==0){
                     A_flip[[i]][posx[[i]][k],posy[[i]][k]] = 0;
                     A_flip[[i]][posx2[[i]][k],posy2[[i]][k]] = 1;
                  
                  p=p+1;
               }
        }   
  } # end while loop

} # end loop species

save(A_flip, file="flips.RData")

# the script works but resample positions already shuffled
  
posx= vector("list",m)
posx2 = vector("list",m)
posy = vector("list",m)
posy2 = vector("list",m)

total=10000

ab = rep(0,m);
ab25 = rep(0,m)

tag= vector("list",m);
tag2 = vector("list",m);
L=51;

for(i in 1:m){
  tag[[i]]= matrix(0,L,L);
  tag2[[i]]= matrix(0,L,L);
  posx[[i]]= floor(runif(total,min=1,max=52));
  posx2[[i]]= floor(runif(total,min=1,max=52));
  posy[[i]]= floor(runif(total,min=1,max=52));
  posy2[[i]]= floor(runif(total,min=1,max=52));
  ab[i]=sum(A_flip[[i]]>0);
  ab25[i] = floor(ab[i]*0.15);
 p=0; flips=ab25[i];
 while(p < ab25[i]){   
       for(k in 1:total){   
            # Randomly change presence of AMF species    
              if(A_flip[[i]][posx[[i]][k],posy[[i]][k]] == 1 && A_flip[[i]][posx2[[i]][k],posy2[[i]][k]]==0 && tag[[i]][posx2[[i]][k],posy2[[i]][k]]==0 && tag2[[i]][posx[[i]][k],posy[[i]][k]] == 0){
                     A_flip[[i]][posx[[i]][k],posy[[i]][k]] = 0;
                     A_flip[[i]][posx2[[i]][k],posy2[[i]][k]] = 1;
                     tag[[i]][posx2[[i]][k],posy2[[i]][k]] = 2;
                     tag2[[i]][posx[[i]][k],posy[[i]][k]] = 2;
                  p=p+1;
               }
        }   
  } # end while loop

} # end loop species

