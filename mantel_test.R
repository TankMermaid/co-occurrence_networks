######  ----  ---- Mantel tests and correlograms ---- ---- ###### 

 
 
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
 
 
 n = 18; 
 m = 15;
 L=51;
 sims=1000;
 coords= expand.grid(seq(1,51),seq(1,51))
 
 a_1=matrix(a1, L*L, 1, byrow=TRUE)
 a_2=matrix(a2, L*L, 1, byrow=TRUE)
 a_3=matrix(a3, L*L, 1, byrow=TRUE)
 a_4=matrix(a4, L*L, 1, byrow=TRUE)
 a_5=matrix(a5, L*L, 1, byrow=TRUE)
 a_6=matrix(a6, L*L, 1, byrow=TRUE)
 a_7=matrix(a7, L*L, 1, byrow=TRUE)
 a_8=matrix(a8, L*L, 1, byrow=TRUE)
 a_9=matrix(a9, L*L, 1, byrow=TRUE)
 a_10=matrix(a10, L*L, 1, byrow=TRUE)
 a_11=matrix(a11, L*L, 1, byrow=TRUE)
 a_12=matrix(a12, L*L, 1, byrow=TRUE)
 a_13=matrix(a13, L*L, 1, byrow=TRUE)
 a_14=matrix(a14, L*L, 1, byrow=TRUE)
 a_15=matrix(a15, L*L, 1, byrow=TRUE)
 
 
 A=cbind(a_1,a_2,a_3,a_4,a_5,a_6,a_7,a_8,a_9,a_10,a_11,a_12,a_13,a_14,a_15)
 
 amf.dist = vegdist(A, method="bray", binary=TRUE, diag=FALSE, upper=TRUE, na.rm = TRUE) 
 
 ## Eliminate rows with zeros only
 dim(B)
 flag=0;
 #flagA=rep(0,L*L);
 for(l in 1:sims){
    if(sum(A[l,])==0){
	flag=c(flag,l);
    }
 }
 flag=flag[-1]
 
 A=A[-flag,]
 coords=coords[-flag,]
    
 amf.correlog2 <- mantel.correlog(amf.dist, XY=coords, cutoff=FALSE, r.type="spearman", nperm=99)

 
 dim(A)
  a22=a1;
  
  for(l in 1:L){
     for(ll in 1:L){
        if(a2[l,ll]==1) a22[l,ll]=2;
	if(a3[l,ll]==1) a22[l,ll]=3;
	if(a4[l,ll]==1) a22[l,ll]=4;
	if(a5[l,ll]==1) a22[l,ll]=5;
     }
  }
      
 
