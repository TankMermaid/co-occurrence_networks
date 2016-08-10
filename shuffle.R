### NULL MODEL 5: SHUFFLING PER COLUMN PLANT-AMF DATA ####

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
 
 Pc= vector("list", n);
 
 for(i in 1:n){
    Pc[[i]] = matrix(P[[i]],L*L,1,byrow=TRUE);
 }
 Pmat=Pc[[1]]
 for(i in 2:n){
    
    Pmat = cbind(Pmat,Pc[[i]])
 }
 
 