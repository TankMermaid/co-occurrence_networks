#### Script to make figure of correlograms for three different null models:
#### CSR , ENV and SWAP

## Load data and variables

load(file="csr_corr.RData") # CSR moran's I values

load(file="env_corr.RData") # ENV moran's I values

load(file="swap_corr.RData") # Swap moran's I values

# observed moran I for Bromus inermis

valbinobs<- matrix(P[[1]],51*51,1, byrow=TRUE)

nmcorobs<- sp.correlogram(coordsfij,valbinobs,order=10, method="I",style = "W")

dev.new(width=14, height=7)
par(mfrow=c(1,3))

## observed versus CSR plant 1 correlogram
plot(moran_csr[[1]],type="l",col="grey", xlab="Lag", ylab="Moran's I",main="CSR",xlim=c(0,11), ylim=c(-0.05,0.15))
for(l in 1:sims){
 lines(moran_csr[[l]],type="l",col="grey")
}
lines(nmcorobs$res[,1],type="l",col="red",lwd=1.8)
lines(seq(0,11),rep(0,12))

## observed versus ENV plant 1 correlogram
plot(moran_env[[1]],type="l",col="grey", xlab="Lag", ylab="Moran's I",main="Env",xlim=c(0,11), ylim=c(-0.05,0.15) )
for(l in 1:sims){
 lines(moran_env[[l]],type="l",col="grey")
}
lines(nmcorobs$res[,1],type="l",col="red",lwd=1.8)
lines(seq(0,11),rep(0,12))

## observed versus swap plant 1 correlogram
plot(moran[[1]],type="l",col="grey", xlab="Lag", ylab="Moran's I",main="Swap",xlim=c(0,11), ylim=c(-0.05,0.15))
for(l in 1:sims){
 lines(moran[[l]],type="l",col="grey")
}
lines(nmcorobs$res[,1],type="l",col="red",lwd=1.8)
lines(seq(0,11),rep(0,12))

# environmental values versus observed
sims=1000;
plot(moran_env[[1]],type="l",col="grey", xlab="lags", ylab="Moran's I",main="Swap",xlim=c(0,11), ylim=c(-0.05,0.15))
for(l in 1:sims){
    lines(moran_env[[l]], type="l",col="grey")
}
lines(seq(0,11),rep(0,12))
lines(nmcorobs$res[,1],type="l",col="red",lwd=1.8)



