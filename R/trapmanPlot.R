library(igraph)
library(Matrix)

# track<-"/media/ptrap/suhome/inst/Rfiles/"
edges0<-read.table("../data/CA-CondMat.txt",header=FALSE)
edges0<-edges0[order(edges0[,1]), ] ###give edges right order
graph0 <- graph.data.frame(edges0,directed=TRUE)
degs=degree(graph0,v=V(graph0))/2 ####/2 because of summed in and oud degrees

mu=mean(degs)
sigma=sd(degs)
kappa=sigma^2/mu+mu-1

#############Parameters to set#############
lambda<-0.1 ###Infection rate
beta<-1 ###recovery rate
gamma<-0.5 ###rate of leaving exposed period

load("../data/trapman_sim.rda")

hist(estimates[,1],xlab="alpha",prob=TRUE,breaks=10,main="Histogram of alpha",xlim=c(0,4))
hist(estimates[,2],xlab="R0",prob=TRUE,breaks=10,main="Histogram of R0 based on alpha",xlim=c(1,5))
hist(estimates[,3],xlab="R0",prob=TRUE,breaks=10,main="Histogram of R0 based on infectiontree",xlim=c(1,4))
hist(estimates[,2]-estimates[,3])
hist(estimates[,2]/estimates[,3])
summary(estimates[,2])
summary(estimates[,3])
summary(estimates[,2]/estimates[,3])
plot(density(estimates[,3]),type="l",xlab="R0",ylim=c(0,1.4),main="Density of R0")
#lines(density(slopes[bool,7]),type="l",lty=2,col="blue")
#lines(density(estimates[,3]),type="l",lty=3,col="red")
lines(density(estimates[,2]),type="l",lty=2,col="red")
#abline(v=R0.mixing,col="red",lwd=2)
#segments(R0.network,0,R0.network,1.3,col="orange",lwd=2)
legend(3.2,1.45,c("homogeneous mixing", "first generations"), col=c("black", "red"), lty=c(1,2))
plot(estimates[,2],estimates[,3],abline(a=0,b=1,untf = FALSE))
regs<-lm(estimates[,3]~estimates[,2])
abline(regs,col="red")
estimates1<-estimates[order(estimates[,2]),]
estimates2<-estimates1[6:95,]
summary(estimates2[,2]/estimates2[,3])
plot(estimates2[,2],estimates2[,3],abline(a=0,b=1,untf = FALSE), main = paste(""),
     xlab = "R0 estimate based on real time growth rate", ylab = "R0 estimate based on real infection process")
hist(estimates2[,2]/estimates2[,3],
     main = paste(""),
     xlab = "Ratio of estimates of R0", cex.lab=1.5, cex.axis=1.5)
hist(estimates2[,2]-estimates2[,3],
     main = paste(""),
     xlab = "Difference of estimates of R0", cex.lab=1.5, cex.axis=1.5)
netest <- kappa*(beta+estimates2[,1])*(gamma+estimates2[,1])/(kappa*beta*gamma +(gamma+ estimates2[,1])*estimates2[,1])
plot(density(estimates2[,3]),type="l",xlab="Transmission potential",ylim=c(0,2.8),main="")
lines(density(netest),type="l",lty=3,col="blue")
lines(density(estimates2[,2]),type="l",lty=2,col="red")
#abline(v=R0.mixing,col="red",lwd=2)
segments(kappa*lambda/(lambda+beta),0,kappa*lambda/(lambda+beta),2.7,col="orange",lwd=2)
#legend(3.2,1.45,c("homogeneous mixing", "first generations"), col=c("black", "red"), lty=c(1,2))

plot(estimates[,2],estimates[,3],abline(a=0,b=1,untf = FALSE), main = paste(""),
     xlab = "R0 estimate based on real time growth rate", ylab = "R0 estimate based on real infection process")

netest<-kappa*(beta+estimates[,1])*(gamma+estimates[,1])/(kappa*beta*gamma +(gamma+ estimates[,1])*estimates[,1])
plot(density(estimates[,3]),type="l",xlab=expression("R"["0"]),ylim=c(0,2.8),main="")
lines(density(netest),type="l",lty=3,col="blue")
lines(density(estimates[,2]),type="l",lty=2,col="red")
#abline(v=R0.mixing,col="red",lwd=2)
segments(kappa*lambda/(lambda+beta),0,kappa*lambda/(lambda+beta),2.7,col="orange",lwd=2)
#legend(3.2,1.45,c("homogeneous mixing", "first generations"), col=c("black", "red"), lty=c(1,2))
boxplot(estimates[,2]/estimates[,3], names=c(""),boxwex=.3)
