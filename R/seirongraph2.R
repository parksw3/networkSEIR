############################################
# Percolation on graph                     #
############################################

#install.packages("igraph")
rm(list=ls(all=TRUE)) 
library(igraph)
library(Matrix)

# track<-"/media/ptrap/suhome/inst/Rfiles/"
edges0<-read.table("../data/CA-CondMat.txt",header=FALSE)
edges0<-edges0[order(edges0[,1]), ] ###give edges right order
graph0 <- graph.data.frame(edges0,directed=TRUE)
degs=degree(graph0,v=V(graph0))/2 ####/2 because of summed in and oud degrees
N=length(V(graph0))

components <- decompose.graph(graph0,mode="weak")
length(components)
comp.sizes <- unlist(lapply(components,vcount))
comp.sizes.d <- table(comp.sizes)
comp.sizes.d

mu=mean(degs)
sigma=sd(degs)
kappa=sigma^2/mu+mu-1

#############Parameters to set#############
nsim<-100
lambda<-0.1 ###Infection rate
beta<-1 ###recovery rate
gamma<-0.5 ###rate of leaving exposed period

I0<-10 ###Initial number of infecteds
###S0=N-I0
startgen<-75 ### For estimating R0 from infection tree, we take the second generation at which the total number of infecteds is larger than startgen divided by the first
uppbound<-log(400)
lowbound<-log(200) ### we estimate alpha based on first time the number of infecteds is larger than e^lowbound and first time it is larger than e^uppbound
stepfreq<-50 ### the number of timesteps per time unit
maxtime<-40


######Simulation#############

simulationepidemic<-function()
{
  
infperiodvert<-rexp(N,beta) ### vector of infectious periods of length the number of vertices
expperiodvert<-rexp(N,gamma) ### vector of exposed periods of length the number of vertices
infperiodedg<-rep(infperiodvert,degs) ### the elements of infperiodvert get multiplicity according to degrees of vertices
expperiodedg<-rep(expperiodvert,degs) ### the elements of expperiodvert get multiplicity according to degrees of vertices
weightsir<-rexp(sum(degs),lambda) ###the time until first contact on an edge since start of infectious period of the starting vertex (if infected)
edges0a<-edges0[order(edges0[,2]), ]
edges0b<-cbind(edges0a,expperiodedg)
edges0c<-edges0b[order(edges0b[,1]), ]
expperiodedg1<-edges0c[,3]
weight<-weightsir+expperiodedg1 ###the weight of an edge is the time of first contact along edge + latent period of receiving vertex
edgesseir1<-cbind(edges0,weight,weightsir,infperiodedg)
edgesseir1<-subset(edgesseir1, edgesseir1[,4]<edgesseir1[,5]) ###we throw away all edges for which the infectious period of the starting vertex is shorter than the time until first contact along edge
edgesseir<-edgesseir1[,c(1:3)]  
graphseir <- graph.data.frame(edgesseir,directed=TRUE) ## we created a weighted infection graph, with weights the serial interval along the edges



Nseir=length(V(graphseir)) ##Number of vertices in graphseir (we do not consider the vertices which will have total degree 0 in graphseir)

v=sample.int(Nseir,I0,replace=FALSE,prob=NULL)  ###Choose I0 vertices uniformly at random

realsum<-function(x){
  realtime<-colSums(shortest.paths(graphseir,v=1:I0,to=1:length(V(graphseir)),mode = c("out"),weights = NULL)<=x)>0
  timesiz<-sum(realtime)
  return(timesiz)}
##Number of vertices within distance x of initial infectives (using edge weights)
genersum<-function(x){
  gen<-colSums(shortest.paths(graphseir,v=1:I0,to=1:length(V(graphseir)),mode = c("out"),weights = NA)<=x)>0
  gensiz<-sum(gen)
  return(gensiz)}
##Number of vertices within distance x of initial infectives (using number of edges)

vecreal<-c(I0)

repeat
{ lvecreal<-length(vecreal);
  vecreal<- cbind(vecreal,realsum(lvecreal/stepfreq));
  if ( (log(realsum(lvecreal/stepfreq)) > uppbound+0.5) | (lvecreal+1 > stepfreq*maxtime) ) 
    break;
}
##We compute the number of vertices infected up to time t until 0.5 + the first time we have e^uppbound infecteds 


lowreal<-vecreal[log(vecreal)<lowbound]
indreallow<-length(lowreal)

uppreal<-vecreal[log(vecreal)<uppbound]
indrealupp<-length(uppreal)

##give times (times stepfreq) when we reach e^lowbound and e^uppbound infecteds

obsdates=c(indreallow:indrealupp)
logvecreal=log(vecreal[obsdates])
obsdates2=obsdates/stepfreq

##reg<-lm(logvecreal ~ obsdates2)
##alphaest=coef(reg)[2]
ifelse(length(obsdates2)>1, alphaest<-coef(lm(logvecreal ~ obsdates2))[2], alphaest<-0)
##alphaest<-coef(lm(logvecreal ~ obsdates2))[2]



#alphaest<-stepfreq*(uppbound-lowbound)/(indrealupp-indreallow)
R0estreal <- (1 + alphaest/beta)*(1+alphaest/gamma) 
#R0estreal
##estimate of alpha and  R0 using alpha

vecgen<-c(I0)
repeat
{ lvecgen<-length(vecgen);
  vecgen<- cbind(vecgen,genersum(lvecgen));
  if ( (genersum(lvecgen)-genersum(lvecgen-1) > startgen) | (lvecgen > 25 ) )
    break;
}

##We compute the number of infecteds per generation up to 2 + the first generation we reach startgen infected individuals 

##lowgenvec<-vecgen[vecgen<startgen]
lowgen<-length(vecgen)
##R0estgen<- (genersum(lowgen+1)-genersum(lowgen))/(max((genersum(lowgen)-genersum(lowgen-1)),1))
R0estgen<- (genersum(lowgen+1)-genersum(1))/(max((genersum(lowgen)-I0),1))

#R0estgen
##estimate of R0 using epidemic network
return(array(c(alphaest, R0estreal, R0estgen)))
}

estimates<-rep(NA,3)
counter<-0

while(counter<nsim)
{
print(counter)
newrow<-simulationepidemic()
print(newrow)
ifelse(newrow[1]>0, counter<-counter+1, counter<-counter)
ifelse(newrow[1]>0, estimates<-rbind(estimates, newrow), estimates<-estimates) 
#print(estimates)
}

estimates <- estimates[-1,]

fn <- "trapman_sim.rda"

save("estimates", file = fn)
