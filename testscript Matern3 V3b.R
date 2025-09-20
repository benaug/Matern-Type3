#Softcore point process with probabilistic exclusion within radius r following halfnormal thinning kernel
#Here, I am not actually estimating r, but setting it to some large value beyond which the probability of 
#thinning from the thinning kernel is effectively 0. Basically used to trim likelihood calculations
#can try to estimate it, but multimodality/identifiability can be very bad.

library(nimble)
library(coda)
source("sim.Matern3 V3.R")
source("NimModel Matern3 V3b.R")
source("Nimble Functions Matern3 V3b.R")
source("mask.check.R")

#If using Nimble version 0.13.1 and you must run this line 
nimbleOptions(determinePredictiveNodesInModel = FALSE)
# #If using Nimble before version 0.13.1, run this line instead
# nimble:::setNimbleOption('MCMCjointlySamplePredictiveBranches', FALSE)

#get some colors
library(RColorBrewer)
cols1 <- brewer.pal(9,"Greens")

#simulate some data
sigma.thin <- 0.75 #halfnormal interaction scale parameter
#probability point i is thinned by point j is p.thin[i,j] <- exp(-(D[i,j])^2/sigma.thin^2)
#r is the maximum interaction radius beyond which we ignore interactions.
#I am assuming it is known and use it to speed up MCMC.
r <- 3*sigma.thin

### Habitat Covariate stuff###
#get x and y extent by buffering state space
xlim <- c(0,15)
ylim <- c(0,15)
#shift X, xlim, ylim, so lower left side of state space is (0,0)
#this is required to use efficient look-up table to find the cell number
#of a continuous location
x.shift <- xlim[1]
y.shift <- ylim[1]
xlim <- xlim-x.shift
ylim <- ylim-y.shift

res <- 0.25 #habitat grid resolution, length of 1 cell side
cellArea <- res^2 #area of one cell
x.vals <- seq(xlim[1]+res/2,xlim[2]-res/2,res) #x cell centroids
y.vals <- seq(ylim[1]+res/2,ylim[2]-res/2,res) #y cell centroids
dSS <- as.matrix(cbind(expand.grid(x.vals,y.vals)))
cells <- matrix(1:nrow(dSS),nrow=length(x.vals),ncol=length(y.vals))
n.cells <- nrow(dSS)
n.cells.x <- length(x.vals)
n.cells.y <- length(y.vals)

#create a density covariate
library(fields)
set.seed(1)
grid <- list(x=x.vals,y=y.vals) 
obj <- Exp.image.cov(grid=grid,aRange=2,setup=TRUE)
D.cov <- sim.rf(obj)
D.cov <- as.numeric(scale(D.cov)) #scale

par(mfrow=c(1,1),ask=FALSE)
image(x.vals,y.vals,matrix(D.cov,n.cells.x,n.cells.y),main="D.cov",xlab="X",ylab="Y",col=cols1)

#Density covariates
D.beta0 <- -0.5
D.beta1 <- 1
#can use a habitat mask
InSS <- rep(1,n.cells) #not excluding any cells from mask here
InSS[1:10] <- 0 #excluding these 10 cells for demonstration
image(x.vals,y.vals,matrix(InSS,n.cells.x,n.cells.y),main="Habitat Mask",col=cols1) #visualize habitat mask
#what is implied expected N in state space?
lambda.cell <- InSS*exp(D.beta0 + D.beta1*D.cov)*cellArea
sum(lambda.cell) #expected N in state space

image(x.vals,y.vals,matrix(lambda.cell,n.cells.x,n.cells.y),main="Expected Density",col=cols1)

#Simulate some data
set.seed(19)
data <- sim.Matern3V3(D.beta0=D.beta0,D.beta1=D.beta1,D.cov=D.cov,InSS=InSS,
                     r=r,sigma.thin=sigma.thin,xlim=xlim,ylim=ylim,res=res)
data$N.primary
mean(data$truth$retain) #percent of points retained

#what is the observed data?
# head(data$s.obs) #the points that were not thinned
data$N.secondary #number of these secondary points

#plot primary and secondary points
par(mfrow=c(2,1),ask=FALSE)
image(x.vals,y.vals,matrix(D.cov,n.cells.x,n.cells.y),main="Primary Points",
      xlab="X",ylab="Y",col=cols1)
points(data$truth$s,pch=16)
image(x.vals,y.vals,matrix(D.cov,n.cells.x,n.cells.y),main="Secondary Points",
      xlab="X",ylab="Y",col=cols1)
points(data$s.obs,pch=16)
par(mfrow=c(1,1))

#function to test for errors in mask set up. 
mask.check(dSS=data$dSS,cells=data$cells,n.cells=data$n.cells,n.cells.x=data$n.cells.x,
           n.cells.y=data$n.cells.y,res=data$res,xlim=data$xlim,ylim=data$ylim,
           x.vals=data$x.vals,y.vals=data$y.vals)

#Augment and initialize
M <- 600 #data augmentation level
N.secondary <- data$N.secondary
xlim <- data$xlim
ylim <- data$ylim
#primary points are observed
s1 <- data$s.obs
s1.cell <- rep(0,N.secondary)
for(i in 1:N.secondary){
  s1.cell[i] <- getCell(s=s1[i,],res=data$res,cells=data$cells,xlim=data$xlim,ylim=data$ylim)
}
age1.init <- 1:N.secondary/(N.secondary+M+1)
#thinned points
s2.init <- matrix(NA,nrow=M,ncol=2)
s2.cell.init <- sample(1:data$n.cells,M,replace=TRUE)
for(i in 1:M){
  tmp <- which(data$cells==s2.cell.init[i],arr.ind=TRUE) #x and y number
  s2.init[i,1] <- runif(1,data$x.vals[tmp[1]]-data$res/2,data$x.vals[tmp[1]+data$res/2])
  s2.init[i,2] <- runif(1,data$y.vals[tmp[2]]-data$res/2,data$y.vals[tmp[2]+data$res/2])
}
#If using a habitat mask, move any s's initialized in non-habitat above to closest habitat
e2dist  <-  function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}
alldists <- e2dist(s2.init,data$dSS)
alldists[,data$InSS==0] <- Inf
for(i in 1:M){
  this.cell <- data$cells[trunc(s2.init[i,1]/data$res)+1,trunc(s2.init[i,2]/data$res)+1]
  if(data$InSS[this.cell]==0){
    cands <- alldists[i,]
    new.cell <- which(alldists[i,]==min(alldists[i,]))
    s2.init[i,] <- data$dSS[new.cell,]
  }
}

#plot to make sure initialized activity centers are inside habitat mask
image(data$x.vals,data$y.vals,matrix(data$InSS,data$n.cells.x,data$n.cells.y))
points(s2.init,pch=16)

age2.init <- (N.secondary+1):(M+N.secondary)/(N.secondary+M+1)

D0.init <- N.secondary/(sum(data$InSS)*data$res^2)
Niminits <- list(N.primary=N.secondary,N.thin=0,D0=D0.init,sigma.thin=sigma.thin,r=data$r,
                 D.beta1=0,s2=s2.init,s2.cell=s2.cell.init,age1=age1.init,age2=age2.init)

#constants for Nimble
constants <- list(M=M,D.cov=data$D.cov,cellArea=data$cellArea,n.cells=data$n.cells,
                  xlim=data$xlim,ylim=data$ylim,res=data$res,N.secondary=N.secondary)

#supply data to nimble
dummy.data1 <- rep(1,N.secondary)
dummy.data2 <- rep(1,M)
dummy.data3 <- rep(1,N.secondary)
Nimdata <- list(s1=s1,s1.cell=s1.cell,InSS=data$InSS,cells=data$cells
                ,dummy.data1=dummy.data1,dummy.data2=dummy.data2,dummy.data3=dummy.data3)

# set parameters to monitor
parameters <- c('D0',"D.beta1","lambda","N.primary","sigma.thin","r")
nt <- 1 #thinning rate
parameters2 <- c("lambda.cell",'D0') #record D0 here for plotting
nt2 <- 5

#Build the model, configure the mcmc, and compile
start.time <- Sys.time()
# deregisterDistributions("dThin") #if you switch model version without restarting, deregister dThin
Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,inits=Niminits)
#if not using density covariate
#if using density covariate, use this and add block sampler for D0 and D.beta1 below
config.nodes <- c("sigma.thin","age1")
conf <- configureMCMC(Rmodel,monitors=parameters, thin=nt,
                      monitors2=parameters2, thin2=nt2,
                      useConjugacy = FALSE, nodes=config.nodes) 

conf$addSampler(target = paste("s2[1:",M,",1:2]"),
                type = 'thinnedSampler',control = list(M=M,N.secondary=data$N.secondary,n.cells.x=data$n.cells.x,n.cells.y=data$n.cells.y,
                                                       res=data$res),
                silent = TRUE)

conf$addSampler(target = c("D0","D.beta1"),
                type = 'AF_slice',control=list(adaptive=TRUE),silent = TRUE)

# Build and compile
Rmcmc <- buildMCMC(conf)
# runMCMC(Rmcmc,niter=10) #this will run in R, used for better debugging
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

# Run the model.
start.time2 <- Sys.time()
Cmcmc$run(2500,reset=FALSE) #can keep running this line to extend sampler
end.time <- Sys.time()
end.time-start.time  # total time for compilation, replacing samplers, and fitting
end.time-start.time2 # post-compilation run time

mvSamples <- as.matrix(Cmcmc$mvSamples)
plot(mcmc(mvSamples[500:nrow(mvSamples),]))

#some targets
data$N.primary
sigma.thin
r

#compare simulated data to final iteration
par(mfrow=c(2,1),ask=FALSE)
image(x.vals,y.vals,matrix(Cmodel$lambda.cell,n.cells.x,n.cells.y),main="True Data",col=cols1)
points(data$truth$s,pch=16)
points(data$truth$s[data$truth$retain==1,],pch=16,col="darkred")
image(x.vals,y.vals,matrix(Cmodel$lambda.cell,n.cells.x,n.cells.y),main="Current Iteration",col=cols1)
points(Cmodel$s[Cmodel$z==1,],pch=16)
points(Cmodel$s[1:data$N.secondary,],pch=16,col="darkred")


#posterior correlation
tmp <- cor(mvSamples[250:nrow(mvSamples),])
round(tmp,2)

mvSamples2  <-  as.matrix(Cmcmc$mvSamples2)
lambda.cell.idx <- grep("lambda.cell",colnames(mvSamples2))
D0.idx <- grep("D0",colnames(mvSamples2))
burnin2 <- 10

#compare expected D plot to truth
#image will show posterior means
lambda.cell.post <- cellArea*mvSamples2[burnin2:nrow(mvSamples2),D0.idx]*mvSamples2[burnin2:nrow(mvSamples2),lambda.cell.idx]
lambda.cell.ests <- colMeans(lambda.cell.post)
#remove non-habitat
lambda.cell.ests[InSS==0] <- NA
lambda.cell[InSS==0] <- NA

par(mfrow=c(1,1),ask=FALSE)
zlim <- range(c(lambda.cell,lambda.cell.ests),na.rm=TRUE) #use same zlim for plots below
#truth
image(x.vals,y.vals,matrix(lambda.cell,n.cells.x,n.cells.y),main="Expected Density",zlim=zlim)
#estimate, posterior means
image(x.vals,y.vals,matrix(lambda.cell.ests,n.cells.x,n.cells.y),main="Expected Density",zlim=zlim)
