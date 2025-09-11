e2dist <-  function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

sim.Matern3V3 <-
  function(D.beta0=NA,D.beta1=NA,D.cov=NA,InSS=NA,r=NA,p.thin=NA,xlim=NA,ylim=NA,res=NA){
    #get expected number of primary points
    cellArea <- res^2
    lambda.cell <- InSS*exp(D.beta0 + D.beta1*D.cov)*cellArea
    lambda.N <- sum(lambda.cell)
    #simulate realized primary points
    N.primary <- rpois(1,lambda.N)
    
    #recreate some Dcov things so we can pass fewer arguments into this function
    x.vals <- seq(xlim[1]+res/2,xlim[2]-res/2,res) #x cell centroids
    y.vals <- seq(ylim[1]+res/2,ylim[2]-res/2,res) #y cell centroids
    dSS <- as.matrix(cbind(expand.grid(x.vals,y.vals)))
    cells <- matrix(1:nrow(dSS),nrow=length(x.vals),ncol=length(y.vals))
    n.cells <- nrow(dSS)
    n.cells.x <- length(x.vals)
    n.cells.y <- length(y.vals)
    
    #simulate a population of activity centers
    pi.cell <- lambda.cell/sum(lambda.cell)
    s.cell <- sample(1:n.cells,N.primary,prob=pi.cell,replace=TRUE)
    #distribute activity centers uniformly inside cells
    s <- matrix(NA,nrow=N.primary,ncol=2)
    for(i in 1:N.primary){
      tmp <- which(cells==s.cell[i],arr.ind=TRUE) #x and y number
      s[i,1] <- runif(1,x.vals[tmp[1]]-res/2,x.vals[tmp[1]+res/2])
      s[i,2] <- runif(1,y.vals[tmp[2]]-res/2,y.vals[tmp[2]+res/2])
    }
    
    #Matern3 thinning process
    age <- runif(N.primary,0,1)
    order <- order(age)
    #put points in age order. Don't really need to reorder for data simulation
    age <- age[order]
    s.cell <- s.cell[order]
    s <- s[order,]
    D <- e2dist(s,s)
    retain <- rep(0,N.primary)
    for(i in 1:N.primary){
      these.points <- which(D[i,]<r&retain) #don't need age since in age order
      n.points <- length(these.points)
      if(n.points>0){
        p.retain <- (1-p.thin)^n.points
        retain[i] <- rbinom(1,1,prob=p.retain)
      }else{
        retain[i] <- 1
      }
    }

    #discard unobserved
    observed <- which(retain==1)
    N.secondary <- length(observed)
    s.obs <- s[observed,]
    
    truth <- list(s=s,s.cell=s.cell,retain=retain,order=order,age=age)
    
    out <- list(N.primary=N.primary,N.secondary=N.secondary,s.obs=s.obs,r=r,
              xlim=xlim,ylim=ylim,x.vals=x.vals,y.vals=y.vals,dSS=dSS,cells=cells,
              n.cells=n.cells,n.cells.x=n.cells.x,n.cells.y=n.cells.y,s.cell=s.cell,
              D.cov=D.cov,InSS=InSS,res=res,cellArea=cellArea,lambda.N=lambda.N,
              truth=truth)
    return(out)
  }
