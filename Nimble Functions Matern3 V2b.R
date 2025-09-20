getCell <- nimbleFunction( #cell 0 not allowed in this model, but leaving in as an error check
  run = function(s=double(1),res=double(0),cells=integer(2),xlim=double(1),ylim=double(1)) {
    returnType(double(0))
    inout <- 1*(s[1]>xlim[1]&s[1]<xlim[2]&s[2]>ylim[1]&s[2]<ylim[2])
    if(inout==1){
      s.cell <- cells[trunc(s[1]/res)+1,trunc(s[2]/res)+1]
    }else{
      s.cell <- 0
    }
    return(s.cell)
  }
)

dCell <- nimbleFunction(
  run = function(x = double(0), pi.cell = double(0),log = integer(0)) {
    returnType(double(0))
    logProb <- log(pi.cell)
    return(logProb)
  }
)

#make dummy random number generator to make nimble happy
rCell <- nimbleFunction(
  run = function(n = integer(0),pi.cell = double(0)) {
    returnType(double(0))
    return(0)
  }
)

getD1 <- nimbleFunction(
  run = function(s=double(2)){
    returnType(double(2))
    M <- nimDim(s)[1]
    D <- matrix(0,M,M)
    for(i in 1:(M-1)){
      for(j in (i+1):M){
        D[i,j] <- sqrt((s[i,1]-s[j,1])^2+(s[i,2]-s[j,2])^2)
        D[j,i] <- D[i,j] #symmetric
      }
    }
    return(D)
  }
)

#distances between primary and secondary points
getD.cross <- nimbleFunction(
  run = function(s1=double(2),s2=double(2),N.thin=double(0)){
    returnType(double(2))
    N.secondary <- nimDim(s1)[1] #max secondary points
    M <- nimDim(s2)[1] #max thinned points
    D <- matrix(0,N.secondary,M)
    for(i in 1:N.secondary){
      for(j in 1:N.thin){
        D[i,j] <- sqrt((s1[i,1]-s2[j,1])^2+(s1[i,2]-s2[j,2])^2)
      }
    }
    return(D)
  }
)

getD1_i <- nimbleFunction(
  run = function(i=double(0), s1=double(2)){
    returnType(double(1))
    N.secondary <- nimDim(s1)[1]
    D <- rep(0,N.secondary)
    for(j in 1:N.secondary){
      if(j!=i){
        D[j] <-  sqrt((s1[i,1]-s1[j,1])^2+(s1[i,2]-s1[j,2])^2)
      }
    }
    return(D)
  }
)

getD.cross_i <- nimbleFunction(
  run = function(s1=double(1),s2=double(2),N.thin=double(0)){
    returnType(double(1))
    M <- nimDim(s2)[1] #max thinned points
    D.cross <- rep(0,M)
    for(j in 1:N.thin){
      D.cross[j] <- sqrt((s1[1]-s2[j,1])^2+(s1[2]-s2[j,2])^2)
    }
    return(D.cross)
  }
)

dThin <- nimbleFunction(
  run = function(x=double(1),D1=double(2),D.cross=double(2),
                 age1=double(1),age2=double(1),r=double(0),p.thin=double(0),N.thin=double(0),log=integer(0)){
    returnType(double(0))
    N.secondary <- nimDim(D1)[1]
    logProb <- 0
    #retained-retained pairs - for hardcore, make sure no n.points>0 or n.points==0
    these.points <- rep(0,N.secondary)
    for(i in 1:N.secondary){
      n.points <- 0
      for(j in 1:N.secondary){
        #secondary points within r where secondary age j is before secondary age i
        if(D1[i,j]<r&age1[j]<age1[i]){
          n.points <- n.points + 1
          these.points[n.points] <- j
        }
      }
      if(n.points>0){
        these.points2 <- these.points[1:n.points] #can't resize these.points
        log.p.retain <- log((1-p.thin)^n.points)
        logProb <- logProb + log.p.retain
      }#else, n.points=0, always retained
    }
    #thinned-retained pairs
    if(N.thin>0 & logProb != -Inf){#don't need to check thinned points if logProb is -Inf from secondary points
      these.points <- rep(0,N.secondary)
      for(j in 1:N.thin){
        n.points <- 0
        for(i in 1:N.secondary){
          #secondary points within r where secondary age is before thinned age
          if(D.cross[i,j]<r&age1[i]<age2[j]){
            n.points <- n.points + 1
            these.points[n.points] <- i
          }
        }
        if(n.points>0){ 
          these.points2 <- these.points[1:n.points] #can't resize these.points
          log.p.retain <- log((1-p.thin)^n.points)
          logProb <- logProb + log(1 - exp(log.p.retain))
        }else{ #all thinned points must have prior secondary points
          logProb <- -Inf
        } 
      }
    }
    return(logProb)
  }
)

# dummy RNG for nimble
rThin <- nimbleFunction(
  run = function(n=integer(0),D1=double(2),D.cross=double(2),
                 age1=double(1),age2=double(1),r=double(0),p.thin=double(0),N.thin=double(0)){
    returnType(double(1))
    N.secondary <- nimDim(D1)[1]
    return(rep(0,N.secondary))
  }
)


#sampler for thinned points and ages
thinnedSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    N.secondary <- control$N.secondary
    M <- control$M
    n.cells.x <- control$n.cells.x
    n.cells.y <- control$n.cells.y
    res <- control$res
    N.primary.node <- model$expandNodeNames("N.primary")
    age2.nodes <- model$expandNodeNames("age2")
    s2.nodes <- model$expandNodeNames("s2")
    s2.cell.nodes <- model$expandNodeNames("s2.cell")
    D1.nodes <- model$expandNodeNames("D1")
    D.cross.nodes <- model$expandNodeNames("D.cross")
    dummy.data2.nodes <- model$expandNodeNames("dummy.data2")
    dummy.data3.nodes <- model$expandNodeNames("dummy.data3")
    # calcNodes <- c(N.primary.node,D.nodes,s.nodes,s.cell.nodes,dummy.data.nodes,age.nodes,z.nodes)
    calcNodes <- c(N.primary.node,s2.nodes,s2.cell.nodes,dummy.data2.nodes,age2.nodes,D1.nodes,D.cross.nodes,dummy.data3.nodes)
  },
  run = function(){
    #simulate a new number of total points
    N.prop <- rpois(1,model$lambda[1])
    #simulate new locations and ages
    s2.cell.prop <- rep(0,N.prop)
    s2.prop <- matrix(0,N.prop,2)
    age2.prop <- rep(0,N.prop)
    for(i in 1:N.prop){
      s2.cell.prop[i] <- rcat(1,model$pi.cell)
      #propose x and y in new cell
      s.cell.x <- s2.cell.prop[i]%%n.cells.x
      s.cell.y <- floor(s2.cell.prop[i]/n.cells.x)+1
      if(s.cell.x==0){
        s.cell.x <- n.cells.x
        s.cell.y <- s.cell.y-1
      }
      xlim.cell <- c(s.cell.x-1,s.cell.x)*res
      ylim.cell <- c(s.cell.y-1,s.cell.y)*res
      s2.prop[i,] <- c(runif(1, xlim.cell[1], xlim.cell[2]), runif(1, ylim.cell[1], ylim.cell[2]))
      #simulate new ages
      age2.prop[i] <- runif(1,0,1)
    }
    
    #get the shadow
    D.prop <- getD.cross(s1=model$s1[1:N.secondary,1:2],s2=s2.prop,N.thin=N.prop)
    #shadow of secondary points: which points in G do you cast shadow on (so including primary and secondary points)
    shadow2D <- matrix(0,N.secondary,N.prop) #shadow of i secondary event cast on j thinned event
    for(i in 1:N.secondary){
      for(j in 1:N.prop){
        #points within r and later than focal secondary point times probability of deletion
        shadow2D[i,j] <- 1*(D.prop[i,j]<model$r[1]&model$age1[i]<age2.prop[j])
      }
    }
    #shadow of all points. prob of thinning
    shadow1D <- rep(0,N.prop)
    for(j in 1:N.prop){
      shadow1D[j] <- 1 - (1-model$p.thin[1])^sum(shadow2D[,j]) #prob of thinning
    }
    
    #which points are thinned?
    keep <- rbinom(N.prop,size=1,prob=shadow1D)
    N.thinned.prop <- sum(keep)
    #if we proposed new thinned points, fill those in
    if(N.thinned.prop>0){
      idx.keep <- which(keep==1) #which simulated points were thinned?
      #if we max out data augmentation, fill up to M (raise M and rerun if this happens, though!)
      if(N.thinned.prop>M){
        N.thinned.prop <- M
        idx.keep <- idx.keep[1:N.thinned.prop]
      }
      #fill in new points and ages
      idx.fill <- 1
      for(i in 1:N.thinned.prop){
        model$s2.cell[idx.fill] <<- s2.cell.prop[idx.keep[i]]
        model$s2[idx.fill,] <<- s2.prop[idx.keep[i],]
        model$age2[idx.fill] <<- age2.prop[idx.keep[i]]
        idx.fill <- idx.fill + 1
      }
    }else{
      idx.fill <- 1
    }
    
    #fill in augmented guys. need to do this to estimate density covariates without bias.
    if(idx.fill<=M){
      for(i in idx.fill:M){
        model$s2.cell[i] <<- rcat(1,model$pi.cell)
        #propose x and y in new cell
        s.cell.x <- model$s2.cell[i]%%n.cells.x
        s.cell.y <- floor(model$s2.cell[i]/n.cells.x)+1
        if(s.cell.x==0){
          s.cell.x <- n.cells.x
          s.cell.y <- s.cell.y-1
        }
        xlim.cell <- c(s.cell.x-1,s.cell.x)*res
        ylim.cell <- c(s.cell.y-1,s.cell.y)*res
        model$s2[i,] <<- c(runif(1, xlim.cell[1], xlim.cell[2]), runif(1, ylim.cell[1], ylim.cell[2]))
        model$age2[i] <<- runif(1,0,1) #dont really need to fill in new ages.
      }
    }
    
    model$N.thin[1] <<- N.thinned.prop
    model$N.primary[1] <<-  N.secondary + N.thinned.prop
    model$D.cross <<- getD.cross(model$s1,s2=model$s2,N.thin=model$N.thin[1])
    model$calculate(calcNodes)
    
    #copy back to mySaved to update logProbs which was not done above
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list( reset = function () {} )
)