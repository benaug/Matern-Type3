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

getD <- nimbleFunction(
  run = function(s=double(2),z=double(1)){
    returnType(double(2))
    M <- nimDim(s)[1]
    D <- matrix(0,M,M)
    for(i in 1:(M-1)){
      if(z[i]==1){
        for(j in (i+1):M){
          if(z[j]==1){
            D[i,j] <- sqrt((s[i,1]-s[j,1])^2+(s[i,2]-s[j,2])^2)
            D[j,i] <- D[i,j] #symmetric
          }
        }
      }
    }
    return(D)
  }
)

#used in thinned point update
getD2 <- nimbleFunction(
  run = function(s1=double(2),s2=double(2)){
    returnType(double(2))
    n1 <- nimDim(s1)[1] #secondary points
    n2 <- nimDim(s2)[1] #thinned points
    D <- matrix(0,n1,n2)
    for(i in 1:n1){
      for(j in 1:n2){
        D[i,j] <- sqrt((s1[i,1]-s2[j,1])^2+(s1[i,2]-s2[j,2])^2)
      }
    }
  return(D)
  }
)

dThin <- nimbleFunction(
  run = function(x=double(1),D=double(2),age=double(1),r=double(0),
                 z=double(1),log=integer(0)){
    returnType(double(0))
    M <- nimDim(D)[1]
    logProb <- 0
    these.points <- rep(0,M)
    for(i in 1:M){
      if(z[i]==1){
        n.points <- 0
        for(j in 1:M){
          if(D[i,j]<r&x[j]==1&age[j]<age[i]){ #only previously retained points cast a shadow
            n.points <- n.points + 1
            these.points[n.points] <- j
          }
        }
        if(n.points>0){
          if(x[i]==1){
            logProb <- logProb - Inf
          }else{
            logProb <- logProb + 0
          }
        }else{
          if(x[i]==0){ #must be retained if z[i]==1 and n.points=0
            logProb <- -Inf
          }
        }
      }
    }
    return(logProb)
  }
)
# dummy RNG for nimble
rThin <- nimbleFunction(
  run = function(n=integer(0),D=double(2),age=double(1),r=double(0),z=double(1)){
    returnType(double(1))
    M <- nimDim(D)[1]
    return(rep(0,M))
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
    retain.nodes <- model$expandNodeNames("retain")
    age.nodes <- model$expandNodeNames("age")
    N.primary.node <- model$expandNodeNames("N.primary")
    D.nodes <- model$expandNodeNames("D")
    s.nodes <- model$expandNodeNames("s")
    s.cell.nodes <- model$expandNodeNames("s.cell")
    dummy.data.nodes <- model$expandNodeNames("dummy.data")
    calcNodes <- c(N.primary.node,D.nodes,s.nodes,s.cell.nodes,dummy.data.nodes,age.nodes,retain.nodes)
  },
  run = function(){
    #simulate a new number of total points
    N.prop <- rpois(1,model$lambda[1])
    #simulate new locations and ages
    s.cell.prop <- rep(0,N.prop)
    s.prop <- matrix(0,N.prop,2)
    age.prop <- rep(0,N.prop)
    for(i in 1:N.prop){
      s.cell.prop[i] <- rcat(1,model$pi.cell)
      #propose x and y in new cell
      s.cell.x <- s.cell.prop[i]%%n.cells.x
      s.cell.y <- floor(s.cell.prop[i]/n.cells.x)+1
      if(s.cell.x==0){
        s.cell.x <- n.cells.x
        s.cell.y <- s.cell.y-1
      }
      xlim.cell <- c(s.cell.x-1,s.cell.x)*res
      ylim.cell <- c(s.cell.y-1,s.cell.y)*res
      s.prop[i,] <- c(runif(1, xlim.cell[1], xlim.cell[2]), runif(1, ylim.cell[1], ylim.cell[2]))
      #simulate new times
      age.prop[i] <- runif(1,0,1)
    }
    
    #get the shadow
    D.prop <- getD2(s1=model$s[1:N.secondary,],s2=s.prop)
    #shadow of secondary points: which points in G do you cast shadow on (so including primary and secondary points)
    shadow2D <- matrix(0,N.secondary,N.prop) #shadow of i secondary event cast on j thinned event
    for(i in 1:N.secondary){
      for(j in 1:N.prop){
        #points within r and later than focal secondary point times probability of deletion
        shadow2D[i,j] <- 1*(D.prop[i,j]<model$r[1]&model$age[i]<age.prop[j])
      }
    }
    #shadow of all points. prob of thinning
    shadow1D <- rep(0,N.prop)
    for(j in 1:N.prop){
      shadow1D[j] <- 1-prod(1-shadow2D[,j])
    }
  
    #which points are thinned?
    keep <- rbinom(N.prop,size=1,prob=shadow1D)
    N.thinned.prop <- sum(keep)
    #if we proposed new thinned points, fill those in
    if(N.thinned.prop>0){
      idx.keep <- which(keep==1) #which simulated points were thinned?
      #if we max out data augmentation, fill up to M (raise M and rerun if this happens, though!)
      if(N.thinned.prop>(M-N.secondary)){
        N.thinned.prop <- M - N.secondary
        idx.keep <- idx.keep[1:N.thinned.prop]
      }
      #fill in new points and ages, set z=1
      idx.fill <- N.secondary + 1
      for(i in 1:length(idx.keep)){
        model$z[idx.fill] <<- 1
        model$s.cell[idx.fill] <<- s.cell.prop[idx.keep[i]]
        model$s[idx.fill,] <<- s.prop[idx.keep[i],]
        model$age[idx.fill] <<- age.prop[idx.keep[i]]
        idx.fill <- idx.fill + 1
      }
    }else{
      idx.fill <- N.secondary + 1
    }
    
    #fill in augmented guys. need to do this to estimate density covariates without bias.
    #set these z=0
    if(idx.fill<=M){
      for(i in idx.fill:M){
        model$z[i] <<- 0
        model$s.cell[i] <<- rcat(1,model$pi.cell)
        #propose x and y in new cell
        s.cell.x <- model$s.cell[i]%%n.cells.x
        s.cell.y <- floor(model$s.cell[i]/n.cells.x)+1
        if(s.cell.x==0){
          s.cell.x <- n.cells.x
          s.cell.y <- s.cell.y-1
        }
        xlim.cell <- c(s.cell.x-1,s.cell.x)*res
        ylim.cell <- c(s.cell.y-1,s.cell.y)*res
        model$s[i,] <<- c(runif(1, xlim.cell[1], xlim.cell[2]), runif(1, ylim.cell[1], ylim.cell[2]))
        model$age[i] <<- runif(1,0,1) #dont really need to fill in new ages.
      }
    }
    model$D <<- getD(s=model$s,z=model$z)
    model$N.primary[1] <<-  N.secondary + N.thinned.prop
    model$calculate(calcNodes)
    
    #copy back to mySaved to update logProbs which was not done above
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list( reset = function () {} )
)