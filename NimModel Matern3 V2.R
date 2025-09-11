NimModel <- nimbleCode({
  #priors
  D0 ~ dunif(0,100)
  D.beta1 ~ dnorm(0,sd=10)
  r ~ dunif(0,10) #interaction radius
  p.thin ~ dunif(0,1) #thinning probability inside interaction radius
  
  #Inhomogenous point process density model
  D.intercept <- D0*cellArea
  lambda.cell[1:n.cells] <- InSS[1:n.cells]*exp(D.beta1*D.cov[1:n.cells])
  pi.cell[1:n.cells] <- lambda.cell[1:n.cells]/pi.denom #expected proportion of total N in cell c
  pi.denom <- sum(lambda.cell[1:n.cells])
  lambda <- D.intercept*pi.denom #Expected N
  N.primary ~ dpois(lambda) #Realized N
  for(i in 1:M){ #M prethinned points
    age[i] ~ dunif(0,1) #continuous RV ordered to create integer order
    #dunif() here implies uniform distribution within a grid cell
    #also tells nimble s's are in continuous space, not discrete
    s[i,1] ~ dunif(xlim[1],xlim[2])
    s[i,2] ~ dunif(ylim[1],ylim[2])
    #get cell s_i lives in using look-up table
    s.cell[i] <- cells[trunc(s[i,1]/res)+1,trunc(s[i,2]/res)+1]
    #categorical likelihood for this cell, equivalent to zero's trick
    #also disallowing s's in non-habitat
    dummy.data[i] ~ dCell(pi.cell[s.cell[i]])
  }
  D[1:M,1:M] <- getD(s[1:M,1:2],z=z[1:M]) #updated with s, sigma, z
  #Matern Type III thinning process
  #retain=1 if secondary point, 0 for thinned points
  retain[1:M] ~ dThin(D=D[1:M,1:M],age=age[1:M],r=r,p.thin=p.thin,z=z[1:M])
})

