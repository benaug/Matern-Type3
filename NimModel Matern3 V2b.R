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
  lambda <- D.intercept*pi.denom #Expected primary points
  N.primary ~ dpois(lambda) #Realized primary points
  for(i in 1:N.secondary){ #N.secondary secondary points
    age1[i] ~ dunif(0,1)
    #dunif() here implies uniform distribution within a grid cell
    #also tells nimble s's are in continuous space, not discrete
    s1[i,1] ~ dunif(xlim[1],xlim[2])
    s1[i,2] ~ dunif(ylim[1],ylim[2])
    #get cell s_i lives in using look-up table
    s1.cell[i] <- cells[trunc(s1[i,1]/res)+1,trunc(s1[i,2]/res)+1]
    #categorical likelihood for this cell, equivalent to zero's trick
    #also disallowing s's in non-habitat
    dummy.data1[i] ~ dCell(pi.cell[s1.cell[i]])
  }
  for(i in 1:M){ #M thinned points
    age2[i] ~ dunif(0,1)
    #dunif() here implies uniform distribution within a grid cell
    #also tells nimble s's are in continuous space, not discrete
    s2[i,1] ~ dunif(xlim[1],xlim[2])
    s2[i,2] ~ dunif(ylim[1],ylim[2])
    #get cell s_i lives in using look-up table
    s2.cell[i] <- cells[trunc(s2[i,1]/res)+1,trunc(s2[i,2]/res)+1]
    #categorical likelihood for this cell, equivalent to zero's trick
    #also disallowing s's in non-habitat
    dummy.data2[i] ~ dCell(pi.cell[s2.cell[i]])
  }
  D1[1:N.secondary,1:N.secondary] <- getD1(s1[1:N.secondary,1:2]) #distances between secondary points
  D.cross[1:N.secondary,1:M] <- getD.cross(s1[1:N.secondary,1:2],s2[1:M,1:2],N.thin=N.thin) #distances between secondary and thinned points
  #Matern Type III thinning process
  dummy.data3[1:N.secondary] ~ dThin(D1=D1[1:N.secondary,1:N.secondary],D.cross=D.cross[1:N.secondary,1:M],
                                     age1=age1[1:N.secondary],age2=age2[1:M],r=r,p.thin=p.thin,N.thin=N.thin)
  N.thin <- N.primary - N.secondary
})

