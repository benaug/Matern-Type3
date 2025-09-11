# Matern-Type3
Matérn Type III repulsive point processes in NIMBLE

This repository contains generalized Matérn Type III repulsive point process models as described by Adams 2009, Rao et al. 2017 and others. See
pages 36-42 of Adams 2009 for a good introduction of Matérn repulsive point processes.

Adams, Ryan Prescott. Kernel methods for nonparametric Bayesian inference of probability densities and point processes. Diss. University of Cambridge, 2009.

Rao, Vinayak, Ryan P. Adams, and David D. Dunson. "Bayesian inference for Matérn repulsive processes." Journal of the Royal Statistical Society Series B: Statistical Methodology 79.3 (2017): 877-897.

Generally, a Matérn Type III process is simulated by simulating a primary point pattern, e.g., a homogenous or inhomogenous point process, simulating a uniform(0,1) age for each point,
sorting the points by age, and applying a thinning process to obtain the secondary points which are observed.

I use the result from Rao et al. showing the thinned points can be deleted and resimulated on each iteration instead of using a birth-death process. 
This approach is much faster and mixes much better. 

I consider three versions of this model.

Version 1: Hardcore point process. Points within a radius r of another, younger point are thinned with probability 1. We then estimate r.

Version 2: Softcore point process where points within a radius r of another, younger point are thinned with probability p.thin between 0 and 1. 
We estimate both r and p.thin. This version often leads to a multimodal posterior and can take a long time for r to converge, at least if p.thin is low.
Specifically, r can get stuck at low values and may take a long time to eventually converge upwards. Perhaps a better initialization method may help, or smarter
priors for r. Should run multiple chains to see if they converge to the same place.

Version 3: Softcore point process where younger points are thinned by a halfnormal thinning kernel with scale sigma.thin. I also include the thinning
radius r, but assume r is in the tail of the thinning kernel and is not estimated. It is used to speed up likelihood calculations. You can estimate it
in principle, but there can be major multimodality issues.

This seems like a good base model to use for spatial capture recapture. But I expect you will need extraordinary data sets with many traps and higher densities
so that many activity centers are localized accurately.