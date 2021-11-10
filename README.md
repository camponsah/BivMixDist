# BivMixDist
The work concerns models with stochastic representation (X, N) where N is a positive integer-valued random variable and X is the sum of N independent/dependent sequence {Xi}.  

In the above context, several models have been developed, which include the bivariate distribution with exponential and geometric marginals (BEG) by Kozubowski and Panorska (2005), the bivariate gamma-geometric law (BGG) by Barreto-Souza W. (2012), the bivariate distribution with Lomax and geometric margins (BLG) by Arendarczyk et al (2018), the bivariate distribution with gamma mixture and discrete Pareto margin (GMDP) Amponsah (2017 and   etc.

We implement various algorithms useful for parameter estimation of  BEG, BGG, GMDP, BLG, and among other. EM-algorithm is implemented for discrete Pareto distribution (see Buddana and Kozubowski, 2014) and testing for exponentiality versus Pareto distribution is also implemented (see Arendarczyk et al, 2018). Others include random generation, density, quantiles and distribution functions. The link to the various articles is listed below.

A mixed bivariate distribution with exponential and geometric marginals (BEG ): https://doi.org/10.1016/j.jspi.2004.04.010 

Bivariate gamma-geometric law and its induced Lévy process (BGG ): https://doi.org/10.1016/j.jmva.2012.03.004

Bivariate gamma mixture discrete Pareto distribution(GMDP): https://scholarworks.unr.edu/bitstream/handle/11714/2065/Amponsah_unr_0139M_12378.pdf?sequence=1&isAllowed=y 

A bivariate distribution with Lomax and geometric margins (BLG): https://doi.org/10.1016/j.jkss.2018.04.006

Testing Exponentiality versus Pareto distribution via likelihood ratio: https://doi.org/10.1080/03610910802439121 

Discrete Pareto distribution (DP):  https://doi.org/10.1515/eqc-2014-0014

A general stochastic model for bivariate episodes driven by a gamma sequence: https://doi.org/10.1186/s40488-021-00120-5

 
# Install from github using 'devtools'
install.packages("devtools")

devtools::install_github("camponsah/BivMixDist")

library(BivMixDist)
