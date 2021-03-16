\name{compas}
\alias{compas}
\title{Simulation of species community data along gradients}

\description{The function generates unimodal curves of different shapes according to the specification of different parameter values. For different combinations of parameters, the response function may resemble the Gaussian curve or its skewed and platykurtic variations. Simulated observations are obtained by generating values for parameters and solving the function for defined x-values. In an ecological context, the response curve may be interpreted as densities of a species along an environmental gradient. Solving the function for a given x-value would be equivalent to sampling the species at the coordinate x. If many different curves are generated, solving the function for a given x-value would be equivalent to sampling a community at the coordinate x of the gradient. The idea is easily expanded to include two or more gradients or dimensions. For the case of two gradients, two sets of coordinates, one for each dimension, are used to obtain a community observation. A theory justifying the use of response curves along gradients and a discussion of the shapes of these curves is found in McGill & Collins (2003). The code is mostly based on Minchin (1987a, 1987b).}

\usage{
compas(S, dims, am, clump=1, beta.R, coords, n.quanti=5, n.quali=0.1, add1)}

\arguments{
  \item{S }{The number of species occurring in the simulated gradients. This IS NOT necessarily the number of species that will appear in the resulting dataset as some species may not be "sampled".}
  \item{dims }{Number of gradients (dimensions).} 
  \item{am }{A vector of abundance of species in its modal point (log scale) for each gradient. This(ese) value(s) is(are) used to sample a lognormal distribution with mean(s) \code{am} and sd=1. The number of supplied means should be the same of the number of \code{dims}.}
  \item{clump }{A non-zero positive integer indicating how strong the species richness gradient is. For \code{clump}=1, modal points for species (see Details) are randomly and uniformly distributed in the space and thus species richness are mostly homogeneous. For higher values, more modal points of species (and thus more species) will be clumped on higher values of the gradient(s) (the upper right-hand corner if \code{dims}=2. See example.}
  \item{beta.R }{Beta diversity (turnover) parameter.}
  \item{coords }{A matrix-like (or vector if \code{dims}=1) object with coordinates of sampling sites. The number of columns (axes) should be the same of the number of \code{dims}.}
  \item{n.quanti }{A value higher than 0 indicating how much noise or variability should be added to abudance values sampled from the response curves. Each abundance value is substituted by a random value obtained from a Negative Binomial distribution with mean equal to the respective abundance value and variance proportional to \code{n.quanti}.} 
  \item{n.quali }{Qualitative Noise. Each specie has probability "1-n.quali" of occurring in a site within its range. The argument control the replacement of \code{n.quali}*100\% of the abundance values by "0".}
  \item{add1 }{A value between 0 and 1. Add (\code{add1}*100)\% of "marginal/vagrant species" occurring randomly with 1 individual in the entire dataset.}
}

\details{
This implementation is based on the software Compas (Minchin 1987b), described in detail in Minchin (1987a). Simulated parameters are random values obtained from distributions such as the normal and the uniform. Some of these parameters can be modified by users. However, some them are fixed and not modified unless you are able to edit the code. The option to fix some of the parameters should simplify the use of the function.

The number of species in the simulated matrix may be smaller than \code{S} because some species may occur outside the gradient. The gradiente is -50 up to 150, but 'sampling' occurs only at the range 0-100. This allows species to have their mode outside the 0-100 gradient. Also, species occurring in the gradient may not be sampled as a result of the quantitative (\code{n.quanti}) and qualitative (\code{n.quali}) noises added.

Parameters alpha and gamma (not available as arguments), which together determine curve symmetry and kurtosis, are obtained from a uniform distribution bounded by 0.1 and 5. 

The abundance at the mode of the curve is determined by random values otained from log-normal distribution with log(mean) \code{am} and  sd=1 (the last one not available as argument). 

The position of the modal point in the gradient (parameter m, not available as argument) will depend on the value of \code{clump}. If \code{clump}=1 coordinates of the modal points are obtained from a random uniform distribution. For higher values of \code{clump}, species modal points will tend to be concentrated on high values of the gradient. In all cases, coordinates are bounded by -50 and +150. For a studied gradient bounded by 0 and 100, the specification of a larger interval allows the position of the modal point to be located outside the gradient. 

The parameter \code{beta.R} determines the range of occurrence of species in the gradient. In terms of diversity, it determines the turnover (beta diversity) along the gradient, and is expressed as \code{beta.R} = 100/mean(r), and r are the ranges of species. Values of r are obtained from a normal distribution with mean 100/\code{beta.R} and standard deviation of 0.3*100/\code{beta.R} (the last one not available as argument). A \code{beta.R} = 2 determines that range of species are, on average, 50 units of the gradient (and 15 SD units). High \code{beta.R} determines restricted ranges and, thus, high beta diversity. Old versions of this function included the argument "beta" and expressed, contrary to the name of the argument, as range in the gradient. The code of the new function is equivalent to the old ones (beta=0.5 was translated as 0.5*100 and determined ranges around 50; this can be attained now as \code{beta.R}=2 which determines ranges 100/2=50). 

Species may be (i) present in densities different from the expected or (ii) absent due to sampling error, historical factors, or the influence of a multitude of environmental factors. Simulated communities are subject to these two types of error or noise to mimic the two phenomena. In order to obtain scattered values of abundances for each species along their gradients, each non-zero value is replaced by a random value obtained from a Negative Binomial distribution (\code{\link[stats]{rnbinom}}) with mean equal to the value to be replaced and variance proportional to \code{n.quanti}. The parametrization of the Negative Binomial is done using mu (mean) and size (dispersion parameter), the later obtained as mu^2 / (variance-mu). The \code{n.quanti} is used to define the variance as mu*(1+\code{n.quanti}). For \code{n.quanti} tending to zero the distribution is similar to the Poisson. For values up to 1 it is still similar to the Poisson. For \code{n.quanti} >5  it is quite distinct with long right-hand tails and many zeros if the mean is low (e.g. <10). Additional absences are achieved by randomly choosing \code{n.quali}*100\% of non-zero values and replacing them by zeros. 

A third characteristic commonly seen in field datasets is the occurrence of a species outside its regular range or in different habitats. These species are termed vagrant or marginal, and usually appear in the dataset with 1 individual. In order to account for this common finding, the function attaches to the simulated community an additional set of species with 1 individual, each species occurring in one sample unit only. The number of these additional species is set as \code{add1}*100\% of the number of species sampled in the simulated communities after the inclusion of the two types of errors cited above.
}

\value{A dataframe of sites (rows) by species (cols) abundances. In case of \code{dims}=1, a plot of response species curves is produced. These curves do not include quantitative noise (\code{n.quanti}), qualitative noise (\code{n.quali}) (see comments above) and 'marginalor vagrant species' (\code{add1}) additions. Notice that the number of species in the result dataframe may vary as (i) some species may not be sampled by \code{coords}, or (ii) species range is located outside the sampling gradient (0-100).
}

\references{
McGill, B. & C. Collins. 2003. A unified theory for macroecology based on spatial patterns of abundance. Evolutionary Ecology Research 5: 469-492.

Minchin, P.R. 1987a. Simulation of multidimensional community patterns: towards a comprehensive model. Vegetatio 71: 145-156.

Minchin, P.R. 1987b. COMPAS: a program for the simulation of multidimensional community patterns based on generalized beta functions. Dep. of Biological Sciences, Southern Illinois University Edwardsville, USA.
}

\author{Adriano Sanches Melo}

\seealso{\code{\link[vegan]{vegdist}}}

\note{As the function includes many parameters, many simulated communities will not mimic the real ones. Users should try different set of options to approxiamte real datasets. A starting point may be the examples provided below.
}

\examples{
# 1 dimension.
coo <- seq(10, 90, 10)
compas(S=30, dims=1, am=2, beta.R=2, coords=coo, 
       n.quanti=5,  n.quali=0.1, add1=0.1)
	
# 2 dimensions.
coo2 <- cbind(coo, coo)
compas(S=50, dims=2, am=c(2,2), beta.R=c(1,1), coords=coo2, 
       n.quanti=5, n.quali=0.1, add1=0.1)

# 2 dimensions. Homogeneous and clumped species distributions.
# Try a few times.
coo.grid <- expand.grid(seq(10,90,10), seq(10,90,10))
plot(coo.grid, xlim=c(0,100), ylim=c(0,100))

mat1 <- compas(S=60, dims=2, am=c(2,2), clump=1, beta.R=c(2,2), coords=coo.grid,
               n.quanti=5, n.quali=0.1, add1=0.1)
S1 <- rowSums(ifelse(mat1>0,1,0))

mat2 <- compas(S=60, dims=2, am=c(2,2), clump=3, 
               beta.R=c(2,2), coords=coo.grid, n.quanti=5, n.quali=0.1, add1=0.1)
S2 <- rowSums(ifelse(mat2>0,1,0))

ind <- coo.grid/10
S1.mat <- matrix(NA, 9, 9)
S2.mat <- matrix(NA, 9, 9)
for(i in 1:length(S1)){
      S1.mat[ind[i,1], ind[i,2]] <- S1[i] 
      S2.mat[ind[i,1], ind[i,2]] <- S2[i] 
   }

plot(coo.grid, cex=S1/4)
plot(coo.grid, cex=S2/4)

image(x=coo, y=coo, z=S1.mat, col=gray(50:1/50))
image(x=coo, y=coo, z=S2.mat, col=gray(50:1/50))

filled.contour(x=coo, y=coo, z=S1.mat)
filled.contour(x=coo, y=coo, z=S2.mat)
}
