\name{dis.chao}
\alias{dis.chao}

\title{Chao et al. dissimilarity indices}

\description{Dissimilarity indices proposed by Chao et al. (2005) that takes into account (rare) unseen shared species. This is done mostly by heavy-weighting shared rare species. The function handle abundance or frequency data and are available for the Jaccard and Sorensen family of indices. Probability versions of the index, that does not heavy-weight rare species, are also included.
}

\usage{
  dis.chao(comm, index="jaccard", version="rare", freq=NULL)
 }

\arguments{
  \item{comm}{Dataframe or matrix with samples in rows and species in columns.}
  \item{index}{The index formula to be used in Chao dissimilarity. Partial match to "jaccard" or "sorensen".}
  \item{version}{The uncorrected probability version of the Chao index or the corrected version that takes into account unseen rare species. Partial match to "probability" or "rare".}
  \item{freq}{A numeric vector indicating total number of sampling units composing each sample (row) of the community data for computing incidence-based Chao indices. Length of \code{freq} must be the same as the number of samples in \code{comm}.}
}

\details{
Communities usually are composed of many rare and only a few common/frequent species. Although rare, a species may provide a valuable information in the estimation of resemblance between samples. However, the use of raw abundances makes the contribution of rare species to be negligible to the resulting index value. For instance, dropping a common species from the community data table usually will make much larger differences in the dissimilarity matrix than dropping a rare species.

There are a few indices that give more importance to an individual belonging to a rare than to a common/frequent species (see \code{\link{dis.goodall}} and \code{\link{dis.nness}}). Notice, however, that this differential weight of species can also be obtained, for instance, by log-transforming or standardizing data (e.g. dividing by maximum within each species) (Melo, in preparation). In this sense, an extreme case in which rare and common species have the same weight is in the use of presence/absence data.

The "probability" version of the abundance-based Chao-Jaccard and Chao-Sorensen indices are distinct from those traditional formulae available, for instance, in \code{\link[vegan]{vegdist}} (Chao et al. 2005). These probability indices do not overweight rare species and are included here for comparability. The version that takes into account unseeen rare species simplifies to the probability version if no rare species (singleton or doubleton for abundance data; unique or duplicate for incidence data) is present (see example below).  
}

\value{
A "dist" object.
} 


\references{
Chao, A., R.L. Chazdon, R.K. Colwell & T. Shen. 2005. A new statistical approach for assessing similarity of species composition with incidence and abundance data. Ecology Letters 8: 148-159.

Melo, A.S. in preparation. Is it possible to improve recovery of multivariate groups by weighting rare species in similarity indices?
}


\author{Adriano Sanches Melo}

\seealso{\code{\link[vegan]{vegdist}}, \code{\link{dis.goodall}}, \code{\link{dis.nness}}.}

\examples{
library(vegan)
### Example 1: Rare species are heavier:
aa <- c(1, 2, 4, 5)
bb <- c(1, 2, 0, 5)
cc <- c(0, 2, 3, 3)
dat3 <- rbind(aa, bb, cc) 
colnames(dat3) <- c("sp1", "sp2", "sp3", "sp4")
dat3

vegdist(dat3, method='jaccard', binary=FALSE) 
# Notice dissimilarity between the pair aa-bb is the same of the pair aa-cc. 
# In fact, bb and cc differ from aa in the same way (one species, and 
# 4 exclusive individuals).

dis.chao(dat3, index="jaccard", version="prob") 
# The probability version of the Chao index, however, produce different 
# dissimilarities for the pairs aa-bb and aa-cc 
# (aa-cc is less dissimilar than aa-bb).

dis.chao(dat3, index="jaccard", version="rare")
# The dissimilarity for the pair aa-cc is the same as that obtained using the 
# probability version. However, the dissimilarity for the pair aa-bb decreased.
# The reason is that aa-bb shares two rare species (sp1, sp2), 
# whereas the pair aa-cc shares a single rare species (sp2). 


### Example 2: "rare" version of the Chao index simplifies to the 
# "probability" version if no rare species (with 1 or 2 individuals) is present.
data(japi)
dim(japi) # 75 sampling units (stones) and 66 morphospecies of 
          # stream macroinvertebrates.
japi.m <- as.matrix(japi)
japi.2 <- ifelse(japi.m==1, 3, japi.m) # no singletons.
japi.3 <- ifelse(japi.2==2, 3, japi.2) # no doubletons.

sort(
     dis.chao(japi.3, index='jac', version='rare')  
   - dis.chao(japi.3, index='jac', version='prob'))


### Example 3: frequency data
# Stones in the japi dataset were sampled from downstream to upstream direction.
# Consecutive stones are spaced 1-6 m. The set of the first 25 stones should be
# more dissimilar to the last set of 25 stones than the middle set 
# (simply because of spatial autocorrelation).
japi.pa <- ifelse(japi.m > 0, 1, 0) # presence-absence data
japi.1st <- japi.pa[ 1:25, ]
japi.2nd <- japi.pa[26:50, ]
japi.3rd <- japi.pa[51:75, ]

japi.inc <- rbind(
                 colSums(japi.1st),
                 colSums(japi.2nd),
                 colSums(japi.3rd)
) # species frequency of occurrence in the three sets of stones.

dis.chao(japi.inc, index="jaccard", freq=c(25, 25, 25) )
}