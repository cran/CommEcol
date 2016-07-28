\name{dis.goodall}
\alias{dis.goodall}
\title{Goodall dissimilarity index}

\description{This is a probability index in which rare species (or descriptors) are overweighted. Accordingly, sharing rare species makes pairs of samples (or objects) more similar than sharing common/frequent species (or descriptors). As Legendre & Legendre (1998) put it when describing this index: "... it is less likely for two sites to both contain the same rare species than a more frequent species. In this sense, agreement for a rare species should be given more importance than for a frequent species, when estimating the similarity between sites.".
}

\usage{
  dis.goodall(comm, p.simi="steinhaus", approach="proportion")
}

\arguments{
  \item{comm}{Dataframe or matrix with samples in rows and species in columns.}
  \item{p.simi}{The partial similarity index to be used in the Goodall index. Partial match to "steinhaus" (raw abundance data), or "gower" (normalized abundance data).}
  \item{approach}{The two approachs to compute Goodall index. Partial match to "proportion" or "chisquare".}
  }

\details{
Communities usually are composed of many rare and only a few common/frequent species. Although rare, a species may provide a valuable information in the estimation of resemblance between samples. However, the use of raw abundances makes the contribution of rare species to be negligible to the resulting index value. For instance, dropping a common species from the community data table usually will make much larger differences in the dissimilarity matrix than dropping a rare species.

There are a few indices that give more importance to an individual belonging to a rare than to a common/frequent species (see \code{\link{dis.chao}} and \code{\link{dis.nness}}). Notice, however, that this differential weight of species can also be obtained, for instance, by log-transforming or standardizing data (e.g. dividing by maximum within each species) (Melo in preparation). In this sense, an extreme case in which rare and common species have the same weight is in the use of presence/absence data.

The implementation of Goodall index (Goodall 1966) follows Legendre & Legendre (1998, pp. 269-273). They suggest to use the "steinhaus" index for raw species abundance data and the "gower" partial similarity index (\code{p.simi}) for normalized data. The index can be calculated in two ways or approaches ("proportion" or "chisquare") and values produced by them should be very different, although this is mostly due to scale; they are monotonically correlated. Notice Legendre & Legendre (1998) present a similarity version of this index. The one produced by this function is a dissimilarity, computed simply as 1-similarity. See the examples below for some behaviors of the index.
}

\value{
A "dist" object.
} 


\references{
Goodall, D.W. 1966. A new similarity index based on probability. Biometrics 22: 882-907.

Legendre, P & L. Legendre. 1998. Numerical Ecology. 2nd ed. Elsevier. 

Melo, A.S. (in preparation) Is it possible to improve recovery of multivariate groups by weighting rare species in similarity indices? 
}


\author{Adriano Sanches Melo}

\seealso{\code{\link[vegan]{vegdist}}, \code{\link{dis.chao}}, \code{\link{dis.nness}}.}

\examples{
library(vegan)
a <- c(1, 1, 0)
b <- c(2, 1, 0)
c <- c(0, 1, 1)
d <- c(0, 1, 2)
e <- c(0, 1, 3)
dat5 <- rbind(a,b,c,d,e)
colnames(dat5) <- c("sp1","sp2","sp3")
dat5

# Notice the samples in the pair a-b differ from each other exactly in the same
# way as samples in the pair c-d. However, a-b shares a rare species (sp1), 
# whereas c-d shares a frequent species (sp3, which is also present in e). Thus,
# the dissimilarity a-b is the same of c-d using Bray-Curtis, but not using 
# Goodall index:
vegdist(dat5, "bray")
dis.goodall(dat5)

# As the importance of a species for the Goodall index depends on its overall
# frequency in the community data, the deletion of a sample changes results:
dis.goodall(dat5[-5,])
}