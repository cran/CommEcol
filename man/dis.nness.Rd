\name{dis.nness}
\alias{dis.nness}
\alias{dis.nness.find.m}

\title{NNESS and NESS dissimilarity indices}

\description{These dissimilarities indices overweight rare species to a desired degree. If rare species are heavy-weighthed, individuals of a rare species are more important than individuals of a common species. In the extreme, a species represented by a single individual will have the same weight of a common species, which is equivalent to the use of presence-absence data.
}

\usage{
  dis.nness(comm, m=NULL, ness=FALSE)
  dis.nness.find.m(comm, ness=FALSE)
  
}

\arguments{
  \item{comm}{Dataframe or matrix with samples in rows and species in columns.}
  \item{m}{The number of individuals to be sampled in the estimation of species shared in NNESS and NESS indices. \code{m} must be an integer and lower than half (NESS) or the total (NNESS) abundance of the smallest sample in the pair under comparison. Higher \code{m} gives more importance to rare species. If no value is provided, the function will use \code{dis.nness.find.m} to find a value that is sensitive to both rare and abundant species (see Details below).}
  \item{ness}{Compute NESS index.}
}

\details{
Communities usually are composed of many rare and only a few common/frequent species. Although rare, a species may provide a valuable information in the estimation of resemblance between samples. However, the use of raw abundances makes the contribution of rare species to be negligible to the resulting index value. For instance, dropping a common species from the community data table usually will make much larger differences in the dissimilarity matrix than dropping a rare species.

There are a few indices that give more importance to an individual belonging to a rare than to a common/frequent species (see \code{\link{dis.chao}} and \code{\link{dis.goodall}}). Notice, however, that this differential weight of species can also be obtained, for instance, by log-transforming or standardizing data (e.g. dividing by maximum within each species) (Melo, in preparation). In this sense, an extreme case in which rare and common species have the same weight is in the use of presence/absence data.

NNESS is a modified or New version of the Normalized Expected Species Shared (NESS). NESS was proposed by Grassle & Smith (1976) and estimates similarity based on the number of species shared between random samples of size \code{m} individuals. Higher \code{m} gives more weight to rare species. For \code{m} = 1, NESS is the same as Morisita and NNESS is the same as Morisita-Horn dissimilarities. As higher values of \code{m} are used, dissimilarities tend to converge to presence-absence Sorensen index (when a rare and a common species have the same weight). NNESS was proposed by Trueblood et al. (1994) and circumvents the inability of NESS to handle samples composed exclusively by singletons (species with 1 individual).

Trueblood et al. (1994) also suggested to use an \code{m} value that is sensitive to both rare and abundant species. The NESS or NNESS index is calculated using values of \code{m} ranging from 1 up to half (NESS) or the total (NNESS) number of individuals in the smallest sample. The Kendall correlation is then calculated for each pair of triangular dissimilarity matrices. The selected \code{m} value is the one which produces a correlation with the matrix obtained with \code{m} = 1 that is most similar to the correlation with \code{m} = max, where max is the maximum value \code{m} may assume (abundance of the smallest sample for NNESS or half of that abundance for NESS). As this procedure may be slow for large datasets, up to 30 \code{m} values are used to compute dissimilarity matrices. NESS and NNESS are most suitable for raw abundance data and, thus, to data expressed as integers. This implementation, however, will work on non-integer data. NESS and NNESS formulae are for similarities and are computed here simply as 1-similarity. 

The calculation of NESS and NNESS involves the use of binomial coefficients. For samples with many individuals (e.g. 1100) and high \code{m} values (e.g. 490), the resulting value is too large to be stored as double precision. Accordingly, this function uses the \code{\link[gmp]{chooseZ}} of the Multiple Precision Arithmetic package when total abundance of at least one sample contains more than 1000 individuals. If not sample contains more than 1000 individuals, computations are performed using the simpler \code{\link[base]{choose}}.

NESS and NNESS are most suitable to abundance data. However, \code{\link[base]{choose}} is able to handle positive non-integers and calculations will be done accordingly. However, \code{\link[gmp]{chooseZ}} is not able to handle non-integers and, thus, for datasets in which at least one sample contains more than 1000 individuals, non-integer abundances will be rounded up to next upper integer.   
}

\value{
A "dist" object for \code{dis.nness} and a integer for \code{dis.nness.find.m}.
} 


\references{
Grassle, J.F. & W. Smith. 1976. A similarity measure sensitive to the contribution of rare species and its use in investigation of variation in marine benthic communities. Oecologia 25: 13-22.

Melo, A.S. Submitted. Is it possible to improve recovery of multivariate groups by weighting rare species in similarity indices?

Trueblood, D.D., E.D. Gallagher & D.M.Gould. Three stages of seasonal succession on the Savin Hill Cove mudflat, Boston Harbor. Limnology and Oceanography 39: 1440-1454.
}


\author{Adriano Sanches Melo}

\seealso{\code{\link[vegan]{vegdist}}, \code{\link{dis.chao}}, \code{\link{dis.goodall}}.}

\examples{
library(vegan)
aa <- c(1, 2, 4, 5)
bb <- c(1, 2, 0, 5)
cc <- c(0, 2, 3, 3)
dat3 <- rbind(aa, bb, cc) 
colnames(dat3) <- c("sp1", "sp2", "sp3", "sp4")
dat3

#  NESS using m=1 is the same as Morisita and
# NNESS using m=1 is the same as Morisita-Horn:
dis.nness(dat3, m=1, ness=TRUE)
vegdist(dat3, method="morisita")

dis.nness(dat3, m=1, ness=FALSE)
vegdist(dat3, method="horn")

# The dissimilarity for the pair aa-bb is reduced if more weight is given to 
# rare species (higher m). The reason is that aa-bb shares two rare 
# species (sp1, sp2), whereas the pair aa-cc shares a single rare species (sp2). 
dis.nness(dat3, m=1, ness=FALSE)
dis.nness(dat3, m=8, ness=FALSE)
}