\name{part.m.tree}
\alias{part.m.tree}
\title{Partition of multi-sample dissimilarity indices using phylogenetic/functional data}

\description{This function computes the partition of multi-sample Sorensen and Jaccard dissimilarity indices described by Baselga (2012), but adapted for phylogenetic/functional data. The two resulting components are dissimilarities due to turnover and nestednes.}

\usage{part.m.tree(comm, tree, index.family="sorensen")}

\arguments{
  \item{comm}{Dataframe or matrix with samples in rows and species in columns.}
  \item{tree}{A phylogenetic/functional tree containing all species listed in \code{comm}}.
  \item{index.family }{The family of dissimilarity indices to be partitionated. Partial match to "sorensen" or "jaccard".}
}

\value{A list including the number of samples and three dissimilarity values. The nomenclature of the dissimilarities follows Baselga (2012). 
   For the jaccard:
  \item{JAC }{Multi-sample Jaccard dissimilarity index.}
  \item{JTU }{The turnover component of the Jaccard index.}
  \item{JNE }{The nestedness component of the Jaccard index.} 

  For the sorensen:
  \item{SOR }{Multi-sample Sorensen dissimilarity index.}
  \item{SIM }{The turnover component of the Sorensen index (Simpson index).}
  \item{SNE }{The nestedness component of the Sorensen index.}
}

\references{
Baselga, A. 2010. Partitioning the turnover and nestdness components of beta diversity. Global Ecology and Biogeography 19, 134-143.

Baselga, A. 2012. The relationship between species replacement, dissimilarity derived from nestedness and nestedness. Global Ecology and Biogeography 21, 1223-1232.

Leprieur, F., C. Albouy, J.D. Bortoli, P.F. Cowman, D.R. Bellwood and D. Mouillot. 2012. Quantifying phylogenetic beta diversity: distinguishing between 'true' turnover of lineages and phylogenetic diversity gradients. PLoS ONE 7(8), e42760. doi:10.1371/journal.pone.0042760
}

\author{Adriano Sanches Melo, but part of the code borrowed from \code{\link[betapart]{beta.multi}} (package betapart) and \code{\link[picante]{phylosor}} (package picante)}

\seealso{\code{\link{part.p.tree}} for pair-wise partitions.}

\examples{
library(picante)
data(sites4.6)
data(tree6)
part.m.tree(comm=sites4.6, tree=tree6, index.family="sorensen")
}
