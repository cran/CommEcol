\name{autosimi}
\alias{autosimi}
\title{Autosimilarity curve of species community data}

\description{Similarities among two subsamples, each one obtained randomly from the same community dataset. Curves are otained for all subsample sizes from 1 up to half the number of sample units in the dataset. Autosimilarity curves can be used to evaluate sample suficiency when sample size is expressed as number of sampling units such as traps or quadrats. This is particularly suitable when the study involves similarities or dissimilarities among samples such as agglomerative clustering, ordination, Mantel test.}

\usage{
autosimi(comm, method="bray", binary=FALSE, log.transf=FALSE, simi=TRUE, permutations=50)}

\arguments{
  \item{comm }{Dataframe or matrix with samples in rows and species in columns.}
  \item{method }{Similarity index obtained from \code{\link[vegan]{vegdist}}. Similarities are obtained simply as 1-dissimilarity. Accordingly, it only makes sense for indices bounded at 0-1. For indices not bounded at 0-1 (distances), use \code{simi=FALSE} to obtain autodissimilarity curves.}
  \item{binary }{Should data be transformed to presence/absence?}
  \item{log.transf }{Transformation \code{log(x+1)} before calculation of similarities (or dissimilarities)}
  \item{simi }{Similarity or dissimilarity curve.}
  \item{permutations }{Number of curves randomly calculated and from which the mean autosimilairty curve is obtained.}
}

\details{The function selects randomly and without replacement an even number of sampling units (n = 2,4,6, ...) from the total sampling units. Next, the first n/2 sampling units are pooled (summed) to create a subsample, and the other n/2 sampling units to create another subsample. If \code{binary=TRUE}, the two subsamples are transformed to presence/absence. If \code{log.transf=TRUE}, \code{log(x+1)} of each value is obtained. The similarity between the two subsamples is calculated and stored. The procedure is repeated for larger subsample sizes until half the size of the full dataset (or up to the integer quotient in the case of odd numbers of sample units). The procedure is then repeated for the requested number of curves. The output is the average curve. This function is a different implementation of the procedures described in Schneck & Melo (2010).
}

\value{A dataframe containing subsample sizes and average similarity values. 
}

\references{
Cao, Y., D.P. Larsen & R.M. Hughes. 2001. Evaluating sampling sufficiency in fish assemblage surveys: a
similarity-based approach. Canadian Jounal of Fisheries and Aquatic Sciences 58: 1782-1793.

Schneck, F. & A.S. Melo. 2010. Reliable sample sizes for estimating similarity among macroinvertebrate assemblages in tropical streams. Annales de Limnologie 46: 93-100.

Weinberg, S. 1978. Minimal area problem in invertebrate communities of Mediterranean rocky substrata. Marine Biology 49: 33-40.
}

\author{Adriano Sanches Melo}

\seealso{\code{\link[vegan]{vegdist}}}

\examples{
x<-matrix(0,4,4)
diag(x)<-1
x4<-rbind(x,x,x,x)
x4
autosimi(x4, binary=TRUE)
plot(autosimi(x4, binary=TRUE))

data(BCI)
simi<-autosimi(BCI, binary=TRUE, permutations=5)
simi
plot(simi, ylim=c(0.5,1)) # maintain the plot window open for the next curve
simi.log<-autosimi(BCI, binary=FALSE, log.transf=TRUE, permutations=5)
points(simi.log, col="red")
}
