\name{tree6}
\alias{tree6}
\docType{data}
\title{Example of phylogenetic tree containing 6 species.}
\description{Example of phylogenetic tree containing 6 species}
\usage{data(tree6)}
\format{
  The format is:
List of 4
 $ edge       : int [1:10, 1:2] 7 8 9 9 8 10 10 7 11 11 ...
 $ Nnode      : int 5
 $ tip.label  : chr [1:6] "taxon_1" "taxon_2" "taxon_4" "taxon_3" ...
 $ edge.length: num [1:10] 1 1 1 1 1 1 1 1 2 2
 - attr(*, "class")= chr "phylo"
}

\examples{
data(tree6)
plot(tree6)
}
\keyword{datasets}
