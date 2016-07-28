\name{stairs6}
\alias{stairs6}
\docType{data}
\title{Example of phylogenetic tree containing 6 species.}
\description{Example of phylogenetic tree containing 6 species}
\usage{data(stairs6)}
\format{
  The format is:
List of 5
 $ edge       : int [1:10, 1:2] 7 8 9 10 11 11 10 9 8 7 ...
 $ Nnode      : int 5
 $ tip.label  : chr [1:6] "taxon_1" "taxon_2" "taxon_3" "taxon_4" ...
 $ edge.length: num [1:10] 1 1 1 1 1 1 2 3 4 5
 $ root.edge  : num 1
 - attr(*, "class")= chr "phylo"
 - attr(*, "order")= chr "cladewise"
}

\examples{
data(stairs6)
plot(stairs6)
}
\keyword{datasets}
