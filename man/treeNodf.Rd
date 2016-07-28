\name{treeNodf}
\alias{treeNodf}
\alias{treeNodfTest}
\title{Tree-like Nestedness of Ecological Metacommunities}

\description{This function implements a generalization of the NODF (Nested Overlap and Decreasing Fill; Almeida-Neto et al. 2008) to quantify nestedness in metacommunities that takes into account relatedness among objects expressed as a tree-like object (Melo, Cianciaruso and Almeida-Neto, submitted).}

\usage{
  treeNodf    (comm, col.tree, order.rows=FALSE, row.tree, order.cols=FALSE)
  treeNodfTest(comm, col.tree, order.rows=FALSE, row.tree, order.cols=FALSE,
                                      null.model="perm.rows", permutations=999)
}

\arguments{
  \item{comm}{Dataframe or matrix with samples in rows and species in columns.}
  \item{col.tree}{A tree-like object containing all species listed in \code{comm}. This tree will be used to quantify (and test) nestedness among objects in rows (e.g. sites).}
  \item{order.rows}{Should rows of \code{comm} be ordered by decreasing Branch-Length (BL) Diversity. See details below.}
  \item{row.tree}{A tree-like object containing all row objects (e.g. sites) listed in \code{comm}. This tree will be used to quantify (and test) nestedness among objects in columns (e.g. species).}
  \item{order.cols}{Should columns of \code{comm} be ordered by decreasing Branch-Length (BL) Diversity. See details below.}
  \item{null.model}{Seven null models are currently implemented: "perm.rows", "perm.cols", "perm.rc", "perm.tip.cols", "perm.tip.rows", "perm.tip.rc" and "ff". See details.}
  \item{permutations}{Number of permutations of the \code{null.model} to assess significance.}
}

\details{This is a direct extension of the NODF metric used to quantify nestedness in metacommunities using presence-absence data. NODF measures the proportion of the species richness present in a species-poor community which is present in a species-rich community. The tree-like version of this metric uses Branch-Length (BL) as a measure of diversity (instead of species richness used in NODF). For phylogenetic trees (a cladogram), BL is the Phylogenetic Diversity (PD, sensu Faith 1992). In this case, treeNODF (or phyloNODF) measures the proportion of the BL (or PD) in a BL-poor (or PD-poor) community that is present in a BL-rich (or PD-rich) community. The same reasoning applies to other tree-like objects such as functional dendrograms (Petchey & Gaston 2006).

The treeNODF follows the same approach of NODF and are calculated for each pair of objects in a given matrix dimension. For instance, for a dataset including 5 sites in its rows, 10 (5*4/2) pairs of values are calculated. The treeNODF for rows is simply the average of these 10 values. The same can be applied for columns. For instance, if columns represent species, one may assess whether there are nestedness in the environmental conditions where they occur (envNODF). For envNODF, users should supply a tree-like object that depicts resemblance among sites (i.e. a dendrogram where sites are classified according to environmental variables of interest). Finally, the treeNODF can be calculated for both rows and columns. In this case, two tree-like objects must be provided. The option to test rows (sites), columns (species) or both will depend exclusively on the hypothesis the user have raised regarding the nestedness structure of his study system.

The treeNODF can be partitionated into two components. The first one, called S.fraction, is the proportion of the species richness (or species incidence in the case of columns) of the BL-poor community that is shared with the BL-rich community. The second one is obtained as topoNODF = treeNODF - S.fraction and measures the effect of tree topology to treeNODF (see Melo et al. submitted). 

NODF and its implementation in R \code{\link[vegan]{nestednodf}} allows one to order rows and columns by frequencies (species richness for sites and species incidences for species). However, evidence of nestedness in a matrix automatically ordered by frequencies does not allows a proper inference of the mechanism generating nestedness (see Almeida-Neto et al 2008, Ulrich et al. 2009). Accordingly, better inferences of mechanisms generating nestedness and tree-nestedness should be done by ordering communities (rows) or species (columns) accordingly to a hypothesis raised \emph{a priori} (and not by their S or BLs).

The "perm.rows" null model permute rows (sites) and is useful to test a site by species matrix in which rows (sites) were ordered by an \emph{a priori} hypothesis (e.g. island area, Lomolino 1996). The "perm.cols" null model permute columns (species) and is useful to test a site by species matrix in which columns (species) were ordered by an \emph{a priori} hypothesis (e.g. species body size). These two models must not be used if diversity (S or BL) is used to order communtiy data. This is because simulated values of the statistic will lower or equal to the observed one, but never higher. Although termed 'null models' here, notice these procedures constitutes a permutation test (rows or columns are permuted, not elements).

The "perm.tip.cols" null model shuffles species (columns) labels across tips of the tree-like object. This null model should be useful to test the effect of resemblance among species (columns) to the treeNODF result. Similarly, the "perm.tip.rows" model shuffles tip labels of the row.tree (usually, sites). The "perm.tip.rc" model shuffles tip labels of both trees.

The "ff" null model is a popular choice in the literature on species co-occurrence analysis. It maintains sum of rows (species richness of sites) and sum of columns (species frequencies) fixed. Notice, however, that treeNODF makes inferences about BL and this type of "diversity" is affected not only by species richness but also by the resemblance among species (or sites) present in communities. The function \code{\link[vegan]{commsimulator}} is used to shuffle matrices according to the "ff" null model using the "quasiswap" option.

}

\value{
   For the \code{treeNodf}, a list including:
  \item{rows }{The average of \code{tree.nodf.rows.dist} and the averages of its two components for rows.}
  \item{cols }{The average of \code{tree.nodf.cols.dist} and the averages of its two components for columns.} 
  \item{mat }{The average of \code{tree.nodf.rows.dist} and \code{tree.nodf.cols.dist} when treeNODF for both rows and columns are calculated. This is the treeNODF for the entire dataset. Averages for the two component metrics are also provided.}
  \item{tree.nodf.rows.dist }{A triangular matrix of treeNODFs for each pair of rows (sites). Each value represents the treeNODF of the site showed in its row in relation to the site showed in its column. Thus, values in the first column represents the treeNODFs of each community in relation to the community supposed to have the highest BL.}
  \item{tree.nodf.cols.dist }{A triangular matrix of treeNODFs for each pair of columns (species). Each value represents the treeNODF of the species showed in its column in relation to the site showed in its row. Thus, values in the first column represents the treeNODFs of each species in relation to the species supposed to have the highest BL (that is, the species supposedly to occurr in the widest set of conditions).}
  \item{s.fraction.rows.dist }{A triangular matrix of s.fration values for each pair of rows (sites).}
  \item{s.fraction.cols.dist }{A triangular matrix of s.fraction values for each pair of columns (species).}
  \item{topo.nodf.rows.dist }{A triangular matrix of topo.nodf values for each pair of rows (sites).}
  \item{topo.nodf.cols.dist }{A triangular matrix of topo.nodf values for each pair of columns (species).}  


  For the \code{treeNodfTest}, a list containing the same list provided by \code{treeNODF} for the observed values and:
  \item{rows.aleats }{treeNODF, S.fraction and topoNODF values for rows obtained using the chosen null model.}
  \item{cols.aleats }{treeNODF, S.fraction and topoNODF values for colums obtained using the chosen null model.}
  \item{mat.aleats }{treeNODF, S.fraction and topoNODF values for the entire matrix obtained using the chosen null model.}
  \item{permutations }{Number of permutations used in the test.}
  } 


\references{
Almeida-Neto M., P. Guimaraes, P.R. Guimaraes, R.D. Loyola & W. Ulrich. 2008. A consistent metric for nestedness analysis in ecological systems: reconciling concept and measurement. Oikos 117:1227-1239.

Faith, D.P. 1992. Conservation evaluation and phylogenetic diversity. Biological Conservation 61:1-10.

Lomolino, M.V. 1996. Investigatin causality of nestedness of insular communities: selective immigration or extinctions? Journal of Biogeography 23:699-703.  

Melo, A.S., M.V. Cianciaruso & M. Almeida-Neto. 2014. treeNODF: nestedness to phylogenetic, functional and other tree-based diversity metrics. Methods in Ecology and Evolution 5: 563-572.

Petchey, O.L. & K.J. Gaston. 2006. Functional diversity: back to basics and looking forward. Ecology Letters 9:741-758.

Ulrich, W., M. Almeida-Neto & N.J. Gotelli. 2009. A comsumer's guide to nestedness analysis. Oikos 118:3-17.}

\author{Adriano Sanches Melo}

\seealso{\code{\link[picante]{pd}} (package picante), \code{\link[betapart]{beta.multi}} (package betapart)}

\examples{
library(picante)
data(sites5.6)
data(tree6)
treeNodf    (sites5.6, col.tree=tree6)
# You likely will need more runs in your permutation test.
treeNodfTest(sites5.6, col.tree=tree6, null.model="perm.rows", permutations=9)

alt5 <- vegdist(1:5, method="euclidean")
alt5 <- hclust(alt5)
alt5 <- as.phylo(alt5)
alt5$tip.label<-rownames(sites5.6)
treeNodf    (sites5.6, row.tree=alt5)
#treeNodfTest(sites5.6, row.tree=alt5, null.model="perm.cols",permutations=99)

treeNodf    (sites5.6, col.tree=tree6, row.tree=alt5)
#treeNodfTest(sites5.6, col.tree=tree6, row.tree=alt5, null.model="perm.rc", permutations=99)
#treeNodfTest(sites5.6, col.tree=tree6, row.tree=alt5, null.model="ff", permutations=99)

# In the example below treeNodf is zero because PD of the first sample is 
# lower than that of the second sample. Notice samples have the same species 
# richness and higher PD for the first sample is due to the presence of a
# distinctive species (e).
tree <- read.tree(text="(((a:1,b:1):1,(c:1,d:1):1):2, e:4);")
plot(tree)
mat  <- matrix(c(1,1,1,1,0, 1,1,1,0,1), byrow=TRUE, nrow=2, dimnames=list(1:2, letters[1:5]))
mat
pd(mat, tree)
treeNodf(mat, tree) 

# Here treeNodf is non-zero... but you would only do that if you have
# an 'a priori' hypothesis to order rows of your matrix this way!
treeNodf(mat[2:1, ], tree) 
}

