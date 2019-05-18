\name{betaRegDisp}
\alias{betaRegDisp}
\title{Beta diversity metrics between sites in a moving window along environmental gradients}

\description{The function computes eight metrics of beta diversity according to an informed environmental gradient. It selects a given number of environmentally-neighborhood sites in a moving window to obtain beta diversity.}

\usage{
betaRegDisp(y, x, xy.coords = NULL, ws = 3, 
            method.1 = "jaccard", method.2 = "ruzicka", 
            method.3 = "ruzicka", 
            independent.data = FALSE, illust.plot = FALSE)}

\arguments{
  \item{y }{Response matrix, where rows are sites and columns are species.}
  \item{x }{Predictor vector. A vector of the environmental gradient under study with the same number of sites as in matrix \code{y}. Make sure the order of samples in \code{y} and \code{x} is the same.}
  \item{xy.coords }{Geographical coordinates. A matrix with two columns of XY decimal degree geographical coordinates, which are used to compute euclidean distance among sites. Rows must be sites in the same order as in \code{y} and \code{x}.}
  \item{ws }{Window size or number of sites to be used in the computation of the distinct beta-diversity metrics or between-site dissimilarities.  It must be a positive integer higher than 2.}
  \item{method.1 }{For beta-diversity metrics 1 to 3 (see details). A dissimilarity index available in the \code{\link[vegan]{vegdist}}. The options are: "euclidean", "jaccard", "bray", "manhattan", "canberra", "kulczynski", "gower", "morisita", "horn", "mountford",  "raup", "binomial", "chao", "altGower", "cao", "mahalanobis".}
  \item{method.2 }{For beta-diversity metrics 4 and 5 (see details). A dissimilarity index available in the \code{\link[adespatial]{beta.div}}. The options are: "hellinger", "chord", "log.chord", "chisquare", "profiles", "percentdiff", "ruzicka", "divergence", "canberra", "whittaker", "wishart", "kulczynski", "jaccard", "sorensen", "ochiai", "ab.jaccard", "ab.sorensen", "ab.ochiai", "ab.simpson", "euclidean".}
  \item{method.3 }{For beta-diversity metrics 6, 7, and 8 (see details). A multisample dissimilarity index available in the \code{\link[betapart]{beta.multi.abund}}. The options are: "ruzicka" and "bray".}
  \item{independent.data }{Should windows not superpose each other? If \code{independent.data=TRUE}, sites do not enter more than once in the dissimilarity calculations, that is, sites are included in a single window.}
  \item{illust.plot }{Should a window plot be open and illustrate how the window moves along the gradient?}
}

\details{The function computes eight beta-diversity metrics among sites included in a set (window) of length \code{ws}. See details in Dala-Corte et al. (2019).

Metrics 1-3 uses dissimilarity indices available in \code{\link[vegan]{vegdist}}:

1. Mean pair-wise dissimilarity between sites in a window;

2. Mean dissimilarity between focal site and the other sites in a window. If an odd number is informed in \code{ws}, the focal site is the central site in relation to its neighbours in the window. If an even number is informed in \code{ws}, the focal site is the first site in the window;

3. Mean distance of sites to their group centroid in a Principal Coordinate (PCoA) space computed using \code{\link[vegan]{betadisper}};

Metrics 4-5 uses dissimilarity indices available in \code{\link[adespatial]{beta.div}}:

4. Total sum of squares (SS) of the window sites (Legendre and De Caceres, 2013);

5. Local contributions to beta diversity (LCBD; Legendre and De Caceres, 2013);

Metrics 6-8 uses dissimilarity indices available in \code{\link[betapart]{beta.multi.abund}}:

6. Total multiple-site dissimilarities for a selected window of sites;

7. Nestedness component of multiple-site dissimilarities for a selected window of sites;

8. Turnover component of multiple-site dissimilarities for a selected window of sites. 
}

\value{A matrix with 10 columns (or 12 if \code{xy.coords} is informed). Values in columns are sorted according to the enviromental gradient, from the lowest to the highest value. Columns correspond to:

1. \code{grad} - The environmental gradient (predictor vector, \code{x});

2. \code{mean.grad} - Mean value of the environmental gradient of sites selected in each window;

3. \code{mean.diss.pairs} - Mean pair-wise dissimilarity between sites in a selected window (metric 1);

4. \code{diss.focal} - Mean dissimilarity between focal site and the other sites (metric 2);

5. \code{mean.dist.cent } - Mean distance of sites to their group centroid in a Principal Coordinate (PCoA) space (metric 3);

6. \code{SS.group} - Total sum of squares (SS) of the sites in a window (metric 4); 

7. \code{SS.focal} - Local contributions to beta diversity (LCBD), which represents how much a focal site contributed to the total window SS;

8. \code{beta.TOT} - Total multiple-site dissimilarity;

9. \code{beta.NES} - Nestedness component of multiple-site dissimilarity;

10. \code{beta.TUR} - Turnover component of multiple-site dissimilarity;

11. \code{mean.geodist} - If \code{xy.coords} is provided, the mean linear euclidean distance between sites in the a window is returned.

12. \code{focal.geodist} - If \code{xy.coords} is provided, the mean linear euclidean distance of the focal site in relation to its neighbours in the window is returned. 
}

\references{
Anderson, M.J., K.E. Ellingsen and B.H. McArdle. 2006. Multivariate dispersion as a measure of beta diversity. Ecology Letters 9: 683-693.

Baselga, A. 2010. Partitioning the turnover and nestedness components of beta diversity. Global Ecology and Biogeography 19: 134-143.

Baselga, A. 2017. Partitioning abundance-based multiple-site dissimilarity into components: balanced variation in abundance and abundance gradients. Methods in Ecology and Evolution 8: 799-808.

Dala-Corte, R.B., L.F. Sgarbi, F.G. Becker and A.S. Melo. 2019. Beta diversity of stream fish communities along anthropogenic environmental gradients at multiple spatial scales. Environmental Monitoring and Assessment 191:288.

Legendre, P. and M. De Caceres. 2013. Beta diversity as the variance of community data: dissimilarity coefficients and partitioning. Ecology Letters 16: 951-963.
}

\author{Luciano F. Sgarbi, Renato B. Dala-Corte and Adriano S. Melo}

\seealso{\code{\link[vegan]{vegdist}}, \code{\link[vegan]{betadisper}}, \code{\link[adespatial]{beta.div}}, \code{\link[betapart]{beta.multi}}}

\examples{
## Example 1. A simmulated community matrix with a known structure of increasing
##  beta diversity by turnover
# n is the total sample sites
# LocS is the number of spp per site
# MaxS is the total number of spp in the matrix

# All samples will contain LocS species. The first sample will contain presences
# for the first LocS species. The subsequent samples will contain LocS presences
# spread over a increasing set of species. The assignment of presences for the
# second sample to the last sample is done randomly. The last sample will
# contain LocS presences assigned randomly to the MaxS species. Thus, for a 
# window size of 3 (ws=3) and a dataset of 10 samples, beta diversity for the
# samples 1-3 will be much lower than for samples 8-10.   

SimComm <- function(n = 21, MaxS = 30, LocS = 10){
    s <- seq (LocS, MaxS, length.out = n)
    mat <- matrix(0, n, MaxS, dimnames = 
                      list(paste("site", 1:n, sep = "_"), 
                           paste("sp", 1:MaxS, sep = "_")))
    for(i in 1:n){
        mat[i, sample(1:s[i], LocS)] <- 1
    }
    mat <- mat[, colSums(mat)!=0]
    return(mat)
}

mat <- SimComm(n = 21, MaxS = 30, LocS = 10)

#Creating an environmental gradient:
grad <- 1:nrow(mat)

b.resu <- betaRegDisp(y = mat, x = grad, xy.coord = NULL, ws = 3, 
                      method.1 = "jaccard", 
                      method.2 = "ruzicka",
                      method.3 = "ruzicka", 
                      independent.data = FALSE, illust.plot = FALSE)

##Ploting all the output of the object for the simmulated community
op <- par(no.readonly = TRUE)
par(mfrow = c(5, 2), oma = c(1, 0, 1, 0.1), mar = c(1.5, 3, .1, .1), cex = 1, las = 0)
for(i in 1:ncol(b.resu)){
  plot(b.resu[, 1], b.resu[, i], ylab = colnames(b.resu)[i], cex.lab = .9, 
       cex.axis = 0.9, tcl = -0.2, mgp = c(1.5, .2, 0), pch = 15, col = "grey")
}
mtext("Environmental gradient", cex = 1.3, 1, -0.1, outer = TRUE)
par(op)


\dontrun{
## Not run: 
##Example 2
data(varespec)
data(varechem)
grad <- varechem[, "Baresoil"]
resu <- betaRegDisp(y = varespec, x = grad, ws = 3, method.1 = "jaccard", 
            method.2 = "ruzicka", method.3 = "ruzicka", 
            independent.data = FALSE, illust.plot = FALSE)

#Plotting all the outputs of the function:
op <- par(no.readonly = TRUE)
par(mfrow = c(5, 2), oma = c(1, 0, 1, 0.1), mar = c(1.5, 3, .1, .1), cex = 1, las = 0)
for(i in 1:ncol(resu)){
  plot(resu[, 1], resu[, i], ylab = colnames(resu)[i], cex.lab = .9, 
       cex.axis = 0.9, tcl = -0.2, mgp = c(1.5, .2, 0), pch = 15, col = "grey")
}
mtext("Environmental gradient", cex = 1.3, 1, 0, outer = TRUE)
par(op)
## End(Not run)
}
}