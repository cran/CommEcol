\name{standExtent}
\alias{standExtent}
\title{Standardization of spatial extent for two metacommunities by subsampling}

\description{Communities far from each other tend to differ more than those nearby. Accordingly, spatial extent usually affects species richness and other metacommunity properties. Comparison between two metacommunities may be biased if spatial extent differs. One potential solution is to standardize spatial extent by subsampling of the metacommunity more spread in space.}

\usage{
standExtent(d.large, d.small, ini.large=NULL, ini.small=NULL)}

\arguments{
  \item{d.large }{A square (symmetrical) or triangular distance matrix of the metacommunity with large spatial extent.}
  \item{d.small }{A square (symmetrical) or triangular distance matrix of the metacommunity with small spatial extent.}
  \item{ini.large }{The first community to be used in the agregation process of the metacommunity with large spatial extent. Make sure it is present in \code{d.large}. If absent, a random  site will be used.}
  \item{ini.small }{The first community to be used in the agregation process of the metacommunity with small spatial extent. Make sure it is present in \code{d.small}. If absent, a random  site will be used.}
}

\details{The standardized of spatial extent is done by selecting the same number of sites in both metacommunities using the criterion that within-metacommunity average distance among sites is similar. 

The procedure consists of: 1) selection of a focal site \code{ini.large} in the  metacommunity more spread in space, 2) selection of a second site that is nearest to the focal one, 3) selection of a third site using the criterion of minimal average distance among them, 4) repeat step 3 until all sites are included. The procedure then 5) selects a focal site \code{ini.small} in the metacommunity less spread in space, 6) selects a second site in the metacommunity less spread in space using the criterion of the distance between the focal site and this second one is the most similar to the distance in step 2 above, 7) selects a third site using the criterion that the distance among the three sites is most similar to the average distance in step 3 above, and 8) repeats step 7 until all sites in the less spread metacommunity are included.

Notice that after some sites are accumulated, the difference between the average distances in the two sets of metacommunities will increase abruptly. This will occur when most sites in the periphery of the less spread metacommunity are included.

Users may specify initial sites in the accumulation process (steps 1 and 5 above). If not specified, randomly selected sites will be used.

 The percentage of difference of the two spatial extents is provided in the form: (d.mean.large - d.mean.small)/d.mean.large, where d.mean.large and d.mean.small are, respectively, the average within-metacommunity distances of the more and less spread metacommunities. This difference may be helpful to decide for the number and identity of sites to be used.

}

\value{A dataframe containing the names of the accumulated sites in both metacommunities, their within-metacommunity average distances and the percentage of difference in spatial extents.
}

\references{
Heino, J., A.S. Melo, J. Jyrkankallio-Mikkola, D.K. Petsch, V.Sa. Saito, K.T. Tolonen, L.M. Bini, T.S.F. Silva, V. Pajunen, J. Soininen & T. Siqueira. (in preparation) Higher richness but lower abundance of stream insects in tropical than boreal regions hold at regional, basin and local scales.
}

\author{Adriano Sanches Melo}

\examples{
### Data
large.lat <- seq(2, 16, 2)
large.lon <- seq(2, 16, 2)
large.coo <- expand.grid(large.lon, large.lat)
large.coo[, 1] <- large.coo[, 1] + rnorm(64, sd=0.4)
large.coo[, 2] <- large.coo[, 2] + rnorm(64, sd=0.4)
rownames(large.coo) <- as.character(paste("large", 1:64, sep=""))
large.distances <- dist(large.coo) 

small.lat <- 1:8 # notice the spatial extent here is much
                  # smaller than that of the "large" set above.
small.lon <- 1:8
small.coo <- expand.grid(small.lon, small.lat)
small.coo[, 1] <- small.coo[, 1] + rnorm(64, sd=0.4)
small.coo[, 2] <- small.coo[, 2] + rnorm(64, sd=0.4)
rownames(small.coo) <- as.character(paste("small", 1:64, sep=""))
small.distances <- dist(small.coo) 

### Example 1 - Graphical demonstration of the subsampling process 
resu <- standExtent(d.large = large.distances, d.small=small.distances,
                    ini.large="large52", ini.small="small45")
                    
plot(2:64, resu[2:64, "percent.dif"]) 
# Notice the subit increase of the curve. This happens because the average distance 
# in the small metacommunity can not follow the increase in its large metacommunity.

par(mfrow=c(2,1))
par(mar=c(1.7, 2, 1.2, 1.2))
plot(large.coo, xlim=c(0, 17), ylim=c(0, 17))
text(10, 16, "Large")
points(large.coo[resu[1, "site.large"],], col="blue", pch=15)
count <- 1
threshold <- 0
while(threshold <= 5 | count <= 7){ # using a threshold of 5 percent or at least 7 
                                      # sites are accumulated. The later is due to
                                      # large variations when number of sites is low.
   count <- count+1
   points(large.coo[resu[count,"site.large"],], col="red", pch=15)
   threshold <- abs(resu[count, "percent.dif"])
   Sys.sleep(0.05)
} # Wait a little bit and watch the plot until the plotting of red dot stops.

par(mar=c(1.7, 2, 1.2, 1.2)) 
plot(small.coo, xlim=c(0, 17), ylim=c(0, 17)) # notice the reduced spatil extent.
text(10, 16, "Small")
points(small.coo[resu[1,"site.small"],], col="blue", pch=15)
count <- 1
threshold <- 0
while(threshold <= 5 | count <= 7){
   count <- count+1
   points(small.coo[resu[count,"site.small"],], col="red", pch=15)
   threshold <- abs(resu[count, "percent.dif"])
   Sys.sleep(0.05)
}

### Example 2 - Generation of 3 paired standardized spatial extents
### using each of the three first sites in the larger spatial extent
### as ini.large.
names.large <- rownames(large.coo)[1:3]
resu.list <- vector(mode="list", length=3)
names(resu.list) <- names.large
for(i in names.large){
   print(i)
   resu.list[[i]] <- standExtent(d.large = large.distances, 
                                 d.small=small.distances,
                                 ini.large=i)
}
resu.list
}
