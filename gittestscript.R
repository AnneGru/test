# test script
# 18 march '21 Anne Grundlehner
setwd("/Users/Anne/Documents/Internship WMR zeehonden/Data Analyses")


#https://www.r-graph-gallery.com/325-background-map-from-geojson-format-in-r.html
#import and plot annotations testfile geojson
library(sp)
library(rgdal)
library(raster)

#### Basics ####
#plot
test.json = readOGR("copy-of-greyseals2_10-02-21_1228.geojson")
plot(test.json) #they are polygons
axis(1) #add x-axis
axis(2) #add y-axis
box()
grid()

#number of seals in full image
length(test.json)
#list of surface area's occupied per seal
test.json$area_m2
#Total occupied surface area in full image and some statistics
sum(test.json$area_m2)
summary(test.json@data)
hist(test.json$area_m2) #distribution of seal sizes (surface area m2)

#### Distances ####
# https://rpubs.com/dgolicher/9458
library(sp)
library(rgeos) #edge to edge distances)
install.packages("sf")
library(sf)
library(spdep) #for distance matrix

#convert CRS to OS grid
gArea(test.json)
test.json1  = spTransform(test.json, CRS("+init=epsg:4326"))
test.json2 <- spTransform(test.json, CRS("+init=epsg:27700"))
#Surface Area m2
gArea(test.json2) # should be similar to sum(test.json$area_m2) > it's somewhat smaller
nf = test.json2

#This matrix contains the distances between each patch and all the others
dist.mat <- gDistance(nf, byid = TRUE)
# obtain the minimum values
# remove self-self distances (0)
mat=dist.mat
mat[mat==0] <- NA
result <- t(sapply(seq(nrow(mat)), function(i) {
  j <- which.min(mat[i,])
  c(paste(rownames(mat)[i], colnames(mat)[j], sep='/'), mat[i,j])
}))
print(result)
nnd = as.numeric[result[,2]]
nnd.50m <- nnd[nnd <= 50] #remove outliers
hist(nnd, breaks=5000, freq=F, xlim=c(0,5))


#### Adults ####
hist(nf$area_m2, breaks=500, freq=F) 
adults = subset(nf, nf$area_m2 >= 0.5) #assuming pups <0.5 m2
length(adults) #number of observations left
mat2 <- gDistance(adults, byid = TRUE) #distance matrix
mat2[mat2==0] <- NA #remove distances to self
result2 <- t(sapply(seq(nrow(mat2)), function(i) {
  j <- which.min(mat2[i,])
  c(paste(rownames(mat2)[i], colnames(mat2)[j], sep='/'), mat2[i,j])
}))
nnd2 = as.numeric(result2[,2])
nnd2.50m <- nnd[nnd <= 50] #remove outliers above 50m distance
hist(nnd2, breaks=5000, freq=F, xlim=c(0,5))


#### Radius ####
# To analyze close to the seal
# Use the distance matrix, remove all values above a certain threshold distance
# in this way, data outside this radius is not taken into account

mat.radius = mat #make copy of matrix to work with
# Determine radius size
r = 20
mat.radius[mat.radius >= r] <- NA #remove all values above 20m become NA
#Maybe make a for-loop

## (1) Number of seals within radius
mat.radius.number = mat.radius
mat.radius.number[mat.radius.number <= r] <- 1 # set all values to 1
hist(mat.radius.number) #check, should only be 1 and NA
#number of individuals around each seal within this radius
#take the sum of each row to get # of seals within [r] m distance of focal seal
n.in.radius = as.data.frame(rowSums(mat.radius.number, na.rm=T))
hist(as.matrix(n.in.radius))

## (2) Occupied surface area within radius
# Use the 1/NA matrix to get a matrix with only the areas of the seals within the radius
m2.seals.within = (mat.radius.number * (test.json$area_m2)) #sizes = m2 occupied per indiv
m2.in.radius = as.data.frame(rowSums(m2.seals.within, na.rm=T)) #take the sum of each row
hist(as.matrix(m2.in.radius))
# can be expressed as fraction of the total area of the radius
areafraction.of.radius = (m2.in.radius) / (pi*(r^2))
head(areafraction.of.radius)
hist(as.matrix(areafraction.of.radius))

#### Randomisation ####
install.packages("spatstat")
library(spatstat)

## Generating a bounding box around our polygons
# In order to have a field for randomisations
#library(sf)
#plot(test.json)
#st_bbox(test.json) #returns coordinates for bounding box
#plot(st_make_grid(test.bbox, n=1), add=TRUE) ## de test is redelijk hoog afgesneden door foutieve annotatie in zee
## OR make a convex hull instead
library(sp)
library(rgdal)
set.seed(1)
dat <- as.data.frame(cbind(test.json$lat, test.json$lon))
ch <- chull(dat)
coords <- dat[c(ch, ch[1]), ]  # closed polygon
plot(dat,cex=0.1)
lines(coords, col="red")


# create single ring polygon, based on the convex hull coords
c1 = coords
r1 = rbind(c1, c1[1, ])  # join (close the polygon)
P1 = Polygon(r1)
Ps1 = Polygons(list(P1), ID = "a")
#simulate points:
rpoints <-spsample(Ps1,n=length(test.json),type="random") #generate random points
plot(dat, cex=0.1) #plot original data
points(rpoints, pch=19, cex=0.1, col="blue") #plot simulated points
axis(1)
axis(2)

# Add simulated points(x,y) to the original polygons
simulations = test.json
sp = SpatialPoints(test.json) #make polygondataframe into spatialpoints object
class(rpoints)
draw.ellipse()
## Shift the original polygons to the new locations
## change their centroid points?
#or create ellipse around the points
### best idea, model seal data shape, based on this create ellipses randomly

#### Random:Circles ####
#creating circles around points, works
library(sampSurf)
x=coords[1,1]
y=coords[1,2]
circle = spCircle(radius=0.0001, spUnits=crs(test.json), centerPoint=c(x=x,y=y))
#plot(test.json)
plot(circle$spCircle, col="red")
points(rpoints, cex=1)

#### Random: shift? ####

# Generate random points using spsample
# Calculate the difference between this and the centroids of the polygons of our actual data
# Use elide() to shift the original polygons with the calculated dx dy

#rpoints are our random points within region
#test.json is the (test) annotation file
simulated.coords = as.matrix(rpoints@coords) 
annotated.coords = as.matrix(cbind(test.json$lat, test.json$lon))
dx = (simulated.coords[,1] - annotated.coords[,1])
dy = (simulated.coords[,2] - annotated.coords[,2])
change.in.coords = as.data.frame(cbind(dx, dy))
#now we want to shift our original polygons with this difference
library(maptools)
polygons.shifted = elide(test.json[1,], shift=c(dx[1], dy[1])) ## dit werkt
#test if changed, WORKS
plot(test.json[1,], axes=TRUE)
plot(polygons.shifted[1], axes=T)
# create a for-loop

polygons.shifted = test.json[1,]

for (i in 1:length(test.json)){
  polygons.shifted[1,] = elide(test.json[i,], shift=c(dx[i], dy[i]))
}

all.polys = test.json[1,]
for (i in 1:length(test.json[i,])){
  polygons.shifted.looped[i,] = elide(test.json@polygons[i], shift=c(dx[i], dy[i]))
  #all.polys = SpatialPolygons(list(all.polys, polygons.shifted.looped[i,]))
}
plot(all.polys)
plot(polygons.shifted.looped)



for (i in 1:length(polygons.shifted)){
  ps1[i,] = elide(test.json[i,], shift=c(dx[i], dy[i]))
}

for(i in 1:length(test.json)){
  poly = elide(test.json[i,], shift=c(dx[i], dy[i]))
  polygons.shifted[i,] <- poly 
}

plot(test.json, axes=T)
plot(polygons.shifted, axes=T, border="red") # does not work yet


#### Random: Ellipse ####

library(spatstat)
#create ellipse coordinates (5 points) based on centre coords
ellipse.coords = ellipse(a=0.001, b=0.002, centre=ellipse.centre[1,], phi=0)
# create single ring feature POLYGON
c2 = as.data.frame(ellipse.coords) #cbind(x1, y1)
r2 = rbind(c2, c2[1:5, ])  # join
P2 = Polygon(r2)
Ps2 = Polygons(list(P2), ID = "a")
plot(ellipse.coords, border="red") #plot ellipse 
points(rpoints[1:10]) #add some points from the random sampling
#make the ellipse polygon a spatial polygon:
Ps2 = SpatialPolygons(list(Ps2))
plot(Ps1, col="grey", main="Artificial Seal")

## loopje dat nog NIET werkt
many.ellipses = test.json
for (i in 1:length(rpoints)){
  ellipse.5p = ellipse(a=0.01, b=0.02, centre=simulated.coords[i,], phi=0)
  ci = as.data.frame(ellipse.5p)
  ri = rbind(ci, ci[1:5,])
  poly.i = Polygon(ci)
  Ps.i = Polygons(list(poly.i), ID = "a")
  many.ellipses@polygons[i] <- Ps.i
}

plot(many.ellipses, axes=T, border="red")
plot(test.json)

#### Visualization ####

library(tripack)
library(ggplot2)
library(ggvoronoi)
plot(voronoi.mosaic(x=test.json[,1], y=test.json[,2], duplicate="remove"))
set.seed(45056)
x <- test.json$lat
y <- test.json$lon
points <- data.frame(x, y,
                     distance = nnd)
# scatter with distances colored
ggplot(points) +
  geom_point(aes(x,y,color=distance), size=0.1)
# Voronoi, colored by distance
ggplot(points) +
  geom_voronoi(aes(x,y,fill=distance))
# similar but outlines by chull
# Ps1 is the convex hull
Ps1 = SpatialPolygons(list(Ps1)) #make ps a spatial polygon version
ggplot(data=points, aes(x=x, y=y, fill=distance)) + 
  geom_voronoi(outline = Ps1) +
  scale_fill_gradient(low="#4dffb8",high="black",guide=F) #ander kleurje
#Voronoi
ggplot(points,aes(x,y)) +
  stat_voronoi(geom="path") +
  geom_point(size=.0125, color="coral")


#force-directed graph drawing algorithms to visualize a distance matrix
dist_m <- dist.mat[1:100, 1:100] ## !! ONLY the first 100 taken as test
dist_mi <- 1/dist_m # one over, as qgraph takes similarity matrices as input
library(qgraph)
jpeg('example_forcedraw.jpg', width=1000, height=1000, unit='px')
qgraph(dist_mi, layout='spring', vsize=3)
dev.off()

#### Polygon to Points ####
# obtain centroids from the seal polygons
#using this ppp object we can ake a distance plot
# and optionaly use it for comparison with the randomly generated point pattern

centers = cbind(test.json$lat, test.json$lon)
plot(centers)
centers = SpatialPoints(centers)
plot(centers)
plot(rpoints)

#### TO DO #####

# Gstat "predicted" concentration

#Spatial patterns: distance map
# Convert to ppp point object
library(spatstat)
ppp(rpoints)
distmap(rpoints[,1], rpoints[,2])
