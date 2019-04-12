---
layout: post
title: Plot a map with rivers in R
date: 2019-03-22 15:47:36 +01:00
categories: Programming
tags: illustration, maps, R
image: /assets/posts/2019-03-22-rivers-map-in-r/Cobitis_invasive_map_ed.jpeg
alt: map with rivers in R
description: R is a great tool to make maps including maps with rivers. I would like to share how I created a map with rivers in R using the simplest code possible.
---

<p>I recently needed to create a map for a publication where we study fish species. So, showing sampling locations on a map with rivers was a requirement. An obvious solution for that was of course to use R. And although making a map with rivers in R turned out to be easy, I spent half the day searching for a solution.</p>

<p>I use R to create maps for <a href="{{ site.baseurl }}/publications/">all my publications</a>. It is a free, simple, and precise way to plot points on a nice-looking map. If I need to add more things to the map, I save it in SVG format and edit it in <a href="https://inkscape.org/" target="_blank">Inkscape</a>, an open-source editor for vector graphics.</p>

<!--more-->

<h2>Make a map plot in R</h2>

<p>In my view, the best way to create maps and plot sampling locations in R is to use the <a href="https://cran.r-project.org/web/packages/maps/index.html" target="_blank">maps</a> library. It produces clean and informative maps.</p>

<p>Usually, I would use the code like the one below.</p>

```r
library(maps)
library(rworldmap) # to plot axes

newmap <- getMap(resolution = "hight")

lat1 <- runif(5,  35, 55)
long1 <- runif(5,  0, 50)
lat2 <- runif(5,  35, 55)
long2 <- runif(5,  0, 50)

svg("map_Europe.svg", height=4, width=6)
par(mar=c(3, 3, 2, 2))
waterColor <- 'cadetblue1'
plot(newmap, xlim = c(0, 50), ylim = c(35, 70),
     asp = 1,lty=2, lwd=1,
     bg=waterColor, col="#ffffff", main = "Europe")
map.axes()
points(x=long1, y=lat1, pch= 15, col="black", cex=1.2)
points(x=long2, y=lat2, pch= 17, col="red", cex=1.2)
legend("topright", bg="white", pt.cex = 1.2,
       legend = c('group1', 'group2'), col=c('black', 'red'), pch = c(15, 17))
dev.off()
```

<p>It produces a map of Europe with random points representing two tentative groups.</p>
<div class="wp-block-image">
<figure class="caption"><img src="{{ site.baseurl }}/assets/posts/2019-03-22-rivers-map-in-r/map_Europe.jpeg" alt="A map of Europe plotted in R" class="wp-image-1598" /></figure>
</div>

<p>You can modify this code to make any map you need by providing your coordinates of the points and changing colors, lines, etc. For example, I used this code with some post-editing in Inkscape to create a map for <a href="https://doi.org/10.1371/journal.pgen.1007949" target="_blank">our recent publication in PLOS Genetics</a>.</p>
<figure class="caption"><img src="{{ site.baseurl }}/assets/posts/2019-03-22-rivers-map-in-r/journal.pgen.1007949.g001.jpeg" alt="The map of sampling points plotted in R, with the distribution ranges added in Inkscape." class="wp-image-1602" />
<figcaption class="caption">The map and sampling points are plotted in R, and the distribution ranges are added in Inkscape.</figcaption>
</figure>

<h2>Add rivers to a map in R</h2>

<p>I liked this map style and I did not want to change it. I only needed to add rivers to this map. However, the online search often provided totally different ways to make maps including ggmap and ggplot2 libraries. This was not what I wanted.</p>

<p>I will not list all the options I tried but believe me I tried many. In the end, <strong>the most optimal solution</strong> was to find a <a href="https://en.wikipedia.org/wiki/Shapefile" target="_blank">shapefile</a> of needed rivers and add them to the map.</p>

<p>I needed a river map of Europe. There are many shapefiles of European water systems online. I tried the top 10 files I found in Google and I liked the most the one published on <a href="https://tapiquen-sig.jimdo.com/english-version/free-downloads/europe/" target="_blank">this website</a>.</p>

<p>I download the <em>Europe_Hydrography.rar</em> file, extracted it, loaded it to R and added rivers to my previous map with this code:</p>

```r
library(rgdal)
riversData <- readOGR("Europe_Hydrography") # load the shapefile
plot(riversData, col=waterColor, add=T) # plot rivers
```

<p>This is the map of Europe with rivers I obtained.</p>
<div class="wp-block-image">
<figure class="caption"><img src="{{ site.baseurl }}/assets/posts/2019-03-22-rivers-map-in-r/map_Europe_rivers.jpeg" alt="A map of Europe with major rivers plotted in R" class="wp-image-1606" />
<figcaption class="caption">Rivers are added using a shapefile and the <strong>rgdal</strong> R library</figcaption>
</figure>
</div>

<h2>Code to make a map with rivers in R</h2>

<p>Here is the whole code to make this rivers map of Europe.</p>

```r
library(maps)
library(rworldmap) # to plot axes
library(rgdal) # to load the shapefile

newmap <- getMap(resolution = "hight")
riversData <- readOGR("Europe_Hydrography") # load the shapefile

lat1 <- runif(5,  35, 55)
long1 <- runif(5,  0, 50)
lat2 <- runif(5,  35, 55)
long2 <- runif(5,  0, 50)

svg("map_Europe.svg", height=4, width=6)
par(mar=c(3, 3, 2, 2))
waterColor <- 'cadetblue1'
plot(newmap, xlim = c(0, 50), ylim = c(35, 70),
     asp = 1,lty=2, lwd=1,
     bg=waterColor, col="#ffffff",
     main = "Europe")
plot(riversData, col=waterColor, add=T) # plot rivers
map.axes()
points(x=long1, y=lat1, pch= 15, col="black", cex=1.2)
points(x=long2, y=lat2, pch= 17, col="red", cex=1.2)
legend("topright", bg="white", pt.cex = 1.2,
       legend = c('group1', 'group2'), col=c('black', 'red'), pch = c(15, 17))
dev.off()
```

<p>Actually,  this is the code I used with some modification to make the rivers map with sampling locations for our upcoming publication. I also added the names of major rivers and migration arrows using Inkscape. You have seen the final image at the top of this post.</p>
<div class="wp-block-image">
<figure class="caption"><img src="{{ site.baseurl }}/assets/posts/2019-03-22-rivers-map-in-r/map_Ukraine.jpeg" alt="A rivers map plot in R" class="wp-image-1608" />
<figcaption class="caption">Sampling locations on the rivers map of Eastern Europe produced with R.</figcaption>
</figure>
</div>

<h2> Conclusion</h2>

<p>R is proven to be a universal tool for making great scientific illustrations for free. Making a map with rivers is R is also simple and effective. You only need to find the right shapefile of rivers you want to plot. Moreover, you can use shapefiles of any other geospatial objects (roads, mountains, climate data, etc) to add them to the map.</p>
