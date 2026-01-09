#Key areas for threatened Neotropical glassfrogs

library(readr)
library(raster)
library(sf)
library(usdm)
library(sdm)
library(ecospat)

occ <- read.csv("Glassfrogs.csv")
occ_sf <- st_as_sf(
  occ,
  coords = c("lon", "lat"),
  crs = 4326)
bio5 <- rast("bio5.tif")
bio5
occ_sf <- st_transform(occ_sf, crs(bio5))
r_occ <- rast(bio5)
values(r_occ) <- NA
occ_vect <- vect(occ_sf)
r_occ <- rasterize(
  occ_vect,
  r_occ,
  field = 1,
  fun = "max")

r_species <- lapply(unique(occ$species), function(sp) {
  
  v <- occ_vect[occ_vect$species == sp]
  
  r <- rast(bio5)
  values(r) <- NA
  
  rasterize(v, r, field = 1, fun = "max")
})
names(r_species) <- unique(occ$species)

#SDM
xy <- read.csv("Glassfrogs_species.csv")
bios <- stack(list.files(path = "D:/Desktop/Glassfrogs/SDM//Bios", pattern = "tif", full.names = T))
Bio2 <- bios$bio_2
plot(Bio2)
abs_obs <- 1000
pseudo_abs_train <- ecospat.rand.pseudoabsences(nbabsences=abs_obs[1],
                                                glob=bios.all,colxyglob=1:2, colvar=1:2, presence=presences,
                                                colxypresence=1:2, mindist=0.04166667)
DataSpecies <- SpatialPointsDataFrame(coords = xy, data = DataSpeciesTrain.spp,
                                      proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
DataSpecies <- DataSpecies[,c(-1,-2)]

d1 <- sdmData(formula=Rep1~., train=DataSpecies, predictors=v1)
m1 <- sdm(Rep1~, data=d1, methods= c("brt", "rf", "svm"), test.percent=30, n=10)
m1
bios_ps <- stack(list.files(path = "D:/Desktop/Glassfrogs/SDM//Bios", pattern = "tif", full.names = T))
names(bios_ps) <- names(bios)
plot(bios_ps[[2]])
e2a <- ensemble(m1,newdata=bios_ps,filename='ps.tif', setting=list(method='weighted',stat='TSS'))
#Binary threshold of 0.9 occurrences
gl1 <- ecospat.mpa(Model, Coordinates, perc = 0.9)
gl1
bin.e1 <- ecospat.binary.model(Model, gl1)
library(letsR)
PAM_Glassfrogs <- lets.presab(
  bin.e1,
  xmn = -82, xmx = -75,
  ymn = -7,  ymx = 3,
  resol = 0.1,
  show.matrix = TRUE,
  crs = CRS("+proj=longlat +datum=WGS84"),
  cover = 0
)
PAM_Glassfrogs

#Future 
bios_future <- stack(list.files(path = "D:/Bios_GISS_370", pattern = "tif", full.names = T))
names(bios_future) <- names(bios)
plot(bios_future[[1]])
e2a <- ensemble(m1,newdata=bios_future,filename='GISS_370.tif', setting=list(method='weighted',stat='TSS'))
plot(e2a)

#TD - PD Glassfrogs 
library(picante)
library(ggplot2)
Phy <- read.nexus("Phy_Glassfrogs")
pd <- pd(PAM_Glassfrogs, Phy, include.root = TRUE)
head(pd)
modelPD <- loess(PD ~ TD, data = pd)
summary(modelPD)
str(modelPD)
pd_res$resPD <- residuals(modelPD)
head(pd_res)

ggplot(pd, aes(x = TD, y = PD)) +
  geom_point(
    color = "steelblue",
    size = 2,
    alpha = 0.8) + geom_smooth(
    method = "loess",
    se = TRUE,
    color = "black",
    linewidth = 1) + labs(
    title = "PD and TD",
    x = "TD",
    y = "PD") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"))

#ED
occ <- read.csv("Glassfrogs.csv")
head(occ)
sp_csv <- unique(occ$species)
length(sp_csv)
Phy <- read.nexus("Phy_Glassfrogs")
Phy
sp_tree <- Phy$tip.label
sp1_v <- intersect(sp_csv, sp_tree)
sp_missing <- setdiff(sp_csv, sp_tree)
length(sp1_v)
length(sp_missing)
sp_missing
Phy_pruned <- drop.tip(Phy, setdiff(sp_tree, sp1_v))
ED <- evol.distinct(Phy_pruned, type = "fair.proportion")
ED_table <- ED[order(-ED$ED), ]
row.names(ED_table) <- NULL
head(ED_table)






