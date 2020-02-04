setwd("H:/Desktop/Spatial/Homework 2")    # Work 
setwd("C:/Users/keete/Documents/Fall 2019/Spatial Stats/Homework 2")   # Not work 
### Homework 3 

# Libraries 

library(geoR)
library(tidyverse)
library(readxl)
library(sf)
library(sp)
library(mvtnorm)
library(gstat)
library(ggfortify)
library(maps)
library(maptools)
library(sf)
library(tmap)
library(geostatsp)
library(RColorBrewer)


# Question 1 

# Prelim 

colo <- read_excel("co_precip.xls")             # Read in 

ggplot(
  data = colo,
  aes(
    x = Elevation,
    y = Precip
  )
) + 
  geom_point(
    
  ) + 
  geom_smooth(
    method = "lm"
  ) + 
  theme_classic(
    
  ) + 
  labs(
    title = "Colorado Precipitation",
    x = "Elevation (meters)",
    y = "Precipitation (mm)"
  )

summary(colo)

colo_sf <- st_as_sf(
  colo,
  coords = c("Longitude", "Latitude"),
  crs = CRS("+proj=utm +zone=13")
)

head(colo_sf)

# Variogram for Nuggest estimation 

cutoff <- .7*max(
  dist(
    cbind(
      colo$Longitude,
      colo$Latitude
    )
  )
)

bins <- 10

variogram <- variogram(
  logPrecip ~ 1, 
  locations = ~Longitude + Latitude,
  data = colo,
  cutoff = cutoff,
  width = cutoff/bins
)

plot(
  variogram,
  cex = 2, 
  pch = 19,  
  ylab = expression(
    paste(
      "Average", (0.7*(Y(s[i]) - Y(s[j])^2))
    )
  ),
  xlab = "Euclidean Distance (km)"
)

# Maximum likelihood by hand 

H <- as.matrix(
  dist(
    st_coordinates(
      colo_sf
    )
  )
)

X <- model.matrix(
  ~ scale(Zelevation),
  data = colo_sf
)

Y <- as.matrix(
  colo_sf$logPrecip
)

MVN_fun <- function(theta, H, X, Y){
  
  part_sill <- theta[1]
  phi <- theta[2]
  p <- theta[3]
  nugget <- theta[4]
  
  Sigma <- part_sill*exp(-(phi*H)^p) + nugget*diag(length(Y))
  
  Sigma_inverse <- solve(Sigma)
  
  Beta <- solve(t(X) %*% Sigma_inverse %*% X) %*% t(X) %*% Sigma_inverse %*% Y 
  
  log_L <- mvtnorm::dmvnorm(
    x = t(Y),
    mean = X %*% Beta,
    sigma = Sigma,
    log = TRUE
    )
  
  return(-1*log_L)
  
}


MLE <- optim(
  c(0.1, 1, 0.5, 0.1),
  fn = MVN_fun,
  H = H,
  X = X,
  Y = Y,
  hessian = TRUE
  )

MLE

# Maximum likelihood by geoR

d <- data.frame(
  st_coordinates(
    colo_sf
  ),
  colo[, -1:-2]
)

d_geo <- as.geodata(
  d,
  coords.col = 1:2,
  data.col = 4,
  covar.col = 6
)

fit_geoR <- geoR::likfit(
  d_geo,
  ini.cov.pars = c(0.09, 1),
  nugget = .1,
  fix.nugget = FALSE,
  trend = ~ scale(Zelevation),
  cov.model = "powered.exponential",
  kappa = 1,
  fix.kappa = FALSE,
  lik.method = "ML"
)

fit_geoR

summary(fit_geoR)

MLE

q1_model_comp <- data.frame(
  model = c("By Hand", "geoR"),
  log.likelihood = c(MLE$value, fit_geoR$loglik),
  nugget = c(MLE$par[4], fit_geoR$nugget),
  partial.sill = c(MLE$par[1], fit_geoR$sigmasq),
  range_parm = c(MLE$par[2], 1/fit_geoR$phi),
  power = c(MLE$par[3], fit_geoR$kappa)
  )
######### Question 2

rad_df <- read_excel("Fukushima_30km_2016.xlsx")

# Splorin' 

summary(rad_df)

rad_df <- rad_df %>%                 # Change Rain to 0/1
  mutate(
    Rain_code = ifelse(
      test = Rain == "Rain",
      yes = 1,
      no = 0
      ),
    interaction = Dist_km*Rain_code
    )


model <- lm(
  log_Dose ~ Dist_km + Dist_km_sq,
  data = rad_df
)

summary.lm(model)

autoplot(
  model,
  smooth.colour = NA
) + 
  theme_classic(
    
  )

rad_df <- rad_df %>% 
  mutate(resid = resid(model))

rad_df_points <- st_as_sf(
  x = rad_df,
  coords = c(
    "LONG",
    "LAT"
  ),
  crs = CRS(
    "+proj=utm +zone=54"
  )
)

data(world.cities)
data(World, metro, rivers, land)

world.cities <- world.cities %>% 
  filter(
    country.etc == "Japan"
    )

japan_cities <- st_as_sf(
  x = world.cities,
  coords = c(
    "long",
    "lat"
  ),
  crs = CRS(
    "+proj=utm +zone=54"
  )
)

summary(japan_cities)
summary(rad_df_points) 

japan_map <- map(
  database = "world",
  region = "Japan",
  plot = F,
  fill = T
)
  
summary(japan_map)

japan_fill <- map2SpatialPolygons(
  japan_map,
  IDs = japan_map$names,
  proj4string = CRS(
    "+proj=utm +zone=54"
  ) 
)

max_lon <- max(rad_df$LONG)
min_lon <- min(rad_df$LONG)

max_lat <- max(rad_df$LAT)
min_lat <- min(rad_df$LAT)

correction <- .12

tm_shape(
  japan_fill,
  proj4string = CRS(
    "+proj=utm +zone=54"
  ),
  ylim = c(
    min_lat - correction, 
    max_lat + correction
  ),
  xlim = c(
    min_lon - correction, 
    max_lon + correction
    )
  ) + 
  tm_fill(
    col = "lightgreen"
  ) + 
  tm_borders(
    lwd = 2,
    col = "black"
  ) + 
  tm_shape(
    rad_df_points
  ) +
  tm_symbols(
    col = "resid",
    palette = "-Spectral",
    n = 7,
    style = "jenks",
    border.lwd = 0.1,
    border.col = 'gray',
    alpha = 0.9,
    scale = .75,
    title.col = "Residuals"
  ) + 
  tm_shape(
    japan_cities
  ) + 
  tm_text(
    "name",
    textNA = "",
    remove.overlap = T,
    shadow = T
  ) +
  tm_legend(
    position = c(
      "left", 
      "bottom"
    ),
    legend.outside = TRUE,
    frame = F,
    main.title = 'Nonspatial Residuals (OLS)'
  ) +
  tm_layout(
    bg.color = "skyblue",
    saturation = .85
  ) + 
  tm_compass(
    north = 0,
    type = "arrow"
  ) 

cutoff_rad <- .65*max(
  dist(
    cbind(
      rad_df$LONG,
      rad_df$LAT
    )
  )
)

bins_rad <- 10

variogram <- variogram(
  log_Dose ~ 1, 
  locations = ~LONG + LAT,
  data = rad_df,
  cutoff = cutoff_rad,
  width = cutoff_rad/bins_rad
)

plot(
  variogram,
  pch = 16, 
  cex = 1.5
  )

######################### Isotropic Models 

### geostatsp

rad_d <- data.frame(
  st_coordinates(
    rad_df_points
  ),
  rad_df[, -1:-2]
)

rad_sp_iso  <- SpatialPointsDataFrame(
  coords = rad_d[, 1:2],
  data = rad_d,
  proj4string = CRS(
    "+proj=utm +zone=54"
  )
)

rad_geostat_iso <- lgm(
  log_Dose ~ Dist_km + Dist_km_sq,
  data = rad_sp_iso,
  grid = 500,
  shape = 0.5,
  fixShape = TRUE,
  fixNugget = FALSE,
  aniso = FALSE,
  reml = TRUE
)

summary(rad_geostat_iso)

################################### Anisotropic Models 

# Geostatsp

rad_fit_aniso <- lgm(
  log_Dose ~ Dist_km + Dist_km_sq,      
  data = rad_sp_iso,                
  grid = 500,                    
  shape = 0.5,
  fixShape = FALSE,   
  nugget = 0.05,
  fixNugget = FALSE,  
  aniso = TRUE,
  reml = TRUE
  ) 

summary(rad_fit_aniso)


# Plots 

reactor <- c(141.0298, 37.4245)

par(mfcol=c(2,2), mai=c(0.5,0.5,0.5,0.5))

plot(rad_fit_aniso$predict[["krigeSd"]],main='Anisotropic Prediction SDs')
points(reactor[1],reactor[2],pch="R")

plot(rad_geostat_iso$predict[["krigeSd"]],main='Isotropic Prediction SDs')
points(reactor[1],reactor[2],pch="R")

plot(rad_fit_aniso$predict[["predict"]],main='Anisotropic Predicted log(Dose)')
points(reactor[1],reactor[2],pch="R")

plot(rad_geostat_iso$predict[["predict"]],main='Isotropic Predicted log(Dose)')
points(reactor[1],reactor[2],pch="R")

model_comp <- data.frame(
  model = c("OLS", "Anisotropic", "Isotropic"),
  AIC = c(AIC(model), AIC(rad_fit_aniso),AIC(rad_geostat_iso)),
  BIC = c(BIC(model), BIC(rad_fit_aniso), BIC(rad_geostat_iso))
  )

model_comp
