## Title: Geographically weighted accuracy for hard and soft land cover classifications: 4 approaches with coded illustrations
# Alexis Comber and Naru Tsutsumida
# If you have any problems with data / code / versions etc please contact Lex Comber a.comber@leeds.ac.uk

### Part 1. Introduction

### Part 2. Data
## 1.1 load packages
library(sf)
library(sp)
library(tidyverse)
library(spgwr)
library(ggspatial)
library(cowplot)
library(cols4all)
library(scales)
# library(devtools)# install_github('chrisbrunsdon/gwxtab')library(gwxtab)
library(robCompositions)

## Load the data
data <- read.csv("validation_data.csv")
# Convert to spatial data - sf format
data_sf = st_as_sf(data, coords = c("East", "North")) 
# confirm the geographic projection 
st_crs(data_sf) = 32633
# Convert to spatial data - sp format
# Load the contextual layer of the coastline
lib <- st_read("northlib_outline.gpkg", quiet = T)
lib = st_transform(lib, 32633)
# Create the Tripoli layer for context
df = data.frame(name = "Tripoli", x = 329266.6, y = 3640653)
Tripoli = st_as_sf(df, coords = c("x","y")) 
st_crs(Tripoli) = 32633
# define 1 km GWR grid
grid = st_make_grid(lib, 1000, what = "centers")
grid.xy = st_coordinates(grid)

## Define a masking function - only used once
poly.outer_sf = function (exo.object, input.poly, extend = 0) 
{
    box <- function(obj, extend) {
    		xy = st_bbox(obj)
    		x = xy[c(1,3)] + c(-extend, extend)
    	y = xy[c(2,4)] + c(-extend, extend)	
    	poly <- data.frame(x,y) %>% st_as_sf(coords = c("x", "y")) %>%
    		st_bbox() %>% st_as_sfc()
    	st_crs(poly) = st_crs(obj)
    	poly 
    }
    choppo <- function(obj1, obj2, extend) {
        res = box(obj1, extend)
        res = st_difference(res, st_union(obj2))
    }
    res = choppo(exo.object, input.poly, extend)
    st_crs(res) = st_crs(input.poly)
    res
}
# create the mask with a 4km buffer
mask = poly.outer_sf(grid, lib, extend = 2000)

## Figure 1
# The locations of the validation points and a 1 km interpolation grid, with Tripoli highlighted in red."
ggplot() + 
 geom_sf(data = lib, col = "black", fill = "NA") +
	geom_sf(data = grid, pch = 3, col = "lightgrey") +
	geom_sf(data = data_sf, col = "black") +
	geom_sf(data = Tripoli, pch = 20, cex = 6, col = "red") +
	theme_bw() + coord_sf() +
	#annotation_north_arrow(pad_x = unit(13.4, "cm"), pad_y = unit(6.3, "cm")) +
  annotation_scale(pad_x = unit(8.1, "cm"), pad_y = unit(5, "cm")) 

## Define some CM accuracy functions
# Overall accuracy
overall.acc <- function(ctmatrix){
  sum(diag(ctmatrix))/sum(ctmatrix)
}
# User accuracies 
user.acc <- function(ctmatrix, index = 2) {
  ctmatrix[index,index] / sum(ctmatrix[index,])
}  
# Producer accuracy
prod.acc <- function(ctmatrix, index = 2){
  ctmatrix[index,index] / sum(ctmatrix[, index])
}
# Kappa coefficient
kp <- function(x) {
	part.1 <- sum(diag(x)) * sum(x)
	part.2 <- sum(colSums(x) * rowSums(x))
	part.3 <- sum(x)^2
	k <- (part.1 - part.2) / (part.3 - part.2)
	return(data.frame(kappa = k))
}

## Define some class names and Table 1
class.names <- c("1. Bare Ground",
                 "2. Grazing Land", 
                 "3. Urban", 
                 "4. Vegetation", 
                 "5. Woodland") 
class.names2 = c("Bare Ground", "Grazing Land", "Urban", "Vegetation", "Woodland")

## Table 1. The confusion matrix between predicted (rows) and observed (columns) land cover class, with user's and producer's accuracies
tab1 <- tab1_orig <- table(Predicted = data$Boolean_Pred, Observed = data$Boolean_Obs)
# apply in a loop
users = producers = vector()
for (i in 1:5) {
  users = append(users, user.acc(tab1, i))
  producers = append(producers, prod.acc(tab1,i))
}
rownames(tab1) = class.names
oa = unlist(overall.acc(tab1))
ka = unlist(kp(tab1))
tab1 = cbind(tab1, round(users,3), round(producers, 3))
colnames(tab1)[6:7] = c("User's", "Producer's")
tab1

## Fuzzy classes, Figures 2 and 3
# define a plot function
plot_fuzzy = function(df_subset, i, Obs = TRUE){
	coords = data[, 2:3]
  if(Obs) tit = paste0(class.names2[i], " Observed")
	if(!Obs) tit = paste0(class.names2[i], " Predicted")
  tmp = cbind(coords,df_subset)
	tmp %>% filter(df_subset > 0) -> tmp
	tmp = st_as_sf(tmp, coords = c("East", "North")) 
	st_crs(tmp) = 32633
	ggplot() + 
		geom_sf(data = lib, fill = "white", col = "black") +
		geom_sf(data = tmp, aes(size = as.numeric(df_subset)), alpha = 0.5) +
		scale_size(breaks = seq(0.25,1,0.25), range = c(0,4)) +
	  ggtitle(label = NULL,subtitle = tit) +
	  theme_bw() +
	  theme(axis.title=element_blank(),
	        axis.text=element_blank(),
          axis.ticks=element_blank(), 
         	legend.title = element_blank(),
          legend.text = element_text(size=10))
       
}
## Fig2: The observed fuzzy class memberships
pfUo = plot_fuzzy(data[,4], i=3) 
# get the legend and then remove 
legFO = get_legend(pfUo)
pfUo = pfUo + theme(legend.position = "none")
pfVo = plot_fuzzy(data[,5], i=4) + theme(legend.position = "none")
pfWo = plot_fuzzy(data[,6], i=5) + theme(legend.position = "none")
pfGo = plot_fuzzy(data[,7], i=2) + theme(legend.position = "none")
pfBo = plot_fuzzy(data[,8], i=1) + theme(legend.position = "none")
cowplot::plot_grid(pfUo, pfVo, pfWo, pfGo, pfBo, legFO, nrow = 2, ncol = 3)


## Fig3: The predicted fuzzy class memberships.
pfUb = plot_fuzzy(data[,11], i=3, Obs = F) 
# get the legend and then remove 
legFB = get_legend(pfUb)
pfUb = pfUb + theme(legend.position = "none")
pfVb = plot_fuzzy(data[,12], i=4, Obs = F) + theme(legend.position = "none")
pfWb = plot_fuzzy(data[,13], i=5, Obs = F) + theme(legend.position = "none")
pfGb = plot_fuzzy(data[,14], i=2, Obs = F) + theme(legend.position = "none")
pfBb = plot_fuzzy(data[,15], i=1, Obs = F) + theme(legend.position = "none")
cowplot::plot_grid(pfUb, pfVb, pfWb, pfGb, pfBb, legFB, nrow = 2, ncol = 3)

### Part 3. Working with hard classes

### 3.1 Accuracy as GLM and GWGLM probabailities

## GLMs

# logit function 
alogit <- function(x){exp(x)/(1+exp(x))}

# create a binary variable: Observed matches Predicted 
res = (data$Boolean_Pred == data$Boolean_Obs) + 0
# create GLM model
mod0 <- glm(res~1,family= binomial) 
# extract the model coefficients
# summary(mod0) 
mod.coefs <- mod0$coefficients
mod.coefs[2] <-sum(mod.coefs) 
# model the probabilities that the 
# probability of x equals 1 given that y equals 1
mod.p.x.eq.1.g.y.eq.1 <- alogit(mod.coefs[2]) 
# this is the overall accuracy
mod.ovac <- mod.p.x.eq.1.g.y.eq.1 
# cat("overall accuracy:", mod.ovac)

## GLM User's and Producer's accuracies
# Manually: define a class list
class.list <- unique(data$Boolean_Obs)[order(unique(data$Boolean_Obs))]
i = 2
class <- class.list[i]	
# define binary variables where the class is present
obs.class <- (data$Boolean_Obs == class) * 1 	
pred.class <- (data$Boolean_Pred == class) * 1	
# join together for use in the GLMs
obs_pred <- data.frame(cbind(obs.class,pred.class)) 
# User's Accuracy
mod1 <- glm(obs.class~pred.class,data = obs_pred,family= binomial) 
mod.coefs <- mod1$coefficients
mod.coefs[2] <-sum(mod.coefs) 
mod.p.x.eq.1.g.y.eq.1 <- alogit(mod.coefs[2])
mod.user <- mod.p.x.eq.1.g.y.eq.1 
#cat("user's accuracy:", round(mod.user, 3))
# Producer - invert terms
mod2 <- glm(pred.class~obs.class,data = obs_pred,family= binomial) 
mod.coefs <- mod2$coefficients
mod.coefs[2] <-sum(mod.coefs) #logit ea+c
mod.p.x.eq.1.g.y.eq.1 <- alogit(mod.coefs[2]) #n1
mod.prod <- mod.p.x.eq.1.g.y.eq.1 
# cat("producer's accuracy:", round(mod.prod, 3))

# GLM User's and Producer's with a function
do_glm_accuracy = function(pred, obs) {
  # create a binary variable: Observed matches Predicted 
  res = (pred == obs) + 0
  # 1. create the Overall Accuracy model
  mod0 <- glm(res~1,family= binomial) 
  # extract the model coefficients
  mod.coefs <- mod0$coefficients
  mod.coefs[2] <-sum(mod.coefs) 
  # model the probabilities that the 
  # probability of x equals 1 given that y equals 1
  mod.p.x.eq.1.g.y.eq.1 <- alogit(mod.coefs[2]) 
  # this is the overall accuracy
  mod.ov <- mod.p.x.eq.1.g.y.eq.1 
  # 2. User and Producer models 
  # define a class list and some empty results vectors
  class.list <- unique(obs)[order(unique(obs))]
  # create empty vectors for the results
  user_result = prod_result = vector()
  # loop through each class 
  for (i in 1:length(class.list) ) {
    class <- class.list[i]	
    # define binary variables where the class is present
    obs.class <- (obs == class) * 1 	# y in the paper
    pred.class <- (pred == class) * 1	# x in the paper
    # join together for use in the GLMs
    obs_pred <- data.frame(cbind(obs.class,pred.class)) 
    # User's Accuracy
    mod1 <- glm(obs.class~pred.class,data = obs_pred,family= binomial) 
    mod.coefs <- mod1$coefficients
    mod.coefs[2] <-sum(mod.coefs) 
    mod.p.x.eq.1.g.y.eq.1 <- alogit(mod.coefs[2])
    mod.user <- mod.p.x.eq.1.g.y.eq.1 
    # add to the result
    user_result = append(user_result, mod.user)
    # Producer - invert terms
    mod2 <- glm(pred.class~obs.class,data = obs_pred,family= binomial) 
    mod.coefs <- mod2$coefficients
    mod.coefs[2] <-sum(mod.coefs) #logit ea+c
    mod.p.x.eq.1.g.y.eq.1 <- alogit(mod.coefs[2]) #n1
    mod.prod <- mod.p.x.eq.1.g.y.eq.1 
    # add to the result
    prod_result = append(prod_result, mod.prod)
  }
  # combine the outputs and return
  out = data.frame(Producer = prod_result, User = user_result)
  rownames(out) = class.list
  return(list(Overall = mod.ov, UserProd = out))
}

# The function can applied to the data as below, resulting in just 2 vectors of the same length of the predicted and observed classes. It returns a list a object with 2 elements,: the overall accuracy measure and a table of the per class user's and producer's accuracies.
# do_glm_accuracy(data$Boolean_Pred, data$Boolean_Obs)

## GWGLM Overall accuracy

# function for gw overall accuracy
gw_overall_accuracy = function(pred, obs, X , Y, grid,
                         # GW parameters
                         bw_defined = FALSE, bw = 0.3, 
                         kernel = gwr.bisquare,
                         verbose = FALSE)
{
  # define logit function
  alogit <- function(x){exp(x)/(1+exp(x))}
  # create a binary variable: Observed matches Predicted 
  # res = (pred == obs) + 0
  # define the coordinates matrix and data.frame for the data
  coords = cbind(X,Y)
  df = data.frame(pred, obs)
  # define a class list
  class.list <- unique(obs)[order(unique(obs))]
  # optimal kernel bandwidth
  if(!bw_defined) {
    NN <- ggwr.sel(res~1, data = df, coords = coords, adapt = TRUE, gweight = kernel,
                   family=binomial, verbose = verbose)
  } else NN = bw	
  gwr.model_overall <- ggwr(res~1, data = df, coords = coords, adapt = NN, 
                            gweight = kernel,
                            fit.points = grid, family= binomial) 
  # extract the probabilities from the SDF of the model 
  gwr.ovac <- alogit(data.frame(gwr.model_overall$SDF)[,2])
  return(list(Bandwidth = NN, GWOverall = gwr.ovac))
}
gw_oa = gw_overall_accuracy(pred = data$Boolean_Pred, obs = data$Boolean_Obs, 
                    X = data[,2], Y = data[,3], grid = grid.xy)
# gw_oa$Bandwidth
# summary(gw_oa$GWOverall)

# function for gw user's and producer'saccuracy
gw_userprod_accuracy = function(pred, obs, X, Y, grid,
                         # GW parameters
                         bw_defined = FALSE, bw = 0.3, 
                         kernel = gwr.bisquare,
                         verbose = FALSE)
{
  # define logit function
  alogit <- function(x){exp(x)/(1+exp(x))}
  # define the coordinates matrix
  coords = cbind(X,Y)
  # define a class list and some empty vectors for results
  class.list <- unique(obs)[order(unique(obs))]
  user_result = prod_result = vector()
  user_bw = prod_bw = vector()
  # then loop through each class to calculate User and Producer measures  
  for (i in 1: length(class.list)) {
    class <- class.list[i]	
    # define binary variables where the class is present
    obs.class <- (obs == class) * 1 	
    pred.class <- (pred == class) * 1	
    # join together for use in the GLMs
    obs_pred <- data.frame(cbind(obs.class,pred.class)) 
    # define the data.frame for the data
    df = data.frame(pred.class, obs.class)
    # User's Accuracy
    if(!bw_defined) {
      NN <- ggwr.sel(obs.class~pred.class, data = df, coords = coords, adapt = TRUE, 
                     gweight = kernel, family=binomial, verbose = verbose)
    } else NN = bw
    gwr.model_user <- ggwr(obs.class~pred.class, data = df, coords = coords, adapt = NN,
                           fit.points=grid, gweight = kernel, family= binomial) 
    coefs <- data.frame(gwr.model_user$SDF)[,2:3]
    coefs[,2] <- rowSums(coefs) #logit ea+c
    p.x.eq.1.g.y.eq.1 <- alogit(coefs[,2]) #n1 (sum coeffs)
    gwr.user <-  p.x.eq.1.g.y.eq.1 
    user_result = cbind(user_result, gwr.user)
    user_bw = append(user_bw, NN)
    # Producers Accuracy 
    if(!bw_defined) {
      NN <- ggwr.sel(pred.class~obs.class, data = df, coords = coords, adapt = TRUE, 
                     gweight = kernel, family=binomial, verbose = verbose)
    } else NN = bw
    gwr.model_user <- ggwr(pred.class~obs.class, data = df, coords = coords, adapt = NN,
                           fit.points=grid, gweight = kernel, family= binomial) 
    coefs <- data.frame(gwr.model_user$SDF)[,2:3]
    coefs[,2] <- rowSums(coefs) #logit ea+c
    p.x.eq.1.g.y.eq.1 <- alogit(coefs[,2]) #n1 (sum coeffs)
    gwr.prod <- p.x.eq.1.g.y.eq.1
    prod_result = cbind(prod_result, gwr.prod)
    prod_bw = append(prod_bw, NN)
    # progress!
    if(verbose) cat(i, "\t")
  }
  colnames(user_result) = paste0("U_",class.list) 
  colnames(prod_result) = paste0("P_",class.list)
  names(user_bw) = class.list
  names(prod_bw) =  class.list
  return(list(UserBW = user_bw, ProdBW = prod_bw,
              GWUser = user_result,
              GWProd = prod_result))
}
gw_up = gw_userprod_accuracy(pred = data$Boolean_Pred, obs = data$Boolean_Obs, 
                         X = data[,2], Y = data[,3], grid = grid.xy)
#gw_up$UserBW
#gw_up$ProdBW
#t(round(apply(gw_up$GWUser, 2, summary),3))
#t(round(apply(gw_up$GWProd, 2, summary),3))


## Table 2: GW overall accuracy
tab2 = c(summary(gw_oa$GWOverall), Bandwidth = gw_oa$Bandwidth, Global = oa)
tab2 = t(as.matrix(tab2))
tab2

## Table 3: GW user's accuracy
tab3 = t(round(apply(gw_up$GWUser, 2, summary),3))
rownames(tab3) <- class.names
tab3 = cbind(tab3, Bandwidths = gw_up$UserBW, Global = users)
tab3

## Table 4: GW producer's accuracy
tab4 = t(round(apply(gw_up$GWProd, 2, summary),3))
rownames(tab4) <- class.names
tab4 = cbind(tab4, Bandwidths = gw_up$ProdBW, Global = producers)
tab4

## Fig4: Spatially distributed user's and producer's accuracies for Grazing Land from a GWGLM.
allgwglm_data = data.frame(grid.xy, overall = gw_oa$GWOverall, 
                           gw_up$GWUser, gw_up$GWProd)
# make spatial if you want! 
# gwr_sf = st_as_sf(gwr_all_part1, coords = c("X", "Y"))
# add geographic projection  
# st_crs(gwr_sf) = 32633
# map with the df object for raster 
p1 = 
  ggplot() + 
		geom_raster(data = allgwglm_data, aes(x = X, y = Y, fill = U_G)) +
		scale_fill_continuous_c4a_seq(palette="kovesi.linear_wh_rd_bk", name = "User's") +
    geom_sf(data = mask, col = NA, fill = "white") +
    theme_bw() +
    coord_sf() +
    theme(legend.position = "bottom",
          axis.title=element_blank(),
	        axis.text=element_blank(),
	        axis.ticks=element_blank())
p2 = 
  ggplot() + 
		geom_raster(data = allgwglm_data, aes(x = X, y = Y, fill = P_G)) +
		scale_fill_continuous_c4a_seq(palette="kovesi.linear_wh_rd_bk", name = "Producer's") +
    geom_sf(data = mask, col = NA, fill = "white") +
    theme_bw() +
    coord_sf() +
    theme(legend.position = "bottom",
          axis.title=element_blank(),
	        axis.text=element_blank(),
	        axis.ticks=element_blank())
plot_grid(p1, p2, nrow = 1)


### 3.2 GW confusion matrices

## Fig5 Two example locations for which geographically weighted grids are determined.
## create the `sp` objects
# convert data - sp format
data.spdf <- SpatialPointsDataFrame(coords = data[,2:3], data = data) 
# convert the outline and the interpolation grid to sf format 
lib.spdf = as(lib, "Spatial")
grid.spdf = as(grid, "Spatial")
# define example points
x <- c(454, 1900 )
x <- st_as_sf(grid.spdf[x,])
x$ID = 1:2
# make the plot
p1 = ggplot() +
  geom_sf(data = st_as_sf(grid.spdf),pch=16,col="#A50F1580",cex=1) +
  geom_sf(data = mask, col = NA, fill = "white") +
  geom_sf(data = x, cex = 2) +
  geom_sf_text(data = x, aes(label = ID), cex = 8, nudge_y = -3000, nudge_x = 2000)+
  theme_bw() +
  coord_sf() +
  theme(axis.title=element_blank(),
	      axis.text=element_blank(),
	      axis.ticks=element_blank())
p1

# set up the gwxt framework with the data
dummy_xtab <- new_spxt(data.spdf,'Boolean_Pred', 'Boolean_Obs')
# define a bandwidth
bw = round(nrow(data.spdf)*0.3, 0)
# get the gwxt function
gwxt <- gwxtab_probe(dummy_xtab,adapt(bw))
# apply to locations

# Table 5: point 1
coords.i = as.vector(st_coordinates(x[1,]))
gwcm1 <- round(gwxt(x = coords.i[1], y = coords.i[2] ), 2)
rownames(gwcm1) <- class.names
gwcm1
# Table 5: point 2
coords.i = as.vector(st_coordinates(x[2,]))
gwcm2 <- round(gwxt(x = coords.i[1], y = coords.i[2] ), 2)
rownames(gwcm2) <- class.names
gwcm2

## Table 7: overall and kappa accuracies 
# kappa
bw = round(nrow(data.spdf)*gw_oa$Bandwidth, 0)
kappa.spdf <- gwxtab_sample(grid.spdf,dummy_xtab,adapt(bw),melt=kp)
# overall 
overall.spdf <- gwxtab_sample(grid.spdf,dummy_xtab,adapt(bw),melt=overall.acc)
names(overall.spdf) = "overall"
# tidy the kappa
negkap2min <- function (x) {
  index <- x < 0
  x[index] <- min(x[!index])
  return(x)}
kappa.spdf$kappa <- negkap2min(kappa.spdf$kappa)
# convert to sf and map
# re-create global CM

tab7 = rbind(summary(overall.spdf$overall), summary(kappa.spdf$kappa))
tab7 = cbind(tab7, Global = c(oa, ka))
tab7

## Table 8 and 9: GWCM user's and producer's
# set up bws and empty result vectors
bws.u = gw_up$UserBW
bws.p = gw_up$ProdBW
user.tab = producer.tab = matrix(nrow = length(grid.spdf), ncol = 5)
for (i in 1:5){
	# define functions with a different index each time
  user.acc <- function(ctmatrix, index = i) {
		ctmatrix[index,index] / sum(ctmatrix[index,])
	}  
	prod.acc <- function(ctmatrix, index = i){
		ctmatrix[index,index] / sum(ctmatrix[, index])
	}
	# calculate user's and producer's accuracies 
	bw = round(nrow(data.spdf)*bws.u[i], 0)
	if(bw == 210) bw = 209
	user.i <- as.vector(unlist(gwxtab_sample(grid.spdf,dummy_xtab,adapt(bw),
	                                  melt=user.acc)@data))
	bw = round(nrow(data.spdf)*bws.p[i], 0)
	if(bw == 210) bw = 209
	producer.i <- as.vector(unlist(gwxtab_sample(grid.spdf,dummy_xtab,adapt(bw),
	                                      melt=prod.acc)@data))
  user.tab[,i] = user.i
  producer.tab[,i] = producer.i
}
colnames(user.tab) = paste0("U_", class.list)
colnames(producer.tab) = paste0("P_",class.list)
# generate tables
tab8 = t(apply(user.tab, 2, summary))
rownames(tab8) <- class.names
tab8 = cbind(tab8, Bandwidths = gw_up$UserBW, Global = users)
tab8

tab9 = t(apply(producer.tab, 2, summary))
rownames(tab9) <- class.names
tab9 = cbind(tab9, Bandwidths = gw_up$ProdBW, Global = producers)
tab9

## Fig6: Two examples of kappa estimates from GW confusion matrices with different bandwidths.
# create a different kappa
bw = round(nrow(data.spdf)*0.3, 0)
kappa2.spdf <- gwxtab_sample(grid.spdf,dummy_xtab,adapt(bw),melt=kp)
# tidy the kappa
kappa2.spdf$kappa <- negkap2min(kappa2.spdf$kappa)
df = data.frame(grid.xy, kappa1 = kappa.spdf$kappa, kappa2 = kappa2.spdf$kappa)
# head(df)
p1 = 
  ggplot() + 
		geom_raster(data = df, aes(x = X, y = Y, fill = kappa1)) +
		scale_fill_continuous_c4a_seq(palette="kovesi.linear_wh_rd_bk", name = "Bandwidth = 0.82") +
    geom_sf(data = mask, col = NA, fill = "white") +
    theme_bw() +
    coord_sf() +
    theme(legend.position = "bottom",
          axis.title=element_blank(),
	        axis.text=element_blank(),
	        axis.ticks=element_blank())
p2 = 
  ggplot() + 
		geom_raster(data = df, aes(x = X, y = Y, fill = kappa2)) +
		scale_fill_continuous_c4a_seq(palette="kovesi.linear_wh_rd_bk", name = "Bandwidth = 0.3") +
    geom_sf(data = mask, col = NA, fill = "white") +
    theme_bw() +
    coord_sf() +
    theme(legend.position = "bottom",
          axis.title=element_blank(),
	        axis.text=element_blank(),
	        axis.ticks=element_blank())
plot_grid(p1, p2, ncol = 2)

### 4. Working with soft classes

## 4.1 GW fuzzy certainty I

# GW Analysis of Fuzzy Difference
gw_fuzzy_accuracy = function(pred, obs, X, Y, grid,
                         # GW parameters
                         bw_defined = F, bw = 0.3, kernel = gwr.Gauss,
                         verbose = F)
{
  # define the coordinates matrix
  coords = cbind(X,Y)
	d <- abs(pred-obs) 		#absolute difference between reference and classified
	df = data.frame(d = d)
	if(!bw_defined) {
      NN <- gwr.sel(d~1, data = df, coords = coords, adapt = TRUE, 
                     gweight = kernel, verbose = verbose)
    } else NN = bw
	# GW difference
	gwr.model<- gwr(d~1, data = df, coords = coords, adapt = NN, fit.points = grid)
	# extract result NOTE that the result is 1-intercept 
	gwr.res = 1- as.vector(unlist(data.frame(gwr.model$SDF$`(Intercept)`)))
	return(list(BW = NN, FuzzyAcc = gwr.res))
}
# define X and Y
X = data[,2]
Y = data[,3]
# gw fuzzy accuracy
urban_gw_acc =      gw_fuzzy_accuracy(pred = data$Urban_Pred, 
                                      obs = data$Urban_Obs,
                                      X, Y, grid.xy) 
vegetation_gw_acc = gw_fuzzy_accuracy(pred = data$Vegetation_Pred, 
                                      obs = data$Vegetation_Obs, 
                                      X, Y, grid.xy) 
woody_gw_acc =      gw_fuzzy_accuracy(pred = data$Woody_Pred, 
                                      obs = data$Woody_Obs,
                                      X, Y, grid.xy) 
grazing_gw_acc =    gw_fuzzy_accuracy(pred = data$Grazing_Pred, 
                                      obs = data$Grazing_Obs,
                                      X, Y, grid.xy) 
bare_gw_acc =       gw_fuzzy_accuracy(pred = data$Urban_Pred,
                                      obs = data$Urban_Obs,
                                      X, Y, grid.xy) 
## Table 10: summarise the results
tab10 = rbind(summary(bare_gw_acc$FuzzyAcc),
              summary(grazing_gw_acc$FuzzyAcc),
              summary(urban_gw_acc$FuzzyAcc),
              summary(vegetation_gw_acc$FuzzyAcc),
              summary(woody_gw_acc$FuzzyAcc))
bws = c(bare_gw_acc$BW, grazing_gw_acc$BW, urban_gw_acc$BW, 
        vegetation_gw_acc$BW, woody_gw_acc$BW)
tab10 = cbind(tab10, Bandwidths = bws)
rownames(tab10) = class.names
tab10

## Fig7: The fuzzy correspondences of the Grazing and Vegetation classes computed from fuzzy difference
df = data.frame(grazing = grazing_gw_acc$FuzzyAcc, 
                vegetation = vegetation_gw_acc$FuzzyAcc, 
                grid.xy)
# map with the df object for raster 
p1 = 
  ggplot() + 
		geom_raster(data = df, aes(x = X, y = Y, fill = grazing)) +
		scale_fill_continuous_c4a_seq(palette="kovesi.linear_wh_yl_gn_bu", 
		                              name = "Grazing",
		                              breaks = pretty_breaks(n = 3)) +
    geom_sf(data = mask, col = NA, fill = "white") +
    theme_bw() +
    coord_sf() +
    theme(legend.position = "bottom",
          axis.title=element_blank(),
	        axis.text=element_blank(),
	        axis.ticks=element_blank())
p2 = 
  ggplot() + 
		geom_raster(data = df, aes(x = X, y = Y, fill = vegetation)) +
		scale_fill_continuous_c4a_seq(palette="kovesi.linear_wh_yl_gn_bu", 
		                              name = "Vegetation",
		                              breaks = pretty_breaks(n = 3)) +
    geom_sf(data = mask, col = NA, fill = "white") +    
    theme_bw() +
    coord_sf() +
    theme(legend.position = "bottom",
          axis.title=element_blank(),
	        axis.text=element_blank(),
	        axis.ticks=element_blank())
plot_grid(p1, p2, ncol = 2)

## 4.2 GW fuzzy certainty II

# fuzzy change function
fuzzy_change_analysis = function(predicted = pred_df, observed = obs_df) {
  c.ns = sort(names(observed))
  pred_df = predicted[,c.ns]
  obs_df = observed[, c.ns]
  # get the "not" fuzzy values
  notpred = 1- pred_df
  notobs  = 1- obs_df
  # set up some output objects
  f.cm = f.cm2 = matrix(0,nrow = 5, ncol = 5) 
  gb_vec = lb_vec = same_vec = NULL
  for (i in 1:length(c.ns)){
    # define row class (predicted!)
    class.i = c.ns[i]
    # empt vectors for gain and loss boundaries
    gb.vec = lb.vec = NULL
    for(j in 1:length(c.ns)){
      class.j = c.ns[j]
      obs.j = as.vector(obs_df[, class.j])
      pred.j = as.vector(pred_df[, class.j])
      if (i == j) { # minimum interval
        same.i = sapply(1:nrow(pred_df), function(x) min(pred.j[x], obs.j[x]))
        # sum and add to fcm
        f.cm[i, j] = sum(same.i)
        # bind to output
        same_vec  = cbind(same_vec, same.i) 
        tit = paste0(class.i, "_", class.j)
        colnames(same_vec)[ncol(same_vec)] <- tit
      }  
      if(i < j ) { # fuzzy Type II errors: 
        # how column wise observed class different to predicted
        # parallel to fuzzy gain
        gb = sapply(1:nrow(pred_df), 
                    function(x) max(0, notpred[x, class.i] + 
                                      min(obs_df[x, class.i], obs.j[x]) -1 ))
        # sum and add to fcm
        f.cm[i,j] = sum(gb)
        # bind to output
        gb.vec = cbind(gb.vec, gb)
        tit = paste0("row_", class.i, "_col_", class.j)
        colnames(gb.vec)[ncol(gb.vec)] <- tit
      }
      if(i > j ) { # fuzzy Type II errors: 
        # how row wise observed is different into predicted
        # parallel to fuzzy loss
        lb = sapply(1:nrow(pred_df), 
                    function(x) max(0, min(pred_df[x, class.i], pred.j[x]) + 
                                      notobs[x, class.i] -1))
        # sum and add to fcm
        f.cm[i,j] = sum(lb)
        # bind to output
        lb.vec = cbind(lb.vec, lb)
        tit = paste0("row_", class.i, "_col_", class.j)
        colnames(lb.vec)[ncol(lb.vec)] <- tit
      }  
    }
    # bind the results together
    if(length(gb.vec) > 0) gb_vec = cbind(gb_vec, gb.vec)  
    if(length(lb.vec) > 0) lb_vec = cbind(lb_vec, lb.vec)  
  }
  # do commission (losses) and omission (gains)
  fom = sapply(1:length(c.ns), function(x) sum(f.cm[-x, x]) )
  fcom = sapply(1:length(c.ns), function(x) sum(f.cm[x, -x]) )
  Fom_Fcom_tab = cbind(fom, fcom)
  rownames(Fom_Fcom_tab) = c.ns
  rownames(f.cm) = paste0("Pred_",c.ns)
  colnames(f.cm) = paste0("Obs_",c.ns)
  #return(list(FCM = f.cm, FUnc_tab = FUnc_tab , 
  #            FType1ConfG = gb_vec, FType2ConfL = lb_vec, FCorresp = same_vec))
  return(list(FCM = f.cm, Fom_Fcom_tab = Fom_Fcom_tab , 
              F_om = gb_vec, F_com = lb_vec, F_corr = same_vec))
}
## Fuzzy data set up - create observed and predicted data
obs_df = data %>% select(ends_with("_Obs"),-Boolean_Obs)
names(obs_df) = gsub("_Obs", "", names(obs_df))
pred_df = data %>% select(ends_with("_Pred"),-Boolean_Pred)
names(pred_df) = gsub("_Pred", "", names(pred_df))
# do some sorting ordering or variables
c.ns = sort(names(obs_df))
pred_df = pred_df[,c.ns]
obs_df = obs_df[, c.ns]

## Table 10: the fuzzy CM
fuzzy_res = fuzzy_change_analysis(predicted = pred_df, observed = obs_df)
fcm = fuzzy_res$FCM
colnames(fcm) = class.list
rownames(fcm) = class.names
fcm

## Table 11: the fuzzy Omission / Commission
colnames(fuzzy_res$Fom_Fcom_tab) = c("Omission", "Commission")
fuzzy_res$Fom_Fcom_tab

# GW Fuzzy Uncertainty
gw_fuzzy_uncertainty = function(class.name = "Grazing",fuzzy_res, X, Y, grid,
                         # GW parameters
                         bw_defined = FALSE, bw = 0.3, kernel = gwr.bisquare,
                         verbose = FALSE)
{
  # get the data
  # fuzzy correspondence
  fuzzy_res$F_cor %>% 
    data.frame() %>% 
    select(contains(class.name)) -> df1
  fcorres_df = data.frame(X, Y, val = as.vector(unlist(df1)))
  
  # fuzzy omission
  col_tit = paste0("col_", class.name)
  fuzzy_res$F_om %>% 
    data.frame() %>%
    select(contains(col_tit)) -> df1
  fuzzy_res$F_com %>% 
    data.frame() %>%
    select(contains(col_tit)) -> df2
  fom_df = data.frame(X, Y, val = rowSums(df1) + rowSums(df2))
  
  # fuzzy commission
  row_tit = paste0("row_", class.name)
  fuzzy_res$F_com %>% 
    data.frame() %>%
    select(contains(row_tit)) -> df1
  fuzzy_res$F_com %>% 
    data.frame() %>%
    select(contains(row_tit)) -> df2
  fcom_df = data.frame(X, Y, val = rowSums(df1) + rowSums(df2))
  # define the coordinates matrix
  coords = cbind(X,Y)
  
  # undertake GW bandwidth
  if(!bw_defined) {
      NN1 <- gwr.sel(val~1, data = fcorres_df, coords = coords, adapt = TRUE, 
                     gweight = kernel, verbose = verbose)
      NN2 <- gwr.sel(val~1, data = fom_df, coords = coords, adapt = TRUE, 
                     gweight = kernel, verbose = verbose)
      NN3 <- gwr.sel(val~1, data = fcom_df, coords = coords, adapt = TRUE, 
                     gweight = kernel, verbose = verbose)
    } else NN1 = NN2 = NN3 = bw
	# GW interpolation
	gwr.model<- gwr(val~1, data = fcorres_df, coords = coords, adapt = NN1, fit.points = grid)
	# extract result NOTE that the result is 1-intercept in each case
	f_corr = 1- as.vector(unlist(data.frame(gwr.model$SDF$`(Intercept)`)))
  # and do for others
	gwr.model<- gwr(val~1, data = fom_df, coords = coords, adapt = NN2, fit.points = grid)
	f_prod = 1- as.vector(unlist(data.frame(gwr.model$SDF$`(Intercept)`)))
	gwr.model<- gwr(val~1, data = fcom_df, coords = coords, adapt = NN3, fit.points = grid)
	f_user = 1-as.vector(unlist(data.frame(gwr.model$SDF$`(Intercept)`)))
	
	return(list(BWs = c(bw_overall = NN1, bw_prod = NN2, bw_user = NN3), 
	            fuzzy_res = cbind(f_corr, f_prod, f_user)))
}

# apply the function
fuzzy_veg = gw_fuzzy_uncertainty("Grazing",fuzzy_res, X, Y, grid.xy)
# have a look at the result
# fuzzy_veg$BWs
# head(fuzzy_veg$fuzzy_res)
# combine with grid coordinates
df = data.frame(fuzzy_veg$fuzzy_res, grid.xy)

## Fig8 The fuzzy correspondence, user's and producers' certainties for the Grazing class, computed from bounded difference
# map with the df object for raster 
tit = paste0("Correspondence, BW: ", round(fuzzy_veg$BWs[1], 2))
p1 = 
  ggplot() + 
		geom_raster(data = df, aes(x = X, y = Y, fill = f_corr)) +
		scale_fill_continuous_c4a_seq(palette="kovesi.linear_wh_yl_gn_bu", 
		                              name = "",
		                              breaks = pretty_breaks(n = 2)) +
    geom_sf(data = mask, col = NA, fill = "white") +
    labs(title = "", subtitle = tit) +
    theme_bw() +
    coord_sf() +
    theme(legend.position = "bottom",
          axis.title=element_blank(),
	        axis.text=element_blank(),
	        axis.ticks=element_blank())
tit = paste0("Producer's, BW: ", round(fuzzy_veg$BWs[2], 2))
p2 = 
  ggplot() + 
		geom_raster(data = df, aes(x = X, y = Y, fill = f_prod)) +
		scale_fill_continuous_c4a_seq(palette="kovesi.linear_wh_yl_gn_bu", 
		                              name = "",
		                              breaks = pretty_breaks(n = 2)) +
    geom_sf(data = mask, col = NA, fill = "white") +    
    labs(title = "", subtitle = tit) +
    theme_bw() +
    coord_sf() +
    theme(legend.position = "bottom",
          axis.title=element_blank(),
	        axis.text=element_blank(),
	        axis.ticks=element_blank())
tit = paste0("User's, BW: ", round(fuzzy_veg$BWs[3], 2))
p3 = 
  ggplot() + 
		geom_raster(data = df, aes(x = X, y = Y, fill = f_user)) +
		scale_fill_continuous_c4a_seq(palette="kovesi.linear_wh_yl_gn_bu", 
		                              name = "",
		                              breaks = pretty_breaks(n = 2)) +
    geom_sf(data = mask, col = NA, fill = "white") +
    labs(title = "", subtitle = tit) +
    theme_bw() +
    coord_sf() +
    theme(legend.position = "bottom",
          axis.title=element_blank(),
	        axis.text=element_blank(),
	        axis.ticks=element_blank())
plot_grid(p1, p2, p3, ncol = 3)

## 4.3 GW compostional error
# define data
Pred = data[,11:15]
Obs = data[,4:8]
# apply aDist function
oAD = aDist(Pred+0.001, Obs+0.001)
AS.all = sapply(1:nrow(data), 
                   function(x) aDist(Pred[x,]+0.001,Obs[x,]+0.001 ))

## Table 13: Aitchison distance
tab13 = c(summary(AS.all), Overall = oAD)
tab13

## Fig9: The Geographically Weighted Aitchison distance (I)
# create data
X = data[,2]
Y = data[,3]
coords = cbind(X,Y)
df = data.frame(ADist = AS.all, coords)
# GW Aich dist
NN <- gwr.sel(ADist~1, data = df, coords = coords, adapt = TRUE, 
              gweight = gwr.bisquare, verbose = F)
gwr.model<- gwr(ADist~1, data = df, coords = coords, adapt = NN,  
                gweight = gwr.bisquare, fit.points = grid.xy)
gwr.res = as.vector(unlist(data.frame(gwr.model$SDF$`(Intercept)`)))
# create data frame and map in the usual way
df = data.frame(AichDist = gwr.res, X = grid.xy[,"X"], Y = grid.xy[,"Y"])
tit = paste0("BW: ", round(NN, 3))
ggplot() + 
  geom_raster(data = df, aes(x = X, y = Y, fill = AichDist)) +
  scale_fill_continuous_c4a_seq(palette="kovesi.linear_wh_yl_gn_bu", 
                                name = "Aitchison Distance",
                                breaks = trans_breaks(identity, identity, n = 4)) +
  geom_sf(data = mask, col = NA, fill = "white") +
  theme_bw() +
  labs(title = "", subtitle = tit) +
  coord_sf() +
  theme(legend.position = "bottom",
        axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank())

# constructing a Geographically Weighted Aitchison distance by hand
# extract location coordinates
X = data[,2]
Y = data[,3]
coords = cbind(X,Y)
# create distance matrix
dist_mat = as.matrix(dist(coords, upper = T, diag = T))
# examine the first few row / columns
# dist_mat[1:6, 1:6]

# function to get index of nearby observations from a distance matrix
gw_get_nearby = function(obs_index, adaptive, bw, dist_mat){
  if( adaptive) {
    return(as.vector(sort(order(dist_mat[obs_index, ])[1:bw])))
  } 
  if(!adaptive) {
    return(as.vector(which(dist_mat[obs_index,] < bw)))
  } 
}
# specify an adaptive bandwidth and convert to integer  
bwa = NN
bw = round(nrow(dist_mat)*bwa,0)
# gw_get_nearby(1, adaptive = T, bw, dist_mat)
# for fixed distance of 10km
# gw_get_nearby(1, adaptive = F, bw = 10000, dist_mat)
get_nearby_weights = function(obs_index, nearby_locs, bw, 
                              dist_mat, kernel = "bisquare", adaptive = T){
  dists = as.vector(dist_mat[obs_index, nearby_locs])
  if(adaptive) bw = max(dists)
  if(kernel == "bisquare") w = (1-(dists/bw)^2)^2 
  if(kernel == "gaussian") w = exp(-.5*(dists/bw)^2)
  if(kernel == "exponential") w = exp(-dists/bw)
  if(kernel == "tricube") w = (1-(dists/bw)^3)^3
  if(kernel == "boxcar") w = rep(1, length(dists))
  w.vec = rep(0, nrow(dist_mat))
  w.vec[nearby_locs] = w
  return(w.vec)
}
# get nearby obs
nearby_index = gw_get_nearby(1, adaptive = T, bw, dist_mat)
# weight observations
wts = get_nearby_weights(1, nearby_index, bw, dist_mat, 
                         kernel = "bisquare", adaptive = T)
# examine the weights
# round(wts, 3)
wts.all <- NULL
# create weights matrix
for(i in 1:nrow(dist_mat)) {
  nearby_index = gw_get_nearby(i, adaptive = T, bw, dist_mat)
  wts.i = get_nearby_weights(i, nearby_index, bw, dist_mat, kernel = "gaussian", adaptive = T)
  wts.all = cbind(wts.all, as.vector(wts.i))
}
# dim(wts.all)
# wts.all[1:6, 1:6]
# define empty matrix
res_ad <-NULL
for(i in 1:nrow(dist_mat)){
  # weight the data 
  Pred.i = Pred * wts.all[i,] 
  Obs.i = Obs * wts.all[i,]   
  # remove observations outside of the bandwidth
  index = rowSums(Pred.i) == 0 | rowSums(Obs.i) == 0
  Pred.i = Pred.i[!index, ]+0.001
  Obs.i = Obs.i[!index, ]+0.001
  ad.i = aDist(Pred.i, Obs.i)
  res_ad = c(res_ad, ad.i)
}
# length(res_ad)
# summary(res_ad)
# now interpolate over the GW grid
df = data.frame(res_ad)
# GW interpolation over grid
gwr.model<- gwr(res_ad~1, data = df, coords = coords, adapt = NN, fit.points = grid.xy)
gwr.res = as.vector(unlist(data.frame(gwr.model$SDF$`(Intercept)`)))

## Fig10. The Geographically Weighted Aitchison distance (II)
# create the data frame for mapping
df = data.frame(AichDist = gwr.res, X = grid.xy[,"X"], Y = grid.xy[,"Y"])
# dim(df)
ggplot() + 
  geom_raster(data = df, aes(x = X, y = Y, fill = AichDist)) +
  scale_fill_continuous_c4a_seq(palette="kovesi.linear_wh_yl_gn_bu", 
                                name = "Aitchison Distance",
                                breaks = scales::trans_breaks(identity, identity, n = 4)) +
  geom_sf(data = mask, col = NA, fill = "white") +
  theme_bw() +
  labs(title = "", subtitle = tit) +
  coord_sf() +
  theme(legend.position = "bottom",
        axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank())

#### END #### 