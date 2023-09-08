## Title: Geographically weighted accuracy for hard and soft land cover classifications: 5 approaches with coded illustrations
# Alexis Comber and Naru Tsutsumida
# If you have any problems with data / code / versions etc please contact Lex Comber a.comber@leeds.ac.uk


## packages and data prep
library(sp)
library(spgwr)
library(ggspatial)
library(cowplot)
library(cols4all)
library(scales)
library(tidyverse)
library(sf)
library(robCompositions)
library(diffeR)
library(knitr)
library(gwxtab)


# load the data
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
# define a masking function
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
ggplot() + 
 geom_sf(data = lib, col = "black", fill = "NA") +
	geom_sf(data = grid, pch = 3, col = "lightgrey") +
	geom_sf(data = data_sf, col = "black") +
	geom_sf(data = Tripoli, pch = 20, cex = 6, col = "red") +
	theme_bw() + coord_sf() +
	#annotation_north_arrow(pad_x = unit(13.4, "cm"), pad_y = unit(6.3, "cm")) +
  annotation_scale(pad_x = unit(8.1, "cm"), pad_y = unit(5, "cm")) 


## Table 1
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
class.names <- c("1. Bare Ground",
                 "2. Grazing Land", 
                 "3. Urban", 
                 "4. Vegetation", 
                 "5. Woodland") 
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
knitr::kable(tab1, 
               caption = paste0("\\label{tab:tab1}The correspondence matrix between predicted (rows) and observed (columns) land cover class, with user's and producer's accuracies. Overall accuracy: ",round(oa,3), "; Kappa estimate: ", round(ka,3), "."))


## Figure 2
# define a vector of names
class.names2 = c("Bare Ground", "Grazing Land", "Urban", "Vegetation", "Woodland")
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
		scale_size(breaks = c(0.1, 0.25, 0.5, 0.75, 0.9, 1), range = c(0,5), limits = c(0,1)) +
	  ggtitle(label = NULL,subtitle = tit) +
	  theme_bw() +
	  theme(axis.title=element_blank(),
	        axis.text=element_blank(),
          axis.ticks=element_blank(), 
         	legend.title = element_blank(),
          legend.text = element_text(size=8))
       
}
pfUo = plot_fuzzy(data[,4], i=3) 
# get the legend and then remove 
legFO = get_legend(pfUo)
pfUo = pfUo + theme(legend.position = "none")
pfVo = plot_fuzzy(data[,5], i=4) + theme(legend.position = "none")
pfWo = plot_fuzzy(data[,6], i=5) + theme(legend.position = "none")
pfGo = plot_fuzzy(data[,7], i=2) + theme(legend.position = "none")
pfBo = plot_fuzzy(data[,8], i=1) + theme(legend.position = "none")
cowplot::plot_grid(pfUo, pfVo, pfWo, pfGo, pfBo, legFO, nrow = 2, ncol = 3)


## Figure 3
pfUb = plot_fuzzy(data[,11], i=3, Obs = F) 
# get the legend and then remove 
legFB = get_legend(pfUb)
pfUb = pfUb + theme(legend.position = "none")
pfVb = plot_fuzzy(data[,12], i=4, Obs = F) + theme(legend.position = "none")
pfWb = plot_fuzzy(data[,13], i=5, Obs = F) + theme(legend.position = "none")
pfGb = plot_fuzzy(data[,14], i=2, Obs = F) + theme(legend.position = "none")
pfBb = plot_fuzzy(data[,15], i=1, Obs = F) + theme(legend.position = "none")
cowplot::plot_grid(pfUb, pfVb, pfWb, pfGb, pfBb, legFB, nrow = 2, ncol = 3)


## Overall Accuracy
# create 1/0 vector of True/False values
vec = (data$Boolean_Pred == data$Boolean_Obs) + 0
# create a binomial model
mod <- glm(vec~1,family= binomial) 
# extract the coefficient estimates
mod.coefs <- as.vector(mod$coefficients)
OA = exp(mod.coefs)/(1+exp(mod.coefs))
# do this with a logit function
alogit <- function(x){exp(x)/(1+exp(x))}
OA = alogit(mod.coefs) 


## Users and Producers Accuracy
class.list <- unique(data$Boolean_Obs)[order(unique(data$Boolean_Obs))]
# Grazing land
i = 2
class <- class.list[i]	
# define binary variables where the class is present
obs.class <- (data$Boolean_Obs == class) * 1 	
pred.class <- (data$Boolean_Pred == class) * 1	
# join together for use in the GLMs
obs_pred <- data.frame(cbind(obs.class,pred.class)) 
# Users Accuracy
mod1 <- glm(obs.class~pred.class,data = obs_pred,family= binomial) 
mod.coefs <- mod1$coefficients
# sum the coefficients
mod.coefs.sum <-sum(mod.coefs) 
# apply the logit function
mod.p.x.eq.1.g.y.eq.1 <- alogit(mod.coefs.sum)
mod.user <- mod.p.x.eq.1.g.y.eq.1 
# Producers - invert terms and repeat
mod2 <- glm(pred.class~obs.class,data = obs_pred,family= binomial) 
mod.coefs <- mod2$coefficients
mod.coefs.sum <-sum(mod.coefs) 
mod.p.x.eq.1.g.y.eq.1 <- alogit(mod.coefs.sum) #n1
mod.prod <- mod.p.x.eq.1.g.y.eq.1 


## A function using GLMs to return OA and all UA / PA values
do_glm_accuracy = function(pred, obs) {
  # 1. create the Overall Accuracy model
  #  binary variable: Observed matches Predicted 
  res = (pred == obs) + 0
  mod0 <- glm(res~1,family= binomial) 
  # extract the model coefficients
  mod.coefs <- mod0$coefficients
  mod.coefs.sum <-sum(mod.coefs) 
  # model the probabilities that the 
  # probability of x equals 1 given that y equals 1
  mod.p.x.eq.1.g.y.eq.1 <- alogit(mod.coefs.sum) 
  # this is the overall accuracy
  mod.OA <- mod.p.x.eq.1.g.y.eq.1 
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
    mod.coefs.sum <-sum(mod.coefs) 
    mod.p.x.eq.1.g.y.eq.1 <- alogit(mod.coefs.sum)
    mod.user <- mod.p.x.eq.1.g.y.eq.1 
    # add to the result
    user_result = append(user_result, mod.user)
    # Producer - invert terms
    mod2 <- glm(pred.class~obs.class,data = obs_pred,family= binomial) 
    mod.coefs <- mod2$coefficients
    mod.coefs.sum <-sum(mod.coefs) 
    mod.p.x.eq.1.g.y.eq.1 <- alogit(mod.coefs.sum) 
    mod.prod <- mod.p.x.eq.1.g.y.eq.1 
    # add to the result
    prod_result = append(prod_result, mod.prod)
  }
  # combine the outputs and return
  out = data.frame(Producer = prod_result, User = user_result)
  rownames(out) = class.list
  return(list(Overall = mod.OA, ProdUser = out))
}

# The function can applied to the data as below, resulting in just 2 vectors of the same length of the predicted and observed classes. It returns a list a object with 2 elements,: the overall accuracy measure and a table of the per class user's and producer's accuracies.
# uncomment to examine!
# do_glm_accuracy(data$Boolean_Pred, data$Boolean_Obs)


## Functions to calculate GW OA, UA and PA
# OA
gw_overall_accuracy = function(pred, obs, X , Y, grid,
                         # GW parameters
                         bw_defined = FALSE, bw = 0.3, 
                         kernel = gwr.bisquare,
                         verbose = FALSE)
{
  # define logit function
  alogit <- function(x){exp(x)/(1+exp(x))}
  # create 1/0 vector of True/False values
  vec = (pred == obs) + 0
  # define the coordinates matrix and data.frame for the data
  coords = cbind(X,Y)
  df = data.frame(pred, obs)
  # optimal kernel bandwidth
  if(!bw_defined) {
    NN <- ggwr.sel(vec~1, data = df, coords = coords, adapt = TRUE, gweight = kernel,
                   family=binomial, verbose = verbose)
  } else NN = bw	
  gwr.model_overall <- ggwr(vec~1, data = df, coords = coords, adapt = NN, 
                            gweight = kernel,
                            fit.points = grid, family= binomial) 
  # extract the probabilities from the SDF of the model 
  gwr.ovac <- alogit(data.frame(gwr.model_overall$SDF)[,2])
  return(list(Bandwidth = NN, GWOverall = gwr.ovac))
}
gw_oa = gw_overall_accuracy(pred = data$Boolean_Pred, obs = data$Boolean_Obs, 
                    X = data[,2], Y = data[,3], grid = grid.xy)
# uncomment to examine!
# gw_oa$Bandwidth
# summary(gw_oa$GWOverall)

# UA & PA
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
    coefs.sum <- rowSums(coefs) #logit ea+c
    p.x.eq.1.g.y.eq.1 <- alogit(coefs.sum) #n1 (sum coeffs)
    gwr.user <-  p.x.eq.1.g.y.eq.1 
    user_result = cbind(user_result, gwr.user)
    user_bw = append(user_bw, NN)
    # Producers Accuracy 
    if(!bw_defined) {
      NN <- ggwr.sel(pred.class~obs.class, data = df, coords = coords, adapt = TRUE, 
                     gweight = kernel, family=binomial, verbose = verbose)
    } else NN = bw
    gwr.model_prod <- ggwr(pred.class~obs.class, data = df, coords = coords, adapt = NN,
                           fit.points=grid, gweight = kernel, family= binomial) 
    coefs <- data.frame(gwr.model_prod$SDF)[,2:3]
    coefs.sum <- rowSums(coefs) #logit ea+c
    p.x.eq.1.g.y.eq.1 <- alogit(coefs.sum) #n1 (sum coeffs)
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
## apply the function: the results are used in the next section 
gw_up = gw_userprod_accuracy(pred = data$Boolean_Pred, obs = data$Boolean_Obs, 
                         X = data[,2], Y = data[,3], grid = grid.xy)
# uncomment to examine!
#gw_up$UserBW
#gw_up$ProdBW
#t(round(apply(gw_up$GWUser, 2, summary),3))
#t(round(apply(gw_up$GWProd, 2, summary),3))


## Table 2
## a little initial explore! 
# overall accuracy
tab2 = c(summary(gw_oa$GWOverall), Bandwidth = gw_oa$Bandwidth, Global = oa)
tab2 = t(as.matrix(tab2))
knitr::kable(tab2, digits = 3, row.names = F, 
             caption = "\\label{tab:tab2}A summary of the local measures of overall accuracy generated by the GW-GLM models, with the adaptive bandwidth and the global overall accuracy measure.")

## users
## Table 3
tab3 = t(round(apply(gw_up$GWUser, 2, summary),3))
rownames(tab3) <- class.names
tab3 = cbind(tab3, Bandwidths = gw_up$UserBW, Global = users)
knitr::kable(tab3, digits = 3, row.names = T, 
             caption = "\\label{tab:tab3}A summary of the local measures of user's accuracy generated by the GW-GLM models, with the adaptive bandwidths and the global user's accuracy measures.")

## producers
## Table 4
tab4 = t(round(apply(gw_up$GWProd, 2, summary),3))
rownames(tab4) <- class.names
tab4 = cbind(tab4, Bandwidths = gw_up$ProdBW, Global = producers)
knitr::kable(tab4, digits = 3, row.names = T, 
             caption = "\\label{tab:tab4}A summary of the local measures of producer's accuracy generated by the GW-GLM models, with the adaptive bandwidths and the global producer's measures.")


## Figure 4
allgwglm_data = data.frame(grid.xy, overall = gw_oa$GWOverall, 
                           gw_up$GWUser, gw_up$GWProd)
# make spatial if you want! 
# gwr_sf = st_as_sf(gwr_all_part1, coords = c("X", "Y"))
# add geographic projection  
# st_crs(gwr_sf) = 32633
# map with the df object for raster 

tit = paste0("Bandwidth = ", round(gw_oa$Bandwidth,2))
p1 = 
  ggplot() + 
		geom_raster(data = allgwglm_data, aes(x = X, y = Y, fill = overall)) +
		scale_fill_continuous_c4a_seq(palette="kovesi.linear_wh_rd_bk", name = "Overall") +
    geom_sf(data = mask, col = NA, fill = "white") +
    theme_bw() +
    coord_sf() +
    labs(subtitle = tit) +
    theme(legend.position = "bottom",
          axis.title=element_blank(),
	        axis.text=element_blank(),
	        axis.ticks=element_blank(),
          legend.text = element_text(size=8))
tit = paste0("Bandwidth = ", round(gw_up$UserBW["G"],2))
p2 = 
  ggplot() + 
		geom_raster(data = allgwglm_data, aes(x = X, y = Y, fill = U_G)) +
		scale_fill_continuous_c4a_seq(palette="kovesi.linear_wh_rd_bk", name = "User's") +
    geom_sf(data = mask, col = NA, fill = "white") +
    theme_bw() +
    coord_sf() +
    labs(subtitle = tit) +
    theme(legend.position = "bottom",
          axis.title=element_blank(),
	        axis.text=element_blank(),
	        axis.ticks=element_blank(),
          legend.text = element_text(size=8))
tit = paste0("Bandwidth = ", round(gw_up$ProdBW["G"],2))
p3 = 
  ggplot() + 
		geom_raster(data = allgwglm_data, aes(x = X, y = Y, fill = P_G)) +
		scale_fill_continuous_c4a_seq(palette="kovesi.linear_wh_rd_bk", name = "Producer's") +
    geom_sf(data = mask, col = NA, fill = "white") +
    theme_bw() +
    coord_sf() +
    labs(subtitle = tit) +
    theme(legend.position = "bottom",
          axis.title=element_blank(),
	        axis.text=element_blank(),
	        axis.ticks=element_blank(),
          legend.text = element_text(size=8))

plot_grid(p1, p2, p3, nrow = 1)


## Figure 5
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

## Tables 5 and 6 
# set up the gwxt framework with the data
dummy_xtab <- new_spxt(data.spdf,'Boolean_Pred', 'Boolean_Obs')
# define a bandwidth
bw = round(nrow(data.spdf)*0.3, 0)
# get the gwxt function
gwxt <- gwxtab_probe(dummy_xtab,adapt(bw))
# apply to locations
# point 1
coords.i = as.vector(st_coordinates(x[1,]))
gwcm1 <- round(gwxt(x = coords.i[1], y = coords.i[2] ), 2)
# point 2
coords.i = as.vector(st_coordinates(x[2,]))
gwcm2 <- round(gwxt(x = coords.i[1], y = coords.i[2] ), 2)
rownames(gwcm1) <- rownames(gwcm2) <- class.names
knitr::kable(gwcm1, digits = 2, caption = "\\label{tab:tab5}The GWCM for Point 1, with predicted (rows) and observed fuzzy class (columns). The values indicate GW pixel counts.", linesep = "") 
knitr::kable(gwcm2, digits = 2, caption = "\\label{tab:tab6}The GWCM for Point 2, with predicted (rows) and observed fuzzy class (columns).The values indicate GW pixel counts.", linesep = "") 


## Tables 7, 8 and 9 
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
rownames(tab7)  <- c("Overall", "Kappa")
knitr::kable(tab7, digits = 3, caption = "\\label{tab:tab7}The local overall and kappa accuracies from a GWCM under an adaptive bandwidth of 0.82, with the  global statistic.", linesep = "") 

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
knitr::kable(tab8, digits = 3, row.names = T, 
             caption = "\\label{tab:tab8}A summary of the local measures of user's accuracy generated by the GW correspondence matrix approach, with the adaptive bandwidths and the global user's accuracy measures.")

tab9 = t(apply(producer.tab, 2, summary))
rownames(tab9) <- class.names
tab9 = cbind(tab9, Bandwidths = gw_up$ProdBW, Global = producers)
knitr::kable(tab4, digits = 3, row.names = T, 
             caption = "\\label{tab:tab9}A summary of the local measures of producer's accuracy generated by the GW correspondence matrix approach, with the adaptive bandwidths and the global producers's accuracy measures.")


## Figure 6
# overallAllocD
bw1 = round(nrow(data.spdf)*gw_oa$Bandwidth, 0)
oAD1.spdf <- gwxtab_sample(grid.spdf,dummy_xtab,adapt(bw1),melt=overallAllocD)
bw2 = round(nrow(data.spdf)*0.3, 0)
oAD2.spdf <- gwxtab_sample(grid.spdf,dummy_xtab,adapt(bw2),melt=overallAllocD)
# join together
df = data.frame(grid.xy, oAD1 = oAD1.spdf$do.call.rbind..res., oAD2 = oAD2.spdf$do.call.rbind..res.)
# head(df)
# summary(df)
p1 = 
  ggplot() + 
		geom_raster(data = df, aes(x = X, y = Y, fill = oAD1)) +
		scale_fill_continuous_c4a_seq(palette="kovesi.linear_wh_rd_bk", name = "Overall Alloc.\nDifference", limits = c(0,15.5)) +
    geom_sf(data = mask, col = NA, fill = "white") +
    labs(subtitle = "Bandwidth = 0.82") +
    theme_bw() +
    coord_sf() +
    theme(
          axis.title=element_blank(),
	        axis.text=element_blank(),
	        axis.ticks=element_blank())
p2 = 
  ggplot() + 
		geom_raster(data = df, aes(x = X, y = Y, fill = oAD2)) +
		scale_fill_continuous_c4a_seq(palette="kovesi.linear_wh_rd_bk", name = "Overall Alloc.\nDifference", limits = c(0,15.5)) +
    geom_sf(data = mask, col = NA, fill = "white") +
    labs(subtitle = "Bandwidth = 0.30") +
    theme_bw() +
    coord_sf() +
    theme(
          axis.title=element_blank(),
	        axis.text=element_blank(),
	        axis.ticks=element_blank())
# remove the legends
legpx = get_legend((p1)); p1 = p1 + theme(legend.position = "none")
p2 = p2 + theme(legend.position = "none")
# plot together
plot_grid(p1,p2, legpx, ncol = 3, rel_widths = c(.8,.8,0.3))


## Table 10
## GW Analysis of Fuzzy Difference
# define a function to do this
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
woodland_gw_acc =      gw_fuzzy_accuracy(pred = data$Woodland_Pred, 
                                      obs = data$Woodland_Obs,
                                      X, Y, grid.xy) 
grazing_gw_acc =    gw_fuzzy_accuracy(pred = data$Grazing_Pred, 
                                      obs = data$Grazing_Obs,
                                      X, Y, grid.xy) 
bare_gw_acc =       gw_fuzzy_accuracy(pred = data$Urban_Pred,
                                      obs = data$Urban_Obs,
                                      X, Y, grid.xy) 
# summarise the results
tab10 = rbind(summary(bare_gw_acc$FuzzyAcc),
              summary(grazing_gw_acc$FuzzyAcc),
              summary(urban_gw_acc$FuzzyAcc),
              summary(vegetation_gw_acc$FuzzyAcc),
              summary(woodland_gw_acc$FuzzyAcc))
bws = c(bare_gw_acc$BW, grazing_gw_acc$BW, urban_gw_acc$BW, 
        vegetation_gw_acc$BW, woodland_gw_acc$BW)
tab10 = cbind(tab10, Bandwidths = bws)
rownames(tab10) = class.names
knitr::kable(tab10, digits = 3, caption = "\\label{tab:tab10}Summaries of the distribution of the GW fuzzy certainty for each class, with optimised kernel bandwidths.", linesep = "") 


## Figure 7
# create a dataset for the plot
df = data.frame(grazing = grazing_gw_acc$FuzzyAcc, 
                vegetation = vegetation_gw_acc$FuzzyAcc, 
                grid.xy)
# map with the df object for raster 
tit = paste0("Grazing Land: Bandwidth = ", round(grazing_gw_acc$BW,2))
p1 = 
  ggplot() + 
		geom_raster(data = df, aes(x = X, y = Y, fill = grazing)) +
		scale_fill_continuous_c4a_seq(palette="kovesi.linear_wh_yl_gn_bu", 
		                              name = "Certainty",
		                              limits = c(0.7, 0.95)) +
    geom_sf(data = mask, col = NA, fill = "white") +
    theme_bw() +
    coord_sf() +
    labs(subtitle = tit) +
    theme(
          axis.title=element_blank(),
	        axis.text=element_blank(),
	        axis.ticks=element_blank())
tit = paste0("Vegetation: Bandwidth = ", round(vegetation_gw_acc$BW,2))
p2 = 
  ggplot() + 
		geom_raster(data = df, aes(x = X, y = Y, fill = vegetation)) +
		scale_fill_continuous_c4a_seq(palette="kovesi.linear_wh_yl_gn_bu", 
		                              name = "",
		                              limits = c(0.7, 0.95)) +
    geom_sf(data = mask, col = NA, fill = "white") +    
    theme_bw() +
    coord_sf() +
    labs(subtitle = tit) +
    theme(
          axis.title=element_blank(),
	        axis.text=element_blank(),
	        axis.ticks=element_blank())
# remove the legends
legpx = get_legend((p1)); p1 = p1 + theme(legend.position = "none")
p2 = p2 + theme(legend.position = "none")
# plot together
plot_grid(p1,p2, legpx, ncol = 3, rel_widths = c(.8,.8,0.2))


## Table 10
## The fuzzy change function - this generates all fuzzy information 
# It takes as input data.frames of the observed and predicted fuzzy memberships from the validation data
# These need to be similarly named
# For each observation in the validation data, it returns 
# 1. the F fuzzy CM, "FCM" in the output
# 2. the fuzzy Omission and Commission table, "Fom_Fcom_tab"
# 3. the per class fuzzy correspondence (diagonal), "F_corr"
# 4. the per class fuzzy omission (lower off-diagonal), "F_om"
# 5. the per class fuzzy commission (upper off-diagonal), "F_com"
fuzzy_change_analysis = function(predicted = pred_df, observed = obs_df) {
  # extract class names
  c.ns = sort(names(observed))
  # make sure the inputs are ordered the same way
  pred_df = predicted[,c.ns]
  obs_df = observed[, c.ns]
  # get the "not" fuzzy values for use in the bounded difference
  notpred = 1- pred_df
  notobs  = 1- obs_df
  # set up some output objects
  f.cm = f.cm2 = matrix(0,nrow = 5, ncol = 5) 
  gb_vec = lb_vec = same_vec = NULL
  for (i in 1:length(c.ns)){
    # define row class (predicted!)
    class.i = c.ns[i]
    # empty vectors for gain and loss boundaries
    gb.vec = lb.vec = NULL
    for(j in 1:length(c.ns)){
      class.j = c.ns[j]
      obs.j = as.vector(obs_df[, class.j])
      pred.j = as.vector(pred_df[, class.j])
      if (i == j) { # minimum interval for diagonal
        same.i = sapply(1:nrow(pred_df), function(x) min(pred.j[x], obs.j[x]))
        # sum and add to fcm
        f.cm[i, j] = sum(same.i)
        # bind to output
        same_vec  = cbind(same_vec, same.i) 
        tit = paste0(class.i, "_", class.j)
        colnames(same_vec)[ncol(same_vec)] <- tit
      }  
      if(i < j ) { # fuzzy Producers / Omission errors:  
        # how column wise observed class is different to predicted
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
      if(i > j ) { # fuzzy User's Commission errors: 
        # how row wise observed class is different to predicted
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
  # do fuzzy commission (losses) and omission (gains)
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
obs_df = data %>% dplyr::select(ends_with("_Obs"),-Boolean_Obs)
names(obs_df) = gsub("_Obs", "", names(obs_df))
pred_df = data %>% dplyr::select(ends_with("_Pred"),-Boolean_Pred)
names(pred_df) = gsub("_Pred", "", names(pred_df))
# do some sorting ordering or variables
c.ns = sort(names(obs_df))
pred_df = pred_df[,c.ns]
obs_df = obs_df[, c.ns]
# apply the function
fuzzy_res = fuzzy_change_analysis(predicted = pred_df, observed = obs_df)
fcm = fuzzy_res$FCM
colnames(fcm) = class.list
rownames(fcm) = class.names
colnames(fuzzy_res$Fom_Fcom_tab) = c("Omission", "Commission")
tab11 = cbind(fcm, fuzzy_res$Fom_Fcom_tab)

knitr::kable(tab11, digits = 3, caption = "\\label{tab:tab11}The fuzzy correspondence matrix of fuzzy predicted (rows) against fuzzy observed (columns) class memberships, with the fuzzy omission and commission uncertainties.", linesep = "") 


## GW fuzzy uncertainty function
# takes the output of the fuzzy_change_analysis function (fuzzy_res)
gw_fuzzy_uncertainty = function(class.name = "Grazing",fuzzy_res, X, Y, grid,
                         # GW parameters
                         bw_defined = FALSE, bw = 0.3, kernel = gwr.bisquare,
                         verbose = FALSE)
{
  # extract the data from the fuzzy_res input
  # fuzzy correspondence
  df1 <- fuzzy_res$F_cor %>% data.frame() %>% select(contains(class.name))
  fcorres_df = data.frame(X, Y, val = as.vector(unlist(df1)))
  # fuzzy omission
  col_tit = paste0("col_", class.name)
  df1 <- fuzzy_res$F_om %>% data.frame() %>% select(contains(col_tit))
  df2 <- fuzzy_res$F_com %>% data.frame() %>% select(contains(col_tit))
  fom_df = data.frame(X, Y, val = rowSums(df1) + rowSums(df2))
  # fuzzy commission
  row_tit = paste0("row_", class.name)
  df1 <- fuzzy_res$F_com %>% data.frame() %>% select(contains(row_tit))
  df2 <- fuzzy_res$F_com %>% data.frame() %>% select(contains(row_tit))
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
fuzzy_bare = gw_fuzzy_uncertainty("Bare",fuzzy_res, X, Y, grid.xy)
fuzzy_grazing = gw_fuzzy_uncertainty("Grazing",fuzzy_res, X, Y, grid.xy)
fuzzy_urban = gw_fuzzy_uncertainty("Urban",fuzzy_res, X, Y, grid.xy)
fuzzy_vegetation = gw_fuzzy_uncertainty("Vegetation",fuzzy_res, X, Y, grid.xy)
fuzzy_woodland = gw_fuzzy_uncertainty("Woodland",fuzzy_res, X, Y, grid.xy)
# have a look at the results
# names(fuzzy_grazing)
# fuzzy_grazing$BWs
# summary(fuzzy_grazing$fuzzy_res)

# summary(fuzzy_bare$fuzzy_res)
# summary(fuzzy_grazing$fuzzy_res)
# summary(fuzzy_urban$fuzzy_res)
# summary(fuzzy_vegetation$fuzzy_res)
# summary(fuzzy_woodland$fuzzy_res)
# higher variation in
# grazing f_corr
# woodland f_user


## Figure 8
# grazing f_corr and woodland f_user

df = data.frame(fuzzy_grazing$fuzzy_res, grid.xy)
tit = paste0("Grazing fuzzy correspondence: BW = ", round(fuzzy_grazing$BWs[1], 2))
p1 = 
  ggplot() + 
		geom_raster(data = df, aes(x = X, y = Y, fill = f_corr)) +
		scale_fill_continuous_c4a_seq(palette="kovesi.linear_wh_yl_gn_bu", 
		                              name = "",
		                              breaks = pretty_breaks(n = 2)) +
    geom_sf(data = mask, col = NA, fill = "white") +
    labs(subtitle = tit) +
    theme_bw() +
    coord_sf() +
    theme(legend.position = "bottom",
          axis.title=element_blank(),
	        axis.text=element_blank(),
	        axis.ticks=element_blank())

df = data.frame(fuzzy_woodland$fuzzy_res, grid.xy)
tit = paste0("Woodland fuzzy user's certainty: BW = ", round(fuzzy_woodland$BWs[3], 2))
p2 = 
  ggplot() + 
		geom_raster(data = df, aes(x = X, y = Y, fill = f_user)) +
		scale_fill_continuous_c4a_seq(palette="kovesi.linear_wh_yl_gn_bu", 
		                              name = "",
		                              breaks = pretty_breaks(n = 2)) +
    geom_sf(data = mask, col = NA, fill = "white") +
    labs(subtitle = tit) +
    theme_bw() +
    coord_sf() +
    theme(legend.position = "bottom",
          axis.title=element_blank(),
	        axis.text=element_blank(),
	        axis.ticks=element_blank())
plot_grid(p1, p2, ncol = 2)


## Table 12

Pred = data[,11:15]
Obs = data[,4:8]
oAD = aDist(Pred+0.001, Obs+0.001)
AS.all = sapply(1:nrow(data), 
                   function(x) aDist(Pred[x,]+0.001,Obs[x,]+0.001 ))
tab12 = c(summary(AS.all), Overall = oAD)

knitr::kable(t(tab12), digits = 2, caption = "\\label{tab:tab12}A summary of the Aitchison distance distribution when calculated for each observation in the validation data, and the overall Aitchison distance value.", linesep = "") 


## Figure 9
# create data# create data# create data
X = data[,2]
Y = data[,3]
coords = cbind(X,Y)
df_point = data.frame(ADist = AS.all, coords)
# GW Aich dist
NN <- gwr.sel(ADist~1, data = df_point, coords = coords, adapt = TRUE, 
              gweight = gwr.bisquare, verbose = F)
gwr.model<- gwr(ADist~1, data = df_point, coords = coords, adapt = NN,  
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
  geom_point(data = df_point, aes(x = X, y = Y, size = ADist), alpha = 0.5) + 
  scale_size_area("Point", breaks = c(1,3,6)) +
  geom_sf(data = mask, col = NA, fill = "white") +
  theme_bw() +
  labs(subtitle = tit) +
  coord_sf() +
  theme(legend.position = "bottom",
        axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank())


## Figure 10

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
## check
## get nearby obs & weight observations
# nearby_index = gw_get_nearby(1, adaptive = T, bw, dist_mat)
# nearby_index
# wts = get_nearby_weights(1, nearby_index, bw, dist_mat, kernel = "bisquare", adaptive = T)
# round(wts, 4)
# define empty matrix
res_ad <-NULL
for(i in 1:nrow(dist_mat)){
  nearby_index = gw_get_nearby(i, adaptive = T, bw, dist_mat)
  wts.i = get_nearby_weights(i, nearby_index, bw, dist_mat, kernel = "gaussian", adaptive = T)
  # weight the data 
  Pred.i = Pred * wts.i 
  Obs.i = Obs * wts.i   
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
# create the data frame for mapping
df = data.frame(AichDist = gwr.res, X = grid.xy[,"X"], Y = grid.xy[,"Y"])
# dim(df)
ggplot() + 
  geom_raster(data = df, aes(x = X, y = Y, fill = AichDist)) +
  scale_fill_continuous_c4a_seq(palette="kovesi.linear_wh_yl_gn_bu", 
                                name = "Aitchison Distance",
                                breaks = scales::trans_breaks(identity, identity, n = 4)) +
  geom_sf(data = mask, col = NA, fill = "white") +
  #geom_point(data = df_point, aes(x = X, y = Y, size = ADist), alpha = 0.5) + 
  #scale_size_area("Point", breaks = c(1,3,6)) +
  theme_bw() +
  labs(subtitle = tit) +
  coord_sf() +
  theme(legend.position = "bottom",
        axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank())


## Table 13

tab13 = tibble(`GW accuracy model` = c("GW-GLM", "GWCM", "GW Fuzzy certainty I", "GW Fuzzy certainty II", 
                                        "GW Aitchison Distance"),
                `Underling GWR model` = c("Generalised GWR", "None", "GWR", "GWR", "GWR"),
                `Bandwidth optimisation` = c("Yes", "No", "Yes", "Yes", "Yes"))

knitr::kable(tab13, row.names = FALSE , caption = "\\label{tab:tab13}A summary of the GW accuracy approaches, their underlying models and their ability to optimise GW kernel bandwidths.", linesep = "") 



