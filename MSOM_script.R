
# Publication details                                                                                                          
# Title: # Anthropogenic and environmental correlates of spatial patterns of co-occurrence of small felids in a montane landscape                                    #
# Authors: Karma Choki, Egil Dröge, Claudio Sillero-Zubiri, David W. Macdonald, Ugyen Penjor                 
# Journal: Global Ecology and Conservation, Volume -, 2024                                                 
# DOI:                                                                                                       

# Multispecies occupancy model for two or more interacting species

library(unmarked)
library(AICcmodavg)
library(dplyr)
library(ggplot2)

# Load data
ch_asiaticgoldencat <- read.csv("Asiangoldencat.csv", header=T, row.names=1, na.strings="NA")
ch_marbledcat <- read.csv("Marbledcat.csv", header=T, row.names=1, na.strings="NA")
ch_leopardcat <- read.csv("Leopardcat.csv", header=T, row.names=1, na.strings="NA")

# Check NAs and convert detection history into a matrix
# This will throw you an error if something is not right in your data.
goodMat <- function(x){
  xm <- as.matrix(x)
  n <- rowSums(!is.na(xm)) 
  noNA <- n > 0
  out <- xm[noNA, ]
}

str(goldencat <- goodMat(ch_asiaticgoldencat))
str(marbledcat <- goodMat(ch_marbledcat))
str(leopardcat <- goodMat(ch_leopardcat))

class(goldencat)

# Import site covariates
cov <- read.csv("sitecovs.csv")
head(cov)

# Select required covariates
trail <- cov$pTrail
effort <- scale(cov$effort)
HDisturb <- scale(cov$pDisturb)
muntjac <- scale(cov$bd_nmix)
HDens500m <- scale(cov$setden_500M)
HDens1km <- scale(cov$setden_1KM)
HDens2km <- scale(cov$setden_2KM)
HDens4km <- scale(cov$setden_4KM)
Tcover500m <- scale(cov$TC500)
Tcover1km <- scale(cov$TC1)
Tcover2km <- scale(cov$TC2)
Tcover4km <- scale(cov$TC4)
RivDen500m <- scale(cov$rivFM_500M)
RivDen1km <- scale(cov$rivFM_1KM)
RivDen2km <- scale(cov$rivFM_2KM)
RivDen4km <- scale(cov$rivFM_4KM)
slope500m <- scale(cov$sloposFM_500M)
slope1km <- scale(cov$sloposFM_1KM)
slope2km <- scale(cov$sloposFM_2KM)
slope4km <- scale(cov$sloposFM_4KM)
ele500m <- scale(cov$ele500)

######################################################################################################
######################################################################################################

# Check for spatial autocorrelation
library(ncf)

# Extract covariates (just to build a glm)
elev <- scale(cov$ele500)
settle <- scale(cov$setden_2KM)
river <- scale(RivDen2km)
forest <- scale(cov$TC2)
slope <- scale(cov$sloposFM_2KM)
prey <- scale(cov$bd_nmix)

# We don't need temporal replicates to check spatial autocorrelation
# So we convert detection history into a matrix and just get one column with det/non-det data

# Convert into a matrix
DHagc <- as.matrix(ch_asiaticgoldencat)

# Get a single column detection/non-detection data for this analysis only
det_agc <- apply(DHagc, 1, max, na.rm=T) #max means max value of index(?)
View(det_agc)

# Fit a logistic regression for naive occupancy on global model
glm_agc  <- glm(det_agc ~ elev + I(elev^2) + settle + river + forest + slope + prey, family=binomial)

# Spline correlogram for naive occupancy model
# AGC
correlog_agc <- spline.correlog(
  x=cov[, 3],
  y=cov[, 4],
  z=residuals(glm_agc, type="pearson"), xmax=50,
  latlon=T)
plot(correlog_agc) 

# Same for other two species
# MC
DHmc <- as.matrix(ch_marbledcat)
det_mc <- apply(DHmc, 1, max, na.rm=T) 
glm_mc  <- glm(det_mc ~ elev + I(elev^2) + settle + river + forest + slope + prey, family=binomial)
correlog_mc <- spline.correlog(
  x=cov[, 3],
  y=cov[, 4],
  z=residuals(glm_mc, type="pearson"), xmax=50,
  latlon=T)
plot(correlog_mc) 

# LC
DHlc <- as.matrix(ch_leopardcat)
det_lc <- apply(DHlc, 1, max, na.rm=T) 
glm_lc  <- glm(det_lc ~ elev + I(elev^2) + settle + river + forest + slope + prey, family=binomial)
correlog_lc <- spline.correlog(
  x=cov[, 3],
  y=cov[, 4],
  z=residuals(glm_lc, type="pearson"), xmax=50,
  latlon=T)
plot(correlog_lc) 

# The CIs for correlation estimates on y-axis includes zero, this means there is no evidence
# of the presence of autocorrelation for all the species

######################################################################################################
######################################################################################################

# Convert to data frame from matrix-array for covariates
cov1 <- as.data.frame(cbind(trail, effort, HDisturb, muntjac, HDens500m,HDens1km,HDens2km,HDens4km, 
                            Tcover500m,Tcover1km, Tcover2km, Tcover4km, 
                            RivDen500m,RivDen1km, RivDen2km, RivDen4km, 
                            slope500m, slope1km, slope2km, slope4km, ele500m))
summary(cov1)

# Rename covaiates 
cov2 <- rename(cov1, 
               "trail"=trail,
               "effort"=V2,
               "HDisturb"=V3,
               "muntjac"=V4,
               "HDens_500m"=V5,
               "HDens_1km"=V6,
               "HDens_2km"=V7,
               "HDens_4km"=V8,
               "Tcover_500m"=V9,
               "Tcover_1km"=V10,
               "Tcover_2km"=V11,
               "Tcover_4km"=V12,
               "RivDen_500m"=V13,
               "RivDen_1km"=V14,
               "RivDen_2km"=V15,
               "RivDen_4km"=V16,
               "slope_500m"=V17,
               "slope_1km"=V18,
               "slope_2km"=V19,
               "slope_4km"=V20,
               "ele500_m"=V21)
summary(cov2)

# We need to bundle the species’ detection history matrix into a list
spsList <- list(goldencat, marbledcat, leopardcat)
names(spsList)
spsNames <- c("Goldencat", "Marbledcat", "Leopardcat")
names(spsList) <- spsNames
str(spsList)

spsList <- c(spsList, list(cov2)) #cov data.frame is converted to list

# Prepare new shorthand names plus covariate
newNames <- c("Goldencat", "Marbledcat", "Leopardcat", "covars")
names(spsList) <- newNames
str(spsList)

# We will combine the detection data for the three species into a named list
ylist <- list(goldencat=spsList$Goldencat, marbledcat=spsList$Marbledcat, leopardcat=spsList$Leopardcat)
str(ylist)
head(ylist$goldencat)

# Convert to unmarked object
umf <- unmarkedFrameOccuMulti(y=ylist, siteCovs=cov2)
summary(umf)

# Only using the top model for prediction
m26 <- occuMulti(
  detformulas=rep('~as.factor(trail) + effort + HDisturb', 3),              
  stateformulas=c('~slope_4km',
                  '~Tcover_4km + RivDen_2km+slope_500m',
                  '~HDens_1km + Tcover_4km+slope_4km',
                  '~ele500_m + HDens_4km',  
                  '~ele500_m + HDens_4km',  
                  '~ele500_m + HDens_4km', 
                  0), 
  data=umf)
summary(m26)

confint(m26, type="det") # confidence interval for detection probability estimates
confint(m26, type="state") # confidence interval for occupancy probability estimates

#######################################################################################################
#CO function of Human Density plot

k <- m26@estimates["state"]
p <- k@estimates[c(13,16,19)]

names(p) <- c("AGC:MC", "AGC:LC", "MC:LC")
str(p)
class(p)

library(tidyverse)
df <- enframe(p, name = "Species", value = "Estimates")
str(df)

z <- confint(m26, type="state") 
t <- z[c(13,16,19), ]
k <- data.frame(t)
j1 <- as.numeric(k[,1])
j2 <- as.numeric(k[,2])
l <- as.data.frame(cbind(j1, j2))
CI <- l%>%rename(
   "Lower"= j1,
  "Upper"= j2
)
CI$Lower <- as.numeric(CI$Lower)
CI$Upper <- as.numeric(CI$Upper)
class(CI$Lower)
CI
v <- cbind(df, CI)
v

#Plot
f <- ggplot(data=v, aes(x=Estimates, y=Species, color=(factor(Species, levels= c("AGC:MC", "AGC:LC", "MC:LC")))))+
  scale_shape_manual(values=c(13, 13, 13))+ 
  scale_color_manual(values=c( "#e8a543","#65b185","#e9493b"), name=NULL) +
  geom_line() +
  geom_vline(xintercept = 0, color="black", linetype = "dotted")+
  geom_point(size=4)+
  geom_errorbarh(aes(xmin=Lower, xmax=Upper), height=0.09)+
    ylab("Pairwise Interacting Species") +
  xlab("β coefficient estimate")+ 
  xlim(c(-0.50, 0.50))+
    theme_bw()+
  theme(text=element_text(size=15, hjust = 0.5))
f

ggsave("β coefficient plot.png", f, width=6.5, height=4.5)


##########################################
########### Plotting marginal occupancy probability#################################

# We can use outputs from predict() to compare marginal occupancy across species with a plot. 
sp1_marg <- predict(m26, type="state", species="goldencat") #predicts the marginal occupancy for all sites #Mention in paper
MO_AGC <- c(mean(sp1_marg[,1]), mean(sp1_marg[,3]), mean(sp1_marg[,4]))
sp1_marg_Det <- predict(m26, type="det", species="goldencat") #predicts the marginal occupancy for all sites
Det_AGC <- c(mean(sp1_marg_Det[,1]), mean(sp1_marg_Det[,3]), mean(sp1_marg_Det[,4]))
#0.2081516 0.1701973 0.2532498

sp2_marg <- predict(m26, type="state", species="marbledcat")
MO_MC <- c(mean(sp2_marg[,1]), mean(sp2_marg[,3]), mean(sp2_marg[,4])) #This line is important
sp2_marg_Det <- predict(m26, type="det", species="marbledcat") #predicts the marginal occupancy for all sites
Det_MC <- c(mean(sp2_marg_Det[,1]), mean(sp2_marg_Det[,3]), mean(sp2_marg_Det[,4]))
#0.2357853 0.1879319 0.2980222

sp3_marg <- predict(m26, type="state", species="leopardcat")
MO_LC <- c(mean(sp3_marg[,1]), mean(sp3_marg[,3]), mean(sp3_marg[,4]))
sp3_marg_Det <- predict(m26, type="det", species="leopardcat") #predicts the marginal occupancy for all sites
Det_LC <- c(mean(sp3_marg_Det[,1]), mean(sp3_marg_Det[,3]), mean(sp3_marg_Det[,4]))
#0.2390136 0.2031544 0.2807317

# Combine all three species
all_marg <- rbind(sp1_marg[1,], sp2_marg[1,], sp3_marg[1,])
all_det <- rbind(sp1_marg_Det[1,], sp2_marg_Det[1,], sp3_marg_Det[1,])
all_marg$Species <- c("Asiatic golden cat", "Marbled cat", "Leopard cat")
all_det$Species <- c("Asiatic golden cat", "Marbled cat", "Leopard cat")
all_marg$Type <- "Occupancy"
all_det$Type <- "Detection"

# Do the plot 
k <- ggplot(data = tidy_est, aes(x = Type, y = Predicted, ymin = lower, ymax = upper, color = Species)) +
  geom_pointrange(position = position_dodge(0.3)) +
  scale_color_manual(values = c("#e9493b","#E7B800","#65b185"), name = NULL) +
  ylab("Predicted probability") +
  xlab("") +
  theme_bw() +
  theme(
    legend.position = c(0.86, 0.9),
    legend.background = element_blank(),
    text = element_text(size = 15),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

k

ggsave("Marginal_occ_det_prob_all_sps_1.png", k, width=6.9, height=4.5)
#######################################################################################################

#Plot the effect of the covariates on marginal occupancy _Rota et al. 2016
#Slope on AGC
nd_marg <- data.frame(
  slope_4km = seq(min(cov2$slope_4km), max(cov2$slope_4km),
                  length.out= 100),
  Tcover_4km = rep(mean(cov2$Tcover_4km), 100),
  HDens_1km  = rep(mean(cov2$HDens_1km), 100),
  RivDen_2km  = rep(mean(cov2$RivDen_2km), 100),
  slope_500m =rep(mean(cov2$slope_500m), 100),
  HDens_4km= rep(mean(cov2$HDens_4km), 100),
  ele500_m = rep(mean(cov2$ele500_m), 100)
)

#Marginal occupancy of AGC using new dataframe (Averages over presence/absence of all other species)
#This is not 1st order natural parameter bc it takes account the presence and absence of other species
AGC_marg <- predict(m26, type= "state", species= "goldencat", newdata= nd_marg) 

#formating data for ggplot
gg_df_marg <- data.frame(
  Slpe= nd_marg$slope_4km,
  occupancy= AGC_marg$Predicted,
  low= AGC_marg$lower,
  high= AGC_marg$upper
)

library(ggplot2)
AGC_marg_fig <- ggplot(gg_df_marg, aes(x=Slpe, y=occupancy))+
  geom_ribbon(aes(ymin=low, ymax=high),fill="#f20574", alpha= 0.2, col="grey93")+
  geom_line(color = "red", size = 0.7) +
  ylim(0, 1)+
  ylab("AGC marginal occupancy")+
  xlab("Slope")+
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank())+
  theme(text= element_text(size= 20))

AGC_marg_fig 
ggsave("AGC_MO_slope.png", last_plot(),width=6, height=5 )


#Slope on MC
nd_marg1 <- data.frame(
  slope_500m = seq(min(cov2$slope_500m), max(cov2$slope_500m),
                   length.out= 100),
  Tcover_4km = rep(mean(cov2$Tcover_4km), 100),
  HDens_1km  =rep(mean(cov2$HDens_1km), 100),
  RivDen_2km  =rep(mean(cov2$RivDen_2km), 100),
  slope_4km   =rep(mean(cov2$slope_4km), 100),
  HDens_4km  = rep(mean(cov2$HDens_4km), 100),
  ele500_m = rep(mean(cov2$ele500_m), 100)
)

#Marginal occupancy of MC using new dataframe (Averages over presence/absence of all other species)
#This is not 1st order natural parameter bc it takes account the presence and absence of other species
MC_marg <- predict(m26, type= "state", species= "marbledcat", newdata= nd_marg1) 

#formating data for ggplot
gg_df_marg1 <- data.frame(
  Slope= nd_marg1$slope_500m,
  occupancy= MC_marg$Predicted,
  low= MC_marg$lower,
  high= MC_marg$upper
)

MC_marg_fig <- ggplot(gg_df_marg1, aes(x=Slope, y=occupancy))+
  geom_ribbon(aes(ymin=low, ymax=high),fill="#f20574", alpha= 0.2, col="grey93")+
  geom_line(color = "red", size = 0.7) + 
  ylab("MC marginal occupancy")+
  xlab("Slope")+ 
  ylim(0, 1)+
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank())+
  theme(text= element_text(size= 20))

MC_marg_fig 
ggsave("MC_MO_slope.png", last_plot(),width=6, height=5 )


#Slope on LC
nd_marg2 <- data.frame(
  slope_4km = seq(min(cov2$slope_4km), max(cov2$slope_4km),
                  length.out= 100),
  Tcover_4km = rep(mean(cov2$Tcover_4km), 100),
  HDens_1km  =rep(mean(cov2$HDens_1km), 100),
  RivDen_2km  =rep(mean(cov2$RivDen_2km), 100),
  slope_500m   =rep(mean(cov2$slope_500m), 100),
  HDens_4km  = rep(mean(cov2$HDens_4km), 100),
  ele500_m = rep(mean(cov2$ele500_m), 100)
)

LC_marg <- predict(m26, type= "state", species= "leopardcat", newdata= nd_marg2) 

#formating data for ggplot
gg_df_marg2 <- data.frame(
  Slope= nd_marg2$slope_4km,
  occupancy= LC_marg$Predicted,
  low= LC_marg$lower,
  high= LC_marg$upper
)

LC_marg_fig <- ggplot(gg_df_marg2, aes(x=Slope, y=occupancy))+
  geom_ribbon(aes(ymin=low, ymax=high),fill="#f20574", alpha= 0.2, col="grey93")+
  geom_line(color = "red", size = 0.7) + 
  ylab("LC marginal occupancy")+
  xlab("Slope")+ 
  ylim(0, 1)+
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank())+
  theme(text= element_text(size= 20))

LC_marg_fig 
ggsave("LC_MO_slope.png", last_plot(),width=6, height=5 )


#Forest cover on MC
nd_margF <- data.frame(
  Tcover_4km = seq(min(cov2$Tcover_4km), max(cov2$Tcover_4km),
                   length.out= 100),
  slope_500m = rep(mean(cov2$slope_500m), 100),
  HDens_1km  =rep(mean(cov2$HDens_1km), 100),
  RivDen_2km  =rep(mean(cov2$RivDen_2km), 100),
  slope_4km   =rep(mean(cov2$slope_4km), 100),
  HDens_4km  = rep(mean(cov2$HDens_4km), 100),
  ele500_m = rep(mean(cov2$ele500_m), 100)
)

MC_marg2 <- predict(m26, type= "state", species= "marbledcat", newdata= nd_margF) 

#formating data for ggplot
gg_df_marg4 <- data.frame(
  TC= nd_margF$Tcover_4km,
  occupancy= MC_marg2$Predicted,
  low= MC_marg2$lower,
  high= MC_marg2$upper
)

library(ggplot2)
MC_marg_fig_FC <- ggplot(gg_df_marg4, aes(x=TC, y=occupancy))+
  geom_ribbon(aes(ymin=low, ymax=high),fill="green", alpha= 0.2, col="grey93")+
  geom_line(color = "dark green", alpha=4, size = 1) + 
  ylab("MC marginal occupancy")+
  xlab("Forest Cover")+ 
  ylim(0,1)+
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank())+
  theme(text= element_text(size= 20))

MC_marg_fig_FC 

ggsave("MC_MO_FC.png", last_plot(),width=6, height=5 )

#River Density on MC
nd_margRiv <- data.frame(
  RivDen_2km = seq(min(cov2$RivDen_2km), max(cov2$RivDen_2km),
                   length.out= 100),
  Tcover_4km = rep(mean(cov2$Tcover_4km), 100),
  slope_500m = rep(mean(cov2$slope_500m), 100),
  HDens_1km  =rep(mean(cov2$HDens_1km), 100),
  slope_4km   =rep(mean(cov2$slope_4km), 100),
  HDens_4km  = rep(mean(cov2$HDens_4km), 100),
  ele500_m = rep(mean(cov2$ele500_m), 100)
)

MC_margRi <- predict(m26, type= "state", species= "marbledcat", newdata= nd_margRiv) 

#formating data for ggplot
gg_df_margRi <- data.frame(
  Ri= nd_margRiv$RivDen_2km,
  occupancy= MC_margRi$Predicted,
  low= MC_margRi$lower,
  high= MC_margRi$upper
)

library(ggplot2)
MC_marg_fig_Ri <- ggplot(gg_df_margRi, aes(x=Ri, y=occupancy))+
  geom_ribbon(aes(ymin=low, ymax=high),fill="blue", alpha= 0.18, col="grey93")+
  geom_line(color = "blue", size = 1) + 
  ylab("MC marginal occupancy")+
  xlab("River Density")+ 
  ylim(0,1)+
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank())+
  theme(text= element_text(size= 20))

MC_marg_fig_Ri 

ggsave("MC_MO_RivDen.png", last_plot(),width=6, height=5 )

#Forest cover on LC
nd_margFC <- data.frame(
  Tcover_4km = seq(min(cov2$Tcover_4km), max(cov2$Tcover_4km),
                   length.out= 100),
  slope_500m = rep(mean(cov2$slope_500m), 100),
  HDens_1km  =rep(mean(cov2$HDens_1km), 100),
  RivDen_2km  =rep(mean(cov2$RivDen_2km), 100),
  slope_4km   =rep(mean(cov2$slope_4km), 100),
  HDens_4km  = rep(mean(cov2$HDens_4km), 100),
  ele500_m = rep(mean(cov2$ele500_m), 100)
)

LC_marg2 <- predict(m26, type= "state", species= "leopardcat", newdata= nd_margFC) 

#formating data for ggplot
gg_df_marg5 <- data.frame(
  TC1= nd_margFC$Tcover_4km,
  occupancy= LC_marg2$Predicted,
  low= LC_marg2$lower,
  high= LC_marg2$upper
)

library(ggplot2)
LC_marg_fig_FC <- ggplot(gg_df_marg5, aes(x=TC1, y=occupancy))+
  geom_ribbon(aes(ymin=low, ymax=high),fill="green", alpha= 0.2, col="grey93")+
  geom_line(color = "dark green", alpha=4, size = 1) + 
  ylab("LC marginal occupancy")+
  xlab("Forest Cover")+ 
  ylim(0,1)+
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank())+
  theme(text= element_text(size= 20))

LC_marg_fig_FC 

ggsave("LC_MO_FC.png", last_plot(),width=6, height=5 )

#housing density on LC
nd_margHDen <- data.frame(
  HDens_1km = seq(min(cov2$HDens_1km),  max(cov2$HDens_1km),
                 length.out= 100),
  Tcover_4km= rep(mean(cov2$Tcover_4km), 100),
  slope_500m = rep(mean(cov2$slope_500m), 100),
  RivDen_2km  =rep(mean(cov2$RivDen_2km), 100),
  slope_4km   =rep(mean(cov2$slope_4km), 100),
  HDens_4km  = rep(mean(cov2$HDens_4km), 100),
  ele500_m = rep(mean(cov2$ele500_m), 100)
)

LC_margHD <- predict(m26, type= "state", species= "leopardcat", newdata= nd_margHDen) 

#formating data for ggplot
gg_df_margHD <- data.frame(
  HD= nd_margHDen$HDens_1km,
  occupancy= LC_margHD$Predicted,
  low= LC_margHD$lower,
  high= LC_margHD$upper
)

library(ggplot2)
LC_marg_fig_HD <- ggplot(gg_df_margHD, aes(x=HD, y=occupancy))+
  geom_ribbon(aes(ymin=low, ymax=high),fill="yellow", alpha= 0.2, col="grey93")+
  geom_line(color = "orange", size = 1) + 
  ylab("LC marginal occupancy")+
  xlab("Housing Density")+ 
  ylim(0,1)+
  xlim(0,6)+
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank())+
  theme(text= element_text(size= 29), axis.text = element_text(size = 24.5))

LC_marg_fig_HD 

ggsave("LC_MO_HD.png", last_plot(),width=6.5, height=5 )

##############################################################################
#######################################################################################################

### Plotting conditional occupancy probability
# LC and MC in presence or absence of AGC

# MC conditional on AGC present
MC_AGCp <- predict(m26, type="state", species="marbledcat", cond="goldencat", nsims=10^3)
# MC conditional on AGC absent
MC_AGCa <- predict(m26, type="state", species="marbledcat", cond="-goldencat", nsims=10^3)

# Plot conditional occupancy
cond_data_mc_agc <- rbind(MC_AGCp[1,], MC_AGCa[1,])
cond_data_mc_agc$cond <- c("AGC Present","AGC Absent")

# LC conditional on AGC present
LC_AGCp <- predict(m26, type="state", species="leopardcat", cond="goldencat", nsims=10^3)
# LC conditional on AGC absent
LC_AGCa <- predict(m26, type="state", species="leopardcat", cond="-goldencat", nsims=10^3)

# Plot conditional occupancy
cond_data_lc_agc <- rbind(LC_AGCp[1,], LC_AGCa[1,])
cond_data_lc_agc$cond <- c("AGC Present","AGC Absent")

# Plot together
cond_data_mc_agc$spec <- "Marbled cat"
cond_data_lc_agc$spec <- "Leopard cat"

tidy2 <- bind_rows(cond_data_mc_agc, cond_data_lc_agc)

e2 <- ggplot(data=tidy2, aes(x=cond, y=Predicted)) +
  geom_pointrange(aes(ymin=lower, ymax=upper, colour=spec), position=position_dodge(0.05), size=1) +
  scale_color_manual(values=c("#E7B800","#65b185"), name=NULL) + # name NULL removes legend title
  xlab("Conditional on...") +
  ylab("Predicted occupancy probability") +
  theme_bw() +
  ylim(c(0, 0.5)) +
  theme(legend.position=c(0.7, 0.9),
    legend.background=element_blank(),
    text=element_text(size=15),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank())
e2

ggsave("CO_MC_and_LC_in_presence_absence_AGC.png", e2, width=4, height=4)


# AGC and MC in presence or absence of LC

# AGC conditional on LC present
AGC_LCp <- predict(m26, type="state", species="goldencat", cond="leopardcat", nsims=10^3)
# AGC conditional on LC absent
AGC_LCa <- predict(m26, type="state", species="goldencat", cond="-leopardcat", nsims=10^3)

# Plot conditional occupancy
cond_data_agc_lc <- rbind(MC_AGCp[1,], MC_AGCa[1,])
cond_data_agc_lc$cond <- c("LC Present","LC Absent")

# MC conditional on LC present
MC_LCp <- predict(m26, type="state", species="marbledcat", cond="leopardcat", nsims=10^3)
# LC conditional on AGC absent
MC_LCa <- predict(m26, type="state", species="marbledcat", cond="-leopardcat", nsims=10^3)

# Plot conditional occupancy
cond_data_mc_lc <- rbind(LC_AGCp[1,], LC_AGCa[1,])
cond_data_mc_lc$cond <- c("LC Present","LC Absent")

# Plot together
cond_data_agc_lc$spec <- "Asiatic golden cat"
cond_data_mc_lc$spec <- "Marbled cat"

tidy2 <- bind_rows(cond_data_agc_lc, cond_data_mc_lc)

e2 <- ggplot(data=tidy2, aes(x=cond, y=Predicted)) +
  geom_pointrange(aes(ymin=lower, ymax=upper, colour=spec), position=position_dodge(0.05), size=1) +
  scale_color_manual(values=c("#e9493b", "#65b185"), name=NULL) + # name NULL removes legend title
  xlab("Conditional on...") +
  ylab("Predicted occupancy probability") +
  theme_bw() +
  ylim(c(0, 0.5)) +
  theme(legend.position=c(0.7, 0.9),
        legend.background=element_blank(),
        text=element_text(size=15),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())
e2

ggsave("CO_AGC_and_MC_in_presence_absence_LC.png", e2, width=4, height=4)


########################################################################################################
############################## CONDITIONAL OCCUPANCY PLOTS  ############################################
########################################################################################################

# Effect of A on B as a function of elevation covariate

nd_cond_ele <- data.frame(
  slope_500m=rep(mean(cov2$slope_500m), 100),
  HDens_1km=rep(mean(cov2$HDens_1km), 100),
  RivDen_2km=rep(mean(cov2$RivDen_2km), 100),
  slope_4km=rep(mean(cov2$slope_4km), 100),
  Tcover_4km=rep(mean(cov2$Tcover_4km), 100),
  HDens_4km=rep(mean(cov2$HDens_4km), 100),
  ele500_m=seq(min(cov2$ele500_m), max(cov2$ele500_m), length.out=100)
)

# AGC conditional on LC
# AGC occurrence when LC are present
AGC_cond_LC_ele <- predict(m26, type="state", species="goldencat", cond="leopardcat", newdata=nd_cond_ele, nsims=10^4)
# AGC occurrence when LC are absent
AGC_cond_noLC_ele <- predict(m26, type="state", species="goldencat", cond="-leopardcat", newdata=nd_cond_ele, nsims=10^4)

# Formating data for plotting in ggplot
gg_df_condele_agc_lc <- data.frame(
  ele500_m=nd_cond_ele$ele500_m,
  occupancy=c(AGC_cond_LC_ele$Predicted,
              AGC_cond_noLC_ele$Predicted),
  low=c(AGC_cond_LC_ele$lower,
        AGC_cond_noLC_ele$lower),
  high=c(AGC_cond_LC_ele$upper,
         AGC_cond_noLC_ele$upper),
  conditional=rep(c("leopard cat present", "leopard cat absent"), each=100)
)

# Plot
AGC_LC_cond_fig_ele <- ggplot(gg_df_condele_agc_lc, aes(x=ele500_m, y=occupancy, group=conditional)) +
  geom_ribbon(aes(ymin=low, ymax=high, fill=conditional), alpha=0.5) +
  geom_line() +
  ylab("AGC conditional occupancy") +
  xlab("Elevation") +
  ylim(0,1)+ xlim(-2,1) +
  labs(fill="") +
  theme_bw() +  
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        text=element_text(size=16),
        legend.position=c(0.35, 0.86),
        legend.background=element_blank()
  )

AGC_LC_cond_fig_ele 

ggsave("AGC_LC_CO_ele1.png", AGC_LC_cond_fig_ele, width=4, height=4)

##########################################################################################################

# Effect of B on A

# LC occurrence when AGC are present
LC_cond_AGC_ele <- predict(m26, type="state", species="leopardcat", cond="goldencat", newdata=nd_cond_ele, nsims=10^4)
# LC occurrence when AGC are absent
LC_cond_noAGC_ele <- predict(m26, type="state", species="leopardcat", cond="-goldencat", newdata=nd_cond_ele, nsims=10^4)

# Formating data for plotting in ggplot
gg_df_condele_lc_agc <- data.frame(
  ele500_m=nd_cond_ele$ele500_m,
  occupancy=c(LC_cond_AGC_ele$Predicted,
              LC_cond_noAGC_ele$Predicted),
  low=c(LC_cond_AGC_ele$lower,
        LC_cond_noAGC_ele$lower),
  high=c(LC_cond_AGC_ele$upper,
         LC_cond_noAGC_ele$upper),
  conditional=rep(c("golden cat present", "golden cat absent"), each=100)
)

# Plot
LC_AGC_cond_fig_ele <- ggplot(gg_df_condele_lc_agc, aes(x=ele500_m, y=occupancy, group=conditional)) +
  geom_ribbon(aes(ymin=low, ymax=high, fill=conditional), alpha=0.5) +
  geom_line() +
  ylab("LC conditional occupancy") +
  xlab("Elevation") +
  ylim(0,1) + xlim(-2,1) +
  labs(fill="") +
  theme_bw() +  
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        text=element_text(size=16),
        legend.position=c(0.35, 0.86),
        legend.background=element_blank()
  )

LC_AGC_cond_fig_ele 

##########################################################################################################
# But this just for one interacting pair. You have three, hence it gets tiring and boring to do the same
# for all other pairs. So, I have automated the plotting codes to produce the figure 
########################################### ELEVATION ##################################################
#------------------------------------------------------------------------------------------------------#
nd_cond_ele <- data.frame(
  slope_500m=rep(mean(cov2$slope_500m), 100),
  HDens_1km=rep(mean(cov2$HDens_1km), 100),
  RivDen_2km=rep(mean(cov2$RivDen_2km), 100),
  slope_4km=rep(mean(cov2$slope_4km), 100),
  Tcover_4km=rep(mean(cov2$Tcover_4km), 100),
  HDens_4km=rep(mean(cov2$HDens_4km), 100),
  ele500_m=seq(min(cov2$ele500_m), max(cov2$ele500_m), length.out=100)
)

MC_cond_AGC_ele <- predict(m26, type="state", species="marbledcat", cond="goldencat", newdata=nd_cond_ele, nsims=10^4)
MC_cond_noAGC_ele <- predict(m26, type="state", species="marbledcat", cond="-goldencat", newdata=nd_cond_ele, nsims=10^4)
head(MC_cond_AGC_ele)
head(MC_cond_noAGC_ele)

LC_cond_AGC_ele <- predict(m26, type="state", species="leopardcat", cond="goldencat", newdata=nd_cond_ele, nsims=10^4)
LC_cond_noAGC_ele <- predict(m26, type="state", species="leopardcat", cond="-goldencat", newdata=nd_cond_ele, nsims=10^4)

AGC_cond_MC_ele <- predict(m26, type="state", species="goldencat", cond="marbledcat", newdata=nd_cond_ele, nsims=10^4)
AGC_cond_noMC_ele <- predict(m26, type="state", species="goldencat", cond="-marbledcat", newdata=nd_cond_ele, nsims=10^4)

LC_cond_MC_ele <- predict(m26, type="state", species="leopardcat", cond="marbledcat", newdata=nd_cond_ele, nsims=10^4)
LC_cond_noMC_ele <- predict(m26, type="state", species="leopardcat", cond="-marbledcat", newdata=nd_cond_ele, nsims=10^4)

AGC_cond_LC_ele <- predict(m26, type="state", species="goldencat", cond="leopardcat", newdata=nd_cond_ele, nsims=10^4)
AGC_cond_noLC_ele <- predict(m26, type="state", species="goldencat", cond="-leopardcat", newdata=nd_cond_ele, nsims=10^4)

MC_cond_LC_ele <- predict(m26, type="state", species="marbledcat", cond="leopardcat", newdata=nd_cond_ele, nsims=10^4)
MC_cond_noLC_ele <- predict(m26, type="state", species="marbledcat", cond="-leopardcat", newdata=nd_cond_ele, nsims=10^4)

# You might want to edit the blank plots (the orange lines - see after plotting) later after you save them.

# Elevation values to plot
ele500m = seq(min(cov2$ele500_m), max(cov2$ele500_m), length.out=100) # standardised values
#ele500m = seq(min(cov$ele500), max(cov$ele500), length.out=100) # original values

# Plot format #dummy plot for dumb girl

#-------#-------#-------#
#  spA  #  spB  #  spC  #         
#-------#-------#-------#-------#
#       #       #       #       #
# blank #  BA   #  CA   # condA #
#       #       #       #       #
#-------#-------#-------#-------#
#       #       #       #       #
#  AB   # blank #  CB   # condB #
#       #       #       #       #
#-------#-------#-------#-------#
#       #       #       #       #
#  AC   #  BC   # blank # condC #
#       #       #       #       #
#-------#-------#-------#-------#

# create a blank plot
blank1 <- MC_cond_AGC_ele
blank1 <- replace(blank1, blank1>0, 0)

# The diagonal plots are all blank - where you will edit and insert animal images
blank2 <- blank3 <- blank1

######### 1st ROW ###########

########### BLANK #####################
blank1$cat <- "pra_blank"             #
blank1$species <- "Asian golden cat"  #
blank1$cat2 <- "condA"                #
blank1$check <- "Absent"              #
blank1$x <- ele500m                   #
#######################################

pm_na <- MC_cond_noAGC_ele
pm_na$cat <- "prm_noagc"
pm_na$species <- "Marbled cat"
pm_na$cat2 <- "condA"
pm_na$check <- "Absent"
pm_na$x <- ele500m

pm_a <- MC_cond_AGC_ele
pm_a$cat <- "prm_agc"
pm_a$species <- "Marbled cat"
pm_a$cat2 <- "condA"
pm_a$check <- "Present"
pm_a$x <- ele500m

pl_na <- LC_cond_noAGC_ele
pl_na$cat <- "prl_noagc"
pl_na$species <- "Leopard cat"
pl_na$cat2 <- "condA"
pl_na$check <- "Absent"
pl_na$x <- ele500m

pl_a <- LC_cond_AGC_ele
pl_a$cat <- "prl_agc"
pl_a$species <- "Leopard cat"
pl_a$cat2 <- "condA"
pl_a$check <- "Present"
pl_a$x <- ele500m

######### 2nd ROW ###########

pa_nm <- AGC_cond_noMC_ele
pa_nm$cat <- "pra_nomc"
pa_nm$species <- "Asian golden cat"
pa_nm$cat2 <- "condM"
pa_nm$check <- "Absent"
pa_nm$x <- ele500m

pa_m <- AGC_cond_MC_ele
pa_m$cat <- "pra_mc"
pa_m$species <- "Asian golden cat"
pa_m$cat2 <- "condM"
pa_m$check <- "Present"
pa_m$x <- ele500m

########### BLANK ################
blank2$cat <- "prm_blank"        #
blank2$species <- "Marbled cat"  #
blank2$cat2 <- "condM"           #
blank2$check <- "Absent"         #
blank2$x <- ele500m              #
##################################

pl_nm <- LC_cond_noMC_ele
pl_nm$cat <- "prl_nmc"
pl_nm$species <- "Leopard cat"
pl_nm$cat2 <- "condM"
pl_nm$check <- "Absent"
pl_nm$x <- ele500m

pl_m <- LC_cond_MC_ele
pl_m$cat <- "prl_mc"
pl_m$species <- "Leopard cat"
pl_m$cat2 <- "condM"
pl_m$check <- "Present"
pl_m$x <- ele500m

######### 3rd ROW ###########

pa_nl <- AGC_cond_noLC_ele
pa_nl$cat <- "pra_nolc"
pa_nl$species <- "Asian golden cat"
pa_nl$cat2 <- "condL"
pa_nl$check <- "Absent"
pa_nl$x <- ele500m

pa_l <- AGC_cond_LC_ele
pa_l$cat <- "pra_lc"
pa_l$species <- "Asian golden cat"
pa_l$cat2 <- "condL"
pa_l$check <- "Present"
pa_l$x <- ele500m

pm_nl <- MC_cond_noLC_ele
pm_nl$cat <- "prm_nolc"
pm_nl$species <- "Marbled cat"
pm_nl$cat2 <- "condL"
pm_nl$check <- "Absent"
pm_nl$x <- ele500m

pm_l <- MC_cond_LC_ele
pm_l$cat <- "prm_lc"
pm_l$species <- "Marbled cat"
pm_l$cat2 <- "condL"
pm_l$check <- "Present"
pm_l$x <- ele500m

########### BLANK #################
blank3$cat <- "prl_blank"         #
blank3$species <- "Leopard cat"   #
blank3$cat2 <- "condL"            #
blank3$check <- "Absent"          #
blank3$x <- ele500m               #
###################################

# Bind all of them into a data frame
cond_psi_elev <- rbind(blank1, pm_na, pm_a, pl_na, pl_a, 
                       pa_nm, pa_m, blank2, pl_nm, pl_m,
                       pa_nl, pa_l, pm_nl, pm_l, blank3)

# Rename for convenience
dat <- cond_psi_elev

# New columns with species names and labels
dat$species2 <- factor(dat$species, levels=c("Asian golden cat", "Marbled cat", "Leopard cat"))
dat$cat2A <- factor(dat$cat2, levels=c('condA', 'condM', 'condL'),
                    labels=c('Conditional on\nAsiatic golden cat', 'Conditional on\nMarbled cat', 'Conditional on\nLeopard cat'))

# Rename for the facet grid
cat2_names <- c(
  'Asian golden cat'='Asiatic golden cat...',
  'Marbled cat'='Marbled cat...',
  'Leopard cat'='Leopard cat...'
)

# Plot
p <- ggplot(dat, aes(x=x, y=Predicted)) +
  ylim(0, 1) +
  geom_ribbon(data=dat, aes(ymin=lower, ymax=upper, group=check, fill=check), alpha=0.4) +
  geom_line(aes(color=check), lwd=0.8, na.rm=T) +
  scale_fill_manual(values=c('#9932c3', '#e8a543')) +
  scale_color_manual(values=c('#9932c3', '#e8a543')) +
  xlab("Elevation (standardised)") +
  ylab('Occupancy probability') +
  facet_grid(cat2A ~ species2, 
             labeller=labeller(species2=as_labeller(cat2_names))) +
  theme_bw() +
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=12),
        strip.text=element_text(size=12),
        legend.title=element_blank(),
        legend.text=element_text(size=12),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_rect(fill=adjustcolor('transparent', 0)),
        strip.background=element_rect(fill=adjustcolor('grey', 0.5)))

p

ggsave('conditional_plot_ELEV_KC.png', p, width=7, height=7)

########################################################################################################

######################################## HOUSING DENSITY ###############################################
#------------------------------------------------------------------------------------------------------#

nd_cond_hden <- data.frame(
  slope_500m=rep(mean(cov2$slope_500m), 100),
  HDens_1km=rep(mean(cov2$HDens_1km), 100),
  RivDen_2km=rep(mean(cov2$RivDen_2km), 100),
  slope_4km=rep(mean(cov2$slope_4km), 100),
  Tcover_4km=rep(mean(cov2$Tcover_4km), 100),
  ele500_m=rep(mean(cov2$ele500_m), 100),
  HDens_4km=seq(min(cov2$HDens_4km), max(cov2$HDens_4km), length.out=100)
)

AGC_cond_LC_hden <- predict(m26, type="state", species="goldencat", cond="leopardcat", newdata=nd_cond_hden, nsims=10^4)
AGC_cond_noLC_hden <- predict(m26, type="state", species="goldencat", cond="-leopardcat", newdata=nd_cond_hden, nsims=10^4)

AGC_cond_MC_hden <- predict(m26, type="state", species="goldencat", cond="marbledcat", newdata=nd_cond_hden, nsims=10^4)
AGC_cond_noMC_hden <- predict(m26, type="state", species="goldencat", cond="-marbledcat", newdata=nd_cond_hden, nsims=10^4)

LC_cond_MC_hden <- predict(m26, type="state", species="leopardcat", cond="marbledcat", newdata=nd_cond_hden, nsims=10^4)
LC_cond_noMC_hden <- predict(m26, type="state", species="leopardcat", cond="-marbledcat", newdata=nd_cond_hden, nsims=10^4)

LC_cond_AGC_hden <- predict(m26, type="state", species="leopardcat", cond="goldencat", newdata=nd_cond_hden, nsims=10^4)
LC_cond_noAGC_hden <- predict(m26, type="state", species="leopardcat", cond="-goldencat", newdata=nd_cond_hden, nsims=10^4)

MC_cond_AGC_hden <- predict(m26, type="state", species="marbledcat", cond="goldencat", newdata=nd_cond_hden, nsims=10^4)
MC_cond_noAGC_hden <- predict(m26, type="state", species="marbledcat", cond="-goldencat", newdata=nd_cond_hden, nsims=10^4)

MC_cond_LC_hden <- predict(m26, type="state", species="marbledcat", cond="leopardcat", newdata=nd_cond_hden, nsims=10^4)
MC_cond_noLC_hden <- predict(m26, type="state", species="marbledcat", cond="-leopardcat", newdata=nd_cond_hden, nsims=10^4)

# Covariate values
hden = seq(min(cov2$HDens_4km), max(cov2$HDens_4km), length.out= 100) # standardised
hden = seq(min(cov$setden_4KM), max(cov$setden_4KM), length.out= 100) 

# create a blank plot
blank1 <- MC_cond_AGC_hden
blank1 <- replace(blank1, blank1>0, 0)

blank2 <- blank3 <- blank1

######### 1st ROW #####################
blank1$cat <- "pra_blank"             #
blank1$species <- "Asian golden cat"  #
blank1$cat2 <- "condA"                #
blank1$check <- "Absent"              #
blank1$x <- hden                      #
#######################################

pm_na <- MC_cond_noAGC_hden
pm_na$cat <- "prm_noagc"
pm_na$species <- "Marbled cat"
pm_na$cat2 <- "condA"
pm_na$check <- "Absent"
pm_na$x <- hden

pm_a <- MC_cond_AGC_hden
pm_a$cat <- "prm_agc"
pm_a$species <- "Marbled cat"
pm_a$cat2 <- "condA"
pm_a$check <- "Present"
pm_a$x <- hden

pl_na <- LC_cond_noAGC_hden
pl_na$cat <- "prl_noagc"
pl_na$species <- "Leopard cat"
pl_na$cat2 <- "condA"
pl_na$check <- "Absent"
pl_na$x <- hden

pl_a <- LC_cond_AGC_hden
pl_a$cat <- "prl_agc"
pl_a$species <- "Leopard cat"
pl_a$cat2 <- "condA"
pl_a$check <- "Present"
pl_a$x <- hden

######### 2nd ROW ###########

pa_nm <- AGC_cond_noMC_hden
pa_nm$cat <- "pra_nomc"
pa_nm$species <- "Asian golden cat"
pa_nm$cat2 <- "condM"
pa_nm$check <- "Absent"
pa_nm$x <- hden

pa_m <- AGC_cond_MC_hden
pa_m$cat <- "pra_mc"
pa_m$species <- "Asian golden cat"
pa_m$cat2 <- "condM"
pa_m$check <- "Present"
pa_m$x <- hden

###################################
blank2$cat <- "prm_blank"         #
blank2$species <- "Marbled cat"   #
blank2$cat2 <- "condM"            #
blank2$check <- "Absent"          #
blank2$x <- hden                  #
###################################

pl_nm <- LC_cond_noMC_hden
pl_nm$cat <- "prl_nmc"
pl_nm$species <- "Leopard cat"
pl_nm$cat2 <- "condM"
pl_nm$check <- "Absent"
pl_nm$x <- hden

pl_m <- LC_cond_MC_hden
pl_m$cat <- "prl_mc"
pl_m$species <- "Leopard cat"
pl_m$cat2 <- "condM"
pl_m$check <- "Present"
pl_m$x <- hden

######### 3rd ROW ###########

pa_nl <- AGC_cond_noLC_hden
pa_nl$cat <- "pra_nolc"
pa_nl$species <- "Asian golden cat"
pa_nl$cat2 <- "condL"
pa_nl$check <- "Absent"
pa_nl$x <- hden

pa_l <- AGC_cond_LC_hden
pa_l$cat <- "pra_lc"
pa_l$species <- "Asian golden cat"
pa_l$cat2 <- "condL"
pa_l$check <- "Present"
pa_l$x <- hden

pm_nl <- MC_cond_noLC_hden
pm_nl$cat <- "prm_nolc"
pm_nl$species <- "Marbled cat"
pm_nl$cat2 <- "condL"
pm_nl$check <- "Absent"
pm_nl$x <- hden

pm_l <- MC_cond_LC_hden
pm_l$cat <- "prm_lc"
pm_l$species <- "Marbled cat"
pm_l$cat2 <- "condL"
pm_l$check <- "Present"
pm_l$x <- hden

###################################
blank3$cat <- "prl_blank"         #
blank3$species <- "Leopard cat"   #
blank3$cat2 <- "condL"            #
blank3$check <- "Absent"          #
blank3$x <- hden                  #
###################################

cond_psi_hden <- rbind(blank1, pm_na, pm_a, pl_na, pl_a, 
                       pa_nm, pa_m, blank2, pl_nm, pl_m,
                       pa_nl, pa_l, pm_nl, pm_l, blank3)

dat_hden <- cond_psi_hden

dat_hden$species2 <- factor(dat_hden$species, levels = c("Asian golden cat", "Marbled cat", "Leopard cat"))
dat_hden$cat2A <- factor(dat_hden$cat2, levels=c('condA', 'condM', 'condL'),
                         labels=c('Conditional on\nAsiatic golden cat', 'Conditional on\nMarbled cat', 'Conditional on\nLeopard cat'))
cat2_names <- c(
  'Asian golden cat' = 'Asiatic golden cat...',
  'Marbled cat' = 'Marbled cat...',
  'Leopard cat' = 'Leopard cat...'
)


p1 <- ggplot(dat_hden, aes(x=x, y=Predicted)) +
  ylim(0, 1) +
  geom_ribbon(data=dat_hden, aes(ymin=lower, ymax=upper, group=check, fill=check), alpha=0.4) +
  geom_line(aes(color=check), lwd=0.8, na.rm=T) +
  scale_fill_manual(values=c( '#038f3c', '#e79612')) +
  scale_color_manual(values=c( '#038f3c', '#e79612')) +
  xlab("Housing density") +
  ylab('Occupancy probability') +
  facet_grid(cat2A ~ species2, 
             labeller=labeller(species2 = as_labeller(cat2_names))) +
  theme_bw() +
  theme(axis.title = element_text(size=18),
        axis.text = element_text(size=12),
        strip.text = element_text(size=12),
        legend.title = element_blank(),
        legend.text = element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill=adjustcolor('transparent', 0)),
        strip.background = element_rect(fill=adjustcolor('grey', 0.5)))

p1

ggsave('conditional_plot_HDENS_KC.png', p1, width=7, height=7)

