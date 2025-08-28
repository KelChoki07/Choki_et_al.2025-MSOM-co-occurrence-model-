
#   _____  _             _     _          _            _      ___    ___  ___   _____ 
# / ____|| |           | |   (_)        | |          | |    |__ \  / _ \|__ \ | ____|
#| |     | |__    ___  | | __ _     ___ | |_    __ _ | |       ) || | | |  ) || |__  
#| |     | '_ \  / _ \ | |/ /| |   / _ \| __|  / _` || |      / / | | | | / / |___ \ 
#| |____ | | | || (_) ||   < | |  |  __/| |_  | (_| || | _   / /_ | |_| |/ /_  ___) |
# \_____||_| |_| \___/ |_|\_\|_|   \___| \__|  \__,_||_|(_) |____| \___/|____||____/ 
                                                                                     
                                                                                     

# Publication details                                                                                                          
# Title: # Anthropogenic and environmental correlates of spatial patterns of co-occurrence of small felids in a montane landscape                                    #
# Authors: Karma Choki, Egil Dröge, Claudio Sillero-Zubiri, David W. Macdonald, Ugyen Penjor                 
# Journal: Global Ecology and Conservation, Volume 58, 2025                                                 
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

##########################################


ggsave("Marginal_occ_det_prob_all_sps_1.png", k, width=6.9, height=4.5)
#######################################################################################################

