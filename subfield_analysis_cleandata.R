# Nichole Bouffard, Aug 2023
# Threshold 0.4 - The happy medium between liberal and conservative
# AC run on ICA-FIX clean data from HCP
#
# Load libraries
library(tidyverse)
library(afex)
library(emmeans)
library(effectsize)


setwd('~/')

######################################
######## PLOTS AND STATS #############
######################################

### Read in all values
master <- read.csv('/data/allsubs_Subfield_regthresh0.4_AC_Values_cleandata.csv')
master <- master %>%  select(-X)

master$participant <- as.factor(master$participant)
#master$roi <- factor(master$roi, levels = c("subiculum", "CA1", "CA2CA3","CA4DG"))
master$roi <- factor(master$roi, levels = c("subiculum", "CA1", "CA2CA3","CA4DG","SRLM"))
master$hemi <- as.factor(master$hemi)

# Gather into long format
fulldata <- master %>% 
  gather('lag', 'ac',4:7)

# Normalize the AC values 
fulldata <- fulldata %>% 
  group_by(participant) %>% 
  mutate(zAC = scale(ac))


# ROI AC means and SD collapsed across lag
fulldata_summary <- fulldata %>% 
  group_by(hemi,roi) %>% 
  summarize(mean = mean(zAC), sd= sd(zAC))

# Just lag1 FOR ANALYSES
fulldata_lag1 <- fulldata %>% 
  filter(lag == 'lag1')

# Models
# Model 1: Effect of ROI, Lag, and Hemi and sig interactions (zAC gives isSingular)
# mdl1 <- mixed(ac ~ roi*lag*hemi + (1|participant), data = fulldata)
# Model 2: Effect of ROI and Hemi (zAC gives isSingular)
mdl2 <- mixed(ac ~ roi*hemi + (1|participant), data = fulldata_lag1)

# Thresh 0.4 - Main effect of ROI: subiculum > CA1 + CA2CA3 + CA4DG + SRLM; CA1 < subiculum + CA2/3 + CA4DG + SRLM; CA2CA3 > CA4DG
ems.condition <- emmeans(mdl2, c('roi'))
pairs(ems.condition)

# Thresh 0.4 - Main effect of Hemisphere: Right > Left
ems.condition <- emmeans(mdl2, c('hemi'))
pairs(ems.condition)

# Thresh 0.4 - Main effect of Lag: Sig diff between each lag
#ems.condition <- emmeans(mdl1, c('lag'))
#pairs(ems.condition)


# Interaction of ROI x Hemi: Sig diff between each lag
ems.condition <- emmeans(mdl1, c('roi','hemi'))
pairs(ems.condition)

# ROI x lag
# lag x Hemi


################# PLOTS FOR POSTER #################
hem.labs <- c("Left hippocampus", "Right hippocampus")
names(hem.labs) <- c("L", "R")


# Plot for Paper Fig1

plotV3<-fulldata_lag1 %>% 
  ggplot(aes(x=roi, y= zAC, color = roi)) +
  stat_summary(fun.data = mean_se, geom = "linerange", position = position_dodge(width = .9), linewidth = 2) +
  stat_summary(fun = "mean", geom = "point", position = position_dodge(width = .9), size = 4) +
  facet_wrap(~hemi, strip.position = 'top',labeller = labeller(hemi = hem.labs), scales = 'free_x') +
  ylab('Average Single Voxel Autocorrelation (Z)') +
  scale_color_manual(breaks = c("subiculum","CA1", "CA2CA3", "CA4DG","SRLM"),
                    values=c("#3E02FF", "#63FEFF", "#4EFF00","#F11B01","#F7C670")) +
  theme_light(base_size = 24) +
  theme(legend.title = element_blank(), legend.position="bottom", legend.box = "horizontal")+
  theme(strip.background = element_blank(),strip.text.x = element_text(colour = "black",size=24),strip.text=element_text(hjust=0), axis.title.x = element_blank(),legend.key = element_rect(colour = "transparent", fill = "transparent"))+
  labs(color='Condition') 

ggsave(plotV3, width = 12, height = 7.5, file='Fig1.png',dpi=300)


