# Aug 2023
# Nichole Bouffard 
# This script reads in the compiled data, separate csv files for X and Y axes. 
# Xaxis = Posterior-anterior, Yaxis = medial-lateral
#
# Load libraries
library(tidyverse)
library(afex)
library(emmeans)
library(ggpubr)
library(ggbeeswarm)
library(R.matlab)

setwd('~/data/')

masterXaxis <- read.table("allsubs_allrois_Xaxis.csv",header = TRUE, sep = ",")
masterYaxis <-read.table("allsubs_allrois_Yaxis.csv",header = TRUE, sep = ",")

masterXaxis$subj <- as.factor(masterXaxis$subj)
masterYaxis$subj <- as.factor(masterYaxis$subj)

xaxis <- masterXaxis %>% 
  group_by(subj) %>% 
  mutate(zAC = scale(avgAC))
yaxis <- masterYaxis %>% 
  group_by(subj) %>% 
  mutate(zAC = scale(avgAC))

# *I need to renumber the medial lateral axis, in one hemisphere the most medial slice is 0 and in the other hemisphere it is the most lateral slice

############################## TRIMMING SLICES ######################################
# Ant/Post Med/Lat CA1 and Subiculum Analysis

# Figure out how to cut off anterior slices that only a few subjects have. 
# Only include slices that at least 75% of subs have (~18 subs)

xaxisSubCA1 <- xaxis %>% 
  filter(roi == 'subiculum' | roi == 'CA1')
yaxisSubCA1 <- yaxis %>% 
  filter(roi == 'subiculum' | roi == 'CA1')


# POSTERIOR-ANTERIOR
xaxisSubCA1 %>% 
  ggplot(aes(x = slice))+
  geom_histogram(stat='count')+
  facet_wrap(~roi)
# Cut off in this plot is 18*2 = 36
# CA1 cut off at slice 21 (L+R), any slice greater than 21 is not being included because not enough subjects have slices
# Subiculum cut off at slice 19 (L+R), any slice greater than 19 is not being included because not enough subjects have slices

xaxisSubCA1short <- xaxisSubCA1 %>%
  filter(ifelse(roi == 'CA1', !slice %in% c('22','23','24','25'),  !slice %in% c('20','21','22','23','24','25')))
# Sanity check
xaxisSubCA1short %>% 
  ggplot(aes(x = slice))+
  geom_histogram(stat='count')+
  facet_wrap(~roi)


# MEDIAL_LATERAL
yaxisSubCA1 %>% 
  ggplot(aes(x = slice))+
  geom_histogram(stat='count')+
  facet_wrap(hemi~roi)
# Cut off in this plot is (18*2 = 36)
# CA1 cut off at slice 13 (L+R), any slice greater than 13 is not being included because not enough subjects have slices
# Subiculum cut off at slice 8 (L+R), any slice greater than 8 is not being included because not enough subjects have slices

yaxisSubCA1short <- yaxisSubCA1 %>% 
  filter(ifelse(roi == 'CA1', slice %in% c('1','2','3','4','5','6','7','8','9','10','11','12','13'), slice %in% c('1','2','3','4','5','6','7','8')))
# Sanity check
yaxisSubCA1short %>% 
  ggplot(aes(x = slice))+
  geom_histogram(stat='count')+
  facet_wrap(~roi)



# CA2CA3, CA4DG, SRLM TRIMMING
# cut off anterior slices that only a few subjects have. 
# Only include slices that at least 75% of subs have (~18 subs)

xaxisOther <- xaxis %>% 
  filter(roi == 'CA2CA3' | roi == 'CA4DG'| roi == 'SRLM')
yaxisOther <- yaxis %>% 
  filter(roi == 'CA2CA3' | roi == 'CA4DG'| roi == 'SRLM')


# POSTERIOR-ANTERIOR
xaxisOther %>% 
  ggplot(aes(x = slice))+
  geom_histogram(stat='count')+
  facet_wrap(~roi)
# Cut off in this plot is (18*2 = 36)
# CA2CA3 cut off at slice 18 (L+R), any slice greater than 18 is not being included because not enough subjects have slices
# CA4DG cut off at slice 16 (L+R), any slice greater than 16 is not being included because not enough subjects have slices
# SRLM cut off at slice 19 (L+R), any slice greater than 19 is not being included because not enough subjects have slices

xaxisOthershort <- xaxisOther %>%
  filter(ifelse(roi == 'CA2CA3', !slice %in% c('19','20','21'),
                ifelse(roi == 'CA4DG', !slice %in% c('17','18'), !slice %in% c('20','21', '22'))))
# Sanity check
xaxisOthershort %>% 
  ggplot(aes(x = slice))+
  geom_histogram(stat='count')+
  facet_wrap(~roi)


# MEDIAL_LATERAL
yaxisOther %>% 
  ggplot(aes(x = slice))+
  geom_histogram(stat='count')+
  facet_wrap(hemi~roi)
# CA2CA3 cut off at slice 12 (L+R), any slice greater than 12 is not being included because not enough subjects have slices
# CA4DG cut off at slice 8 (L+R), any slice greater than 8 is not being included because not enough subjects have slices
# SRLM cut off at slice 11 (L+R), any slice greater than 11 is not being included because not enough subjects have slices

yaxisOthershort <- yaxisOther %>% 
  filter(ifelse(roi == 'CA2CA3', !slice %in% c('13','14','15', '16'),
                ifelse(roi == 'CA4DG', !slice %in% c('9','10','11','12'), !slice %in% c('12','13','14','15'))))
# Sanity check
yaxisOthershort %>% 
  ggplot(aes(x = slice))+
  geom_histogram(stat='count')+
  facet_wrap(~roi)


#### COMBINE ALL TRIMMED SUBFIELDS ####

xaxisTrimmed <- rbind(xaxisSubCA1short,xaxisOthershort)
yaxisTrimmed <- rbind(yaxisSubCA1short,yaxisOthershort)

################################################################
###################### PLOTS FOR PAPER #########################
################################################################
setwd('~/data/')


###################### ALL SUBFIELDS TRIMMED SLICES #########################
# Plot Fig 2A
hem.labs <- c("Left hippocampus", "Right hippocampus")
names(hem.labs) <- c("L", "R")

xaxisPlotShort <- xaxisTrimmed%>% 
  ggplot(aes(x=slice, y=zAC, color=roi, group = roi)) +
  stat_summary(geom="line", fun = "mean", linewidth =1.5)+
  stat_summary(inherit.aes = FALSE, aes(x=slice, y=zAC), geom="line", fun = "mean", color="black", linewidth=1.5, linetype="dashed")+
  facet_wrap(~hemi,labeller = labeller(hemi = hem.labs)) +
  theme_bw(base_size =20)+
  scale_color_manual(breaks = c("subiculum","CA1", "CA2CA3", "CA4DG","SRLM"),
                    values=c("#3E02FF", "#63FEFF", "#4EFF00","#F11B01","#F7C670")) +
  # theme(legend.position = 'none')+
  ylab('Average Single Voxel Autocorrelation (Z)') +
  ggtitle('Posterior-anterior axis') 

ggsave(xaxisPlotShort, width = 11, height = 6, file='Fig2A.png',dpi=300)


####### MEDIAL LATERAL WITH LEFT AND RIGHT FLIPPED #####

# Plot Fig 2B

# Flipping left hemisphere by multiplying by -1
yaxisTrimmedFlipped <- yaxisTrimmed %>% 
  mutate(slice_flip = slice) %>% 
  mutate(slice_flip = ifelse(hemi == 'L', slice*-1, slice))

yaxisPlotShort <-yaxisTrimmedFlipped%>% 
  ggplot(aes(x=slice_flip, y=zAC, color=roi, group = roi)) +
  stat_summary(geom="line", fun = "mean",linewidth =1.5)+
  stat_summary(inherit.aes = FALSE, aes(x=slice_flip, y=zAC), geom="line", fun = "mean", color="black", linewidth=1.5, linetype="dashed")+
  facet_wrap(~hemi,labeller = labeller(hemi = hem.labs), scales = 'free_x') +
  theme_bw(base_size =20)+
  scale_color_manual(breaks = c("subiculum","CA1", "CA2CA3", "CA4DG","SRLM"),
                     values=c("#3E02FF", "#63FEFF", "#4EFF00","#F11B01","#F7C670")) +
  # theme(legend.position = 'none')+
  ylab('Average Single Voxel Autocorrelation (Z)') +
  ggtitle('Medial-lateral axis')

ggsave(yaxisPlotShort, width = 11, height = 6, file='Fig2B.png',dpi=300)
