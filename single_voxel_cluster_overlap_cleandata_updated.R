# Nichole Bouffard, Aug 2023
# This script evaluates the amount of overlap between the functional space subfields and the AC clusters
# Registration threshold 0.4 
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
master <- read.csv('/data/allsubs_single_voxel_cluster_comparison_values_cleandata_updated.csv')
master <- master %>%  select(-X)

master$participant <- as.factor(master$participant)
master$voxel <- as.factor(master$voxel)
master$cluster_code <- as.factor(master$cluster_code)
master$subfield_code <- as.factor(master$subfield_code)
master$hemi <- as.factor(master$hemi)


# Update cluster codes
fulldata <- master %>% 
  mutate(cluster = ifelse(cluster_code == '1', 'posterior_lateral',
                          ifelse(cluster_code == '2', 'intermediate', 'anterior_medial'))) %>% 
  mutate(subfield = ifelse(subfield_code == '1', 'subiculum',
                    ifelse(subfield_code == '2', 'CA1', 
                    ifelse(subfield_code == '3', 'CA2CA3',
                    ifelse(subfield_code == '4', 'CA4DG', 
                    ifelse(subfield_code == '5', 'SRLM',
                    ifelse(subfield_code == '0', 'unlabeled', 'combo')))))))


# Normalize the AC values 
# I think grouping by participant first and z scoring within each person makes sense...but I'm not sure
fulldata <- fulldata %>% 
  group_by(participant) %>% 
  mutate(zAC = scale(ac))


# PLOTS FOR PAPER - FIG 3D
hem.labs <- c("Left hippocampus", "Right hippocampus")
names(hem.labs) <- c("L", "R")
lab.bottom <- c("Cluster 1\n(Anterior-medial)","Cluster 2\n(Intermediate)", "Cluster 3\n(Posterior-lateral)")
names(lab.bottom) <- c("anterior_medial", "intermediate","posterior_lateral")

clusterPlot <- fulldata %>% 
  group_by(participant, cluster, hemi) %>% 
  summarise(mean_zAC = mean(zAC)) %>% 
  ggplot(aes(x=cluster, y=mean_zAC, color = cluster)) +
  geom_jitter(aes(group = cluster),color = 'lightgrey', alpha = .5, position = position_jitter(width = .2), size = 4)+
  #geom_point(aes(group = cluster),color = 'lightgrey', position = position_dodge(width = .9), size = 2.5)+
  stat_summary(fun = "mean", geom = "pointrange", position = position_dodge(width = .9), size = 1.5) +
  stat_summary(fun.data = mean_se, geom = "errorbar", position = position_dodge(width = .9), linewidth = 1.5, width = .3) +
  scale_color_manual(breaks = c("anterior_medial","intermediate", "posterior_lateral"),
                     values=c("#FFD700", "#F4A92A","#F15620")) +
  scale_x_discrete(labels = lab.bottom)+
  theme_bw(base_size = 36)+
  theme(legend.position = 'none')+
  ylab('Average Single Voxel Autocorrelation (Z)')+
  facet_wrap(~hemi, labeller = labeller(hemi = hem.labs))

ggsave(clusterPlot, width = 21, height = 12, file='Fig3D.png',dpi=300)



##########################################################################################
##################### Proportion of overlap with AC clusters #############################
##########################################################################################

# Compute proportions of voxels 
fulldata_cluster_summary <- fulldata %>% 
  group_by(participant, hemi, cluster_code) %>% 
  summarise(total_n_vox = n())
fulldata_subfield_summary <- fulldata %>% 
  group_by(participant, hemi, cluster_code, subfield_code) %>% 
  summarise(n_vox_subfield = n())

# sanity check
test <- fulldata %>% 
  filter(participant == '126426') %>% 
  filter(hemi == 'L') %>% 
  filter(cluster_code == '1') %>% 
  filter(subfield_code == '0')

# Combine summary dfs to compute proportions
fulldata_summary <- left_join(fulldata_subfield_summary, fulldata_cluster_summary, by = c('participant', 'hemi', 'cluster_code'))

fulldata_summary <- fulldata_summary %>% 
  mutate(proportion = n_vox_subfield/total_n_vox)

# sanity check
test <- fulldata_summary %>% 
  filter(participant == '130114') %>% 
  filter(hemi == 'L') %>% 
  filter(cluster_code == '1') %>% 
  summarise(sumVox = sum(n_vox_subfield), sumprop = sum(proportion))


# Update cluster codes
fulldata_summary <- fulldata_summary %>% 
  mutate(cluster = ifelse(cluster_code == '1', 'posterior_lateral',
                   ifelse(cluster_code == '2', 'intermediate', 'anterior_medial'))) %>% 
  mutate(subfield = ifelse(subfield_code == '1', 'subiculum',
                    ifelse(subfield_code == '2', 'CA1', 
                    ifelse(subfield_code == '3', 'CA2CA3',
                    ifelse(subfield_code == '4', 'CA4DG', 
                    ifelse(subfield_code == '5', 'SRLM',
                    ifelse(subfield_code == '0', 'unlabeled', 'combo')))))))



# FOR PAPER - All Subfields - Stacked Barplot - Fig3B
hem.labs <- c("Left hippocampus", "Right hippocampus")
names(hem.labs) <- c("L", "R")
lab.bottom <- c("Anterior- \nmedial","Intermediate", "Posterior-\nlateral")
names(lab.bottom) <- c("anterior_medial", "intermediate","posterior_lateral")

plotStacked <- fulldata_summary %>% 
  filter(subfield != 'combo') %>% 
  filter(subfield != 'unlabeled') %>% 
  #filter(cluster != 'intermediate') %>% 
  ggplot(aes(x=cluster, y= proportion, fill = subfield)) +
  geom_bar(position="fill", stat="identity") +
  ylab('Proportion of overlap with autocorrelation clusters') +
  scale_fill_manual(breaks = c("subiculum","CA1", "CA2CA3", "CA4DG","SRLM"),
                     values=c("#3E02FF", "#63FEFF", "#4EFF00","#F11B01","#F7C670")) +
  scale_x_discrete(labels = lab.bottom)+
  theme_bw(base_size = 60)+
 # theme_minimal()+
  theme(axis.title.x = element_blank())+
  facet_wrap(~hemi, strip.position = 'top', labeller = labeller(hemi = hem.labs))


ggsave(plotStacked, width = 32, height = 19, file='Fig3B.png',dpi=300)


##########################################################################################
################## Stats: Proportion of Subfields in Clusters ############################
##########################################################################################

# Models
# Model 1: Proportion of subfields - Effect of subfield and Hemi with cluster as a nested variable? (gives isSingular)
mdl1 <- mixed(proportion ~ subfield*hemi + (1|participant) + (1|cluster), data = fulldata_summary)
mdl1 <- mixed(proportion ~ cluster*subfield*hemi + (1|participant), data = fulldata_summary)

goodsubfields_only <- fulldata_summary %>% 
  filter(subfield != 'unlabeled') %>% 
  filter(subfield != 'combo')

# Gives isSingular
mdl2 <- mixed(proportion ~ cluster*subfield*hemi + (1|participant), data = goodsubfields_only)

ems.condition <- emmeans(mdl2, c('subfield','hemi','cluster'))
# RIGHT HEM
rightHem <- subset(ems.condition, hemi == 'R')
pairs(rightHem)
# LEFT HEM
leftHem <- subset(ems.condition, hemi == 'L')
pairs(leftHem)


# Breaking up antmed and postlat clusters
# Just antmed clusters
fulldata_antmed_only <- fulldata_summary %>% 
  filter(cluster == 'anterior_medial') %>% 
  filter(subfield != 'unlabeled') %>% 
  filter(subfield != 'combo')
# Model 2: Anterior-medial HPC proportion of subfields (gives isSingular)
mdl2 <- mixed(proportion ~ subfield*hemi + (1|participant), data = fulldata_antmed_only)

ems.condition <- emmeans(mdl2, c('subfield','hemi'))
# RIGHT HEM
rightHem <- subset(ems.condition, hemi == 'R')
pairs(rightHem)
# LEFT HEM
leftHem <- subset(ems.condition, hemi == 'L')
pairs(leftHem)

# Just postlat clusters
fulldata_postlat_only <- fulldata_summary %>% 
  filter(cluster == 'posterior_lateral') %>% 
  filter(subfield != 'unlabeled') %>% 
  filter(subfield != 'combo')
# Model 3: Posterior-lateral HPC proportion of subfields (zAC gives isSingular)
mdl3 <- mixed(proportion ~ subfield*hemi + (1|participant), data = fulldata_postlat_only)


ems.condition <- emmeans(mdl1, c('subfield'))
pairs(ems.condition)

