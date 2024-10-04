######################################################
# Handling the big ciliary beat frequency datasets.  #
# These need a beefy computer.                       #
######################################################

### libraries
library(ggplot2)
library(dplyr)
library(tidyverse)
library(data.table)
library(effsize)

###########################
## Making the dataset... ##
###########################
# load in data, removing extra junk (columns etc)
continual2 <- fread("continual2_concatenated.csv")
longcilia2 <- fread("longcilia2_concatenated.csv")
multitest <- fread("multitest_concatenated.csv")
cilia <- fread("cilia_240524.csv")
cilia <- cilia %>%
  select(-V1)
longcilias <- fread("longcilias_240527.csv")
longcilias <- longcilias %>%
  select(-V1)
ruxcilia <- fread("ruxcilia_240524.csv")
ruxcilia <- ruxcilia %>%
  select(-V1)
continual2 <- subset(continual2,  culturepoint < 9 & culturepoint != 3) #omit #3, which had a yeast infection, and pts 9-12 which are not differentiated & venus n construct virus.
good <- rbind(longcilia2, multitest, longcilias, cilia, ruxcilia, continual2)

rm(longcilia2)
rm(multitest)
rm(cilia)
rm(longcilias)
rm(ruxcilia)
rm(continual2)
#write.csv(good, file = "allovem_novenus.csv")
continual2 <- subset(good,
                     experiment == 'continual2')
good <- subset(good, 
               (condition == 'mock' | condition == 'infected') # removing the ruxolitinib treated cultures
               & (experiment != 'continual2') # this got rinsed all the time and i think that affected things.
               & !(experiment == 'ruxcilia' & hpi == -2)) # this is only mock and kind of useless.
  
good <- good %>%
  mutate(hpi_category = case_when(
    hpi < 24 ~ 'early',
    TRUE ~ as.character(hpi)  # Keep other hpi values as they are
  ))

good$hpi_category <- factor(good$hpi_category, levels = c("early", "24", "48", "72", "96", "120", "168"))

# dfs for 1. fovwise visualization and 2. annotating fovs for gfp +/-
thres <- 1100
gfpfov.onlypost <- good %>%
  group_by(video, hpi_category, vidpoint, donor, condition, experiment) %>%
  mutate(gfppos = ifelse(gfp > thres, 'yes', 'no')) %>%
  group_by(video, hpi_category, vidpoint, donor, condition, experiment, gfppos) %>%
  summarise(total = n(), powerful = sum(power > 2e-4)/n(), plaus = sum(domfreq <15)/n(), medfreq = median(domfreq)) %>%
  filter(all(c('yes', 'no') %in% gfppos))
gfpfov.onlypost <- gfpfov.onlypost %>% mutate(status =
                                                case_when(gfppos == 'yes' ~ 'infected', 
                                                          gfppos == 'no' ~ "near")
)

gfpfov.onlynegt <- good %>%
  group_by(video, hpi_category, vidpoint, donor, condition, experiment) %>%
  mutate(gfppos = ifelse(gfp > thres, 'yes', 'no')) %>%
  group_by(video, hpi_category, vidpoint, donor, condition, experiment, gfppos) %>%
  summarise(total = n(), powerful = sum(power > 2e-4)/n(), plaus = sum(domfreq <15)/n(), medfreq = median(domfreq)) %>%
  filter(all(!c('yes') %in% gfppos))
gfpfov.onlynegt$status <- 'distant'

gfpfov.all <- rbind(gfpfov.onlypost, gfpfov.onlynegt)
gfpfov.4cat <- gfpfov.all %>% mutate(status = case_when(condition == 'mock' ~ 'mock',
                                                        condition == 'infected' ~ status))
clean <- subset(gfpfov.4cat, total > 3)
clean$jointstatus <- paste(clean$hpi_category, clean$status) 
level_order <- c('infected', 'near', 'distant', 'mock') 


# labelling pixels by status.
good <- good %>%
  mutate(
    status = case_when(
      condition == 'mock' ~ 'mock',
      video %in% gfpfov.onlynegt$video ~ 'distant',
      video %in% gfpfov.onlypost$video & gfp > thres ~ 'infected',
      video %in% gfpfov.onlypost$video & gfp <= thres ~ 'near',
      #TRUE ~ NA_character_  # To handle cases when no condition is met
      TRUE ~ 'infected'
    )
  )
tf <- good %>% filter(is.na(status))



###################################
## 5d. bulk violin plot of cbfs. ##
###################################
ggplot(good, aes(y=domfreq, x = hpi_category, fill = condition)) +
  labs(x = 'Hours post infection', y ="Dominant pixel frequency (Hz)", color = "Condition")+
  #ylim(0, 20) +
  scale_fill_discrete(labels = c("Infected", "Mock"), name = NULL)+
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), size = 0.2) +
  theme_bw() +
  #facet_wrap(~experiment) +
  # stat_summary(fun.y = mean,
  #              fun.ymin = function(x) mean(x) - sd(x), 
  #              fun.ymax = function(x) mean(x) + sd(x), 
  #              size = 0.1,
  #              linewidth = 0.1,
  #              geom = "pointrange", position = position_dodge(width = 0.9), color = 'black') +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=6),
        plot.margin = margin(1, 0, 0, 0),
        legend.key.size = unit(4, 'point'),
        legend.box.margin = margin(0, 0, 0, -10),
        legend.text = element_text(size=6), legend.title = element_text(size=6))
ggsave("fig5_big_domfreq.pdf", width = 6, height = 4, units = "in")

g1 <- log10(subset(good, condition == 'mock' & hpi_category == '48')$power)
g2 <- log10(subset(good, condition == 'infected' & hpi_category == '48')$power)
cohen.d(g1, g2)
# LOG10POWER 168 hpi mock vs infected, d estimate: 0.5075015 (medium) 0.506981 0.508022 
# LOG10POWER 120 hpi mock vs infected, d estimate: 0.6816253  (medium) 0.6808114 0.6824392 
# LOG10POWER 96 hpi mock vs infected, d estimate: 0.1877761 (negligible) 0.1872624 0.1882898 
# LOG10POWER 72 hpi mock vs infected, d estimate: -0.07107229 (negligible) -0.07148849 -0.07065609 
# LOG10POWER 48 hpi mock vs infected, d estimate: -0.1467173 (negligible) -0.1472622 -0.1461724 
# LOG10POWER 24 hpi mock vs infected, d estimate: 0.1222385 (negligible) 0.1217517 0.1227253 
# LOG10POWER early mock vs infected, d estimate: -0.001164318 (negligible) -0.0017345500 -0.0005940866 

# DOMFREQ 96 hpi mock vs infected, 
calc_cohens_d <- function(df, variable, group) {
  df <- df %>% filter(!is.na(!!sym(variable)))
  infected <- df %>% filter(!!sym(group) == 'infected') %>% pull(!!sym(variable))
  mock <- df %>% filter(!!sym(group) == 'mock') %>% pull(!!sym(variable))
  d <- cohen.d(infected, mock)
  return(d$estimate)
}
  
bulk_ds <- good %>%
  group_by(hpi_category) %>%
  summarize(cohens_d = calc_cohens_d(cur_data(), 'domfreq', 'condition'))



########################################################################
## calculating cohen's d for the bulk mock vs. infected comparisons   ##
########################################################################
calc_cohens_d <- function(df, variable, group) {
  df <- df %>% filter(!is.na(!!sym(variable)))
  infected <- df %>% filter(!!sym(group) == 'infected') %>% pull(!!sym(variable))
  mock <- df %>% filter(!!sym(group) == 'mock') %>% pull(!!sym(variable))
  
  if(length(infected) == 0 | length(mock) == 0) {
    return(list(
      cohen_d = NA,
      sd = NA,
      ci_lower = NA,
      ci_upper = NA,
      magnitude = NA
    ))
  }
  
  d <- cohen.d(infected, mock)
  
  # Determine the magnitude of Cohen's d
  magnitude <- case_when(
    abs(d$estimate) < 0.2 ~ "small",
    abs(d$estimate) < 0.5 ~ "medium",
    abs(d$estimate) < 0.8 ~ "large",
    TRUE ~ "very large"
  )
  
  return(list(
    cohen_d = d$estimate,
    sd = d$sd,
    ci_lower = d$conf.int[1],
    ci_upper = d$conf.int[2],
    magnitude = magnitude
  ))
}
good <- good %>%
  mutate(log10_power = log10(power))
# Apply the function to each group and combine results
log10power_bulk_ds <- good %>%
  group_by(hpi_category) %>%
  group_modify(~ {
    res <- calc_cohens_d(.x, 'log10_power', 'condition')
    tibble(cohen_d = res$cohen_d, sd = res$sd, ci_lower = res$ci_lower, ci_upper = res$ci_upper, magnitude = res$magnitude)
  })
write.csv(log10power_bulk_ds, 'log10power_bulk_cohends.csv')


domfreq_bulk_ds <- good %>%
  group_by(hpi_category) %>%
  group_modify(~ {
    res <- calc_cohens_d(.x, 'domfreq', 'condition')
    tibble(cohen_d = res$cohen_d, sd = res$sd, ci_lower = res$ci_lower, ci_upper = res$ci_upper, magnitude = res$magnitude)
  })
write.csv(domfreq_bulk_ds, 'domfreq_bulk_cohends.csv')
####################################################
## calculating cohen's d for labelled comparisons ##
####################################################
calc_cohens_d <- function(df, variable, group) {
  df <- df %>% filter(!is.na(!!sym(variable)))
  infected <- df %>% filter(!!sym(group) == 'infected') %>% pull(!!sym(variable))
  near <- df %>% filter(!!sym(group) == 'near') %>% pull(!!sym(variable))
  distant <- df %>% filter(!!sym(group) == 'distant') %>% pull(!!sym(variable))
  mock <- df %>% filter(!!sym(group) == 'mock') %>% pull(!!sym(variable))
  
  d1 <- cohen.d(infected, mock)
  d2 <- cohen.d(infected, near)
  d3 <- cohen.d(infected, distant)
  d4 <- cohen.d(distant, mock)
  
  return(list(
    infmock_cohen_d = d1$estimate,
    infmock_sd = d1$sd,
    infmock_ci_lower = d1$conf.int[1],
    infmock_ci_upper = d1$conf.int[2],
    infnear_cohen_d = d2$estimate,
    infnear_sd = d2$sd,
    infnear_ci_lower = d2$conf.int[1],
    infnear_ci_upper = d2$conf.int[2],
    infdist_cohen_d = d3$estimate,
    infdist_sd = d3$sd,
    infdist_ci_lower = d3$conf.int[1],
    infdist_ci_upper = d3$conf.int[2],
    distmock_cohen_d = d4$estimate,
    distmock_sd = d4$sd,
    distmock_ci_lower = d4$conf.int[1],
    distmock_ci_upper = d4$conf.int[2],
  ))
}
good <- good %>%
  mutate(log10_power = log10(power))
# Apply the function to each group and combine results
#looking at change in domfreq among the plausibles...
testt <- subset(good, domfreq < 15)
plausible_domfreq_labelled_ds <- subset(good, domfreq < 15) %>%
  group_by(hpi_category) %>%
  group_modify(~ {
    res <- calc_cohens_d(.x, 'domfreq', 'status')
    tibble(
      infmock_cohen_d = res$infmock_cohen_d,
      infmock_sd = res$infmock_sd,
      infmock_ci_lower = res$infmock_ci_lower,
      infmock_ci_upper = res$infmock_ci_upper,
      infnear_cohen_d = res$infnear_cohen_d,
      infnear_sd = res$infnear_sd,
      infnear_ci_lower = res$infnear_ci_lower,
      infnear_ci_upper = res$infnear_ci_upper,
      infdist_cohen_d = res$infdist_cohen_d,
      infdist_sd = res$infdist_sd,
      infdist_ci_lower = res$infdist_ci_lower,
      infdist_ci_upper = res$infdist_ci_upper,
      distmock_cohen_d = res$distmock_cohen_d,
      distmock_sd = res$distmock_sd,
      distmock_ci_lower = res$distmock_ci_lower,
      distmock_ci_upper = res$distmock_ci_upper
    )
  })
write.csv(log10power_labelled_ds, 'log10power_labelled_cohends.csv')
####
#this one is the one that works ->

# Define the function to calculate Cohen's d
calc_cohens_d <- function(df, variable, group) {
  df <- df %>% filter(!is.na(!!sym(variable)))
  infected <- df %>% filter(!!sym(group) == 'infected') %>% pull(!!sym(variable))
  near <- df %>% filter(!!sym(group) == 'near') %>% pull(!!sym(variable))
  distant <- df %>% filter(!!sym(group) == 'distant') %>% pull(!!sym(variable))
  mock <- df %>% filter(!!sym(group) == 'mock') %>% pull(!!sym(variable))
  
  if(length(infected) == 0 | length(mock) == 0 | length(near) == 0 | length(distant) == 0) {
    return(list(
      infmock_cohen_d = NA,
      infmock_ci_lower = NA,
      infmock_ci_upper = NA,
      infnear_cohen_d = NA,
      infnear_ci_lower = NA,
      infnear_ci_upper = NA,
      infdist_cohen_d = NA,
      infdist_ci_lower = NA,
      infdist_ci_upper = NA,
      distmock_cohen_d = NA,
      distmock_ci_lower = NA,
      distmock_ci_upper = NA
    ))
  }
  
  d1 <- cohen.d(infected, mock)
  d2 <- cohen.d(infected, near)
  d3 <- cohen.d(infected, distant)
  d4 <- cohen.d(distant, mock)
  
  return(list(
    infmock_cohen_d = d1$estimate,
    infmock_ci_lower = d1$conf.int[1],
    infmock_ci_upper = d1$conf.int[2],
    infnear_cohen_d = d2$estimate,
    infnear_ci_lower = d2$conf.int[1],
    infnear_ci_upper = d2$conf.int[2],
    infdist_cohen_d = d3$estimate,
    infdist_ci_lower = d3$conf.int[1],
    infdist_ci_upper = d3$conf.int[2],
    distmock_cohen_d = d4$estimate,
    distmock_ci_lower = d4$conf.int[1],
    distmock_ci_upper = d4$conf.int[2]
  ))
}

good <- good %>%
  mutate(log10_power = log10(power))

# Apply the function to each group and combine results
plausible_domfreq_labelled_ds <- subset(good, domfreq < 15) %>%
  group_by(hpi_category) %>%
  group_modify(~ {
    res <- calc_cohens_d(.x, 'domfreq', 'status')
    tibble(
      infmock_cohen_d = res$infmock_cohen_d,
      infmock_ci_lower = res$infmock_ci_lower,
      infmock_ci_upper = res$infmock_ci_upper,
      infnear_cohen_d = res$infnear_cohen_d,
      infnear_ci_lower = res$infnear_ci_lower,
      infnear_ci_upper = res$infnear_ci_upper,
      infdist_cohen_d = res$infdist_cohen_d,
      infdist_ci_lower = res$infdist_ci_lower,
      infdist_ci_upper = res$infdist_ci_upper,
      distmock_cohen_d = res$distmock_cohen_d,
      distmock_ci_lower = res$distmock_ci_lower,
      distmock_ci_upper = res$distmock_ci_upper
    )
  })

write.csv(domfreq_labelled_ds, 'domfreq_labelled_cohends.csv')
write.csv(plausible_domfreq_labelled_ds, 'plausonly_domfreq_labelled_cohends.csv')


####################################
## 5e. bulk density plot of power.##
####################################
ggplot(good, aes(x=power, color = condition, fill = condition)) +
  labs(x = 'Hours post infection', y ="Power Density (pixel units^2/Hz)", color = NULL)+
  #ylim(0, 20) +
  theme_bw() +
  facet_wrap(~hpi_category) +
  geom_density(alpha = 0.1) +
  scale_fill_discrete(labels = c("Infected", "Mock"), name = NULL)+
  geom_vline(xintercept = 2e-4, linetype = 'dashed') +
  scale_x_log10() +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=6),
        plot.margin = margin(1, 0, 0, 0),
        legend.key.size = unit(4, 'point'),
        legend.box.margin = margin(0, 0, 0, -10),
        legend.text = element_text(size=6), legend.title = element_text(size=6))
ggsave("fig5_big_density.pdf", width = 6, height = 4, units = "in")

ggplot(good, aes(x=domfreq, color = condition, fill = condition)) +
  labs(x = 'Hours post infection', y ="Power Density (pixel units^2/Hz)", color = NULL)+
  #ylim(0, 20) +
  theme_bw() +
  facet_wrap(~hpi_category) +
  geom_density(alpha = 0.1) +
  scale_fill_discrete(labels = c("Infected", "Mock"), name = NULL)+
  geom_vline(xintercept = 2e-4, linetype = 'dashed') +
  #scale_x_log10() +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=6),
        plot.margin = margin(1, 0, 0, 0),
        legend.key.size = unit(4, 'point'),
        legend.box.margin = margin(0, 0, 0, -10),
        legend.text = element_text(size=6), legend.title = element_text(size=6))
ggsave("fig5_big_densitystyle_domfreq.pdf", width = 6, height = 4, units = "in")




###################################################################################
###  plot of power vs. freq in both, with labels for the sides. bulk & labelled ###
###################################################################################
# Define the threshold values
x_threshold <- 2e-4
y_threshold <- 15

# Calculate fractions
good <- good %>%
  mutate(quadrant = case_when(
    power <= x_threshold & domfreq <= y_threshold ~ "weakbeatin",
    power > x_threshold & domfreq <= y_threshold ~ "strongbeatin",
    power <= x_threshold & domfreq > y_threshold ~ "weaknoise",
    power > x_threshold & domfreq > y_threshold ~ "strongnoise"
  ))

quadrant_counts <- good %>%
  group_by(hpi_category, condition, quadrant) %>%
  summarise(count = n()) %>%
  mutate(total = sum(count), fraction = count / total)

# Subset data to get fraction for each facet
fraction_data <- good %>%
  group_by(hpi_category, condition) %>%
  summarise(
    weakbeatin = sum(power <= x_threshold & domfreq <= y_threshold) / n(),
    strongbeatin = sum(power > x_threshold & domfreq <= y_threshold) / n(),
    weaknoise = sum(power <= x_threshold & domfreq > y_threshold) / n(),
    strongnoise = sum(power > x_threshold & domfreq > y_threshold) / n()
  )

# Bulk
set.seed(70)
good.smol <- sample_frac(good, 0.05)
ggplot(good.smol, aes(y = domfreq, x = power)) +
  geom_density_2d() +  
  theme_bw() +
  scale_x_log10() +
  geom_vline(xintercept = x_threshold, linetype = 'dashed', color = 'deeppink2') +
  geom_hline(yintercept = y_threshold, linetype = 'dashed', color = 'deeppink2') +
  facet_grid(hpi_category ~ condition) +
  theme(axis.text = element_text(size = 6), axis.title = element_text(size = 6),
        plot.margin = margin(1, 0, 0, 0),
        legend.key.size = unit(4, 'point'),
        legend.box.margin = margin(0, 0, 0, -10),
        legend.text = element_text(size = 6), legend.title = element_text(size = 6)) +
  geom_text(data = fraction_data, aes(x = 1e-6, y = 0, label = paste0("Q1: ", round(weakbeatin, 2))), size = 3, hjust = 0, vjust = 0) +
  geom_text(data = fraction_data, aes(x = 1e-1, y = 0, label = paste0("Q2: ", round(strongbeatin, 2))), size = 3, hjust = 1, vjust = 0) +
  geom_text(data = fraction_data, aes(x = 1e-6, y = 60, label = paste0("Q3: ", round(weaknoise, 2))), size = 3, hjust = 0, vjust = 1) +
  geom_text(data = fraction_data, aes(x = 1e-1, y = 60, label = paste0("Q4: ", round(strongnoise, 2))), size = 3, hjust = 1, vjust = 1)
ggsave("bulk_powervsfreq.pdf", width = 4, height = 9, units = "in")


# Labeled
quadrant_counts <- good %>%
  filter(hpi > 10) %>%
  group_by(hpi_category, status, quadrant) %>%
  summarise(count = n()) %>%
  mutate(total = sum(count), fraction = count / total)

# Subset data to get fraction for each facet
fraction_data <- good %>%
  filter(hpi > 10) %>%
  group_by(hpi_category, status) %>%
  summarise(
    weakbeatin = sum(power <= x_threshold & domfreq <= y_threshold) / n(),
    strongbeatin = sum(power > x_threshold & domfreq <= y_threshold) / n(),
    weaknoise = sum(power <= x_threshold & domfreq > y_threshold) / n(),
    strongnoise = sum(power > x_threshold & domfreq > y_threshold) / n()
  )
# Define the order of the status levels
level_order <- c('infected', 'near', 'distant', 'mock')

# Convert the status variable to a factor with the specified order
good.smol <- good.smol %>%
  mutate(status = factor(status, levels = level_order))

ggplot(subset(good.smol, hpi > 10), aes(y = domfreq, x = power)) +
  geom_density_2d() +  
  theme_bw() +
  scale_x_log10() +
  geom_vline(xintercept = x_threshold, linetype = 'dashed', color = 'deeppink2') +
  geom_hline(yintercept = y_threshold, linetype = 'dashed', color = 'deeppink2') +
  facet_grid(hpi_category ~ status) +
  theme(axis.text = element_text(size = 6), axis.title = element_text(size = 6),
        plot.margin = margin(1, 0, 0, 0),
        legend.key.size = unit(4, 'point'),
        legend.box.margin = margin(0, 0, 0, -10),
        legend.text = element_text(size = 6), legend.title = element_text(size = 6)) +
  geom_text(data = fraction_data, aes(x = 1e-6, y = 0, label = paste0("Q1: ", round(weakbeatin, 2))), size = 3, hjust = 0, vjust = 0) +
  geom_text(data = fraction_data, aes(x = 1e-1, y = 0, label = paste0("Q2: ", round(strongbeatin, 2))), size = 3, hjust = 1, vjust = 0) +
  geom_text(data = fraction_data, aes(x = 1e-6, y = 60, label = paste0("Q3: ", round(weaknoise, 2))), size = 3, hjust = 0, vjust = 0) +
  geom_text(data = fraction_data, aes(x = 1e-1, y = 60, label = paste0("Q4: ", round(strongnoise, 2))), size = 3, hjust = 1, vjust = 0)
ggsave("labelled_powervsfreq.pdf", width = 6, height = 9, units = "in")

#############################################
## 5f. violin plot of labelled pixel cbfs. ##
#############################################
level_order <- c('infected', 'near', 'distant', 'mock') 
ggplot(subset(good, hpi >10), aes(y=domfreq, x = factor(hpi_category), fill = factor(status, levels = level_order))) +
  labs(x = 'Hours post infection', y ="Dominant pixel frequency (Hz)", color = NULL)+
  #ylim(0, 20) +
  theme_bw() +
  scale_fill_manual(values=c('chartreuse3', 'darkorchid4', 'darkorchid2', 'deeppink1'), labels = c("GFP+", "Near GFP-", "Far GFP-", "Mock"), name = NULL)+
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), size = 0.2) +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=6),
        plot.margin = margin(1, 0, 0, 0),
        legend.key.size = unit(4, 'point'),
        legend.box.margin = margin(0, 0, 0, -10),
        legend.text = element_text(size=6), legend.title = element_text(size=6))
ggsave("fig5_big_labelled_domfreq.pdf", width = 6, height = 4, units = "in")

################################################
## 5g. density blot of labelled pixel powers. ##
################################################
ggplot(subset(good, hpi >10), aes(x=power, color = factor(status), fill = factor(status))) +
  labs(x = 'Power Density (pixel units^2/Hz)', y ="Density", color = NULL)+
  theme_bw() +
  facet_wrap(~hpi_category) +
  scale_fill_manual(values=c('chartreuse3', 'darkorchid4', 'darkorchid2', 'deeppink1'), labels = c("GFP+", "Near GFP-", "Far GFP-", "Mock"), name = NULL)+
  scale_color_manual(values=c('chartreuse3', 'darkorchid4', 'darkorchid2', 'deeppink1'), labels = c("GFP+", "Near GFP-", "Far GFP-", "Mock"), name = NULL)+
  geom_density(alpha = 0.1) +
  geom_vline(xintercept = 2e-4, linetype = 'dashed') +
  scale_x_log10() +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=6),
        plot.margin = margin(1, 0, 0, 0),
        legend.key.size = unit(4, 'point'),
        legend.box.margin = margin(0, 0, 0, -10),
        legend.text = element_text(size=6), legend.title = element_text(size=6))
ggsave("fig5_big_labelled_power.pdf", width = 6, height = 4, units = "in")


######################
## MAKIN BEE SWARMS ##
######################
fxns_labelled_fov <- good %>%
  group_by(hpi_category, status, experiment, vidpoint) %>%
  summarize(total = n(),
            donor = first(donor),
            hpi = first(hpi),
            powerful = sum(power > 2e-4)/n(), 
            plaus = sum(domfreq < 15)/n(),
            ultraplaus = sum(domfreq < 15 & domfreq > 2)/n(),
            low = sum(domfreq < 2)/n(),
            medpower = median(power),
            medlogpower = median(log10(power)),
            medfreq = median(domfreq),
            medofplaus = median(ifelse(domfreq < 15, domfreq, NA), na.rm = TRUE),
            medofultraplaus = median(ifelse(domfreq > 2 & domfreq < 15, domfreq, NA), na.rm = TRUE),
            powerbeatin = sum(domfreq < 15 & power > 2e-4)/n())
write.csv(fxns_labelled_fov, file ='big_fxns_labelled_byvidpoint.csv')

fxns_labelled_culturepoint <- good %>%
  group_by(hpi_category, status, experiment, culturepoint) %>%
  summarize(total = n(),
            donor = first(donor),
            hpi = first(hpi),
            powerful = sum(power > 2e-4)/n(), 
            plaus = sum(domfreq < 15)/n(),
            ultraplaus = sum(domfreq < 15 & domfreq > 2)/n(),
            low = sum(domfreq < 2)/n(),
            medpower = median(power),
            medlogpower = median(log10(power)),
            medfreq = median(domfreq),
            medofplaus = median(ifelse(domfreq < 15, domfreq, NA), na.rm = TRUE),
            medofultraplaus = median(ifelse(domfreq > 2 & domfreq < 15, domfreq, NA), na.rm = TRUE),
            powerbeatin = sum(domfreq < 15 & power > 2e-4)/n())
write.csv(fxns_labelled_culturepoint, file ='big_fxns_labelled_byculturepoint.csv')

fxns_labelled_bystatus <- good %>%
  group_by(hpi_category, status) %>%
  summarize(total = n(),
            powerful = sum(power > 2e-4)/n(), 
            plaus = sum(domfreq < 15)/n(),
            ultraplaus = sum(domfreq < 15 & domfreq > 2)/n(),
            low = sum(domfreq < 2)/n(),
            medpower = median(power),
            medlogpower = median(log10(power)),
            medfreq = median(domfreq),
            medofplaus = median(ifelse(domfreq < 15, domfreq, NA), na.rm = TRUE),
            medofultraplaus = median(ifelse(domfreq > 2 & domfreq < 15, domfreq, NA), na.rm = TRUE),
            powerbeatin = sum(domfreq < 15 & power > 2e-4)/n())
write.csv(fxns_labelled_bystatus, file ='big_fxns_labelled_bystatusonly.csv')


fxns_bulk_bycondition <- good %>%
  group_by(hpi_category, condition) %>%
  summarize(total = n(),
            powerful = sum(power > 2e-4)/n(), 
            plaus = sum(domfreq < 15)/n(),
            ultraplaus = sum(domfreq < 15 & domfreq > 2)/n(),
            low = sum(domfreq < 2)/n(),
            medpower = median(power),
            medlogpower = median(log10(power)),
            medfreq = median(domfreq),
            medofplaus = median(ifelse(domfreq < 15, domfreq, NA), na.rm = TRUE),
            medofultraplaus = median(ifelse(domfreq > 2 & domfreq < 15, domfreq, NA), na.rm = TRUE),
            powerbeatin = sum(domfreq < 15 & power > 2e-4)/n())
write.csv(fxns_bulk_bycondition, file ='big_fxns_bulk_byconditiononly.csv')

fxns_bulk_fov <- good %>%
  group_by(hpi_category, condition, experiment, vidpoint) %>%
  summarize(total = n(),
            donor = first(donor),
            hpi = first(hpi),
            powerful = sum(power > 2e-4)/n(), 
            plaus = sum(domfreq < 15)/n(),
            ultraplaus = sum(domfreq < 15 & domfreq > 2)/n(),
            low = sum(domfreq < 2)/n(),
            medpower = median(power),
            medlogpower = median(log10(power)),
            medfreq = median(domfreq),
            medofplaus = median(ifelse(domfreq < 15, domfreq, NA), na.rm = TRUE),
            medofultraplaus = median(ifelse(domfreq > 2 & domfreq < 15, domfreq, NA), na.rm = TRUE),
            powerbeatin = sum(domfreq < 15 & power > 2e-4)/n())
write.csv(fxns_bulk_fov, file ='big_fxns_bulk_byvidpoint.csv')

fxns_bulk_culturepoint <- good %>%
  group_by(hpi_category, condition, experiment, culturepoint) %>%
  summarize(total = n(),
            donor = first(donor),
            hpi = first(hpi),
            powerful = sum(power > 2e-4)/n(), 
            plaus = sum(domfreq < 15)/n(),
            ultraplaus = sum(domfreq < 15 & domfreq > 2)/n(),
            low = sum(domfreq < 2)/n(),
            medpower = median(power),
            medlogpower = median(log10(power)),
            medfreq = median(domfreq),
            medofplaus = median(ifelse(domfreq < 15, domfreq, NA), na.rm = TRUE),
            medofultraplaus = median(ifelse(domfreq > 2 & domfreq < 15, domfreq, NA), na.rm = TRUE),
            powerbeatin = sum(domfreq < 15 & power > 2e-4)/n())
write.csv(fxns_bulk_culturepoint, file ='big_fxns_bulk_byculturepoint.csv')


fxns_labelled_fov
fxns_labelled_culturepoint
fxns_labelled_bystatus

fxns_bulk_bycondition
fxns_bulk_fov
fxns_bulk_culturepoint

ggplot(fxns_bulk_fov, aes(x = hpi_category, 
                                          color = factor(condition),
                                          #color = factor(condition, levels = level_order), 
                                          y = ultraplaus)) +
  #facet_wrap(~experiment) +
  geom_jitter(alpha = 0.1)+
  geom_jitter(data = fxns_bulk_culturepoint, aes(x = hpi_category, y = ultraplaus, shape = as.factor(donor)))+
  labs(y = 'Fraction Beating Surface Area', x = 'Hours Post Infection') +
  stat_summary(fun = mean,
               fun.min = function(x) mean(x) - sd(x), 
               fun.max = function(x) mean(x) + sd(x),
               size = 1, 
               alpha = 0.8,
               linewidth = 0.5,
               geom = "pointrange") +
  theme_bw() +
  geom_hline(yintercept = 0.2166, linetype = 'dashed') 
  #scale_color_manual(values=c('chartreuse3', 'darkorchid4', 'darkorchid2', 'deeppink1'), labels = c("GFP+", "Near GFP-", "Far GFP-", "Mock"), name = NULL)
g1 <- subset(fxns_bulk_fov, condition == 'infected' & hpi == '72')$ultraplaus
g2 <- subset(fxns_bulk_fov, condition == 'mock' & hpi == '72')$ultraplaus
t.test(g1, g2)

ggplot(fxns_bulk_fov, aes(x = hpi_category, 
                          color = factor(condition),
                          #color = factor(condition, levels = level_order), 
                          y = medofultraplaus)) +
  #facet_wrap(~experiment) +
  geom_jitter(alpha = 0.1)+
  geom_violin()+
  geom_jitter(data = fxns_bulk_culturepoint, aes(x = hpi_category, y = medofultraplaus, shape = as.factor(donor)))+
  labs(y = 'Beat Frequency of Plausibly Motile Culture Area', x = 'Hours Post Infection') +
  stat_summary(fun = mean,
               fun.min = function(x) mean(x) - sd(x), 
               fun.max = function(x) mean(x) + sd(x),
               size = 1, 
               alpha = 0.8,
               linewidth = 0.5,
               geom = "pointrange") +
  theme_bw()


#scale_color_manual(values=c('chartreuse3', 'darkorchid4', 'darkorchid2', 'deeppink1'), labels = c("GFP+", "Near GFP-", "Far GFP-", "Mock"), name = NULL)
g1 <- subset(fxns_bulk_fov, condition == 'infected' & hpi == '72')$ultraplaus
g2 <- subset(fxns_bulk_fov, condition == 'mock' & hpi == '72')$ultraplaus
t.test(g1, g2)


ggplot(fxns_labelled_fov, aes(x = hpi_category, 
                          #color = factor(status),
                          color = factor(status, levels = level_order), 
                          y = ultraplaus)) +
  #facet_wrap(~experiment) +
  geom_jitter(alpha = 0.1)+
  geom_jitter(data = fxns_labelled_culturepoint, aes(x = hpi_category, y = ultraplaus, shape = as.factor(donor)))+
  labs(y = 'Fraction Beating Surface Area', x = 'Hours Post Infection') +
  stat_summary(fun = mean,
               fun.min = function(x) mean(x) - sd(x), 
               fun.max = function(x) mean(x) + sd(x),
               size = 1, 
               alpha = 0.8,
               linewidth = 0.5,
               geom = "pointrange") +
  theme_bw() +
  geom_hline(yintercept = 0.2166, linetype = 'dashed') +
  scale_color_manual(values=c('chartreuse3', 'darkorchid4', 'darkorchid2', 'deeppink1'), labels = c("GFP+", "Near GFP-", "Far GFP-", "Mock"), name = NULL)
g1 <- subset(fxns_bulk_fov, condition == 'infected' & hpi == '72')$ultraplaus
g2 <- subset(fxns_bulk_fov, condition == 'mock' & hpi == '72')$ultraplaus
t.test(g1, g2)


ggplot(fxns_labelled_fov, aes(x = hpi_category, 
                              #color = factor(status),
                              color = factor(status, levels = level_order), 
                              y = medpower)) +
  #facet_wrap(~experiment) +
  geom_jitter(alpha = 0.1)+
  scale_y_log10() +
  geom_jitter(data = fxns_labelled_culturepoint, aes(x = hpi_category, y = medpower, shape = as.factor(donor)))+
  labs(y = 'Fraction Beating Surface Area', x = 'Hours Post Infection') +
  stat_summary(fun = mean,
               fun.min = function(x) mean(x) - sd(x), 
               fun.max = function(x) mean(x) + sd(x),
               size = 1, 
               alpha = 0.8,
               linewidth = 0.5,
               geom = "pointrange") +
  theme_bw() +
  geom_hline(yintercept = 2e-4, linetype = 'dashed') +
  scale_color_manual(values=c('chartreuse3', 'darkorchid4', 'darkorchid2', 'deeppink1'), labels = c("GFP+", "Near GFP-", "Far GFP-", "Mock"), name = NULL)

g1 <- subset(fxns_labelled_culturepoint, status == 'mock' & hpi == '120')$medofultraplaus
g2 <- subset(fxns_labelled_culturepoint, status == 'distant' & hpi == '120')$medofultraplaus
t.test(g1, g2)

ggplot(subset(fxns, hpi_category != 'early'), aes(x = hpi, color = factor(status, levels = level_order), y = ultraplaus)) +
  #facet_wrap(~experiment) +
  geom_jitter(alpha = 0.1, position = position_dodge2(width = 20))+
  stat_summary(fun = mean,
               fun.min = function(x) mean(x) - sd(x), 
               fun.max = function(x) mean(x) + sd(x),
               size = 1, 
               alpha = 0.8,
               linewidth = 0.5,
               position = position_dodge2(width = 20),
               geom = "pointrange") +
  geom_hline(yintercept = 0.15) +
  scale_color_manual(values=c('chartreuse3', 'darkorchid4', 'darkorchid2', 'deeppink1'), labels = c("GFP+", "Near GFP-", "Far GFP-", "Mock"), name = NULL)

ggplot(subset(fxns, hpi_category != 'early'), aes(x = hpi, color = factor(status, levels = level_order), y = medofplaus)) +
  #facet_wrap(~experiment) +
  geom_jitter(alpha = 0.1, position = position_dodge2(width = 20))+
  stat_summary(fun = mean,
               fun.min = function(x) mean(x) - sd(x), 
               fun.max = function(x) mean(x) + sd(x),
               size = 1, 
               alpha = 0.8,
               linewidth = 0.5,
               position = position_dodge2(width = 20),
               geom = "pointrange") +
  geom_hline(yintercept = 0.15) +
  scale_color_manual(values=c('chartreuse3', 'darkorchid4', 'darkorchid2', 'deeppink1'), labels = c("GFP+", "Near GFP-", "Far GFP-", "Mock"), name = NULL)

ggplot(subset(fxns, hpi_category != 'early'), aes(x = hpi, color = factor(status, levels = level_order), y = medofultraplaus)) +
  #facet_wrap(~experiment) +
  geom_jitter(alpha = 0.1, position = position_dodge2(width = 20))+
  stat_summary(fun = mean,
               fun.min = function(x) mean(x) - sd(x), 
               fun.max = function(x) mean(x) + sd(x),
               size = 1, 
               alpha = 0.8,
               linewidth = 0.5,
               position = position_dodge2(width = 20),
               geom = "pointrange") +
  geom_hline(yintercept = 0.15) +
  scale_color_manual(values=c('chartreuse3', 'darkorchid4', 'darkorchid2', 'deeppink1'), labels = c("GFP+", "Near GFP-", "Far GFP-", "Mock"), name = NULL)
