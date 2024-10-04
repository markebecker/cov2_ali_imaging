# Changes in ciliary motion during infection (Figs. 7 & 8; Supp. Fig 11)
# Figs. 12 & 13 and the more data-intense portions of figure 8 were plotted in '5_cilia_in_infxn_largedataset.R'

# libraries
library(readr)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(effsize)
library(survival)
library(survminer)
library(data.table)

########################################
## Mucus disc spinning cessation time ##
########################################

exptlog_iii <- read_csv("data/exptlog_iii.csv")
spots <- read_csv("data/240702_goodspots.csv")

# take the spots over time csv and extract specific values of interest
daily.spots.summary <- spots %>%
  group_by(csv) %>%
  summarise(good.peak.spots = max(norfreq), 
            good.peak.spots.time = time[which.max(norfreq)], 
            onedayspots = norfreq[time == 23],
            twodayspots = norfreq[time == 47],
            threedayspots = norfreq[time == 71],
            fourdayspots = norfreq[time == 95],
            fivedayspots = norfreq[time == 119],
            aucspots = sum(norfreq)) 
# cross back with the info i want
daily.spots.summary.j <- inner_join(daily.spots.summary, exptlog_iii)

relevant <- subset(daily.spots.summary.j, 
                   (condition == 'vanilla' | condition == 'static'| condition == 'shaken')
                   & (virus != 'venus')
                   & (usable == 1)
                   & (mucus.stop.frame.end == 'yes' | mucus.stop.frame.end == 'no'))
relevant <- relevant %>%
  mutate(mucus.stop.frame.end = as.character(mucus.stop.frame.end))

relevant <- relevant %>%
  mutate(mucus.stop.frame.end = case_when(
    mucus.stop.frame.end == 'no' ~ 1,
    mucus.stop.frame.end == 'yes' ~ 0,
    TRUE ~ as.numeric(mucus.stop.frame.end)  # Convert other numeric values to numeric type
  ))


surv_obj <- Surv(time=2*(relevant$mucus.stop.frame -1)+1, event = relevant$mucus.stop.frame.end)
surv_fit <- survfit(surv_obj ~ virus, data = relevant)

surv_diff <- survdiff(surv_obj ~ virus, data = relevant)
print(surv_diff)
color_palette <- c("mock" = "deeppink2", "infected" = "cyan4")

relevant$mucus.stop.time <- 2*(relevant$mucus.stop.frame -1)+1
f5a.small <- relevant[, c("csv", "donor", "expt.name", "sample.no", "mucus.stop.frame", 'mucus.stop.time', 'virus')]
write.csv(f5a.small,"sourcedata/fig5a.csv", row.names = FALSE)



p <- ggsurvplot(surv_fit, data = relevant, pval = TRUE, 
                conf.int = TRUE,
                risk.table = FALSE,
                palette = c("cyan4", "deeppink2"),
                conf.int.style = 'ribbon',
                ggtheme = theme_bw(),
                pval.coord = c(0.15, 0.4),
                legend = "bottom",
                pval.size = 4,
                #pval.method = TRUE,
                axes.offset = TRUE,
                xlab = "Hours post infection", ylab = "Probability of Mucus Still Moving",
                legend.title = "Condition",
                legend.labs = c("Infected", "Mock"))
p <- ggpar(p, font.x = 14, font.y = 14, font.tickslab = 12, font.legend = 12)
p$plot <- p$plot + theme(legend.key.size = unit(1, "lines"), legend.position = c(0.2, 0.15))
print(p)
ggsave("plots/fig5a.pdf", plot = p$plot, width = 4, height = 4, units = "in")


cox_model <- coxph(surv_obj ~ virus + donor + good.peak.spots, data = relevant)
cox_model <- coxph(surv_obj ~ virus, data = relevant)

summary(cox_model)

###############################
## 5b, Cell survival duration ##
###############################
hunting <- read.csv("data/hunting_ii.csv")
hunting <- hunting %>%
  mutate(time.disappears = as.numeric(time.disappears))
hunting <- hunting %>%
  mutate(cat = case_when(
    when.identified == 'first' ~ 'first',
    when.identified != 'first' ~ 'late',
    TRUE ~ NA_character_
  ))

hunting.filtered <- hunting %>% filter(cat %in% c('first', 'late'))
hunting.filtered$eqtreat <- (hunting.filtered$time.disappears - hunting.filtered$time.appears)

write.csv(hunting.filtered,"sourcedata/fig5c.csv", row.names = FALSE)

# Create the survival object
surv_obj <- Surv(time = hunting.filtered$eqtreat, event = rep(1, nrow(hunting.filtered)))

# Fit the survival model
surv_fit <- survfit(surv_obj ~ 1, data = hunting.filtered)

# Plot the survival curve
p <- ggsurvplot(surv_fit, data = hunting.filtered, 
                conf.int = TRUE,
                risk.table = FALSE,
                conf.int.style = 'ribbon',
                ggtheme = theme_bw(),
                legend = "none",
                #xlim = c(0, 90),
                axes.offset = TRUE,
                xlab = "GFP+ Lifetime (Hours)", ylab = "Survival Probability")
p <- ggpar(p, font.x = 14, font.y = 14, font.tickslab = 12)
print(p)
ggsave("plots/fig5c.pdf", plot = p$plot, width = 3.2, height = 2.4, units = "in")


summary(hunting$fov)
#######################################
# Supplement w/ survival time by strain ##
#######################################
## supplement for orf7a+ survival time


# libraries
library(readr)
library(ggplot2)
library(dplyr)
library(survival)
library(survminer)

###############################
## 5b, Cell survival duration ##
###############################
hunting <- read.csv("data/hunting_ii.csv")
vo <- read.csv("data/venus_omicron_survival.csv")
hunting <- hunting %>%
  mutate(time.disappears = as.numeric(time.disappears))
hunting <- hunting %>%
  mutate(experiment = "d7a-eGFP")

vo <- vo %>%
  filter(experiment != 'omicron')
combined_data <- bind_rows(hunting, vo)

write.csv(combined_data,"sourcedata/s11.csv", row.names = FALSE)

# Create the survival object
surv_obj <- Surv(time = combined_data$gfp.survival.time, event = rep(1, nrow(combined_data)))

# Fit the survival model
surv_fit <- survfit(surv_obj ~ experiment, data = combined_data)
p <- ggsurvplot(surv_fit, data = combined_data, 
                pval = TRUE,
                conf.int = TRUE,
                risk.table = FALSE,
                xlim = c(0, 85),
                palette = c("cyan4", "coral2"),
                conf.int.style = 'ribbon',
                ggtheme = theme_bw(),
                legend = "right",
                xlab = "GFP+ Lifetime (hours)", ylab = "Survival Probability",
                pval.coord = c(0.15, 0.2),
                legend.labs = c("WA-1/d7a-eGFP", "WA-1/Venus-N"))
p <- ggpar(p, font.x = 14, font.y = 14, font.tickslab = 12, font.legend = 12)
p$plot <- p$plot + theme(legend.key.size = unit(1, "lines"))
print(p)
ggsave("plots/fig5_orf7a_supp_nomicron.pdf", plot = p$plot, width = 6, height = 4, units = "in")


### other potential supplement
# Create the survival object
surv_obj <- Surv(time = hunting.filtered$eqtreat, event = rep(1, nrow(hunting.filtered)))

# Fit the survival model
surv_fit <- survfit(surv_obj ~ cat, data = hunting.filtered)
p <- ggsurvplot(surv_fit, data = hunting.filtered, 
                pval = TRUE,
                conf.int = TRUE,
                risk.table = FALSE,
                xlim = c(0, 85),
                conf.int.style = 'ribbon',
                ggtheme = theme_bw(),
                legend = "bottom",
                xlab = "GFP+ Lifetime (hours)", ylab = "Survival Probability",
                legend.title = "Identified:",
                pval.coord = c(0.15, 0.4),
                legend.labs = c("16-24 HPI", "48+ HPI"))
p <- ggpar(p, font.x = 14, font.y = 14, font.tickslab = 12, font.legend = 12)
p$plot <- p$plot + theme(legend.key.size = unit(1, "lines"), legend.position = c(0.2, 0.15))
print(p)
ggsave("C:/Users/Mark/Desktop/figure_images/5/revised/fig5C_SUPPLEMENT.pdf", plot = p$plot, width = 4, height = 4, units = "in")

# About how long does it take to first see gfp in the first round of infected cells?
mk <- subset(hunting.filtered, cat == 'first')
mean(mk$frame.appears)
median(mk$frame.appears)


###############################
## 5d, Cell beating duration ##
###############################
hunting.long <- hunting %>%
  pivot_longer(cols = starts_with("median.beat."),
               names_to = "timepoint",
               names_prefix = "median.beat.",
               values_to = "beat_frequency")

# Convert the timepoint to a numeric variable
hunting.long$timepoint <- as.numeric(hunting.long$timepoint)

# adjust timepoint times to reflect the times that actually occurred in different experiments.
hunting.long <- hunting.long %>%
  mutate(timepoint.real = case_when(
    experiment == "ruxcilia" & timepoint == 24 ~ 25.1,
    experiment == "ruxcilia" & timepoint == 48 ~ 51.34,
    experiment == "ruxcilia" & timepoint == 72 ~ 73.8,
    experiment == "longcilia2" & timepoint == 24 ~ 29.4,
    experiment == "longcilia2" & timepoint == 48 ~ 54.1,
    experiment == "longcilia2" & timepoint == 72 ~ 77.4,
    experiment == "longcilia2" & timepoint == 96 ~ 101.8,
    experiment == "multitest" & timepoint == 24 ~ 26.8,
    experiment == "multitest" & timepoint == 72 ~ 75.5,
    TRUE ~ NA_real_  # Default case if no condition is met
  ))
hunting.long <- hunting.long %>%
  mutate(beat_frequency = as.numeric(beat_frequency),
         befr = case_when(
           beat_frequency > 15 ~ 0,
           TRUE ~ beat_frequency)) %>%
  filter(!is.na(frame.appears))

# Ensure correct column names and data types
hunting.long <- hunting.long %>%
  mutate(timepoint.real = as.numeric(timepoint.real),
         befr = as.numeric(befr))
write.csv(hunting.long,"sourcedata/fig5d.csv", row.names = FALSE)

# Create a dataframe with the last observed 'befr' value for each 'whole.coords'
last_observed <- hunting.long %>%
  group_by(whole.coords) %>%
  filter(timepoint.real == max(timepoint.real, na.rm = TRUE)) %>%
  summarize(time.disappears = first(time.disappears), last_befr = first(befr), .groups = 'drop')


# Create dataframes for the segments
segments <- hunting.long %>%
  group_by(whole.coords) %>%
  filter(timepoint.real == max(timepoint.real, na.rm = TRUE)) %>%
  mutate(time.disappears = first(time.disappears),
         last_befr = first(befr)) %>%
  select(whole.coords, timepoint.real, befr, time.disappears, last_befr)

segments_from_appear <- hunting.long %>%
  group_by(whole.coords) %>%
  filter(timepoint.real == min(timepoint.real, na.rm = TRUE)) %>%
  mutate(time.appears = first(time.appears),
         timepoint.real = first(timepoint.real),
         first_befr = first(befr)) %>%
  select(whole.coords, timepoint.real, time.appears, first_befr)

# Plot the data
ggplot(subset(hunting.long), aes(x = timepoint.real, y = as.numeric(befr), group = whole.coords)) +
  geom_line(aes(color = factor(whole.coords))) +
  geom_point(aes(color = factor(whole.coords)),  size = 0.5) +
  #geom_point(data = last_observed, alpha = 0.4, aes(x = time.disappears, y = last_befr, color = factor(whole.coords))) +
  geom_segment(data = segments, alpha = 0.4, aes(x = timepoint.real, y = as.numeric(befr), xend = time.disappears, yend = as.numeric(last_befr), color = factor(whole.coords))) +
  geom_segment(data = segments_from_appear, alpha = 0.4, aes(x = time.appears, y = as.numeric(first_befr), xend = timepoint.real, yend = as.numeric(first_befr), color = factor(whole.coords))) +
  labs(x = 'Hours post infection', y = 'Ciliary Beat Frequency', color = 'Cell ID') +
  #xlim(24, 105) +
  theme_bw() +
  #guides(color = guide_legend(ncol = 1)) +
  #stat_smooth() +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=6),
        plot.margin = margin(1, 2, 0, 0),
        legend.position = 'none')
        legend.key.size = unit(4, 'point'),
        legend.box.margin = margin(0, 0, 0, -10),
        legend.text = element_text(size=7), legend.title = element_text(size=6))
ggsave("plots/fig5d_test.pdf", width = 1.6, height = 1.4, units = "in")

# how long were things observed beating?
first <- subset(hunting, when.identified != 'third' & !is.na(frame.appears))
mean(first$n.plauscbfs.b4death, na.rm=TRUE)

table(first$n.plauscbfs.b4death)

all <- subset(hunting,!is.na(frame.appears))

# focus initiation.
longsurviving <- subset(first, n.plauscbfs.b4death >=3)
table(longsurviving$focus.initiator)
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

##############################################################
## Fig. 5e. Fxn beating surface area over time in aggregate ##
##############################################################
fxns_bulk_fov <- fread('data/big_fxns_bulk_byvidpoint.csv')
fxns_bulk_culturepoint <- fread('data/big_fxns_bulk_byculturepoint.csv')
fxns_bulk_fov$hpi_category <- factor(fxns_bulk_fov$hpi_category, levels = c("early", "24", "48", "72", "96", "120", "168"))
fxns_bulk_culturepoint$hpi_category <- factor(fxns_bulk_culturepoint$hpi_category, levels = c("early", "24", "48", "72", "96", "120", "168"))

fxns_bulk_fov <- subset(fxns_bulk_fov, hpi_category != 120)
fxns_bulk_culturepoint <- subset(fxns_bulk_culturepoint, hpi_category != 120)

ggplot(fxns_bulk_fov, aes(x = hpi_category, color = factor(condition), y = ultraplaus)) +
  geom_jitter(alpha = 0.05, size = 0.4, position = position_jitterdodge(dodge.width = 0.5))+
  geom_point(data = fxns_bulk_culturepoint, size = 0.4, alpha = 0.4, position = position_jitterdodge(dodge.width = 0.5), aes(shape = as.factor(donor)))+
  labs(y = 'Fraction Beating Surface Area', x = 'Hours Post Infection', shape = 'Donor:', color = 'Condition:') +
  stat_summary(fun = mean,
               fun.min = function(x) mean(x) - sd(x), 
               fun.max = function(x) mean(x) + sd(x),
               size = 0.4, 
               alpha = 0.8,
               linewidth = 0.5,
               position = position_dodge(width = 0.5),
               geom = "pointrange") +
  theme_bw() +
  ylim(0, 1.1) +
  scale_color_manual(values=c('cyan4', 'deeppink1'), labels = c("Infected", "Mock"), name = NULL) +
  #geom_hline(yintercept = 0.2166, linetype = 'dashed')  +
  theme(axis.text = element_text(size = 6), 
        axis.title = element_text(size = 8),
        plot.margin = margin(0, 0, 0, 0),
        legend.margin = margin(0, 0, 0, -8),
        legend.key.size = unit(0.5, 'point'),
        legend.box.margin = margin(-10, -10, 0, 0),
        legend.text = element_text(size = 5), 
        legend.title = element_text(size = 5),
        legend.position = "bottom")
ggsave("plots/fig5e.pdf", width = 54, height = 77, units = "mm")
write.csv(fxns_bulk_fov,"sourcedata/fig5e.csv", row.names = FALSE)

# stats
qqnorm(fxns_bulk_fov$ultraplaus)
qqline(fxns_bulk_fov$ultraplaus, col = "red")

# Histogram
hist(fxns_bulk_fov$ultraplaus, breaks = 30, main = "Histogram of ultraplaus", xlab = "ultraplaus")

time = '168'
g1 <- subset(fxns_bulk_fov, condition == 'infected' & hpi_category == time)$ultraplaus
g2 <- subset(fxns_bulk_fov, condition == 'mock' & hpi_category == time)$ultraplaus
wilcox_test_result <- wilcox.test(g1, g2)
print(wilcox_test_result)

# Manual calculation of adjusted p-value using Bonferroni correction
# Number of comparisons (here, assuming one test per hpi category)
num_comparisons <- length(unique(fxns_bulk_fov$hpi_category))
adjusted_p_value <- p.adjust(wilcox_test_result$p.value, method = "bonferroni", n = num_comparisons)
print(adjusted_p_value)

################################
## Fig. 5f. CBF of Ultraplaus ##
################################
# see figure5_largedataset.R
# best done on a workstation.

#################################################################
## Fig. 5g. Fxn beating surface area by gfp category over time ##
#################################################################

fxns_labelled_fov <- fread('data/big_fxns_labelled_byvidpoint.csv')
fxns_labelled_culturepoint <- fread('data/big_fxns_labelled_byculturepoint.csv')
fxns_labelled_fov$hpi_category <- factor(fxns_labelled_fov$hpi_category, levels = c("early", "24", "48", "72", "96", "120", "168"))
fxns_labelled_culturepoint$hpi_category <- factor(fxns_labelled_culturepoint$hpi_category, levels = c("early", "24", "48", "72", "96", "120", "168"))

level_order <- c('infected', 'near', 'distant', 'mock') 

fxns_labelled_fov <- subset(fxns_labelled_fov, hpi_category != 120)
fxns_labelled_culturepoint <- subset(fxns_labelled_culturepoint, hpi_category != 120)
write.csv(fxns_labelled_fov,"sourcedata/fig5g.csv", row.names = FALSE)

ggplot(subset(fxns_labelled_fov, hpi > 12), aes(x = hpi_category, color = factor(status, levels = level_order), y = ultraplaus)) +
  geom_jitter(alpha = 0.15, size = 0.4, position = position_jitterdodge(dodge.width = 0.7))+
  geom_point(data = subset(fxns_labelled_culturepoint, hpi > 12), size = 0.4, alpha = 0.4, position = position_jitterdodge(dodge.width = 0.7), aes(shape = as.factor(donor)))+
  labs(y = 'Fraction Beating Surface Area', x = 'Hours Post Infection', shape = 'Donor:', color = 'Condition:') +
  stat_summary(fun = mean,
               fun.min = function(x) mean(x) - sd(x), 
               fun.max = function(x) mean(x) + sd(x),
               size = 0.4, 
               alpha = 0.8,
               linewidth = 0.5,
               position = position_dodge(width = 0.7),
               geom = "pointrange") +
  theme_bw() +
  ylim(0, 1.1) +
  #geom_hline(yintercept = 0.2166, linetype = 'dashed')  +
  scale_color_manual(values=c('chartreuse3', 'darkorchid4', 'darkorchid2', 'deeppink1'), labels = c("GFP+", "Near GFP-", "Far GFP-", "Mock"), name = NULL) +
  theme(axis.text = element_text(size = 6), 
        axis.title = element_text(size = 8),
        plot.margin = margin(1, 1, 1, 1),
        legend.margin = margin(0, 0, 0, -8),
        legend.key.size = unit(0.5, 'point'),
        legend.box.margin = margin(-10, -10, 0, 0),
        legend.text = element_text(size = 5), 
        legend.title = element_text(size = 5),
        legend.position = "bottom")
ggsave("plots/fig5g.pdf", width = 90, height = 77, units = "mm")

 # stats
qqnorm(fxns_labelled_fov$ultraplaus)
qqline(fxns_labelled_fov$ultraplaus, col = "red")
# Histogram
hist(fxns_labelled_fov$ultraplaus, breaks = 30, main = "Histogram of ultraplaus", xlab = "ultraplaus")
# again, not normal
time = '72'

g1 <- subset(fxns_labelled_fov, status == 'infected' & hpi_category == time)$ultraplaus
g2 <- subset(fxns_labelled_fov, status == 'near' & hpi_category == time)$ultraplaus
g3 <- subset(fxns_labelled_fov, status == 'distant' & hpi_category == time)$ultraplaus
g4 <- subset(fxns_labelled_fov, status == 'mock' & hpi_category == time)$ultraplaus

time = '168'

#g1 <- subset(fxns_labelled_culturepoint, status == 'infected' & hpi_category == time)$ultraplaus
#g2 <- subset(fxns_labelled_culturepoint, status == 'near' & hpi_category == time)$ultraplaus
#g3 <- subset(fxns_labelled_culturepoint, status == 'distant' & hpi_category == time)$ultraplaus
#g4 <- subset(fxns_labelled_culturepoint, status == 'mock' & hpi_category == time)$ultraplaus

wilcox_test_result <- wilcox.test(g1, g4)
print(wilcox_test_result)
adjusted_p_value <- p.adjust(wilcox_test_result$p.value, method = "bonferroni", n = 24)
print(adjusted_p_value)

wilcox_test_result <- wilcox.test(g1, g3)
print(wilcox_test_result)
adjusted_p_value <- p.adjust(wilcox_test_result$p.value, method = "bonferroni", n = 24)
print(adjusted_p_value)

wilcox_test_result <- wilcox.test(g3, g4)
print(wilcox_test_result)
adjusted_p_value <- p.adjust(wilcox_test_result$p.value, method = "bonferroni", n = 24)
print(adjusted_p_value)

wilcox_test_result <- wilcox.test(g1, g2)
print(wilcox_test_result)
adjusted_p_value <- p.adjust(wilcox_test_result$p.value, method = "bonferroni", n = 24)
print(adjusted_p_value)

glm_model <- glm(log(ultraplaus) ~ hpi_category + status, family = inverse.gaussian(link = "log"), data = fxns_labelled_fov)
summary(glm_model)

ok <- fxns_bulk_fov %>%
  group_by(experiment, culturepoint, condition, hpi, donor) %>%
  summarise(
    count = n()
  )

summary <- data %>%
  group_by(donor, timepoint, fov, experiment, condition) %>%
  summarise(count = n())

