# Plots for ciliary motion KOs (Figs. 5 & 6; Supp. Figs. 9 & 10)

# libraries
library(readr)
library(ggplot2)
library(dplyr)
library(ggbeeswarm)
#library(viridis)
library(effsize)

colors <- c("NT" = "#EF756A", "CYPA" = "#80AE00", "DNAH5" = "#3FBFC5", "DNAI1" = "#C37DFF", 'Shuffled' = 'black')
#scale_color_manual(values = colors, labels = c("nt" = "NT", "cypa" = "CYPA", "dnah5" = "DNAH5", "dnai1" = "DNAI1", "Shuffled" = "Shuffled")) +
###############################################
## 4d. Beat frequency & power of KOs vs. WT. ##
###############################################

# loading in data
ciliako2_concatenated <- read_csv("data/ciliako2_concatenated.csv")
ciliako44_concatenated <- read_csv("data/ciliako44_concatenated.csv")
# i named the kos inconsistently :/
# and dnai1-5 is the culture in ciliako2 where the well flooded before infecting so no useful data
ciliako2_fixed <- ciliako2_concatenated %>%
  mutate(condition = case_when(
    condition == 'dnah5-1' ~ 'dnah5',
    condition == 'dnai1-2' ~ 'dnai1',
    TRUE ~ condition
  )) %>%
  filter(condition != 'dnai1-5')

# join them
all_kos <- rbind(ciliako2_fixed, ciliako44_concatenated)
all_kos$conditionn <- factor(all_kos$condition, levels = c("nt", "cypa", "dnah5", "dnai1"))

# for testing plots. bc there are simply too many points.
subset_data <- all_kos %>% sample_frac(0.01)

# for overlaying FOV dots for statistical tests & vis
kos_resultbyfov <- all_kos %>%
  group_by(experiment, vidpoint) %>%
  summarise(total = n(), 
            medpower = median(power), 
            medfreq = median(domfreq), 
            conditionn = first(conditionn),
            donor = first(donor),
            video = first(video),
            .groups = 'drop')

ggplot(kos_resultbyfov, aes(x = conditionn, y = medpower, color = conditionn)) +
  geom_jitter()


# violin plot
ggplot(subset_data, aes(as.factor(conditionn), domfreq)) +
  xlab("KO") +
  ylab("Dominant pixel frequency (Hz)") +
  scale_x_discrete(labels = c("NT", "CYPA", "DNAH5", "DNAI1")) +
  theme_bw() +
  ylim(0, 80) +
  #ylim(0, 20) +
  #facet_wrap(~donor)+
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), adjust = 1.5, width = 0.8) +
  geom_quasirandom(data = kos_resultbyfov, width = 0.2, alpha = 0.1, aes(color = conditionn, x = conditionn, y = medfreq)) +
  #scale_color_viridis()+
  #stat_summary(fun = median,
  #             fun.min = function(x) quantile(x, 0.25), 
  #             fun.max = function(x) quantile(x, 0.75), 
  #             geom = "pointrange", position = position_dodge(width = 0.9), color = '#05e895', size = 0.1) +
   theme(axis.text = element_text(size=6, face = 'italic'), axis.title = element_text(size=6),
        plot.margin = margin(1, 0, 0, 0),
        legend.position = 'none')
ggsave("plots/fig4d.pdf", width = 1.6, height = 1.5, units = "in")
write.csv(subset_data,"sourcedata/fig4d_pix.csv", row.names = FALSE)
write.csv(kos_resultbyfov, 'sourcedata/fig4d_fov.csv', row.names = FALSE)

g1 <- subset(all_kos, condition == 'nt')$domfreq
g2 <- subset(all_kos, condition == 'cypa')$domfreq
g3 <- subset(all_kos, condition == 'dnah5')$domfreq
g4 <- subset(all_kos, condition == 'dnai1')$domfreq
cohen.d(g1, g2) # -0.129; -0.129 -0.127 negligible
cohen.d(g1, g3) #  -1.162; -1.164 -1.161 large
cohen.d(g1, g4) # -0.939; -0.930 -0.928 large



r##############################
## 4e, Plot o power density ##
##############################

# density plot of power
ggplot(all_kos, aes(color = as.factor(conditionn), fill = as.factor(conditionn), power)) +
  labs(x = 'Power Density\n(pixel units^2/Hz)', y ="Density", color = "KO:")+
  theme_bw() +
  geom_density(alpha = 0.1, show.legend = c(color = TRUE, fill = FALSE)) +
  scale_x_log10() +
  geom_vline(xintercept = 2e-4, linetype = 'dashed') +
  scale_color_discrete(labels = c('nt' ="NT", 'cypa' = "CYPA", 'dnah5' = "DNAH5", 'dnai1' = "DNAI1"))+
  theme(axis.text = element_text(size=6), axis.title = element_text(size=6),
        legend.text = element_text(size = 6, face = 'italic'),
        plot.margin = margin(1, 0, 0, 0),
        legend.title = element_text(size = 6),
        legend.key.size = unit(0.5, 'lines'),
        legend.margin = margin(0, 0, 0, 0),
        legend.spacing.x = unit(0, 'lines'),
        legend.box.margin = margin(-8, 0, 0, 0))
ggsave("plots/fig4e.pdf", width = 2, height = 1.5, units = "in")

g1 <- subset(all_kos, condition == 'nt')$power
g2 <- subset(all_kos, condition == 'cypa')$power
g3 <- subset(all_kos, condition == 'dnah5')$power
g4 <- subset(all_kos, condition == 'dnai1')$power
cohen.d(g1, g2) # 0.14; 0.139-0.141; negligible
cohen.d(g1, g3) #  0.21; 0.205-0.207
cohen.d(g1, g4) # 

 fxns <- all_kos %>%
  group_by(conditionn) %>%
  summarize(npow = sum(power > 2e-4), 
            total = n(), 
            fxnpow = npow/total, 
            powerbeatin = sum(domfreq < 18 & power > 2e-4), 
            fxnpowerbeat = powerbeatin/npow, 
            naivebeatin = sum(domfreq < 18), 
            fxnnaivebeat = naivebeatin/total)


#########################################
## 4h. # spots over time in KOs vs WT. ##
#########################################
spots <- read_csv("data/240702_goodspots.csv")
kos <- subset(spots, expt.name == 'ciliako2' | expt.name == 'ciliako44')
kos <- kos %>%
  mutate(condition = case_when(
    condition == 'nt-3' ~ 'nt',
    condition == 'cypa-2' ~ 'cypa',
    condition == 'dnah5-1' ~ 'dnah5',
    condition == 'dnai1-2' ~ 'dnai1',
    TRUE ~ condition
  )) %>%
  filter(usable != 0)
kos$condition <- factor(kos$condition, levels = c("nt", "cypa", "dnah5", "dnai1"))

kos$combid <- paste(kos$expt.name, kos$sample.no)

ggplot(kos, aes(time, as.numeric(norfreq), group = condition, colour=factor(condition), shape = factor(donor))) +
  #scale_color_met_d("Hiroshige", labels = c("1% agarose", "2% agarose", "Mock", "Plain"))+
  labs(x = 'Hours post infection', y ="# of GFP+ spots", color = "KO:", shape = 'Donor:')+
  scale_color_discrete(labels = c("NT", "CYPA", "DNAH5", "DNAI1"))+
  theme_bw() + 
  guides(colour = guide_legend(override.aes = list(size=1), order = 1),
         shape = guide_legend(override.aes = list(size = 1),
                              order = 2,
                              label.theme = element_text(face = 'plain', size = 5) )) +
  geom_point(alpha = 0.4, size = 0.4) +
  #geom_line(alpha = 0.4, aes(group = combid)) +
  stat_smooth(se = TRUE) +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=8),
        plot.margin = margin(0, 0, 0, 0),
        legend.key.size = unit(0.5, 'point'),
        legend.box.margin = margin(0, 0, -10, -10),
        #legend.background = element_rect(fill = NA, color == NA),
        legend.position = c(0.2, 0.7),
        legend.spacing.y = unit(0.01, 'cm'),
        legend.text = element_text(size=5, face = 'italic'), legend.title = element_text(size=5))
ggsave("plots/fig4g.pdf", width = 2.0, height = 1.6, units = "in")
f4h.small <- kos[, c("csv", "donor", "expt.name", "sample.no", "time", 'norfreq', 'condition')]
write.csv(f4h.small,"sourcedata/fig4h.csv", row.names = FALSE)


# kinetics....
# Fit splines and calculate the rate of change
fit.splines <- function(df, spar.value = 0.1) {
  
  df <- df %>%
    arrange(time)
  
  # Fit a natural cubic spline
  #spline_fit <- smooth.spline(df$time, df$norfreq)
  spline.fit <- smooth.spline(df$time, df$norfreq, spar = spar.value)
  # Predict the fitted values and first derivative (rate of change)
  fitted.values <- predict(spline.fit, df$time)
  first.derivative <- predict(spline.fit, df$time, deriv = 1)
  second.derivative <- predict(spline.fit, df$time, deriv = 2)
  
  df <- df %>%
    mutate(fitted.norfreq = fitted.values$y,
           rate.of.change = first.derivative$y,
           accel = second.derivative$y)
  
  return(df)
}

# Apply the function to each culture
data <- kos %>%
  group_by(csv) %>%
  do(fit.splines(.))

kos <- kos %>%
  left_join(data %>% select(csv, time, fitted.norfreq, rate.of.change, accel), 
            by = c("csv", "time"))

kos$combid <- paste(kos$expt.name, kos$sample.no)
ggplot(kos, aes(time, as.numeric(rate.of.change), colour=factor(condition), group = factor(combid))) +
  labs(x = 'Hours post infection', y =expression(Delta ~ "GFP+ spots per hour"), color = "Donor:")+
  geom_hline(yintercept = 0, color = 'black', linetype = 'dashed') + 
  geom_line(alpha = 0.4) +
  theme_bw() + 
  xlim(0, 120) +
  #stat_summary(aes(group = donor)) +
  stat_smooth(method = "loess", span = 0.1, aes(group = condition)) +
  guides(colour = guide_legend(override.aes = list(size=2))) +
  #geom_point(alpha = 0.6, size = 0.01) +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=8),
        plot.margin = margin(1, 6, 0, 1),
        legend.position = 'none')
        #legend.key.size = unit(0.1, 'point'),
        #legend.box.margin = margin(0, 0, 0, -10),
        #legend.position = c(0.2, 0.8),
        #legend.text = element_text(size=6), legend.title = element_text(size=6))
ggsave("plots/fig4i.pdf", width = 1.6, height = 1.6, units = "in")
f4i.small <- kos[, c("csv", "donor", "expt.name", "sample.no", "time", 'rate.of.change', 'condition')]
write.csv(f4i.small,"sourcedata/fig4i.csv", row.names = FALSE)


#########################
## 4j. N copies vs. KO ##
#########################

# for info about each culture
exptlog_iii <- read_csv("data/exptlog_iii.csv")
# take the spots over time csv and extract specific values of interest
daily.spots.summary <- kos %>%
  group_by(csv) %>%
  summarise(good.peak.spots = max(norfreq), 
            max.rate = max(rate.of.change),
            good.peak.spots.time = time[which.max(norfreq)], 
            onedayspots = norfreq[time == 23],
            twodayspots = norfreq[time == 47],
            threedayspots = norfreq[time == 71],
            fourdayspots = norfreq[time == 95],
            fivedayspots = norfreq[time == 113], # i am aware this is not 119. however this is the latest timepoint for ciliako44.
            aucspots = sum(norfreq)) 
# cross back with the info i want
daily.spots.summary.j <- inner_join(daily.spots.summary, exptlog_iii)
# only take the ciliary motion kos
ko.spots.summary.j <- subset(daily.spots.summary.j, expt.name == 'ciliako2' | expt.name == 'ciliako44')
# add back the couple cultures that are missing spots info (ciliako44, 19 and 20 got got by an x motor issue)
missing <- subset(exptlog_iii, expt.name == 'ciliako44' & (sample.no == 19 | sample.no == 20))
missing <- missing %>%
  mutate(
    good.peak.spots = NA,
    good.peak.spots.time = NA,
    onedayspots = NA,
    twodayspots = NA,
    threedayspots = NA,
    fourdayspots = NA,
    fivedayspots = NA,
    aucspots = NA)
ko.spots.all.j <- bind_rows(ko.spots.summary.j, missing)
# rename these two columns that differ between this frame and the cbf dataset
ko.spots.all.j <- ko.spots.all.j %>%
  rename(
    experiment = expt.name,
    culturepoint = sample.no
  )
# now make summary data about the cbfs
kos_resultbyculture <- all_kos %>%
  group_by(experiment, culturepoint) %>%
  summarise(total = n(), 
            medpower = median(power), 
            medfreq = median(domfreq), 
            fracbeat = sum(domfreq < 15 & power > 2e-4)/n(),
            fracpower = sum(power > 2e-4)/n(),
            fracplaus = sum(domfreq < 15)/n(),
            condition = first(condition),
            donor = first(donor),
            .groups = 'drop')
bigkodata <- inner_join(kos_resultbyculture, ko.spots.all.j, by = c("experiment", "culturepoint"))


# fix the little things that will annoy me forever graphing.
bigkodata$condition <- factor(bigkodata$condition.x, levels = c("nt", "cypa", "dnah5", "dnai1"))
bigkodata$n.mm.2mean <- as.numeric(bigkodata$n.mm.2mean)
# making sure that the numbers aren't super super messed up here...
min(bigkodata$n.mm.2mean)

# the graph.
ggplot(bigkodata, aes(condition, n.mm.2mean, shape = factor(donor.y), color = condition, group = condition)) +
   labs(x = 'KO', y = bquote('N copies per'~mm^2))+
  scale_y_log10(limits = c(1e7, 3e9)) +
  theme_bw() +
  geom_jitter(width = 0.1, alpha = 0.7) +
  stat_summary(fun = mean,
               fun.min = function(x) mean(x) - sd(x), 
               fun.max = function(x) mean(x) + sd(x),
               size = 0.25, 
               alpha = 0.8,
               linewidth = 0.5, 
               geom = "pointrange",
               color = 'black') + 
  scale_x_discrete(labels = c("NT", "CYPA", "DNAH5", "DNAI1")) +
  theme(axis.text = element_text(size=5, face = 'italic'), axis.title = element_text(size=6),
        plot.margin = margin(1, 0, 0, 0),
        #legend.position = c(0.8, 0.8),
        legend.position = 'none')
ggsave("plots/fig4j.pdf", width = 1.5, height = 1.2, units = "in")
f4j.small <- bigkodata[, c("experiment", "culturepoint", "donor.y", "n.mm.2mean", 'condition')]
write.csv(f4j.small,"sourcedata/fig4j.csv", row.names = FALSE)


anova_result <- aov(n.mm.2mean ~ condition, data = bigkodata)
summary(anova_result) # p = 0.0145
tukey_result <- TukeyHSD(anova_result)
print(tukey_result)
# cypa-nt 0.99
# nt-dnah5 0.03
# nt-dnai1 0.25
plot(anova_result, 1)
library(car)
leveneTest(as.numeric(n.mm.2mean) ~ condition, data = bigkodata)
# tests out ok.

plot(anova_result, 2)
# looks ok 2 me
anova_residuals <- residuals(object = anova_result )
# shapiro wilk
shapiro.test(x = anova_residuals )

####################################
##  Supp. Fig. 9: Fxn ciliation   ##
####################################
dataa = subset(daily.spots.summary.j, 
               (expt.name == 'ciliako2' | expt.name == 'ciliako44')) 
dataa$conditionn <- factor(dataa$condition, levels = c("nt-3", "cypa-2", "dnah5-1", "dnai1-2"))

ggplot(dataa, aes(y=as.numeric(t0spy350)/as.numeric(t0spy0), x = conditionn, color = conditionn, shape = factor(donor))) +
  labs(y = 'Fraction SPY650-Tubulin+ Area', x = 'Condition') +
  geom_jitter(alpha = 0.4, width = 0.2) +
  theme_bw() +
  ylim(0.5, 1.1) +
  stat_summary(fun = mean,
               fun.min = function(x) mean(x) - sd(x), 
               fun.max = function(x) mean(x) + sd(x),
               size = 0.25, 
               alpha = 0.8,
               linewidth = 0.5, 
               geom = "pointrange",
               color = 'black', 
               aes(group = conditionn)) + 
  scale_x_discrete(labels = c("NT", "CYPA", "DNAH5", "DNAI1")) +
  theme(axis.text = element_text(size=5, face = 'italic'), axis.title = element_text(size=6),
        plot.margin = margin(1, 0, 0, 0),
        #legend.position = c(0.8, 0.8),
        legend.position = 'none')
ggsave("plots/s9.pdf", width = 2, height = 2, units = "in")
s9a.small <- dataa[, c("expt.name", "sample.no", "donor", "conditionn", 't0spy350', 't0spy0', 'n.mm.2mean', 'good.peak.spots')]
write.csv(s9a.small,"sourcedata/s9a.csv", row.names = FALSE)

median(as.numeric(dataa$t0spy350)/as.numeric(dataa$t0spy0))

anova_result <- aov(as.numeric(t0spy350)/as.numeric(t0spy0) ~ condition, data = dataa)
summary(anova_result) # p = 0.0145
tukey_result <- TukeyHSD(anova_result)
print(tukey_result)


ggplot(dataa, aes(x=as.numeric(t0spy350)/as.numeric(t0spy0), y = as.numeric(n.mm.2mean), color = conditionn, shape = factor(donor))) +
  labs(x = 'Fraction SPY650-Tubulin+ Area', y = 'N per mm^2 culture area') +
  geom_point(alpha = 1, width = 0.2) +
  theme_bw() +
  stat_smooth(method = 'lm', aes(group = 'none')) +
  theme(axis.text = element_text(size=5, face = 'italic'), axis.title = element_text(size=6),
        plot.margin = margin(1, 0, 0, 0),
        #legend.position = c(0.8, 0.8),
        legend.position = 'none')
ggsave("C:/Users/Mark/Desktop/figure_images/4/s9b.pdf", width = 2, height = 2, units = "in")

ggplot(dataa, aes(x=as.numeric(t0spy350)/as.numeric(t0spy0), y = good.peak.spots, color = conditionn, shape = factor(donor))) +
  labs(x = 'Fraction SPY650-Tubulin+ Area', y = 'Peak GFP+ Spots') +
  geom_point(alpha = 1, width = 0.2) +
  theme_bw() +
  stat_smooth(method = 'lm', aes(group = 'none')) +
  theme(axis.text = element_text(size=5, face = 'italic'), axis.title = element_text(size=6),
        plot.margin = margin(1, 0, 0, 0),
        #legend.position = c(0.8, 0.8),
        legend.position = 'none')
ggsave("C:/Users/Mark/Desktop/figure_images/4/s9c.pdf", width = 2, height = 2, units = "in")

g1 <- dataa$t0spy350/dataa$t0spy0
g2 <- dataa$good.peak.spots
g3 <- as.numeric(dataa$n.mm.2mean)
cor.test(g1, g2, method = 'pearson')


#######
# what if we look at peak # spots? instead of n.
anova_result <- aov(good.peak.spots~condition, data = bigkodata)
summary(anova_result)
# p = 0.00117
tukey_result <- TukeyHSD(anova_result)
print(tukey_result)
# dnah5 and dnai1 both significantly differ from nt & cypa.
plot(anova_result, 1)
library(car)
leveneTest(as.numeric(good.peak.spots) ~ condition, data = merged)
# tests out ok.

plot(anova_result, 2)
# looks ok 2 me
aov_residuals <- residuals(object = anova_result )
# Run Shapiro-Wilk test
shapiro.test(x = aov_residuals )
# still looks ok! inchresting. graph it:
ggplot(bigkodata, aes(condition, good.peak.spots, shape = factor(donor.y), color = condition, group = condition)) +
  labs(x = 'KO', y = 'Peak # GFP+ spots')+
  #scale_y_log10(limits = c(1e7, 3e9)) +
  ylim(0, 40000) +
  theme_bw() +
  geom_jitter(width = 0.1, alpha = 0.7) +
  stat_summary(fun = mean,
               fun.min = function(x) mean(x) - sd(x), 
               fun.max = function(x) mean(x) + sd(x),
               size = 0.25, 
               alpha = 0.8,
               linewidth = 0.5, 
               geom = "pointrange",
               color = 'black') + 
  scale_x_discrete(labels = c("NT", "CYPA", "DNAH5", "DNAI1")) +
  theme(axis.text = element_text(size=5, face = 'italic'), axis.title = element_text(size=6),
        plot.margin = margin(1, 0, 0, 0),
        #legend.position = c(0.8, 0.8),
        legend.position = 'none')
ggsave("C:/Users/Mark/Desktop/figure_images/4/fig4j.pdf", width = 1.3, height = 1.2, units = "in")
#coolio! it looks v. different. how does good.peak.spots compare to n.mm.2mean?


###############################
## 4k. spots vs. fxn beating ##
###############################

g1 <- bigkodata$good.peak.spots
g2 <- bigkodata$fracplaus
cor.test(g1, g2, method = 'pearson')

ggplot(bigkodata, aes(fracplaus, as.numeric(fourdayspots), color = condition, shape = factor(donor.y))) +
  labs(x = 'Fraction Beating Area', y ="Peak # GFP+ spots")+
  geom_point(alpha = 0.8) +
  theme_bw() +
  geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, aes(group = FALSE, color = NULL), show.legend = FALSE) +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=6),
        plot.margin = margin(1, 0, 0, 0),
        legend.key.size = unit(0.1, 'point'),
        legend.box.margin = margin(0, 0, 0, -10),
        legend.position = 'none',
        legend.text = element_text(size=6, face = 'italic'), legend.title = element_text(size=6))
ggsave("plots/fig4k.pdf", width = 1.5, height = 1.2, units = "in")
f4k.small <- bigkodata[, c("experiment", "culturepoint", "donor.y", "fourdayspots", "fracplaus", 'condition')]
write.csv(f4k.small,"sourcedata/fig4k.csv", row.names = FALSE)

library(broom)
fit <- lm(good.peak.spots ~ fracplaus, data = bigkodata)
summary(fit)


#######################################
## S10: Spots over time significance ##
#######################################
library(tidyr)
library(ggplot2)

# Reshape your data to a long format
bigkodata_long <- bigkodata %>%
  pivot_longer(cols = c(onedayspots, twodayspots, threedayspots, fourdayspots, fivedayspots),
               names_to = "day",
               values_to = "spots")
bigkodata_long$day <- factor(bigkodata_long$day, 
                             levels = c("onedayspots", "twodayspots", "threedayspots", "fourdayspots", "fivedayspots"),
                             labels = c("1 Day Post Infection", "2 Days Post Infection", "3 Days Post Infection", "4 Days Post Infection", "5 Days Post Infection"))
# set y factor to give me a little room for plotting while keeping the axis free
max_values <- bigkodata_long %>%
  group_by(day) %>%
  summarize(max_spots = max(spots, na.rm = TRUE)) %>%
  mutate(max_spots_with_padding = max_spots * 1.5)

# Plot using facets
ggplot(bigkodata_long, aes(x = condition, y = spots, shape = factor(donor.y), color = condition, group = condition)) +
  labs(x = 'KO', y = 'Peak # GFP+ spots') +
  theme_bw() +
  geom_jitter(width = 0.1, alpha = 0.7) +
  stat_summary(fun = mean,
               fun.min = function(x) mean(x) - sd(x), 
               fun.max = function(x) mean(x) + sd(x),
               size = 0.25, 
               alpha = 0.8,
               linewidth = 0.5, 
               geom = "pointrange",
               color = 'black') + 
  scale_x_discrete(labels = c("NT", "CYPA", "DNAH5", "DNAI1")) +
  facet_wrap(~ day, scales = "free_y", nrow = 1) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.3))) +  # This adds extra space on the upper end
  theme(axis.text = element_text(size = 5, face = 'italic'), 
        axis.title = element_text(size = 6),
        plot.margin = margin(1, 0, 0, 0),
        strip.text = element_text(size = 6),
        legend.position = 'none')
ggsave("plots/s9_alldays.pdf", width = 6.8, height = 2.2, units = "in")
s10.small <- bigkodata_long[, c("experiment", "culturepoint", "donor.y", "condition", "spots", "day")]
write.csv(s10.small,"sourcedata/s10.csv", row.names = FALSE)



######################
#######################r
# things that were interesting but not going in the final version.

# graph o power done just like the cbfs.
ggplot(subset_data, aes(as.factor(conditionn), power)) +
  xlab("KO") +
  labs(x = 'KO', y = bquote('Power Density (unitless^2/Hz)'))+
  #ylab(bquote("Power Density "~(Unitless^2~" /Hz") +
  scale_x_discrete(labels = c("NT", "CYPA", "DNAH5", "DNAI1")) +
  theme_bw() +
  scale_y_log10()+
  #geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  #geom_quasirandom(data = kos_resultbyfov, width = 0.2, alpha = 0.2, aes(color = conditionn, x = conditionn, y = medpower)) +
  #scale_color_viridis()+
  #stat_summary(fun = median,
  #             fun.min = function(x) quantile(x, 0.25), 
  #             fun.max = function(x) quantile(x, 0.75), 
  #             geom = "pointrange", position = position_dodge(width = 0.9), color = '#05e895', size = 0.1) +
  theme(axis.text = element_text(size=6, face = 'italic'), axis.title = element_text(size=6),
        plot.margin = margin(1, 0, 0, 0),
        legend.position = 'none')
ggsave("C:/Users/Mark/Desktop/figure_images/4/fig4etest.pdf", width = 1.6, height = 1.5, units = "in")


####
# Add a line profile for the shuffled denosied ciliako44 dataset. as an example of what no signal looks like.
# i think this would be ideal if it worked but it was turning into a huge time suck to get it formatted nice.

shuffled <- read_csv("D:/shufflin/shuffleset.csv")

all_kos <- all_kos %>% mutate(source = "Raw")
shuffled <- shuffled %>% mutate(source = "Shuffled")
shuffled <- subset(shuffled, process == 'ck.denoised.frameshuf.smal')
shuffled$condition <- 'Shuffled'


colors <- c("NT" = "#EF756A", "CYPA" = "#80AE00", "DNAH5" = "#3FBFC5", "DNAI1" = "#C37DFF", 'Shuffled' = 'black')
fills <- c("NT" = "#EF756A", "CYPA" = "#80AE00", "DNAH5" = "#3FBFC5", "DNAI1" = "#C37DFF", 'Shuffled' = NA)

ggplot(combined_data, aes(x = power, color = as.factor(condition), fill = as.factor(condition))) +
  geom_density(aes(linetype = source), alpha = 0.1) +
  labs(x = 'Power Density\n(pixel units^2/Hz)', y = "Density", color = "KO:") +
  theme_bw() +
  scale_x_log10() +
  scale_color_manual(values = colors, labels = c("nt" = "NT", "cypa" = "CYPA", "dnah5" = "DNAH5", "dnai1" = "DNAI1", "Shuffled" = "Shuffled")) +
  scale_fill_manual(values = fills, labels = c("nt" = "NT", "cypa" = "CYPA", "dnah5" = "DNAH5", "dnai1" = "DNAI1", "Shuffled" = "Shuffled")) +
  geom_vline(xintercept = 2e-4) +
  #scale_color_discrete(labels = c("NT", "CYPA", "DNAH5", "DNAI1")) +
  theme(
    axis.text = element_text(size = 6, face = 'italic'),
    axis.title = element_text(size = 6),
    plot.margin = margin(1, 0, 0, 0),
    legend.title = element_text(size = 6),
    legend.key.size = unit(0.5, 'lines'),
    legend.margin = margin(0, 0, 0, 0),
    legend.spacing.x = unit(0, 'lines'),
    legend.box.margin = margin(-8, 0, 0, 0),
    legend.text = element_text(size = 6)
  )

ggplot(combined_data, aes(x = condition, y = power, color = as.factor(condition), fill = as.factor(condition))) +
  geom_violin() +
  labs(x = 'Power Density\n(pixel units^2/Hz)', y = "Density", color = "KO:") +
  theme_bw() +
  scale_y_log10()+
  #geom_vline(xintercept = 2e-4) +
  #scale_color_discrete(labels = c("NT", "CYPA", "DNAH5", "DNAI1")) +
  theme(
    axis.text = element_text(size = 6, face = 'italic'),
    axis.title = element_text(size = 6),
    plot.margin = margin(1, 0, 0, 0),
    legend.title = element_text(size = 6),
    legend.key.size = unit(0.5, 'lines'),
    legend.margin = margin(0, 0, 0, 0),
    legend.spacing.x = unit(0, 'lines'),
    legend.box.margin = margin(-8, 0, 0, 0),
    legend.text = element_text(size = 6)
  )

relevant <- subset(combined, process == 'ck.denoised.frameshuf.smal')
ggplot(relevant, aes(x = power, color = condition)) +
  #geom_point(alpha = 0.001) +  
  #geom_density_2d(color = 'blue') +
  geom_density() +
  geom_vline(xintercept = 2e-4) +
  scale_x_log10()

###
###
# Playing with the fxn beating area vs. spread numbers.




ggplot(merged, aes(y = as.numeric(good.peak.spots), x = fracbeat, color = condition)) +
  geom_point() +
  ylim(0, 30000) +
  geom_smooth(method = 'lm', formula = y ~ x, se = TRUE, aes(group = FALSE, color = NULL), show.legend = FALSE)
cor.test(as.numeric(merged$good.peak.spots), merged$fracbeat, method = 'pearson')

ggplot(merged, aes(y = as.numeric(n.mm.2mean), x = fracbeat, color = condition)) +
  geom_point() +
  scale_y_log10() +
  geom_smooth(method = 'lm', formula = y ~ x, se = TRUE, aes(group = FALSE, color = NULL), show.legend = FALSE)
cor.test(as.numeric(merged$n.mm.2mean), merged$fracbeat, method = 'pearson')



ggplot(merged, aes(x= as.numeric(n.mm.2mean), y = good.peak.spots, color = condition)) +
  geom_point() +
  geom_smooth(method = 'lm', formula = y ~ x, se = TRUE, aes(group = FALSE, color = NULL), show.legend = FALSE)

cor.test(as.numeric(merged$good.peak.spots), as.numeric(merged$n.mm.2mean), method = 'pearson')


cor.test(as.numeric(merged$good.peak.spots), merged$fracbeat, method = 'pearson')


### looking at all the different days with the goodpeakspots/n.mm.2 differences.
ggplot(bigkodata, aes(y=n.mm.2mean, x = good.peak.spots, color = condition, shape = factor(donor.y))) +
  geom_point()+
  geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, aes(group = FALSE, color = NULL), show.legend = FALSE)
reg <- lm(n.mm.2mean ~ good.peak.spots, data = daytest_merge)
summary(reg)


# how do the plots look at different days?
#dayone
ggplot(bigkodata, aes(y= onedayspots, x= condition, color = condition, group = condition, shape = factor(donor.y))) +
  theme_bw() +
  stat_summary(fun = mean,
               fun.min = function(x) mean(x) - sd(x), 
               fun.max = function(x) mean(x) + sd(x)) +
  geom_jitter(width = 0.1, alpha = 0.7)
anova_result <- aov(onedayspots~condition, data = bigkodata)
summary(anova_result)
# p = 0.63

#day two
ggplot(bigkodata, aes(y= twodayspots, x= condition, color = condition, group = condition, shape = factor(donor.y))) +
  theme_bw() +
  stat_summary(fun = mean,
               fun.min = function(x) mean(x) - sd(x), 
               fun.max = function(x) mean(x) + sd(x)) +
  geom_jitter(width = 0.1, alpha = 0.7)
# i can see a difference~ both kos very low, a couple low nts. cypa hoigh.
anova_result <- aov(twodayspots~condition, data = bigkodata)
summary(anova_result)
# p = 0.00192
tukey_result <- TukeyHSD(anova_result)
print(tukey_result)
# dnai5 and dnai1 sig dif cypa. cypa-nt 0.1.

# day three
ggplot(bigkodata, aes(y= threedayspots, x= condition, color = condition, group = condition, shape = factor(donor.y))) +
  theme_bw() +
  stat_summary(fun = mean,
               fun.min = function(x) mean(x) - sd(x), 
               fun.max = function(x) mean(x) + sd(x)) +
  geom_jitter(width = 0.1, alpha = 0.7)
# DIFFERENT
anova_result <- aov(threedayspots~condition, data = bigkodata)
summary(anova_result)
# p = 0.000171
tukey_result <- TukeyHSD(anova_result)
print(tukey_result)
# cypa-nt 0.11; dh5-nt di1-nt 0.045ish; dh5-cypa di1-cypa 0.0005ish.

# day four
ggplot(bigkodata, aes(y= fourdayspots, x= condition, color = condition, group = condition, shape = factor(donor.y))) +
  theme_bw() +
  stat_summary(fun = mean,
               fun.min = function(x) mean(x) - sd(x), 
               fun.max = function(x) mean(x) + sd(x)) +
  geom_jitter(width = 0.1, alpha = 0.7)
# DIFFERENT
anova_result <- aov(fourdayspots~condition, data = bigkodata)
summary(anova_result)
# p = 0.0008
tukey_result <- TukeyHSD(anova_result)
print(tukey_result)
# cypa-nt 0.68; dh5-nt di1-nt 0.027 ish; dh5-cypa di1-cypa 0.004ish.
# differences declining as WT and KO reach their peak while the kos continue rising.

#day five slash final measurements.
ggplot(bigkodata, aes(y= fivedayspots, x= condition, color = condition, group = condition, shape = factor(donor.y))) +
  theme_bw() +
  stat_summary(fun = mean,
               fun.min = function(x) mean(x) - sd(x), 
               fun.max = function(x) mean(x) + sd(x)) +
  geom_jitter(width = 0.1, alpha = 0.7)
# merging
anova_result <- aov(fivedayspots~condition, data = bigkodata)
summary(anova_result)
# p = 0.0124
tukey_result <- TukeyHSD(anova_result)
print(tukey_result)
# cypa-nt 0.99; dh5-nt di1-nt ~ 0.06; dh5-cypa di1-cypa ~ 0.07


# how does n.mm.2mean relate to different days, or the sum of days?

# sum of days first bc i can see an argument for this mattering
# if virions remain infectious in the mucus for extended periods of time
# which may or may not be true. ¯\_(ツ)_/¯
bigkodata$sumofdays <- bigkodata$onedayspots + bigkodata$twodayspots + bigkodata$threedayspots + bigkodata$fourdayspots + bigkodata$fivedayspots
ggplot(bigkodata, aes(y=n.mm.2mean, x = sumofdays, color = condition, shape = factor(donor.y))) +
  geom_point()+
  geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, aes(group = FALSE, color = NULL), show.legend = FALSE)
reg <- lm(n.mm.2mean ~ sumofdays, data = bigkodata)
summary(reg) # r sq 0.56
cor.test(x=bigkodata$sumofdays, bigkodata$n.mm.2mean, method = 'pearson')
#0.75, p = 0.9e-5

# one day
ggplot(bigkodata, aes(y=n.mm.2mean, x = onedayspots, color = condition, shape = factor(donor.y))) +
  geom_point()+
  geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, aes(group = FALSE, color = NULL), show.legend = FALSE)
reg <- lm(n.mm.2mean ~ onedayspots, data = bigkodata)
summary(reg) #r sq 0.043
cor.test(x=bigkodata$onedayspots, bigkodata$n.mm.2mean, method = 'pearson')
#0.2, p = 0.36

# two days
ggplot(bigkodata, aes(y=n.mm.2mean, x = twodayspots, color = condition, shape = factor(donor.y))) +
  geom_point()+
  geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, aes(group = FALSE, color = NULL), show.legend = FALSE)
# clearly and obviously not linear
reg <- lm(n.mm.2mean ~ twodayspots, data = bigkodata)
summary(reg) #r sq 0.35
cor.test(x=bigkodata$twodayspots, bigkodata$n.mm.2mean, method = 'pearson')
#0.59, p = 0.005

# three days
ggplot(bigkodata, aes(y=n.mm.2mean, x = threedayspots, color = condition, shape = factor(donor.y))) +
  geom_point()+
  geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, aes(group = FALSE, color = NULL), show.legend = FALSE)
reg <- lm(n.mm.2mean ~ threedayspots, data = bigkodata)
summary(reg) #r sq 0.53
cor.test(x=bigkodata$threedayspots, bigkodata$n.mm.2mean, method = 'pearson')
#0.72, p = 0.0002

# four days
ggplot(bigkodata, aes(y=n.mm.2mean, x = fourdayspots, color = condition, shape = factor(donor.y))) +
  geom_point()+
  geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, aes(group = FALSE, color = NULL), show.legend = FALSE)
reg <- lm(n.mm.2mean ~ fourdayspots, data = bigkodata)
summary(reg) #r sq 0.60
cor.test(x=bigkodata$fourdayspots, bigkodata$n.mm.2mean, method = 'pearson')
#0.77, p = 4e-5

# five days
ggplot(bigkodata, aes(y=n.mm.2mean, x = fivedayspots, color = condition, shape = factor(donor.y))) +
  geom_point()+
  geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, aes(group = FALSE, color = NULL), show.legend = FALSE)
reg <- lm(n.mm.2mean ~ fivedayspots, data = bigkodata)
summary(reg) #r sq 0.48
cor.test(x=bigkodata$fivedayspots, bigkodata$n.mm.2mean, method = 'pearson')
#0.70, p = 0.0004
# looking at all the graphs, the donor 7 nt that came out HOTT in the N department was never the spottiest,
# while the donor 3 nt that was always among the spottiest  just didn't have that much N. Interesting.

# peak spots
ggplot(bigkodata, aes(y=n.mm.2mean, x = good.peak.spots, color = condition, shape = factor(donor.y))) +
  geom_point()+
  geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, aes(group = FALSE, color = NULL), show.legend = FALSE)
reg <- lm(n.mm.2mean ~ good.peak.spots, data = bigkodata)
summary(reg) #r sq 0.59
cor.test(x=bigkodata$good.peak.spots, bigkodata$n.mm.2mean, method = 'pearson')
#0.77, p = 4.4e-5

# what if you look at just days 3&4, when most of the virus looks like it's being made?
bigkodata$peakdays <-  bigkodata$threedayspots + bigkodata$fourdayspots
ggplot(bigkodata, aes(y=n.mm.2mean, x = peakdays, color = condition, shape = factor(donor.y))) +
  geom_point()+
  geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, aes(group = FALSE, color = NULL), show.legend = FALSE)
reg <- lm(n.mm.2mean ~ peakdays, data = bigkodata)
summary(reg) #r sq 0.57
cor.test(x=bigkodata$peakdays, bigkodata$n.mm.2mean, method = 'pearson')
#0.76, p = 6e-5

# what if you look at days 4 & 5, for the freshest virus with a little points for freshness?
bigkodata$l8 <- bigkodata$fourdayspots + bigkodata$fivedayspots
ggplot(bigkodata, aes(y=n.mm.2mean, x = l8, color = condition, shape = factor(donor.y))) +
  geom_point()+
  geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, aes(group = FALSE, color = NULL), show.legend = FALSE)
reg <- lm(n.mm.2mean ~ l8, data = bigkodata)
summary(reg) #r sq 0.56
cor.test(x=bigkodata$l8, bigkodata$n.mm.2mean, method = 'pearson')
# 0.75, p = 0.0001

# what about the actual area under the curve?
ggplot(bigkodata, aes(y=n.mm.2mean, x = aucspots, color = condition, shape = factor(donor.y))) +
  geom_point()+
  geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, aes(group = FALSE, color = NULL), show.legend = FALSE)
reg <- lm(n.mm.2mean ~ aucspots, data = bigkodata)
summary(reg) #r sq 0.54
cor.test(x=bigkodata$aucspots, bigkodata$n.mm.2mean, method = 'pearson')
# 0.74, 0.0001
# not it. which kind of makes sense perhaps?

