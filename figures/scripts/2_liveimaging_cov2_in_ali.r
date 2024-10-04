# Live imaging of SARS-CoV-2 in ALI cultures (Fig. 3, Supp. Figs 4-7)

# libraries
library(readr)
library(ggplot2)
library(dplyr)
library(ggbeeswarm)

spots <- read_csv("data/240702_goodspots.csv")
# fix bad frames with way too high spots.
# these are frames where clearly there was an error.
spots <- spots %>%
  mutate(freq = case_when(
    freq > 50000 ~ 0,
    TRUE ~ freq)) %>%
  group_by(csv) %>%
  mutate(norfreq = as.numeric(freq) - min(as.numeric(freq)))

# removing problematic parts
spots <- subset(spots, 
                !(expt.name == 'd5d6' & FRAME > 58) & # erroneously set freq to 0; data does not exist
                !(expt.name == 'd4shake' & FRAME > 59) & #erroneously set freq to 0; data does not exist
                !(expt.name == 'd2shake' & (FRAME > 46 & FRAME < 50)) & # it was out of focus & detected no spots in this time
                !(expt.name == 'agarose_d3' & (FRAME > 35 & FRAME < 48))
                & !((expt.name == 'd7d8') & (FRAME > 60))) # after rinse and thus not apples to apples
# graphing to id problematic frame intervals.        
problem <- subset(spots, expt.name == 'd2d4ag')
ggplot(problem, aes(FRAME, as.numeric(norfreq), colour=factor(donor), group = interaction(expt.name, sample.no))) +
  labs(x = 'Hours post infection', y ="# of GFP+ spots", color = "Donor:")+
  theme_bw() +
  geom_vline(xintercept = 35) +
  geom_vline(xintercept = 48) +
  #geom_label() +
  #stat_smooth(method = "loess", se = TRUE, aes(group = 'none')) +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=8), 
        legend.key.size = unit(1, 'point'),
        legend.box.margin = margin(0, 0, 0, -10),
        legend.text = element_text(size=6), legend.title = element_text(size=8),
        plot.margin = margin(1, 0, 0, 0)) +
  guides(colour = guide_legend(override.aes = list(size=2))) +
  geom_line(alpha = 0.6, linewidth = 0.2)          


# measuring kinetics
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
data <- spots %>%
  group_by(csv) %>%
  do(fit.splines(.))

spots <- spots %>%
  left_join(data %>% select(csv, time, fitted.norfreq, rate.of.change, accel), 
            by = c("csv", "time"))



# for info about each culture
exptlog_iii <- read_csv("data/exptlog_iii.csv")
# take the spots over time csv and extract specific values of interest
daily.spots.summary <- spots %>%
  group_by(csv) %>%
  summarise(good.peak.spots = max(norfreq), 
            good.peak.spots.time = time[which.max(norfreq)],
            aucspots = sum(norfreq)) 
# cross back with the info i want
daily.spots.summary.j <- inner_join(daily.spots.summary, exptlog_iii)
daily.spots.summary.j <- daily.spots.summary.j %>%
  mutate(crypts = case_when(
    crypts == 'pointycrypts' ~ 'clean',
    TRUE ~ crypts
  ))
daily.spots.summary.j <- daily.spots.summary.j %>%
  mutate(mucus.simple = case_when(
    mucus.simple == 'static' ~ 'disorganized',
    TRUE ~ mucus.simple
  ))



# Figure 2c: fxn gfp+ culture area over time
relevant = subset(spots, 
                  (condition == 'vanilla' | condition == 'static' | condition == 'shaken') 
                  & !((expt.name == 'd7d8') & (time > 120)) # this one got rinsed then so not a vanilla infection
                  & !((expt.name == 'ruxcilia') & time > 72) # big gap after 72 on this one and it flooded at that point.
                  #& (fromold == 1) # in infections in the bsl3 core i couldn't rinse so it affected the kinetics.
                  & (!((expt.name == 'd5d6') & (sample.no == 2))) # this one i'm p sure was contaminated with non-reporter virus. Ramon sequenced and it came out funky.
                  & (virus == 'egfp'))
ggplot(relevant, aes(time, as.numeric(norfreq), colour=factor(donor), group = interaction(expt.name, sample.no))) +
  labs(x = 'Hours post infection', y ="# of GFP+ spots", color = "Donor:")+
  theme_bw() +
  #geom_label() +
  #facet_wrap(~expt.name) +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=8), 
        legend.key.size = unit(1, 'point'),
        legend.position = c(0.85, 0.6),
        legend.box.margin = margin(0, 0, 0, -10),
        legend.text = element_text(size=6), legend.title = element_text(size=8),
        plot.margin = margin(1, 0, 0, 0)) +
  guides(colour = guide_legend(override.aes = list(size=2))) +
  geom_line(alpha = 0.6, size = 0.3) +
  stat_smooth(method = "loess", se = TRUE, aes(group = 'none'))
  #geom_point(alpha = 0.6, size = 0.2)
#geom_point(alpha = 0.3) +
#geom_smooth(aes(color = NULL))
ggsave("plots/fig2c.pdf", width = 2.05, height = 1.8, units = "in")
relevant.small <- relevant[, c("csv", "donor", "expt.name", "sample.no", "time", 'norfreq')]

write.csv(relevant.small,"sourcedata/fig2c.csv", row.names = FALSE)

n <- relevant[, c('expt.name', 'sample.no')]
n <- relevant[, c('expt.name')]

uniquen <- unique(n)

# Figure 2d: peak # gfp+ spots vs. n rna copies per mm2
relevant = subset(daily.spots.summary.j, 
                  (condition == 'vanilla' | condition == 'static' | condition == 'shaken') 
                  #& (fromold == 1) # in infections in the bsl3 core i couldn't rinse so it affected the kinetics.
                  & !((expt.name == 'd5d6') & (sample.no == 2)) # this one i'm p sure was contaminated with non-reporter virus. Ramon sequenced and it came out funky.
                  & (virus == 'egfp'))

ggplot(relevant, (aes(y=as.numeric(n.mm.2mean), x= good.peak.spots,  color = factor(donor)))) + 
  labs(y = bquote('N copies per '~mm^2~''), x ="Peak # GFP+ Spots", color = "Donor:")+
  geom_point(show.legend = TRUE) +
  scale_y_log10()+
  theme_bw() + 
  geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, aes(group = FALSE, color = NULL), show.legend = FALSE) +
  #geom_smooth(method = 'lm', formula = y ~ x, se = FALSE) +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=8),
        plot.margin = margin(1, 0, 0, 0),
        legend.position = 'none')
        #legend.key.size = unit(0.1, 'point'),
        #legend.box.margin = margin(0, 0, 0, -10),
        #legend.text = element_text(size=6), legend.title = element_text(size=6))
ggsave("plots/fig2d.pdf", width = 2, height = 1.8, units = "in")
relevant.small <- relevant[, c("csv", "donor", "expt.name", "sample.no", "n.mm.2mean", 'good.peak.spots')]
write.csv(relevant.small,"sourcedata/fig2d.csv", row.names = FALSE)

# interesting to look at colored by rotary or by donor.
cor.test(log10(as.numeric(relevant$n.mm.2mean)), relevant$good.peak.spots, method = 'pearson')
ok <- lm(log10(as.numeric(n.mm.2mean)) ~ good.peak.spots, data = relevant)

### looking at initial spy # vs N titer
dataa = subset(relevant, (condition == 'vanilla' | condition == 'shaken' | condition == 'static') 
               & (dyes != "cmo") & (dyes != 'undyed') 
               & (usable == 1)
               & !((expt.name == 'd5d6') & (sample.no == 2)) # this one i'm p sure was contaminated with non-reporter virus. Ramon sequenced and it came out funky.
               & (virus == 'egfp')
               # these experiments had spytub added like right before imaging instead of >48 hours preinfection
               & (expt.name != 'd2shake') & (expt.name != 'd1d9') & expt.name != 'multitest')

ggplot(dataa, (aes(y=as.numeric(good.peak.spots), x= as.numeric(t0spy500/t0spy0),  color = factor(focus.simple)))) + 
  #labs(y = bquote('N copies per '~mm^2~''), x ="Peak # GFP+ Spots", color = "Donor:")+
  #bquote('N copies per '~mm^2' culture')
  geom_point(show.legend = TRUE) +
  #geom_point(show.legend = FALSE) +
  #scale_y_log10()+
  #xlim(0, 30000) +
  theme_bw() + 
  geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, aes(group = FALSE, color = NULL), show.legend = FALSE) +
  #geom_smooth(method = 'lm', formula = y ~ x, se = FALSE) +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=8),
        plot.margin = margin(1, 0, 0, 0),
        legend.key.size = unit(0.1, 'point'),
        legend.box.margin = margin(0, 0, 0, -10),
        legend.text = element_text(size=6), legend.title = element_text(size=6))
ggsave("C:/Users/Mark/Desktop/figure_images/2/fig2dsup_cryptinfxn.pdf", width = 2, height = 1.8, units = "in")




#Figure 2e: number of initial foci
level_order <- c('mock', 'infected')
f2edata <- subset(daily.spots.summary.j, 
                  (condition == "vanilla" | condition == 'static' | condition == 'shaken')
                  & !((expt.name == 'd5d6') & (sample.no == 2)) # this one i'm p sure was contaminated with non-reporter virus. Ramon sequenced and it came out funky.
                  & (virus == 'mock' | virus == 'egfp')
                  & usable == 1)
ggplot(f2edata, (aes(y=`20hpi.spots`, x=factor(virus, level = level_order), color = factor(donor)))) + 
  labs(x = 'Condition', y ="# GFP+ spots at 20 HPI", color = "Donor:")+
  geom_hline(yintercept = 0, color = 'black') +
  geom_quasirandom(alpha = 0.7, show.legend = FALSE) +
  stat_summary(fun = median,
               fun.min = function(z) { quantile(z,0.25) }, 
               fun.max = function(z) { quantile(z,0.75) }, 
               geom = "pointrange",
               show.legend = FALSE,
               aes(group = factor(virus, level = level_order))) +
  scale_x_discrete(labels = c("Mock", "Infected"))+
  theme_bw() + 
  coord_flip() +
  #scale_y_log10() +
  theme(axis.text = element_text(size=5), axis.title = element_text(size=8),
        plot.margin = margin(1, 0, 0, 0),
        legend.key.size = unit(0.1, 'point'),
        legend.box.margin = margin(0, 0, 0, -10),
        legend.text = element_text(size=6), legend.title = element_text(size=6))
ggsave("plots/fig2e.pdf", width = 2.4, height = 1.5, units = "in")
f2edata.small <- f2edata[, c("csv", "donor", "expt.name", "sample.no", "20hpi.spots", 'virus')]
write.csv(f2edata.small,"sourcedata/fig2e.csv", row.names = FALSE)


ok <- subset(f2edata, virus == 'egfp')
median(ok$`20hpi.spots`, na.rm = 'TRUE')

#Figure 2f: time to peak egfp+ area
level_order <- c('mock', 'egfp') 
f2fdata <- subset(daily.spots.summary.j, 
                  (condition == "vanilla" | condition == 'static' | condition == 'shaken') 
                  & (2*video.end.frame+1 > good.peak.spots.time) 
                  & (expt.name != 'd7d8')
                  & !((expt.name == 'd5d6') & (sample.no == 2)) # this one i'm p sure was contaminated with non-reporter virus. Ramon sequenced and it came out funky.
                  & (virus == 'egfp' | virus == 'mock'))
                  #& !(expt.name == 'agarose_d3' & sample.no == 3)) not sure why i left this one out of the first draft? maybe someday I'll realize
tf <- subset(spots, csv == 'longcilia2_pt10_aligned_cropped-all-spots.csv')
# it seems some of the spots files really, really did not import correctly. I'm just going to omit them.
f2fdata <- subset(f2fdata, good.peak.spots > 10)
#ggplot(f2fdata, (aes(y=good.peak.spots, x=good.peak.spots.time, color = factor(donor)))) + 
#  labs(y = 'Peak GFP+ spots', x ="Time of Peak GFP+ Spots (HPI)", color = "Donor:")+
  #geom_quasirandom(size = 0.3, alpha = 0.8) +
#  geom_point() +
#  scale_y_log10() +
#  theme_bw() +
#  stat_smooth(method = 'lm', aes(group = 'none'))+
#  theme(axis.text = element_text(size=5), axis.title = element_text(size=8),
#        plot.margin = margin(0, -5, 0, 0),
#        legend.key.width = unit(0.1, 'in'),
#        legend.key.height = unit(0.1, 'in'),
#        legend.box.margin = margin(0, 0, 0, -10),
#        legend.text = element_text(size=4), legend.title = element_text(size=5))
ggplot(f2fdata, (aes(y=good.peak.spots.time, x=factor(virus, level = level_order), color = good.peak.spots))) + 
  labs(x = 'Condition', y ="Time of Peak GFP+ Spots (HPI)", color = "# Spots")+
  geom_quasirandom(size = 0.3, alpha = 0.8) +
  stat_summary(fun.y = median,
               fun.ymin = function(z) { quantile(z,0.25) }, 
               fun.ymax = function(z) { quantile(z,0.75) }, 
               color = 'black',
               linewidth = 0.2,
               size = 0.3,
               geom = "crossbar") +
  scale_color_continuous(high = '#05e895', low = '#402375', limits = c(NA, NA)) +
  coord_flip() +
  ylim(0, 200) +
  scale_x_discrete(labels = c("Mock", "Infected"))+
  theme_bw() +
  theme(axis.text = element_text(size=5), axis.title = element_text(size=8),
        plot.margin = margin(0, 1, 0, 0),
        legend.key.width = unit(0.1, 'in'),
        legend.position = c(0.9, 0.35),
        legend.key.height = unit(0.1, 'in'),
        legend.box.margin = margin(0, 0, 0, -10),
        legend.text = element_text(size=4), legend.title = element_text(size=5))
ggsave("plots/fig2f.pdf", width = 2.4, height = 1.5, units = "in")
f2fdata.small <- f2fdata[, c("csv", "donor", "expt.name", "sample.no", "good.peak.spots", 'virus', 'good.peak.spots.time')]
write.csv(f2fdata.small,"sourcedata/fig2f.csv", row.names = FALSE)

cor.test(f2fdata$good.peak.spots.time, f2fdata$good.peak.spots, method = 'pearson')


relevant <- subset(good.spots.summary.j, condition == "vanilla" & (2*video.end.frame+1 > good.peak.spots.time) & !(expt.name == 'agarose_d3' & sample.no == 3) & virus == 'egfp')
f2finf <- subset(f2fdata, virus == 'egfp')
median(f2finf$good.peak.spots.time, na.rm = 'TRUE')
cor.test((f2finf$good.peak.spots.time, log10(f2finf$good.peak.spots), method = 'pearson')


#Figure 2h: correlation of focus type w/ mucus mvmt & peak egfp
level_order <- c('mock', 'egfp') 
level_order <- c('plaque', 'comet', 'diffuse', 'crypt', 'mock') 

relevant <- subset(daily.spots.summary.j, 
                   (condition == 'vanilla' | condition == 'shaken' | condition == 'static') 
                   & (virus == 'egfp') 
                   & (good.peak.spots > 10)
                   & (expt.name != 'd5d6' & sample.no != 2))
level_order <- c('plaque', 'comet', 'diffuse') 
ggplot(relevant, (aes(y=as.numeric(good.peak.spots), x=factor(focus.simple, levels = level_order), color = (mucus.simple == 'rotary')))) + 
  labs(x = 'Focus Type', y ="Peak GFP+ Spots", color = "MCC Type:")+
  scale_color_manual(labels = c("Disorganized", "Rotary"), values = c('#4A4CA0', 'deeppink2'))+
  theme_bw() + 
  geom_quasirandom(alpha = 0.7, width = 0.2) +
  ylim(0, 40000) +
  scale_x_discrete(limits = level_order, labels = c('Plaque', 'Comet', 'Diffuse'))+
  stat_summary(fun = median,
               fun.min = function(x) mean(x) - sd(x), 
               fun.max = function(x) mean(x) + sd(x), 
               geom = "crossbar",
               width = 0.3,
               aes(group = FALSE, color = NULL)) +
  guides(shape = guide_legend(override.aes = list(size = 0.01))) +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=8),
        plot.margin = margin(1, -5, 0, 0),
        legend.key.width = unit(0.1, 'in'),
        legend.key.height = unit(0.1, 'in'),
        legend.box.margin = margin(0, 0, 0, -10),
        legend.text = element_text(size=5), legend.title = element_text(size=6))
ggsave("plots/fig2h.pdf", width = 2, height = 1.35, units = "in")
f2hdata.small <- relevant[, c("csv", "donor", "expt.name", "sample.no", "good.peak.spots", 'virus', 'focus.simple', 'mucus.simple')]
write.csv(f2hdata.small,"sourcedata/fig2h.csv", row.names = FALSE)

wilcox.test(subset(relevant, focus.simple == 'comet')$good.peak.spots, subset(relevant, focus.simple == 'plaque')$good.peak.spots)
wilcox.test(as.numeric(subset(relevant, focus.simple == 'comet')$n.mm.2mean), as.numeric(subset(relevant, focus.simple == 'plaque')$n.mm.2mean))

t.test(as.numeric(subset(relevant, focus.simple == 'diffuse')$good.peak.spots), as.numeric(subset(relevant, focus.simple == 'comet')$good.peak.spots))
t.test(as.numeric(subset(relevant, focus.simple == 'comet')$n.mm.2mean), as.numeric(subset(relevant, focus.simple == 'diffuse')$n.mm.2mean))

g1 <- relevant$good.peak.spots
qqnorm(g1)
qqline(g1, col = "red")

###################################
## S04: CMO vs SPYTUB ~kinetics  ##
###################################
comp <- subset(daily.spots.summary.j,
               (condition == 'vanilla' | condition == 'static' | condition == 'shaken') &
                 (donor == 2 | donor == 3 | donor == 4) &
                 (dyes == 'cmo' | dyes == 'spy + nuc') &
                 virus == 'egfp')

compp <- subset(comp, 2*video.end.frame-1 != good.peak.spots.time)
ggplot(compp, (aes(y=as.numeric(good.peak.spots.time), x=factor(dyes), color = factor(donor)))) + 
  labs(x = 'Dye', y ="Time of Peak GFP+ Spots (HPI)", color = "Donor:")+
  theme_bw() + 
  ylim(0, 200) +
  scale_x_discrete(labels = c('CMO', 'NV + SPYtub'))+
  stat_summary(fun = mean,
               fun.min = function(x) mean(x) - sd(x), 
               fun.max = function(x) mean(x) + sd(x), 
               geom = "crossbar",
               width = 0.7,
               aes(group = FALSE, color = NULL)) +
  geom_jitter(width = 0.3)+
  #guides(shape = guide_legend(override.aes = list(size = 0.01))) +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=8),
        plot.margin = margin(1, -5, 0, 0),
        legend.key.width = unit(0.1, 'in'),
        legend.key.height = unit(0.1, 'in'),
        legend.box.margin = margin(0, 0, 0, -10),
        legend.text = element_text(size=5), legend.title = element_text(size=6))
ggsave("plots/s2a.pdf", width = 2, height = 2, units = "in")

g1 <- subset(compp, dyes == 'cmo')$good.peak.spots.time
g2 <- subset(compp, dyes != 'cmo')$good.peak.spots.time
wilcox.test(g1, g2)
t.test(g1, g2)

ggplot(comp, (aes(y=as.numeric(`20hpi.spots`), x=factor(dyes), color = factor(donor)))) + 
  labs(x = 'Dye', y ="GFP+ Spots at 20 HPI", color = "Donor:")+
  theme_bw() + 
  ylim(0, 350) +
  scale_x_discrete(labels = c('CMO', 'NV + SPYtub'))+
  stat_summary(fun = median,
               fun.min = function(x) quantile(x, 0.75), 
               fun.max = function(x) quantile(x, 0.25), 
               geom = "crossbar",
               width = 0.7,
               aes(group = FALSE, color = NULL)) +
  geom_jitter(width = 0.3)+
  theme(axis.text = element_text(size=6), axis.title = element_text(size=8),
        plot.margin = margin(1, -5, 0, 0),
        legend.key.width = unit(0.1, 'in'),
        legend.key.height = unit(0.1, 'in'),
        legend.box.margin = margin(0, 0, 0, -10),
        legend.text = element_text(size=5), legend.title = element_text(size=6))
ggsave("plots/s2b.pdf", width = 2, height = 2, units = "in")

g1 <- subset(comp, dyes == 'cmo')$`20hpi.spots`
g2 <- subset(comp, dyes != 'cmo')$`20hpi.spots`
wilcox.test(g1, g2)
t.test(g1, g2)


ggplot(compp, (aes(y=as.numeric(good.peak.spots), x=factor(dyes), color = factor(donor)))) + 
  labs(x = 'Dye', y ="Peak # of GFP+ Spots", color = "Donor:")+
  theme_bw() + 
  ylim(0, 38000) +
  scale_x_discrete(labels = c('CMO', 'NV + SPYtub'))+
  stat_summary(fun = mean,
               fun.min = function(x) mean(x) - sd(x), 
               fun.max = function(x) mean(x) + sd(x),
               geom = "crossbar",
               width = 0.7,
               aes(group = FALSE, color = NULL)) +
  geom_jitter(width = 0.3)+
  theme(axis.text = element_text(size=6), axis.title = element_text(size=8),
        plot.margin = margin(1, -5, 0, 0),
        legend.key.width = unit(0.1, 'in'),
        legend.key.height = unit(0.1, 'in'),
        legend.box.margin = margin(0, 0, 0, -10),
        legend.text = element_text(size=5), legend.title = element_text(size=6))
ggsave("plots/s2c.pdf", width = 2, height = 2, units = "in")

g1 <- subset(compp, dyes == 'cmo')$good.peak.spots
g2 <- subset(compp, dyes != 'cmo')$good.peak.spots
wilcox.test(g1, g2)
t.test(g1, g2)


comp$mucus.stop.time <- 2*comp$mucus.stop.frame-1
comppp <- subset(comp, mucus.stop.frame.end == 'no')
ggplot(comppp, (aes(y=as.numeric(mucus.stop.time), x=factor(dyes), color = factor(donor)))) + 
  labs(x = 'Dye', y ="Peak # of GFP+ Spots", color = "Donor:")+
  theme_bw() + 
  ylim(0, 140) +
  scale_x_discrete(labels = c('CMO', 'NV + SPYtub'))+
  stat_summary(fun = mean,
               fun.min = function(x) mean(x) - sd(x), 
               fun.max = function(x) mean(x) + sd(x),
               geom = "crossbar",
               width = 0.7,
               aes(group = FALSE, color = NULL)) +
  geom_jitter(width = 0.3)+
  theme(axis.text = element_text(size=6), axis.title = element_text(size=8),
        plot.margin = margin(1, -5, 0, 0),
        legend.key.width = unit(0.1, 'in'),
        legend.key.height = unit(0.1, 'in'),
        legend.box.margin = margin(0, 0, 0, -10),
        legend.text = element_text(size=5), legend.title = element_text(size=6))
ggsave("plots/s2d.pdf", width = 2, height = 2, units = "in")

g1 <- subset(comppp, dyes == 'cmo')$mucus.stop.time
g2 <- subset(comppp, dyes != 'cmo')$mucus.stop.time
wilcox.test(g1, g2)
t.test(g1, g2)

comp$video.end.time <- 2*comp$video.end.frame-1
s4data.small <- comp[, c("csv", "donor", "expt.name", "sample.no", "good.peak.spots", 'good.peak.spots.time', '20hpi.spots', 'mucus.stop.time', 'video.end.time')]
write.csv(s4data.small,"sourcedata/s4.csv", row.names = FALSE)

##########################################
## S05: N vs. GFP colored by focus type ##
##########################################

# Figure 2d: peak # gfp+ spots vs. n rna copies per mm2
relevant = subset(daily.spots.summary.j, 
                  (condition == 'vanilla' | condition == 'static' | condition == 'shaken') &
                  #& (fromold == 1) # in infections in the bsl3 core i couldn't rinse so it affected the kinetics.
                  & !((expt.name == 'd5d6') & (sample.no == 2) # this one i'm p sure was contaminated with non-reporter virus. Ramon sequenced and it came out funky.
                  & (usable == 1)
                  & (virus == 'egfp'))

ggplot(relevant, (aes(y=as.numeric(n.mm.2mean), x= good.peak.spots,  color = factor(focus.simple)))) + 
  labs(y = bquote('N copies per '~mm^2~''), x ="Peak # GFP+ Spots", color = "Donor:")+
  geom_point(show.legend = TRUE) +
  scale_y_log10()+
  theme_bw() + 
  geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, show.legend = FALSE) +
  #geom_smooth(method = 'lm', formula = y ~ x, se = FALSE) +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=8),
        plot.margin = margin(1, 0, 0, 0),
        legend.position = 'none')
#legend.key.size = unit(0.1, 'point'),
#legend.box.margin = margin(0, 0, 0, -10),
#legend.text = element_text(size=6), legend.title = element_text(size=6))
ggsave("plots/fig2d.pdf", width = 2, height = 1.8, units = "in")
relevant.small <- relevant[, c("csv", "donor", "expt.name", "sample.no", "n.mm.2mean", 'good.peak.spots')]
write.csv(relevant.small,"sourcedata/fig2d.csv", row.names = FALSE)




#rell <- subset(rel, csv != 'mb_longcilia_pt04_aligned_cropped-all-spots.csv') # this one was super out of focus. looking and seeing how it affects the line.
rel <- subset(daily.spots.summary.j, usable == 1 & virus == 'egfp' & !is.na(focus.simple) & focus.simple != 'mock')
# im not doing a notably worse job counting spots for the plaques vs. diffuse
ggplot(rel, (aes(y=as.numeric(good.peak.spots), x= as.numeric(maxgfp1100/maxgfp0),  color = factor(focus.simple)))) + 
  labs(y = 'Peak GFP+ Spots', x = "Fractional GFP+ Culture Area", color = "Focus type:")+
  #bquote('N copies per '~mm^2' culture')
  geom_point(show.legend = TRUE, alpha = 0.6) +
  #geom_point(show.legend = FALSE) +
  #scale_y_log10()+
  #xlim(0, 30000) +
  theme_bw() + 
  geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, aes(group = focus.simple), show.legend = FALSE) +
  #geom_smooth(method = 'lm', formula = y ~ x, se = FALSE) +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=8),
        plot.margin = margin(1, 0, 0, 0),
        legend.key.size = unit(0.1, 'point'),
        legend.box.margin = margin(0, 0, 0, -10),
        legend.text = element_text(size=6), legend.title = element_text(size=6))
ggsave("plots/S5.pdf", width = 2.8, height = 2.5, units = "in")
rel.small <- rel[, c("csv", "donor", "expt.name", "sample.no", "n.mm.2mean", 'good.peak.spots', "focus.simple")]
write.csv(rel.small,"sourcedata/s5.csv", row.names = FALSE)



######
## Supplemental Figure 6a




###################
## S06, kinetics ##
###################

# same subset as fig. 2c
relevant = subset(spots, 
                  (condition == 'vanilla' | condition == 'static' | condition == 'shaken') 
                  & !((expt.name == 'd7d8') & (time > 120)) # this one got rinsed then so not a vanilla infection
                  & !((expt.name == 'ruxcilia') & time > 72) # big gap after 72 on this one and it flooded at that point.
                  #& (fromold == 1) # in infections in the bsl3 core i couldn't rinse so it affected the kinetics.
                  & (!((expt.name == 'd5d6') & (sample.no == 2))) # this one i'm p sure was contaminated with non-reporter virus. Ramon sequenced and it came out funky.
                  & (virus == 'egfp'))

kin <- subset(relevant,
              expt.name != 'ruxcilia'
              & expt.name != 'd2d4ag') #these two experiments had major focus changes that erroneously caused massive spikes/dips
length(unique(kin$csv))
ggplot(kin, aes(time, as.numeric(rate.of.change), colour=factor(donor), group = interaction(expt.name, sample.no))) +
  labs(x = 'Hours post infection', y = expression(Delta ~ "GFP+ spots per hour"), color = "Donor:")+
  theme_bw() +
  #geom_label() +
  #facet_wrap(~expt.name) +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=8), 
        legend.key.size = unit(1, 'point'),
        legend.box.margin = margin(0, 0, 0, -10),
        legend.text = element_text(size=6), legend.title = element_text(size=8),
        plot.margin = margin(1, 0, 0, 0)) +
  guides(colour = guide_legend(override.aes = list(size=2))) +
  geom_line(alpha = 0.6, size = 0.3) +
  stat_smooth(method = "loess", se = TRUE, span = 0.2, aes(group = 'none'))
#geom_point(alpha = 0.6, size = 0.2)
#geom_point(alpha = 0.3) +
#geom_smooth(aes(color = NULL))
ggsave("plots/s06.pdf", width = 5, height = 4, units = "in")
kin.small <- kin[, c("time", "rate.of.change", "donor", "expt.name", "sample.no")]
write.csv(kin.small,"sourcedata/s06.csv", row.names = FALSE)


#######################################################
## S07. Kinetics by focus type & culture abnormality ##
#######################################################
kin2 <- subset(kin, focus.simple == 'comet' | focus.simple == 'plaque' | focus.simple == 'diffuse' | focus.simple == 'crypt')
ggplot(kin2, aes(time, as.numeric(rate.of.change), colour=factor(focus.simple), group = interaction(expt.name, sample.no))) +
  labs(x = 'Hours post infection', y = expression(Delta ~ "GFP+ spots per hour"), color = "Focus type:")+
  theme_bw() +
  geom_hline(yintercept = 0) +
  #geom_label() +
  #facet_wrap(~expt.name) +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=8), 
        legend.key.size = unit(1, 'point'),
        legend.box.margin = margin(0, 0, 0, -10),
        legend.text = element_text(size=6), legend.title = element_text(size=8),
        plot.margin = margin(1, 0, 0, 0)) +
  guides(colour = guide_legend(override.aes = list(size=2))) +
  geom_line(alpha = 0.4, size = 0.3) +
  stat_smooth(method = "loess", se = TRUE, span = 0.1, aes(group = focus.simple))
#geom_point(alpha = 0.6, size = 0.2)
#geom_point(alpha = 0.3) +
#geom_smooth(aes(color = NULL))
ggsave("plots/s07kinetics_byfocustype.pdf", width = 3.3, height = 3, units = "in")

fkin2.small <- kin2[, c("time", "rate.of.change", "donor", "expt.name", "sample.no", "focus.simple", "crypts", "jammed", "cyst", "mucus.simple")]
write.csv(fkin2.small,"sourcedata/s07ab.csv", row.names = FALSE)

kin2 <- kin2 %>%
  mutate(crypts = case_when(
    crypts == 'pointycrypts' ~ 'clean',
    TRUE ~ crypts
  ))
kin2 <- kin2 %>%
  mutate(mucus.simple = case_when(
    mucus.simple == 'static' ~ 'disorganized',
    TRUE ~ mucus.simple
  ))


library(tidyr)
# Combine the different conditions into one long format
kin2_long <- kin2 %>%
  pivot_longer(cols = c(crypts, jammed, cyst, mucus.simple), 
               names_to = "condition_type", 
               values_to = "defect")


library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
# Create individual plots for each condition

# Define color palettes for each condition
color_palette_crypts <- c("clean" = "purple", "edgefurrow" = "orange", "furrows" = "chartreuse4")  # Replace with your actual colors
crypts_labels <- c("clean" = "Clean", "edgefurrow" = "Edge Furrows", "furrows" = "Furrows")  # Replace with your desired labels

color_palette_jammed <- c("jammed" = "darkorchid4", "edge" = "coral", "unjammed" = "cyan4")  # Replace with your actual colors
jammed_labels <- c("jammed" = "Jammed", "edge" = "Edge Unjamming", "unjammed" = "Widespread Unjamming")  # Replace with your desired labels

color_palette_cyst <- c("clean" = "brown1", "cysts" = "deepskyblue", "densecysts" = "forestgreen")  # Replace with your actual colors
cyst_labels <- c("clean" = "Clean", "cysts" = "Few cysts", "densecysts" = "Dense cysts")  # Replace with your desired labels

color_palette_mucus <- c("disorganized" = "#4A4CA0", "rotary" = "deeppink2")  # Replace with your actual colors
mucus_labels <- c("disorganized" = "Disorganized", "rotary" = "Rotary")  # Replace with your desired labels

plot_crypts <- ggplot(filter(kin2_long, condition_type == "crypts"), aes(time, as.numeric(rate.of.change), colour=factor(defect), group = interaction(expt.name, sample.no))) +
  labs(x = 'Hours post infection', y = expression(Delta ~ "GFP+ spots per hour"), color = "Crypts:") +
  theme_bw() +
  xlim(10, 200) +
  geom_hline(yintercept = 0) +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=6), 
        legend.key.size = unit(1, 'point'),
        legend.box.margin = margin(0, 0, 0, -10),
        legend.text = element_text(size=5), legend.title = element_text(size=6),
        plot.margin = margin(1, 0, 0, 0)) +
  guides(colour = guide_legend(override.aes = list(size=2))) +
  geom_line(alpha = 0.4, size = 0.3) +
  stat_smooth(method = "loess", se = TRUE, span = 0.1, aes(group = defect)) +
  scale_color_manual(values = color_palette_crypts, labels = crypts_labels)


plot_jammed <- ggplot(filter(kin2_long, condition_type == "jammed"), aes(time, as.numeric(rate.of.change), colour=factor(defect), group = interaction(expt.name, sample.no))) +
  labs(x = 'Hours post infection', y = expression(Delta ~ "GFP+ spots per hour"), color = "Jammedness:") +
  theme_bw() +
  xlim(10, 200) +
  geom_hline(yintercept = 0) +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=6), 
        legend.key.size = unit(1, 'point'),
        legend.box.margin = margin(0, 0, 0, -10),
        legend.text = element_text(size=5), legend.title = element_text(size=6),
        plot.margin = margin(1, 0, 0, 0)) +
  guides(colour = guide_legend(override.aes = list(size=2))) +
  geom_line(alpha = 0.4, size = 0.3) +
  stat_smooth(method = "loess", se = TRUE, span = 0.1, aes(group = defect)) +
  scale_color_manual(values = color_palette_jammed, labels = jammed_labels)


plot_cyst <- ggplot(filter(kin2_long, condition_type == "cyst"), aes(time, as.numeric(rate.of.change), colour=factor(defect), group = interaction(expt.name, sample.no))) +
  labs(x = 'Hours post infection', y = expression(Delta ~ "GFP+ spots per hour"), color = "Cysts:") +
  theme_bw() +
  xlim(10, 200) +
  geom_hline(yintercept = 0) +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=6), 
        legend.key.size = unit(1, 'point'),
        legend.box.margin = margin(0, 0, 0, -10),
        legend.text = element_text(size=5), legend.title = element_text(size=6),
        plot.margin = margin(1, 0, 0, 0)) +
  guides(colour = guide_legend(override.aes = list(size=2))) +
  geom_line(alpha = 0.4, size = 0.3) +
  stat_smooth(method = "loess", se = TRUE, span = 0.1, aes(group = defect)) +
  scale_color_manual(values = color_palette_cyst, labels = cyst_labels)

plot_mucus <- ggplot(filter(kin2_long, condition_type == "mucus.simple"), aes(time, as.numeric(rate.of.change), colour=factor(defect), group = interaction(expt.name, sample.no))) +
  labs(x = 'Hours post infection', y = expression(Delta ~ "GFP+ spots per hour"), color = "MCT:") +
  theme_bw() +
  xlim(10, 200) +
  geom_hline(yintercept = 0) +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=6), 
        legend.key.size = unit(1, 'point'),
        legend.box.margin = margin(0, 0, 0, -10),
        legend.text = element_text(size=5), legend.title = element_text(size=6),
        plot.margin = margin(1, 0, 0, 0)) +
  guides(colour = guide_legend(override.aes = list(size=2))) +
  geom_line(alpha = 0.4, size = 0.3) +
  stat_smooth(method = "loess", se = TRUE, span = 0.1, aes(group = defect)) +
  scale_color_manual(values = color_palette_mucus, labels = mucus_labels)

# Combine the plots with patchwork
combined_plot <- plot_crypts + plot_jammed + plot_cyst + plot_mucus + plot_layout(ncol = 1, nrow = 4)
print(combined_plot)
# Save the combined plot
ggsave("plots/sfig7kinetics_combined_grid.pdf", combined_plot, width = 3.5, height = 5.5, units = "in")



## Rest of sfig 7c.

###########################
## Jammedness regression ##
###########################

relevant <- subset(daily.spots.summary.j, 
                   (condition == 'vanilla' | condition == 'shaken' | condition == 'static') 
                   & (virus == 'egfp') 
                   & (good.peak.spots > 10)
                   & (usable == 1)
                   & (expt.name != 'd5d6' & sample.no != 2))

relevant <- relevant %>%
  mutate(donor = as.factor(donor),
         crypts = as.factor(crypts),
         cyst = as.factor(cyst),
         n.mm.2mean = as.numeric(n.mm.2mean),
         good.peak.spots = as.numeric(good.peak.spots),
         good.peak.spots.time = as.numeric(good.peak.spots.time),
         jammed = as.factor(jammed),
         infxn.in.crypt = as.factor(infxn.in.crypt),
         dyes = as.factor(dyes),
         mucus.simple = as.factor(mucus.simple),
         focus.simple = as.factor(focus.simple),
         s20hpi.spots = as.numeric(`20hpi.spots`),
         piv.median.speed = as.numeric(piv.median.speed),
         piv.mean.speed = as.numeric(piv.mean.speed),
         ciliation = case_when(
           (!(expt.name %in% c("d2shake", "d1d9", "multitest")) & !(dyes %in% c("cmo", "undyed"))) ~ t0spy350 / t0spy0,
           TRUE ~ NA_real_)
  )
nmodelcrypts <- lm(log10(n.mm.2mean) ~ crypts +cyst + jammed + mucus.simple + as.factor(donor), data = relevant, subset = !(as.factor(donor) %in% c(5, 6)))
# adding in crypts doesn't help; crypts instead of 
nmodel <- lm(log10(n.mm.2mean) ~ infxn.in.crypt + cyst + jammed + mucus.simple + as.factor(donor), data = relevant, subset = !(as.factor(donor) %in% c(5, 6)))
# ^ this is the good one !
nmodelpiv <- lm(log10(n.mm.2mean) ~ infxn.in.crypt + cyst + jammed + mucus.simple + as.factor(donor), data = relevant, subset = !(as.factor(donor) %in% c(5, 6)))
# ^ subbing piv speed (mean or median) for qualitative jamming judgement is not helpful.

spotsmodel_sameasn <- lm(log10(good.peak.spots) ~ infxn.in.crypt + cyst + jammed + mucus.simple + as.factor(donor), data = relevant, subset = !(as.factor(donor) %in% c(5, 6)))
# adding focus.simple to n seems to be the best way to get it to a higher r^2

simplenmodel <- lm(log10(n.mm.2mean) ~ crypts + cyst + jammed + mucus.simple, data = relevant)
#eliminated donor due to identical aic
summary(spotsmodel)
AIC(nmodel,nmodelpiv)

## making a mosaic plot to demonstrate.
relevant <- relevant %>%
  mutate(jammed = case_when(
    jammed == 'jammed' ~ 'Jammed',
    jammed == 'edge' ~ 'Edge Unjamming',
    jammed == 'unjammed' ~ 'Widespread Unjamming',
    jammed == 'hypermobile' ~ 'Hypermobile',
    TRUE ~ jammed
  ))
relevant$jammed <- factor(relevant$jammed, levels = c("Jammed", "Edge Unjamming", "Widespread Unjamming", "Hypermobile"))

mosaic(~ focus.simple +jammed + mucus.simple, data = relevant, shade = TRUE, legend = TRUE, rot_labels = c(50, 20, 20, 90))


par(mfrow=c(2,2))
plot(nmodel)
shapiro.test(residuals(nmodel))
vif(nmodel)
alias(nmodel)

simplenmodel.adjusted <- lm(log10(n.mm.2mean) ~ crypts + cyst + jammed + mucus.simple + as.factor(donor), 
                            data = relevant, 
                            subset = !(as.factor(donor) %in% c(5, 6)))
par(mfrow=c(2,2))
plot(simplenmodel.adjusted)
shapiro.test(residuals(simplenmodel.adjusted))
vif(simplenmodel.adjusted)
summary(simplenmodel.adjusted)
anova(nmodel, spotsmodel_sameasn)
AIC(bigmodel, mucusonly)

# Print the summary of the model
summary(jammingonlyspots)
par(mfrow = c(2, 2))
plot(minmodel)


library(car)
vif(nmodel)
alias(nmodel)
relevant_adjusted <- relevant %>%
  filter(donor != 8, infxn.in.crypt != "yes")
nmodel <- lm(as.numeric(n.mm.2mean) ~ crypts + cyst + jammed + mucus.simple + focus.simple + as.factor(donor) + infxn.in.crypt +ciliation, data = relevant_adjusted)

spotsmodel <- lm(good.peak.spots ~ crypts + cyst + jammed + mucus.simple + focus.simple + as.factor(donor) + infxn.in.crypt + ciliation, data = relevant)
summary(spotsmodel)

AIC(nmodel, spotsmodel)


