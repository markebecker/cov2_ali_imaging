# Basic live imaging of ALIs. (Supp. Figs 2 & 3)



## Supp. Fig. 2
## What is a good power threshold?
## Does the denoising do bad things? How much better is it really?

# libraries
library(ggplot2)
library(ggbeeswarm)
library(data.table)
library(dplyr)
library(readr)
library(tidyverse)
# load csvs with processed data to look at
ciliako44.raw <- read.csv("data/ko_cbfs/ciliako44_raw_unshuf_concatenated.csv")
ciliako44.raw.frameshuf <- read.csv('data/ko_cbfs/ciliako44_raw_shuffled_concatenated.csv')
ciliako44.denoised <- read.csv("data/ko_cbfs/ciliako44_denoised_unshuf_concatenated.csv")
ciliako44.denoised.frameshuf <- read.csv("data/ko_cbfs/ciliako44_denoised_shuffled_concatenated.csv")
# ignoring shuffled pixels!  I did look at these however there is spatial structure (bright and dim patches)
# and you get powerful patches where there is bright spytub even in the frameshufs
# so randomizing pixel positions within each frame changes the power density distribution.
# using frameshuffs only to retain spatial structure while randomizing temporal structure


#they are too big. remove extra columns
ck.raw = subset(ciliako44.raw, select = -c(refimage, spytub, nucview, gfp, pol))
ck.raw.frameshuf = subset(ciliako44.raw.frameshuf, select = -c(spytub, nucview, gfp, pol))
ck.denoised = subset(ciliako44.denoised, select = -c(spytub, nucview, gfp, pol))
ck.denoised.frameshuf = subset(ciliako44.denoised.frameshuf, select = -c(spytub, nucview, gfp, pol))

rm(ciliako44.raw, ciliako44.denoised, ciliako44.denoised.frameshuf, ciliako44.raw.frameshuf)

# still unwieldy. smallen them to like 5%
ck.raw.smal <- ck.raw %>% sample_frac(0.05)
ck.raw.frameshuf.smal <- ck.raw.frameshuf %>% sample_frac(0.05)
ck.denoised.smal <- ck.denoised %>% sample_frac(0.05)
ck.denoised.frameshuf.smal <- ck.denoised.frameshuf %>% sample_frac(0.05)

namecol <- function(df, name) {
  df$process <- name
  return(df)
}

data_list <- c('ck.raw.smal', 'ck.raw.frameshuf.smal','ck.denoised.smal', 'ck.denoised.frameshuf.smal')
combined <- do.call(rbind, Map(namecol, mget(data_list), data_list))
fwrite(combined, 'D:/shufflin/shuffleset.csv')
#so i don't have to load all this data again.
combined <- read_csv("data/shuffleset.csv")

combined.smoller <- combined %>% slice_sample(prop = 0.02)
fwrite(combined.smoller, 'sourcedata/s022.csv')

combined <- subset(combined, condition == 'nt' | condition == 'cypa')
# What power level are 85% of pixels below?
quantiles <- combined %>%
  group_by(process) %>%
  summarize(quantile_value = quantile(power, probs = 0.85))
# Overwrite quantile values for unshuffled processes with shuffled ones. for graphing purposes.
quantiles$quantile_value[quantiles$process %in% c('ck.raw.smal', 'ck.denoised.smal')] <- 
  quantiles$quantile_value[quantiles$process %in% c('ck.raw.frameshuf.smal', 'ck.denoised.frameshuf.smal')]

# What fraction of powerful pixels fall above or below 15 Hz at that power thres?
fxns <- combined %>%
  group_by(process) %>%
  summarize(npow = sum(power > quantiles$quantile_value[quantiles$process %in% c(process)]),
            total = n(), 
            fxnpow = npow/total, 
            powerbeatin = sum((domfreq < 15) & (power > quantiles$quantile_value[quantiles$process %in% c(process)])), 
            fxnpowerbeat = powerbeatin/total, 
            naivebeatin = sum(domfreq < 15), 
            fxnnaivebeat = naivebeatin/total,
            powerlessbeatin = fxnnaivebeat - fxnpowerbeat)
labels <- c(
  "ck.raw.smal" = "Raw",
  "ck.denoised.smal" = "Denoised",
  "ck.raw.frameshuf.smal" = "Raw, Shuffled Frames",
  "ck.denoised.frameshuf.smal" = "Denoised, Shuffled Frames"
)

ggplot(combined, aes(x = power, y = domfreq)) +
  labs(x = 'Power Density (pixel units^2/Hz)', y = 'Frequency (Hz)') +
  #geom_point(alpha = 0.001) +  
  #geom_density_2d(color = 'blue') +
  geom_density_2d(bins = 10) +
  #geom_jitter(alpha = 0.0006) +
  facet_wrap(~process, labeller = labeller(process = labels)) +
  geom_hline(yintercept = 15, color = 'deeppink') +
  geom_vline(data = quantiles, aes(xintercept = quantile_value), color = "deeppink") +  # Add vlines
  #annotate('text', label = fxns$fxnpowerbeat) +
  scale_x_log10() +
  theme_bw() +
  theme(axis.text = element_text(size=8, face = 'italic'), axis.title = element_text(size=8),
        plot.margin = margin(1, 0, 0, 0),
        legend.position = 'none')
ggsave("plots/s3b.pdf", width = 4, height = 2.7, units = "in")

denoised <- subset(combined, process == 'ck.denoised.smal' | process == 'ck.denoised.frameshuf.smal')
raw <- subset(combined, process == 'ck.raw.smal' | process == 'ck.raw.frameshuf.smal')

ggplot(denoised, aes(x = power, y = domfreq)) +
  labs(x = 'Power Density(pixel units^2/Hz)', y = 'Frequency (Hz)') +
  #geom_point(alpha = 0.001) +  
  #geom_density_2d(color = 'blue') +
  geom_density_2d(bins = 10) +
  #geom_jitter(alpha = 0.0006) +
  facet_wrap(~process, labeller = labeller(process = labels)) +
  geom_hline(yintercept = 15, color = 'deeppink') +
  geom_vline(data = quantiles, aes(xintercept = quantile_value), color = "deeppink") +  # Add vlines
  #annotate('text', label = fxns$fxnpowerbeat) +
  scale_x_log10() +
  theme_bw() +
  theme(axis.text = element_text(size=8, face = 'italic'), axis.title = element_text(size=8),
        plot.margin = margin(1, 0, 0, 0),
        legend.position = 'none')
ggsave("C:/Users/Mark/Desktop/figure_images/s3/s3b.pdf", width = 4, height = 2.7, units = "in")


# Overwrite quantile values for unshuffled processes with shuffled ones
quantiles$quantile_value[quantiles$process %in% c('ck.raw.smal', 'ck.denoised.smal')] <- 
  quantiles$quantile_value[quantiles$process %in% c('ck.raw.frameshuf.smal', 'ck.denoised.frameshuf.smal')]


#######################################

levels = seq(from = 0, to = 1, by = 0.01)

leveltestt <- data.frame()
for (i in seq_along(levels)) {
  n <- levels[i]
  # What power level are 85% of pixels below?
  quantiles <- combined %>%
    group_by(process) %>%
    summarize(quantile_value = quantile(power, probs = n))

# What fraction of powerful pixels fall above or below 18 Hz at that power thres?
  fxns <- combined %>%
    group_by(process) %>%
    summarize(npow = sum(power > quantiles$quantile_value[quantiles$process %in% c(process)]),
              powerbeatin = sum(domfreq < 18 & power > quantiles$quantile_value[quantiles$process %in% c(process)]), 
              fxnpowerbeat = powerbeatin/npow)
  leveltestt <- rbind(leveltestt, c(n, t(fxns$fxnpowerbeat)))
}
colnames(leveltestt) <- c("level", "denoised.frameshuf", "denoised", "raw.frameshuf", "raw")

ggplot(leveltestt, aes(x = level)) +
  geom_line(aes(y = denoised.frameshuf, color = "denoised.frameshuf")) +
  geom_line(aes(y = denoised, color = "denoised")) +
  geom_line(aes(y = raw.frameshuf, color = "raw.frameshuf")) +
  geom_line(aes(y = raw, color = "raw")) +
  geom_vline(xintercept = 0.85) +
  scale_color_manual(values = c("denoised.frameshuf" = "palegreen4", "denoised" = 'lightskyblue4', "raw.frameshuf" = "tomato3", "raw" = "orange"))



####################
##  Supp. Fig. 3  ##
####################

exptlog_iii <- read_csv("data/exptlog_iii.csv")
spots <- read_csv("D:/000_rebuttal/240702_goodspots.csv")

# take the spots over time csv and extract specific values of interest
daily.spots.summary <- spots %>%
  group_by(csv) %>%
  summarise(good.peak.spots = max(norfreq), 
            good.peak.spots.time = time[which.max(norfreq)], 
            aucspots = sum(norfreq)) 
# cross back with the info i want
daily.spots.summary.j <- inner_join(daily.spots.summary, exptlog_iii)


daily.spots.summary.j <- daily.spots.summary.j %>%
  mutate(jammed = case_when(
    jammed == 'jammed' ~ 'Jammed',
    jammed == 'edge' ~ 'Edge Migration',
    jammed == 'unjammed' ~ 'Widespread Migration',
    jammed == 'hypermobile' ~ 'Hypermobile',
    TRUE ~ jammed
  ))

# this is the correct dataset- including the KOs, I have inconsistent level names for crypts & cysts.
dataa = subset(daily.spots.summary.j, ((condition == 'vanilla') | (condition == 'shaken') | (condition == 'static') | (condition == 'rinsed') & (usable == 1)))
dataa$jammed <- factor(dataa$jammed, levels = c("Jammed", "Edge Migration", "Widespread Migration", "Hypermobile"))
dataa <- dataa %>%
  mutate(crypts = case_when(
    crypts == 'clean' ~ 'No Furrows',
    crypts == 'edgefurrow' ~ 'Edge Furrows',
    crypts == 'furrows' ~ 'Center Furrows',
    crypts == 'pointycrypts' ~ 'No Furrows'
  )) %>%
  mutate(mucus.simple = case_when(
    mucus.simple == 'rotary' ~ 'Rotary',
    mucus.simple == 'disorganized' ~ 'Disorganized'))
dataa$crypts <- factor(dataa$crypts, levels = c("No Furrows", "Edge Furrows", "Center Furrows"))

ggplot(dataa, aes(x = jammed, y = as.numeric(piv.median.speed), color = as.factor(crypts))) +
  geom_jitter() + 
  geom_hline(yintercept = 3, linetype = 'dashed', color = 'black') +
  labs(x = "Jamming Category (over entire imaging period)",
       y = "Median SPY650-tubulin speed\nat 4 HPI (um/hr)",
       color = "Crypts:") +
  theme_bw() + 
  #geom_hline(yintercept = 118, color = 'red') + 
  theme(axis.text = element_text(size=5), 
        axis.title = element_text(size=8), 
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        legend.key.size = unit(0.5, 'lines'),
        legend.margin = margin(0, 0, 0, 0),
        legend.spacing.x = unit(0, 'lines'),
        legend.box.margin = margin(-8, -8, -8, -8),
        legend.position = c(0.2, 0.8),
        plot.margin = margin(5, 10, 0, 5))
ggsave("plots/s3a.pdf", width = 3.2, height = 3, units = "in")
s3a.small <- dataa[, c("csv", "donor", "expt.name", "sample.no", "piv.median.speed", 'crypts', 'jammed')]
write.csv(s3a.small,"sourcedata/s3a.csv", row.names = FALSE)


largs <- list(set_varnames = list(jammed = "Jammed Axis Title", mucus.simple = "Mucus Simple Axis Title"))


# cross back with the info i want
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
dataa = subset(daily.spots.summary.j, ((condition == 'vanilla') | (condition == 'shaken') | (condition == 'static') | (condition == 'rinsed') & (usable == 1)))
dataa$crypts <- as.factor(dataa$crypts)
dataa$jammed <- as.factor(dataa$jammed)
dataa$mucus.simple <- as.factor(dataa$mucus.simple)
#dataa$crypts <- factor(dataa$crypts, levels = c("clean", "edgefurrow", "furrows"))
#dataa$jammed <- factor(dataa$jammed, levels = c("jammed", "edge", "unjammed", 'hypermobile'))


dataa$crypts <- factor(dataa$crypts, levels = c("No Furrows", "Edge Furrows", "Center Furrows"))
dataa$jammed <- factor(dataa$jammed, levels = c("Jammed", "Edge Migration", "Widespread Migration", 'Hypermobile'))


table_data <- table(dataa$jammed, dataa$mucus.simple)
largs <- list(
  set_varnames = c(jammed = "Migration Extent", mucus.simple = "Mucus Pattern")
)
largs <- list(
  set_varnames = c(jammed = "Migration Extent", mucus.simple = "Mucus Pattern", crypts = 'Furrows'),
  rot_labels = c(0, 20, 10, 0),  # Adjust the rotation angles (top, right, bottom, left)
  label_args = list(
    cex = 0.2,  # Adjust the size of the labels
    font = 2    # Font style (1 = plain, 2 = bold, etc.)
  )
)

mosaic(~ jammed + mucus.simple, data = dataa,  
       labeling_args = largs, 
       shade = TRUE, 
       #gp = shading_max,
       legend = TRUE)
# save as pdf 7.5 x 7.5 in
mosaic(~ jammed + crypts, data = dataa,  
       labeling_args = largs, 
       shade = TRUE,
       #gp = shading_max,
       legend = TRUE)

s3bc.small <- dataa[, c("csv", "donor", "expt.name", "sample.no", "piv.median.speed", 'crypts', 'jammed', 'mucus.simple')]
write.csv(s3bc.small,"sourcedata/s3bc.csv", row.names = FALSE)


mosaic(~ jammed + crypts, data = dataa, shade = TRUE, legend = TRUE, rot_labels = c(0, 20, 0, 90))
mosaic(~ jammed + crypts, data = dataa, shade = TRUE, legend = TRUE, rot_labels = c(0, 20, 0, 90))
mosaic(~ jammed + mucus.simple, data = dataa, shade = TRUE, legend = TRUE, rot_labels = c(0, 20, 90, 90))

mosaic(~ jammed + crypts, data = dataa, shade = TRUE, legend = TRUE, rot_labels = c(0, 20, 0, 90))

mosaic(~ jammed + mucus.simple, data = dataa, shade = TRUE, legend = TRUE, rot_labels = c(0, 20, 90, 90))



