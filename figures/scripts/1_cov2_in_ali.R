# SARS-CoV-2 in ALI cultures (Fig. 2)

#libraries
library(readr)
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)

######################
## 1d. % ciliation. ##
######################
goodspots_j <- read_csv("data/240620_goodspots_summary.csv")
dataa = subset(goodspots_j, (condition == 'vanilla' | condition == 'shaken' | condition == 'static') 
               & (dyes != "cmo" & dyes != 'undyed') 
               & (usable == 1)
               # these experiments had spytub added like right before imaging instead of >48 hours preinfection
               & (expt.name != 'd2shake') & (expt.name != 'd1d9') & expt.name != 'multitest')

ggplot(dataa, aes(x=as.numeric(dataa$t0spy350)/as.numeric(dataa$t0spy0))) +
  geom_histogram(bins = 30) +
  theme_bw() +
  scale_x_continuous(breaks = c(0, 0.5, 1), limits = c(0, 1)) +
  labs(x = "Fraction Tubulin+\n Area",
       y = "Count") +
  theme(axis.text = element_text(size=6), 
        axis.title = element_text(size=6), 
        plot.margin = margin(5, 10, 0, 5))
ggsave("plots/fig1dnew.pdf", width = 0.9, height = 1.1, units = "in")
dataa.small <- dataa[, c("csv", "donor", "expt.name", "sample.no", "t0spy350", 't0spy0')]

write.csv(dataa.small,"sourcedata/fig1d.csv", row.names = FALSE)

median(as.numeric(dataa$t0spy350)/as.numeric(dataa$t0spy0))

#############################
## 1g. N copies over time. ##
#############################

qpcr_j <- read_csv("data/qpcr_p3oldandflippies.csv")
keys <- colnames(qpcr_j)[!grepl('copies.mm.2',colnames(qpcr_j))]
X <- as.data.table(qpcr_j)
X$copies.mm.2 <- as.numeric(X$copies.mm.2)
qpcrcurve = X[,list(mm= mean(copies.mm.2), sd = sd(copies.mm.2), hpi = hpi, type = type, donor = donor),'sample.name']
trimmedq = distinct(qpcrcurve)

#p 0.054 at 48hpi for grouped
t.test(subset(trimmedq, hpi == 1)$mm, subset(trimmedq, hpi == 48)$mm)
#p 0.0035 at 72hpi for grouped
t.test(subset(trimmedq, hpi == 168)$mm, subset(trimmedq, hpi == 72)$mm)
# no significant differences between timepoints with matched inverted and conventional points
t.test(subset(trimmedq, hpi == 24 & type == "inv")$mm, subset(trimmedq, hpi == 24 & type == "conv")$mm)
t.test(subset(trimmedq, hpi == 48 & type == "inv")$mm, subset(trimmedq, hpi == 48 & type == "conv")$mm)
t.test(subset(trimmedq, hpi == 72 & type == "inv")$mm, subset(trimmedq, hpi == 72 & type == "conv")$mm)
t.test(subset(trimmedq, hpi == 144 & type == "inv")$mm, subset(trimmedq, hpi == 144 & type == "conv")$mm)
t.test(subset(trimmedq, hpi == 168 & type == "inv")$mm, subset(trimmedq, hpi == 168 & type == "conv")$mm)

t.test(subset(trimmedq, hpi == 168)$mm, subset(trimmedq, hpi == image)$mm)

shapiro.test(subset(trimmedq, hpi == 120)$mm)

#mean of the final timepoint
mean(subset(trimmedq, hpi == 168)$mm)

#the plot
colors <- c("inv" = "deeppink2", "conv" = "chartreuse3")
ggplot(trimmedq, aes(x=hpi, y=mm, color = factor(type))) + 
  labs(x = 'Hours Post Infection', y = bquote('N copies per'~mm^2), color = "ALI type")+
  scale_y_log10(limits = c(-1, 3e9)) +
  geom_point(alpha = 0.4, size = 1, position = position_dodge(width = 1)) +
  scale_color_manual(values = colors, labels = c("conv" = "Conventional", "inv" = "Inverted")) +
  #geom_smooth(method = "loess", se = FALSE)+
  stat_summary(aes(group = factor(type)), 
               fun = mean, 
               fun.min = function(x) mean(x) - sd(x), 
               fun.max = function(x) mean(x) + sd(x), 
               size = 0.25, 
               alpha = 0.8,
               position = position_dodge(width = 10),
               linewidth = 0.5, 
               geom = "pointrange") +
  theme_bw() + 
  #geom_hline(yintercept = 118, color = 'red') + 
  theme(axis.text = element_text(size=6), 
        axis.title = element_text(size=8), 
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        legend.key.size = unit(0.5, 'lines'),
        legend.margin = margin(0, 0, 0, 0),
        legend.spacing.x = unit(0, 'lines'),
        legend.box.margin = margin(-8, -8, -8, -8),
        #legend.justification = c("right", "bottom"),
        legend.position = c(0.8, 0.15),
        plot.margin = margin(5, 10, 0, 5))
ggsave("plots/fig1g.pdf", width = 3.1, height = 1.85, units = "in")
write.csv(trimmedq,"sourcedata/fig1g.csv", row.names = FALSE)



############################
## 1j & i. line profiles. ##
############################

cell <- read.csv("data/cell_values.csv")
cellslong <- cell %>%
  pivot_longer(cols = c("N", "S", "dsRNA"), 
               names_to = "channel", 
               values_to = "pixel_value")
colors <- c("N" = "deeppink1", "S" = "chartreuse3", "dsRNA" = "blue", "nuclei" = "grey")

ggplot(cellslong, aes(x = distance, y = pixel_value, color = channel)) +
  geom_line() +
  scale_color_manual(values = colors) +
  coord_flip() +
  scale_y_continuous(breaks = c(0, 50000)) +
  labs(x = "Distance (µm)",
       y = "Pixel Value",
       color = "Channel") +
  theme_bw() + 
  #geom_hline(yintercept = 118, color = 'red') + 
  theme(axis.text = element_text(size=6), 
        axis.title = element_text(size=8), 
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        legend.key.size = unit(0.5, 'lines'),
        legend.margin = margin(0, 0, 0, 0),
        legend.spacing.x = unit(0, 'lines'),
        legend.box.margin = margin(-8, -8, -8, -8),
        plot.margin = margin(5, 10, 0, 5))
ggsave("plots/fig1i.pdf", width = 2, height = 1.15, units = "in")
write.csv(cellslong,"sourcedata/fig1i.csv", row.names = FALSE)

point <- read.csv("data/point_values.csv")
pointslong <- point %>%
  pivot_longer(cols = c("N", "S", "dsRNA"), 
               names_to = "channel", 
               values_to = "pixel_value")
colors <- c("N" = "deeppink1", "S" = "chartreuse3", "dsRNA" = "blue", "nuclei" = "grey")

ggplot(pointslong, aes(x = distance, y = pixel_value, color = channel)) +
  geom_line() +
  scale_color_manual(values = colors) +
  coord_flip() +
  xlim(0, 4) +
  labs(x = "Distance (µm)",
       y = "Pixel Value",
       color = "Channel") +
  theme_bw() + 
  scale_y_continuous(breaks = c(0, 3000)) +
  #geom_hline(yintercept = 118, color = 'red') + 
  theme(axis.text = element_text(size=6), 
        axis.title = element_text(size=8), 
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        legend.key.size = unit(0.5, 'lines'),
        legend.margin = margin(0, 0, 0, 0),
        legend.spacing.x = unit(0, 'lines'),
        legend.box.margin = margin(-8, -8, -8, -8),
        plot.margin = margin(5, 10, 0, 5))
ggsave("plots/fig1j.pdf", width = 2, height = 1.15, units = "in")
write.csv(pointslong,"sourcedata/fig1j.csv", row.names = FALSE)
