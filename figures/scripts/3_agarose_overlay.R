# Agarose overlay (Fig. 4)

# libraries
library(readr)
library(ggplot2)
library(dplyr)
library(ggbeeswarm)

spots <- read_csv("data/240702_goodspots.csv")
spots <- subset(spots, 
                !(expt.name == 'd5d6' & FRAME > 58) & # erroneously set freq to 0; data does not exist
                  !(expt.name == 'd4shake' & FRAME > 59) & #erroneously set freq to 0; data does not exist
                  !(expt.name == 'd2shake' & (FRAME > 46 & FRAME < 50)) & # think it was out of focus & detected no spots in this time
                  !(expt.name == 'agarose_d3' & (FRAME > 35 & FRAME < 48)))
relevant <- subset(spots, expt.name == 'agarose' | expt.name == 'agarose_d3' | expt.name == 'd2d4ag' | expt.name == 'd3ag2')

# fig 3d. plot points over time from the agarose experiments
f3data <- subset(relevant, (condition == "vanilla" | condition == "agarose") & (donor == 2 | donor == 3 | donor == 4))
ggplot(subset(relevant, (condition == "vanilla" | condition == "agarose") & (donor == 2 | donor == 3 | donor == 4)), aes(time, as.numeric(norfreq), colour=factor(condition),  group = interaction(expt.name, sample.no))) +
  labs(x = 'Hours post infection', y ="# of GFP+ spots", color = "Condition")+
  theme_bw() +
  #facet_wrap(~expt.name) +
  scale_color_discrete(labels = c("Agarose", "Plain"), name = 'Condition:')+
  theme(axis.text = element_text(size=6), axis.title = element_text(size=8), 
        legend.key.size = unit(1, 'point'),
        legend.box.margin = margin(0, 0, 0, -10),
        legend.text = element_text(size=6), legend.title = element_text(size=8),
        plot.margin = margin(1, 0, 0, 0)) +
  guides(colour = guide_legend(override.aes = list(size=2))) +
  geom_line(alpha = 0.3, size = 0.5) +
  stat_smooth(method = "loess", se = TRUE, aes(group = factor(condition)))
ggsave("plots/fig3d.pdf", width = 2, height = 1.35, units = "in")
f3d.small <- f3data[, c("csv", "donor", "expt.name", "sample.no", "time", 'norfreq', 'condition')]
write.csv(f3d.small,"sourcedata/fig3d.csv", row.names = FALSE)

# plotting peak spots... i'd rather not bc some are still increasing at last time on the scope.
relevantt <- subset(daily.spots.summary.j, (expt.name == 'agarose' | expt.name == 'agarose_d3' | expt.name == 'd2d4ag' | expt.name == 'd3ag2') & virus == 'egfp')
ggplot(subset(relevantt, (condition == "vanilla" | condition == 'agarose') & virus == 'egfp'), (aes(y=good.peak.spots, x=factor(condition)))) + 
  labs(x = 'Condition', y ="Peak # GFP+ spots", color = "Donor:")+
  geom_quasirandom(alpha = 0.8, show.legend = TRUE) +
  stat_summary(fun.y = median,
               fun.ymin = function(z) { quantile(z,0.25) }, 
               fun.ymax = function(z) { quantile(z,0.75) }, 
               geom = "pointrange",
               show.legend = FALSE,
               aes(group = factor(condition))) +
  scale_x_discrete(labels = c("Agarose", "Plain"))+
  theme_bw() + 
  ylim(0, 25000) +
  theme(axis.text = element_text(size=5), axis.title = element_text(size=8),
        plot.margin = margin(1, 0, 0, 0),
        #legend.key.size = unit(0.1, 'point'),
        legend.box.margin = margin(0, 0, 0, -10),
        legend.text = element_text(size=6), legend.title = element_text(size=6))
ggsave("plots/fig3e_NOoo.pdf", width = 1, height = 1.35, units = "in")
wilcox.test(subset(relevantt, condition == 'vanilla')$good.peak.spots, subset(relevantt, condition == 'agarose')$good.peak.spots)

# plotting spots at 4 dpi. this is good enough imo.
relevantt <- subset(daily.spots.summary.j, (expt.name == 'agarose' | expt.name == 'agarose_d3' | expt.name == 'd2d4ag' | expt.name == 'd3ag2') & virus == 'egfp')
f3edata <- subset(relevant, (condition == "vanilla" | condition == 'agarose') & virus == 'egfp' & time == 95)
ggplot(subset(relevant, (condition == "vanilla" | condition == 'agarose') & virus == 'egfp' & time == 95), (aes(y=norfreq, x=factor(condition)))) + 
  labs(x = 'Condition', y ="GFP+ spots at 4 DPI", color = "Donor:")+
  geom_quasirandom(alpha = 0.6, show.legend = TRUE) +
  stat_summary(fun = mean,
               fun.min = function(x) mean(x) - sd(x), 
               fun.max = function(x) mean(x) + sd(x), 
               geom = "pointrange",
               show.legend = FALSE,
               aes(group = factor(condition))) +
  scale_x_discrete(labels = c("Agarose", "Plain"))+
  theme_bw() + 
  ylim(-100, 25000) +
  theme(axis.text = element_text(size=5), axis.title = element_text(size=8),
        plot.margin = margin(1, 0, 0, 0),
        #legend.key.size = unit(0.1, 'point'),
        legend.box.margin = margin(0, 0, 0, -10),
        legend.text = element_text(size=6), legend.title = element_text(size=6))
ggsave("plots/fig3e.pdf", width = 1, height = 1.35, units = "in")
wilcox.test(subset(relevant, condition == 'vanilla' & time == 95 & virus == 'egfp')$norfreq, subset(relevant, condition == 'agarose' & time == 95 & virus == 'egfp')$norfreq)
f3e.small <- f3edata[, c("csv", "donor", "expt.name", "sample.no", "time", 'norfreq', 'condition')]
write.csv(f3e.small,"sourcedata/fig3e.csv", row.names = FALSE)
