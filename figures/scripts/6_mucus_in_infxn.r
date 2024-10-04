# Mucus changes in infection (Fig. 9, Supp. Fig 14)
# libraries
library(ggplot2)
library(data.table)
library(dplyr)
library(tidyr)
library(readr)
library(splines)


######################################
## Fig. 6b. Dot blot quantification ##
######################################
blot_quant <- read_csv("data/blot_quant.csv")
mock_24hpi_mean <- blot_quant %>%
  filter(dpi == 1, condition == 'mock') %>%
  summarize(mean.intden = mean(as.numeric(int.den), na.rm = TRUE)) %>%
  pull(mean.intden)
blot_quant <- blot_quant %>%
  mutate(norintden = as.numeric(int.den) / mock_24hpi_mean)
write.csv(blot_quant,"sourcedata/fig6b.csv", row.names = FALSE)


f6bdata <- subset(blot_quant, condition == 'mock' | condition == 'infected')
ggplot(f6bdata, aes(x = dpi, y = norintden, color = condition, group = condition)) +
  labs(x = 'Days post infection', y = bquote('MUC5AC\n(Relative Density)'))+
  #scale_y_log10(limits = c(1e7, 3e9)) +
  theme_bw() +
  geom_point(alpha = 0.6, size = 0.25) +
  scale_color_manual(values=c('cyan4', 'deeppink1'), labels = c("Infected", "Mock"), name = NULL) +
  stat_summary(fun = mean,
               fun.min = function(x) mean(x) - sd(x), 
               fun.max = function(x) mean(x) + sd(x),
               size = 0.5, 
               alpha = 0.8,
               linewidth = 0.5, 
               geom = "pointrange") + 
  #scale_x_discrete(labels = c("NT", "CYPA", "DNAH5", "DNAI1")) +
  theme(axis.text = element_text(size=5, face = 'italic'), axis.title = element_text(size=6),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.5, 'lines'),
        legend.margin = margin(0, 0, 0, 0),
        legend.spacing.x = unit(0, 'lines'),
        legend.box.margin = margin(-8, -8, -8, -8),
        plot.margin = margin(5, 1, 1, 1),
        legend.position = c(0.3, 0.85))
        #legend.position = 'none')
ggsave("plots/fig6b.pdf", width = 1.3, height =2, units = "in")
timee <- 6
g1 <- subset(f6bdata, condition == 'mock' & dpi == timee)$norintden
g2 <- subset(f6bdata, condition == 'infected' & dpi == timee)$norintden
t.test(g1, g2)

ok <- lm(norintden ~ dpi * condition, data = f6bdata)
summary(ok)
plot(ok)

# Function to perform a permutation test
permutation_test <- function(x, y, num_permutations = 10000) {
  observed_diff <- abs(mean(x) - mean(y))
  combined <- c(x, y)
  n <- length(x)
  permutation_diffs <- replicate(num_permutations, {
    permuted <- sample(combined)
    abs(mean(permuted[1:n]) - mean(permuted[(n+1):length(combined)]))
  })
  p_value <- mean(permutation_diffs >= observed_diff)
  return(p_value)
}

# Run the permutation test
p_value <- permutation_test(g1, g2)
print(paste("Permutation test p-value:", p_value))


library(lme4)
library(lmerTest) # for p-values
library(ggplot2)  # for visualization

model <- lmer(int.den ~ condition * dpi + (1 | donor), data = f6bdata)

summary(model)

########################
## Infection kinetics ##
########################

spots <- read_csv("data/240702_goodspots.csv")
# removing problematic parts
spots <- subset(spots, 
                !(expt.name == 'd5d6' & FRAME > 58) & # erroneously set freq to 0; data does not exist
                  !(expt.name == 'd4shake' & FRAME > 59) & #erroneously set freq to 0; data does not exist
                  !(expt.name == 'd2shake' & (FRAME > 46 & FRAME < 50)) & # think it was out of focus & detected no spots in this time
                  !(expt.name == 'agarose_d3' & (FRAME > 35 & FRAME < 48)))

# Fix frames with erroneously high spots.
spots <- spots %>%
  mutate(freq = case_when(
    freq > 50000 ~ 0,
    TRUE ~ freq
  )) %>%
  group_by(csv) %>%
  mutate(norfreq = as.numeric(freq) - min(as.numeric(freq)))

infected <- subset(spots, virus != 'mock')
rinsies <- subset(infected, (expt.name == 'd7d8' | expt.name == 'multitest' | expt.name == 'continual2'))
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
data <- rinsies %>%
  group_by(csv) %>%
  do(fit.splines(.))

rinsies <- rinsies %>%
  left_join(data %>% select(csv, time, fitted.norfreq, rate.of.change, accel), 
            by = c("csv", "time"))

#########################################
## Fig 6fg. D7D8 rinse quantification. ##
#########################################
d7d8 <- subset(rinsies, expt.name == 'd7d8')

ggplot(d7d8, aes(time, as.numeric(norfreq), colour=factor(donor), group = factor(sample.no))) +
  labs(x = 'Hours post infection', y ="# of GFP+ spots", color = "Donor:")+
  geom_vline(xintercept = 119, color = 'black', linetype = 'dashed') + 
  geom_line(alpha = 0.8) +
  theme_bw() + 
  #stat_summary(aes(group = donor)) +
  #stat_smooth(method = "loess", span = 0.1, aes(group = donor)) +
  guides(colour = guide_legend(override.aes = list(size=2))) +
  #geom_point(alpha = 0.6, size = 0.01) +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=8),
        plot.margin = margin(1, 1, 1, 1),
        legend.key.size = unit(0.1, 'point'),
        legend.box.margin = margin(0, 0, 0, -10),
        legend.position = c(0.2, 0.8),
        legend.text = element_text(size=6), legend.title = element_text(size=6))
ggsave("plots/fig6e.pdf", width = 2.3, height = 2.2, units = "in")

ggplot(d7d8, aes(time+4, as.numeric(rate.of.change), colour=factor(donor), group = factor(sample.no))) +
  labs(x = 'Hours post infection', y =expression(Delta ~ "GFP+ spots per hour"), color = "Donor:")+
  geom_vline(xintercept = 119, color = 'black', linetype = 'dashed') + 
  geom_line(alpha = 0.3) +
  theme_bw() + 
  geom_hline(yintercept = 0)+
  #stat_summary(aes(group = donor)) +
  stat_smooth(method = "loess", alpha = 0.3, span = 0.05, aes(group = donor)) +
  guides(colour = guide_legend(override.aes = list(size=2))) +
  #geom_point(alpha = 0.6, size = 0.01) +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=8),
        plot.margin = margin(1, 1, 1, 1),
        legend.position = c(0.2, 0.8),
        legend.key.size = unit(0.1, 'point'),
        legend.box.margin = margin(0, 0, 0, -10),
        legend.text = element_text(size=6), legend.title = element_text(size=6))
ggsave("C:/Users/Mark/Desktop/figure_images/6/revised/fig6f.pdf", width = 2.2, height = 2.2, units = "in")

f6ef.small <- d7d8[, c("csv", "donor", "expt.name", "sample.no", "time", 'norfreq', 'rate.of.change')]
write.csv(f6ef.small,"sourcedata/fig6ef.csv", row.names = FALSE)


# take the spots over time csv and extract specific values of interest
daily.spots.summary <- forspots %>%
  group_by(csv) %>%
  summarise(good.peak.spots = max(norfreq), 
            good.peak.spots.time = time[which.max(norfreq)],
            aucspots = sum(norfreq)) 
# cross back with the info i want
daily.spots.summary.j <- inner_join(daily.spots.summary, exptlog_iii)


################################################
##    S14. Continual2 rinse quantification.   ##
################################################

continual2 <- subset(rinsies, expt.name == 'continual2' & sample.no < 9)

# time +3 bc there was a 3 hour lag between finishing up the infections & getting it on the scope.
ggplot(continual2, aes(time+3, as.numeric(norfreq), colour=factor(donor), group = factor(sample.no))) +
  labs(x = 'Hours post infection', y ="# of GFP+ spots", color = "Donor:")+
  geom_vline(xintercept = 25, linetype = 'dashed') +
  geom_vline(xintercept = 44, linetype = 'dashed') +
  geom_vline(xintercept = 66, linetype = 'dashed') +
  geom_vline(xintercept = 88, linetype = 'dashed') +
  geom_vline(xintercept = 156, linetype = 'dashed') +
  geom_line() +
  theme_bw() + 
  #stat_summary(aes(group = donor)) +
  #stat_smooth(method = "loess", span = 0.1, aes(group = donor)) +
  guides(colour = guide_legend(override.aes = list(size=2))) +
  #geom_point(alpha = 0.6, size = 0.01) +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=8),
        plot.margin = margin(1, 0, 0, 0),
        legend.key.size = unit(0.1, 'point'),
        legend.box.margin = margin(0, 0, 0, -10),
        legend.text = element_text(size=6), legend.title = element_text(size=6))
ggsave("plots/fig6s_a.pdf", width = 2.5, height = 2.2, units = "in")

ggplot(continual2, aes(time+3, as.numeric(rate.of.change), colour=factor(donor), group = factor(sample.no))) +
  labs(x = 'Hours post infection', y = expression(Delta ~ "GFP+ spots per hour"), color = "Donor:")+
  geom_vline(xintercept = 23, linetype = 'dashed') +
  geom_vline(xintercept = 44, linetype = 'dashed') +
  geom_vline(xintercept = 64, linetype = 'dashed') +
  geom_vline(xintercept = 86, linetype = 'dashed') +
  geom_vline(xintercept = 154, linetype = 'dashed') +
  geom_line(alpha = 0.6) +
  theme_bw() + 
  #stat_summary(aes(group = donor)) +
  #stat_smooth(method = "loess", span = 0.1, aes(group = donor)) +
  guides(colour = guide_legend(override.aes = list(size=2))) +
  #geom_point(alpha = 0.6, size = 0.01) +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=8),
        plot.margin = margin(1, 0, 0, 0),
        legend.key.size = unit(0.1, 'point'),
        legend.box.margin = margin(0, 0, 0, -10),
        legend.text = element_text(size=6), legend.title = element_text(size=6))
ggsave("plots/fig6ds_b.pdf", width = 2.5, height = 2.2, units = "in")

continual2$time_corr <- continual2$time+3
s14.small <- continual2[, c("csv", "donor", "expt.name", "sample.no", "time_corr", 'norfreq', 'rate.of.change')]
write.csv(s14.small,"sourcedata/s14ab.csv", row.names = FALSE)



ggplot(rinsies, aes(time, as.numeric(accel), colour=factor(donor), group = interaction(expt.name, sample.no))) +
  labs(x = 'Hours post infection', y ="# of GFP+ spots", color = "Donor:")+
  theme_bw() +
  xlim(24,210) +
  #geom_label() +
  facet_wrap(~expt.name) +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=8), 
        legend.key.size = unit(1, 'point'),
        legend.box.margin = margin(0, 0, 0, -10),
        legend.text = element_text(size=6), legend.title = element_text(size=8),
        plot.margin = margin(1, 0, 0, 0)) +
  guides(colour = guide_legend(override.aes = list(size=2))) +
  geom_line(alpha = 0.6, size = 0.3) +
  stat_smooth(method = "loess", se = TRUE, aes(group = 'none'))