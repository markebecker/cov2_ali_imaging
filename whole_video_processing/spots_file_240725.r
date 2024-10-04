# making a new big and beautiful spots file for the revisions.

# libraries for data manipulation and vis
library(dplyr)
library(data.table)
library(tidyr)
library(readr)
library(ggplot2)

#function for getting the spots info all there
# between the triple hashes is the good version.
# signed mb 6/20/24.
###########################
#######################
#####################
summarize_spots <- function(path, key) {
  # List all files in the specified path
  print("Listing all files in the specified path...")
  filenames <- list.files(path = path, full.names = TRUE)
  filenames_base <- basename(filenames)
  
  # Filter filenames to include only those present in the key
  print("Filtering filenames to include only those present in the key...")
  valid_filenames <- filenames[filenames_base %in% key$csv]
  
  process_csv <- function(file) {
    print(paste("Processing file:", basename(file)))
    data <- read.csv(file, stringsAsFactors = FALSE)
    trimmed <- data[-(1:3),]
    summarized_data <- trimmed %>%
      group_by(FRAME) %>%
      summarize(csv = basename(file),
                freq = n(), 
                med.max.gfp = median(as.numeric(MAX_INTENSITY_CH3)), 
                med.max.spytub = median(as.numeric(MAX_INTENSITY_CH1)), 
                med.max.nv = median(as.numeric(MAX_INTENSITY_CH2)))
    return(summarized_data)
  }
  
  # Process valid CSVs and combine the results
  print("Processing valid CSV files...")
  result <- lapply(valid_filenames, process_csv)
  print("Combining results...")
  combined_result <- do.call(rbind, result)
  combined_result$FRAME <- as.numeric(combined_result$FRAME)
  combined_result$FRAME <- combined_result$FRAME + 1
  
  # Complete frames with missing data
  print("Completing frames with missing data...")
  
  # Filter out rows in key where video.end.frame is NA
  key <- key %>% filter(!is.na(video.end.frame))
  
  all_frames <- key %>%
    group_by(csv, video.end.frame) %>%
    reframe(FRAME = as.integer(seq_len(first(video.end.frame)))) %>%
    ungroup()
  
  # Prune all_frames to include only files in combined_result
  print("Pruning all_frames to include only files in combined_result...")
  all_frames <- semi_join(all_frames, combined_result, by = "csv")
  
  # Add rows from all_frames that are not present in combined_result
  print("Merging all_frames with combined_result...")
  merged <- left_join(all_frames, combined_result, by = c('csv', 'FRAME'))
  
  # Fill rows with information from key
  print("Filling rows with information from key...")
  filled <- left_join(merged, key, by = "csv")
  filled <- filled %>%
    mutate(freq = replace_na(freq, 0))
  
  print("Processing complete.")
  return(filled)
}

# final experiment log
exptlog_iii <-fread("D:/000_rebuttal/exptlog_iii.csv")
path <- "D:/0_trackmates/"
spotsall <- summarize_spots(path, exptlog_iii)
fwrite(spotsall, 'D:/000_rebuttal/240702_spotsall.csv')

spotsall$norfreq 
norspots <- (spotsall %>% group_by(csv) %>% mutate(norfreq = as.numeric(freq)-min(as.numeric(freq))))

# fix the d5d6 timing issue...
# scope crashed in the middle and left it with a fairly large gap.
forspots <- norspots %>%
  mutate(FRAME = ifelse(expt.name == 'd5d6' & FRAME > 29, FRAME + 8, FRAME))

forspots$time <- (forspots$FRAME-1)*2+1

fwrite(forspots, 'D:/000_rebuttal/240702_goodspots.csv')

good.spots.summary <- forspots %>%
  group_by(csv) %>%
  summarise(good.peak.spots = max(norfreq), good.peak.spots.time = time[which.max(norfreq)])
good.spots.summary.j <- inner_join(good.spots.summary, exptlog_iii)

fwrite(good.spots.summary.j, 'D:/000_rebuttal/240620_goodspots_summary.csv')
##########################
##########################
##########################

# making sure it all looks good.
library(ggplot2)

library(broom)



dataa = subset(good.spots.summary.j, (virus == 'egfp') & (condition == 'vanilla') & (dyes != "cmo"))
xvar = as.numeric(dataa$t0spy500)
yvar = as.numeric(dataa$good.peak.spots)
# Fit a linear model
model <- lm(yvar ~ xvar, data = dataa)

model_summary <- tidy(model)

# Calculate R-squared
r_squared <- summary(model)$r.squared


ggplot(dataa, aes(x=xvar, y = yvar)) +
  annotate("text", x = Inf, y = Inf, 
           label = paste("R-squared: ", round(r_squared, 2), "\n",
                         "p-value: ", format.pval(model_summary$p.value[2], digits = 2)), 
           hjust = 1.1, vjust = 2, size = 3, color = "blue") +
  geom_smooth(method = 'lm', formula = y ~ x, se = TRUE, aes(group = FALSE, color = NULL), show.legend = FALSE) +
  geom_point()


ggplot(subset(forspots, (expt.name == "continual2") & (condition == 'drinse')), aes(time, as.numeric(norfreq), colour=factor(sample.no))) +
  labs(x = 'Hours post infection', y ="# of GFP+ spots", color = "Culturepoint")+
  theme_bw() +
  geom_vline(xintercept = c(20, 44, 68, 100, 168, 216))+
  geom_point() 



ggplot(norspots, aes(time, as.numeric(norfreq), colour=factor(from))) +
  labs(x = 'Hours post infection', y ="# of GFP+ spots", color = "Culturepoint")+
  theme_bw() +
  ylim(0, 25000)+
  geom_point() 


ggplot(norspots, aes(x = freq[which(from == "old")], y = freq[which(from == "new")])) +
  geom_point() +
  labs(x = "Frequency (old)", y = "Frequency (new)",
       title = "Comparison of Frequency (old) vs Frequency (new)")


X240605_goodspots %>%
  filter(file %in% c("multitest_pt7_concat_aligned_cropped-all-spots.csv",
                     "longcilia2_pt10_aligned_cropped-all-spots.csv",
                     "mb_longcilia_pt05_aligned_cropped-all-spots.csv",
                     "mb_longcilia_pt07_aligned_cropped-all-spots.csv",
                     "mb_longcilia_pt06_aligned_cropped-all-spots.csv",
                     "mb_longcilia_pt04_aligned_cropped-all-spots.csv",
                     "mb_220616_ruxcilia_pt2_stc_concat_all_aligned_cropped-all-spots.csv")) %>%
  summarise(max_norfreq = max(as.numeric(norfreq)))
