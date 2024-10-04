##looking at .pts files to figure out which points go to which culture.

#   _____                 _   _                   _                                     
#  / ____|               | | | |                 | |                                    
# | |  __  ___   ___   __| | | |     ___  _ __ __| |                                    
# | | |_ |/ _ \ / _ \ / _` | | |    / _ \| '__/ _` |                                    
# | |__| | (_) | (_) | (_| | | |___| (_) | | | (_| |_                                   
# _\_____|\___/_\___/ \__,_|_|______\___/|_|  \__,_( )_                                 
# \ \        / / |         | |   (_)               |/(_)                                
#  \ \  /\  / /| |__   __ _| |_   _ ___    __ _  ___  _ _ __   __ _    ___  _ __        
#   \ \/  \/ / | '_ \ / _` | __| | / __|  / _` |/ _ \| | '_ \ / _` |  / _ \| '_ \       
#    \  /\  /  | | | | (_| | |_  | \__ \ | (_| | (_) | | | | | (_| | | (_) | | | |      
#     \/  \/   |_| |_|\__,_|\__| |_|___/  \__, |\___/|_|_| |_|\__, |  \___/|_| |_|      
#                                          __/ |               __/ |                    
#                                         |___/ _         _   |___/               _ _ _ 
#                                              (_)       | | | |                 | | | |
#                                               _ _ __   | |_| |__   ___ _ __ ___| | | |
#                                              | | '_ \  | __| '_ \ / _ \ '__/ _ \ | | |
#                                              | | | | | | |_| | | |  __/ | |  __/_|_|_|
#                                              |_|_| |_|  \__|_| |_|\___|_|  \___(_|_|_)
                                                                                       
                                                                                       
library(ggplot2)
library(dplyr)

#file = "//fsmresfiles.fsm.northwestern.edu/fsmresfiles/Basic_Sciences/CDB/HopeLab/Users/mbecker/2_sarscov2/240415_longcilia2/raw/longcilia_tfinal.pts"
#file = 'D:/multicilia/multitest_points.pts'
#file = 'D:/longcilia/longcilia_tfinal.pts'
file = "//fsmresfiles.fsm.northwestern.edu/fsmresfiles/Basic_Sciences/CDB/HopeLab/Users/mbecker/2_sarscov2/240612_omicron/omicron_t2.pts"

interest <- read.table(file, quote="\"", comment.char="")
# name the columns
colnames(interest) <- c('point', 'x', 'y', 'z')
# strip the colons off the point names & convert to numeric
interest$point <- as.numeric(sub(":.*", "", interest$point))


#pick what points I want to have labelled. usually all of them is too dense.
#justafew <- subset(interest, as.numeric(interest$point) >=1 & as.numeric(interest$point) < 11)
#justafew <- sample_n(interest, 100)
justafew <- subset(interest, point >= 0 & point <= 500)

# plot the points in x and y
# the thang flips the axes from how it's shown on the scope
ggplot(interest, aes(x = -x, y = -y)) +
  geom_point(alpha = 0.1) +
  geom_text(data = justafew, aes(label = point))

