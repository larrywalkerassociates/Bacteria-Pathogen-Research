##find data correlations: starting point, overall glimpse of data patterns
#becky

# include non-detects

library(tidyverse)
library(readxl)
library(ggplot2)

#import data:
readRenviron(".Renviron")
data_path <- Sys.getenv("DATA_PATH")

data<-read_csv(file.path(data_path, "pathogen_data.csv")) %>% 
  select(-c("c.coli","c.jejuni","c.lari", "16scampy"))

indicators <- c("HF183", "E. coli", "crassphage", "pmmov")

pathogens <- c("norog1", "norog2", "adv")
#ask Ryan- are pathogens typically also log transformed?
#what is 16scampy?

#set pathogen NDs(originally -88) to equal 1
#also set indicator NDs (originally -88) to equal 1
df_1 <- data
df_1[df_1 == -88] <- 1

#add together pathogens
df_2 <- df_1 %>% 
  dplyr::mutate(pathogens_add = norog1 + norog2 + adv) %>% 
  select(-norog1, -norog2, -adv) %>% 
  select(-SampleDate, -StationCode)


df2_num <- as.matrix(df_2)

#attempt to normalize data for non-normal distribution: (https://www.statology.org/transform-data-in-r/)

# #natural log- not sure we need this
# df_ln <- log(df2_num) %>% 
#   as.data.frame()
# hist(df_ln$HF183)

#log10-defer to this for normalization
df_log10 <- log10(df2_num) 

#test dist
df_test <- df_log10 %>% 
  as.data.frame()
hist(df_test$pathogens_add)


library("Hmisc") #for rcorr
spearman_rcorr <- rcorr(df_log10, type = "spearman")
# pearson_rcorr <- rcorr(df_1_log, type = "pearson")

P_values <- as.data.frame(as.table(spearman_rcorr$P)) %>% 
  rename(Spearman_P_value = Freq)
R_values <- as.data.frame(as.table(spearman_rcorr$r)) %>% 
  rename(Spearman_R_value = Freq)
n_values <- as.data.frame(as.table(spearman_rcorr$n)) %>% 
  rename(Spearman_n_value = Freq)

spearman_eval_withNDs <- P_values %>% 
  left_join(R_values) %>% 
  left_join(n_values)

write.csv(spearman_eval_withNDs, file.path(data_path,"spearman_evaluation_withND_2026-05-21.csv"))




  