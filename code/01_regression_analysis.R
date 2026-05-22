##find data correlations: starting point, overall glimpse of data patterns
#becky

#exclude non-detects

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

df_1 <- data %>% 
  dplyr::mutate(norog1 = ifelse(norog1 == -88, 0, norog1)) %>% 
  dplyr::mutate(norog2 = ifelse(norog2 == -88, 0, norog2)) %>% 
  dplyr::mutate(adv = ifelse(adv == -88, 0, adv)) %>% 
  dplyr::mutate(pathogens_add = norog1 + norog2 + adv) %>% 
  select(-norog1, -norog2, -adv) %>% 
  dplyr::mutate(pathogens_add = ifelse(pathogens_add == 0, NA, pathogens_add)) %>% 
  select(-SampleDate, -StationCode)


df1_num <- as.matrix(df_1)

#log10-defer to this for normalization
df_log10 <- log10(df1_num)
df_log10[!is.finite(df_log10)] <- NA

library("Hmisc") #for rcorr
spearman_rcorr <- rcorr(df_log10, type = "spearman")
# pearson_rcorr <- rcorr(df_1_log, type = "pearson")

P_values <- as.data.frame(as.table(spearman_rcorr$P)) %>% 
  rename(Spearman_P_value = Freq)
R_values <- as.data.frame(as.table(spearman_rcorr$r)) %>% 
  rename(Spearman_R_value = Freq)
n_values <- as.data.frame(as.table(spearman_rcorr$n)) %>% 
  rename(Spearman_n_value = Freq)

spearman_eval <- P_values %>% 
  left_join(R_values) %>% 
  left_join(n_values)

write.csv(spearman_eval, file.path(data_path,"spearman_evaluation_2026-05-21.csv"))





#trying multi-linear regression analyses:

#multiple linear regression (baseline mutlivariable correlation)- most direct extension of correlation
    #r2 = variance explained
    #adjusted r2 = penalized r2 (compare models)-- so will increase if number of variables improves the model, and will decrease if number of variables impacts the model
    # helps to decide between fewer predictors or more predictors 
    #p-values-- which indicators matter
    #F-statistic: overall model significance- confusing

#first analysis: multiple linear regression
df_1_log<- as.data.frame(df_1_log)

m_lm_1 <- lm(df_1_log$pathogens_add ~ df_1_log$HF183 + df_1_log$`E. coli` + df_1_log$pmmov + df_1_log$crassphage)
summary(m_lm_1)
plot(m_lm_1)

#note: p-value of crassphage variable to the pathogens is only significant one identified (star)
# p of 0.0348
  #adjusted R2 = 0.32
#f-statistic = 2.578 on 4 and 9 DF
#p-value of 0.1095

#second analysis: multiple linear regression- different variables
m_lm_2 <- lm(df_1_log$pathogens_add ~ df_1_log$HF183 + df_1_log$crassphage)
summary(m_lm_2)
plot(m_lm_2)

#adjusted r2 = 0.36 (increased slightly)
# F statistic = 4.681 on 2 and 11 DF (increased slightly)
# p-value: 0.033 (more significant)

#thir analysis: single value
lm_3 <- lm(df_1_log$pathogens_add ~ df_1_log$crassphage)
summary(lm_3)
plot(lm_3)
#adjusted r2 = 0.37




  