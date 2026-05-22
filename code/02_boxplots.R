##binary evaluation
#becky

# include non-detects

library(tidyverse)
library(readxl)
library(ggplot2)
library(rstatix) #library of functions for statistical analysis
library(coin)

#prep data:---------------------------------------------------------------------------------------------------------------------------
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
df_1[df_1 == -88] <- 1.1 #not just 1, bc log10 of 1 is 0

#add together pathogens
df_2 <- df_1 %>% 
  dplyr::mutate(pathogens_add = norog1 + norog2 + adv) %>% 
  select(-norog1, -norog2, -adv) %>% 
  select(-SampleDate, -StationCode)

#apply log transform, then add binary code
df2_num <- as.matrix(df_2)
df_log10 <- log10(df2_num) %>%
  as.data.frame() %>% 
  mutate(pathogen_binary = ifelse(pathogens_add <= 0.6, "ND", "Detect"))

#for boxplot: pivot from wide to long format data
df_long <- df_log10 %>% 
  pivot_longer(!pathogen_binary, 
               names_to = "indicator",
               values_to = "result")


#boxplot showing individual indicator results compared to detections of pathogens-----------------------------------------------------

## CRASSPHAGE vs. PATHOGENS------------------------------------
crassphage_test <- df_long %>% 
  filter(indicator == "crassphage") 

#### stats: 
    #mann-whitney test (non-parametric and data < 30 samples):
    #mann-whitney: an independent two-samples test (two-sample rank-sum test)

hist(crassphage_test$result)
crass_sum_stats <- get_summary_stats(crassphage_test, type = "median_iqr")

# using library(rstatix):

#wilcox
stat.test <- crassphage_test %>% 
  rstatix::wilcox_test(result ~ pathogen_binary) %>% 
  add_significance()
stat.test
print(stat.test)

#https://www.rdocumentation.org/packages/rstatix/versions/0.7.3/topics/wilcox_effsize
# library(coin)
# coin::wilcox_test()
effect_sz <- crassphage_test %>% 
  wilcox_effsize(result ~ pathogen_binary) #result is numeric variable, pathogen_binary is the factor with one or multiple levels giving the corresponding groups
print(effect_sz)
#effect size: way to quanitfy the difference between two groups: 
# https://www.statology.org/effect-size/:
# while p-value can tell us whether or not there is a statistically sig. difference between two groups, effect size can tell us HOW LARGE this difference actually is
# in practice, effect sizes are much more interesting and useful to know than p-value
# #The r value varies from 0 to close to 1. The interpretation values for r commonly in published literature and on the internet are: 0.10 - < 0.3 (small effect), 0.30 - < 0.5 (moderate effect) and >= 0.5 (large effect).

# using alternate wilcox function
crass_nonp <- wilcox.test(result~pathogen_binary, data= crassphage_test)
print(crass_nonp)


#### plot:
crass <- ggplot(crassphage_test, aes(x = pathogen_binary, y = result, fill = pathogen_binary))+
  geom_boxplot()+
  geom_jitter()+
  labs(title = get_test_label(stat.test, detailed = TRUE),
       subtitle = paste("effect size: ", effect_sz$magnitude, "(r = < 0.3)"),
       x = "Pathogen Detection",
       y = "log(crassphage conc.)")+
  theme_classic()+
  theme(legend.position="none")
crass
ggsave(file.path(data_path,"crassphage_two-sample_test.pdf"), plot = crass, bg = "white", scale = 1)


## HF183 vs. PATHOGENS------------------------------------
HF183_test <- df_long %>% 
  filter(indicator == "HF183")

#wilcox
stat.test <- HF183_test %>% 
  rstatix::wilcox_test(result ~ pathogen_binary) %>% 
  add_significance()
stat.test
print(stat.test)

#check, aother function
stat2 <- wilcox.test(result~pathogen_binary, data= HF183_test)
print(stat2)

#effect size
effect_sz <- HF183_test %>% 
  wilcox_effsize(result ~ pathogen_binary) 
print(effect_sz)

#### plot:
hf183 <- ggplot(HF183_test, aes(x = pathogen_binary, y = result, fill = pathogen_binary))+
  geom_boxplot()+
  geom_jitter()+
  labs(title = get_test_label(stat.test, detailed = TRUE),
       subtitle = paste("effect size: ", effect_sz$magnitude, "(r = < 0.3)"),
       x = "Pathogen Detection",
       y = "log(HF183 conc.)")+
  theme_classic()+
  theme(legend.position="none")
hf183
ggsave(file.path(data_path,"HF183_two-sample_test.pdf"), plot = hf183, bg = "white", scale = 1)



## E.coli vs. PATHOGENS------------------------------------
e.coli_test <- df_long %>% 
  filter(indicator == "E. coli")

#wilcox
stat.test <- e.coli_test %>% 
  rstatix::wilcox_test(result ~ pathogen_binary) %>% 
  add_significance()
stat.test
print(stat.test)

#check, aother function
stat2 <- wilcox.test(result~pathogen_binary, data= e.coli_test)
print(stat2)

#effect size
effect_sz <- e.coli_test %>% 
  wilcox_effsize(result ~ pathogen_binary) 
print(effect_sz)

#### plot:
e.coli <- ggplot(e.coli_test, aes(x = pathogen_binary, y = result, fill = pathogen_binary))+
  geom_boxplot()+
  geom_jitter()+
  labs(title = get_test_label(stat.test, detailed = TRUE),
       subtitle = paste("effect size: ", effect_sz$magnitude, "(r = < 0.3)"),
       x = "Pathogen Detection",
       y = "log(E.coli conc.)")+
  theme_classic()+
  theme(legend.position="none")
e.coli
ggsave(file.path(data_path,"e.coli_two-sample_test.pdf"), plot = e.coli, bg = "white", scale = 1)

## pmmov vs. PATHOGENS------------------------------------
pmmov_test <- df_long %>% 
  filter(indicator == "pmmov")

#wilcox
stat.test <- pmmov_test %>% 
  rstatix::wilcox_test(result ~ pathogen_binary) %>% 
  add_significance()
stat.test
print(stat.test)

#check, aother function
stat2 <- wilcox.test(result~pathogen_binary, data= pmmov_test)
print(stat2)

#effect size
effect_sz <- pmmov_test %>% 
  wilcox_effsize(result ~ pathogen_binary) 
print(effect_sz)

#### plot:
pmmov <- ggplot(pmmov_test, aes(x = pathogen_binary, y = result, fill = pathogen_binary))+
  geom_boxplot()+
  geom_jitter()+
  labs(title = get_test_label(stat.test, detailed = TRUE),
       subtitle = paste("effect size: ", effect_sz$magnitude, "(r = < 0.3)"),
       x = "Pathogen Detection",
       y = "log(pmmov conc.)")+
  theme_classic()+
  theme(legend.position="none")
pmmov
ggsave(file.path(data_path,"ppmov_two-sample_test.pdf"), plot = pmmov, bg = "white", scale = 1)
