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

detect_comparison <- df_long %>% 
  filter(pathogen_binary == "Detect")


#### plot: no nds
dc_plot <- ggplot(detect_comparison, aes(x = indicator, y = result))+
  geom_boxplot()+
  geom_jitter()+
  theme_classic()+
  theme(legend.position="none")
dc_plot
ggsave(file.path(data_path,"detect_comparison.pdf"), plot = dc_plot, bg = "white", scale = 1)


#### plot: nd comparison
dc_plot2 <- ggplot(df_long, aes(x = indicator, y = result, fill = pathogen_binary))+
  geom_boxplot()+
  #geom_jitter()+
  theme_classic()
  #theme(legend.position="none")
dc_plot2
ggsave(file.path(data_path,"all_comparison.pdf"), plot = dc_plot2, bg = "white", scale = 1)


#no log transform-----------------------------------------------------------------------
df2_nolog <- df_2 %>% 
  mutate(pathogen_binary = ifelse(pathogens_add <= 3.4, "ND", "Detect"))

df_long_nolog <- df2_nolog %>% 
  pivot_longer(!pathogen_binary, 
               names_to = "indicator",
               values_to = "result")

#### plot: nd comparison non-log
test <- ggplot(df_long_nolog, aes(x = indicator, y = result, fill = pathogen_binary))+
  geom_boxplot()+
  #geom_jitter()+
  theme_classic()
#theme(legend.position="none")
test

