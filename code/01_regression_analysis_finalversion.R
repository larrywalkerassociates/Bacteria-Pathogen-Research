##find data correlations: starting point, overall glimpse of data patterns
#becky
rm(list = ls())

#exclude non-detects

library(tidyverse)
library(readxl)
library(ggplot2)
library(GGally)

#import data:
readRenviron(".Renviron")
data_path <- Sys.getenv("DATA_PATH")

#2026 data (import on 2026-09-27)
data<-read_csv(file.path(data_path,"RAS analysis","input", "pathogen_data_pt2_RAS.csv"))%>% 
  dplyr::filter(Replicate != 2) |> 
  filter(SampleTypeCode == c("Water_Grab")) |> 
  select(SampleDate, StationCode, SampleTypeCode, AnalyteName, Result,
         EventType, StationType) 
  
data_wide <- data |> 
  pivot_wider(names_from = AnalyteName, values_from = Result) |> 
  select(-c("c.coli","c.jejuni","c.lari", "16scampy")) |> 
  mutate('E. coli' = "", data_source = "data provided on 2026-09-23")
  

#data analysis:
indicators <- c("HF183", "E. coli", "crassphage", "pmmov")

pathogens <- c("norog1", "norog2", "adv")
      #ask Ryan- are pathogens typically also log transformed?
      #what is 16scampy?


df_1_pathogens_add<- data_wide |> 
  dplyr::mutate(norog1_rev = ifelse(norog1 == -88, 0, 
                                    ifelse(is.na(norog1), NA, norog1)),
                norog2_rev = ifelse(norog2 == -88, 0, 
                                    ifelse(is.na(norog2), NA, norog2)),
                adv_rev = ifelse(adv == -88, 0, 
                                 ifelse(is.na(adv), NA, adv))) |> 
  dplyr::mutate(pathogens_add = rowSums(across(norog1_rev:adv_rev), na.rm = TRUE)) |>
  dplyr::mutate(pathogens_add = ifelse((is.na(norog1_rev) & is.na(norog2_rev) & is.na(adv_rev)), NA, pathogens_add)) |> 
  select(-norog1_rev,-norog2_rev,-adv_rev) |> 
  filter(!is.na(pathogens_add)) |> 
  mutate(path_detection = ifelse(pathogens_add == c("0"), paste(c("ND")), paste(c("Detect"))))
write.csv(df_1_pathogens_add, file.path(data_path,"RAS analysis",
                                        "output", "pathogen_sum_2026-09-27.csv"))


#isolate second dataset, only variables no character columns      
df_2 <- df_1_pathogens_add |> 
  select(-SampleDate, -StationCode, -SampleTypeCode, -StationType, -path_detection, -EventType, -data_source, -'E. coli') |> 
  mutate(Entero = as.numeric(Enterococcus)) |> 
  select(-Enterococcus) |> 
  dplyr::mutate(pathogens_add = ifelse(pathogens_add == 0, 1, pathogens_add))

df_2[df_2 == -88] <- 1
#all NDs  should be set equal to 1 at this point (both indicators and pathogens, pathogen sum)

#remove NA pathogens for concentration regression analysis (including NDs wouldn't make sense)
df_2_noND <- df_2 |> 
  filter(pathogens_add != 1)

df2_num <- as.matrix(df_2_noND)


#log10-defer to this for normalization
df_log10 <- log10(df2_num)

#final dataset doesn't include NDs of pathogens)

#idea: when you have a detected pathogen, is there a relationship to the indicator?
df_3 <- as.data.frame(df_log10)


library("Hmisc") #for rcorr
spearman_rcorr <- rcorr(df_log10, type = "spearman") #non-parametric
# pearson_rcorr <- rcorr(df_1_log, type = "pearson") #parametric

P_values <- as.data.frame(as.table(spearman_rcorr$P)) %>% 
  rename(Spearman_P_value = Freq)
R_values <- as.data.frame(as.table(spearman_rcorr$r)) %>% 
  rename(Spearman_R_value = Freq)
n_values <- as.data.frame(as.table(spearman_rcorr$n)) %>% 
  rename(Spearman_n_value = Freq)

spearman_eval <- P_values %>% 
  left_join(R_values) %>% 
  left_join(n_values)

write.csv(spearman_eval, file.path(data_path, "RAS analysis","output", "correlations",  "spearman_evaluation_2026-09-27_pathogendetects.csv"))
# conclusion as of Aug 2026: 
  # need more detected pathogen data to draw reasonable interpretation on pathogen-indicator or pathogen-pathogen relationship

p1 <- ggplot(df_3, aes(x = crassphage, y = pathogens_add))+
  geom_point(color= "blue", size= 2)+
  geom_smooth(method = "lm", se= TRUE, color = "red", linetype = "solid")+
  labs(x = "Result crass (log)",
       y = "Pathogens (log)")+
  theme_minimal()
print(p1)  

p2 <- ggpairs(df_2)
print(p2)
ggsave(file.path(data_path2, "RAS analysis","correlatinos", "ggpair_correlograms.pdf"), plot = p2, bg = "white", scale = 1)




#trying multi-linear regression analyses:-------------------------------------------------------------------------------------------------------------------------

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




  