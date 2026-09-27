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
data_path2 <- Sys.getenv("DATA_PATH2")

#2026 data:
data<-read_csv(file.path(data_path2,"RAS analysis", "pathogen_data_pt2_RAS.csv"))%>% 
  dplyr::filter(Replicate != 2) |> 
  filter(SampleTypeCode == c("Water_Grab")) |> 
  select(SampleDate, StationCode, SampleTypeCode, AnalyteName, Result, 
         EventType, StationType) 

#add pathogens norog1, norog2 and adv. if ND (-88), convert to 0, if all pathogens are ND, flag as ND
df_1_pathogens_add<- data |> 
  pivot_wider(names_from = AnalyteName, values_from = Result) |> 
  select(-c("c.coli","c.jejuni","c.lari", "16scampy")) |> 
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
#write.csv(df_1_pathogens_add, file.path(data_path2,"RAS analysis", "cleaned_data_wideformat", "pathogen_data_wide_clean_2026-08-17.csv"))

#calculate log
df_2_indicators <- df_1_pathogens_add |> 
  select(-norog1, -norog2, -adv) |> 
  #ask Ryan: there are seven rows with all blank for indicators, but one ND under pathogen (see raw data pivot)- why? removed here
  filter(!is.na(crassphage)) |> 
  select(crassphage, pmmov, Enterococcus, HF183, path_detection, EventType, SampleDate, StationCode, pathogens_add) |> 
  pivot_longer((!path_detection & !EventType & !SampleDate & !StationCode),
               names_to = "analyte", 
               values_to = "result") |> 
  mutate(result = ifelse((result == -88 | result == 0) , 1, result)) |> 
  mutate(result_log = log10(result)) |> 
  select(-result)

#pivot wider
df_3_wide <- df_2_indicators |> 
  pivot_wider(names_from = analyte, values_from = result_log) |> 
  filter(!is.na(HF183)) |> #for predictive model, need vectors to be same length, so remove missing rows
  mutate(pathogen_binary = ifelse(path_detection == "ND", 0, 1))


#1) logistic regression:

#multiple variables
test1<- glm(pathogen_binary ~ crassphage + Enterococcus + HF183 + pmmov,
            data = df_3_wide, family = binomial)
summary(test1)
print(test1)
#p-values: show which indicators help predict detection
#coefficients (odds): stregnth of association

#2) evaluate predictive ability (if coefficients exist, test performance)
library(pROC)
roc_obj <- roc(df_3_wide$pathogen_binary, fitted(test1))
auc(roc_obj) #AUC = 0.7149
# AUC ~0.5 → no predictive ability
# AUC >0.7 → useful
# AUC >0.8 → strong

#compare to single variable
test2 <- glm(pathogen_binary ~ crassphage, family = binomial, data = df_3_wide)
roc_obj_2 <- roc(df_3_wide$pathogen_binary, fitted(test2))
auc(roc_obj_2)


#3) compare models (single vs multiple indicators)
AIC(test1, test2)

#If multi-indicator model improves AUC or lowers AIC → combination helps
#multiple variables: AUC is 0.714 and AIC is 59.3
#single variable:    AUC is 0.639 and AIC is 58.2

#4) test for any signal at all
#likelihood ratio test?:
anova(test1, test = "Chisq")

#5) machine learning check-----------------------------------------------------------

#“Do certain patterns of indicator bacteria reliably show when pathogens are present?”
#“When I look at lots of different rule combinations… do they consistently point to detection or not?”

# Random Forest = many decision trees + majority vote
# 
# Builds lots of decision trees from random subsets of data & variables
# Each tree predicts detect/non-detect
# Final prediction = average (or vote) across trees
# 
# 👉 Idea: reduces noise → improves prediction over a single tree

# Trees mostly agree → indicators are useful ✅
# Trees disagree or guess randomly →
# 👉 indicators don’t predict pathogens well ❌


#example decision tree 1: if HF183 > 5 -- detect
#                 tree 2: if ecoli > 30 -- detect
#                 tree 3: if pmov > 3 -- detect
## new sample: HF1=7, ecoli=20, pmov=2 (detect, non detect, nondetect)

#dataframe = df_3_wide

library(randomForest)
rf <- randomForest(as.factor(pathogen_binary) ~ crassphage + Enterococcus + HF183 + pmmov, data = df_3_wide, ntree = 500)
rf$confusion
tail(rf$err.rate)
#output from machine learning model
#[500,] 0.115942 0.06557377 0.5

#0.115942 = overall out-of-bag (OOB) error rate (11.6%)
#0.06557377 = class 1 error rate (6.6%)
#0.5 = class 2 error rate (50%)


# OOB error: 11.6, model is generally accurate (88.4%)
#pathogen binary classes: 
# class 1 is "0": non-detect
# Class 2 is "1": detect


# Random forest tests many combinations of variables
# If one indicator (e.g., HF1) is truly predictive:
#   Trees will mostly split on that variable
# Other variables get used less or ignored
# Key idea- should I include all indicators in the model?
# 👉 You do not need to pre-choose
# 👉 RF automatically finds the strongest signal

#so to summarize, the model is generally accurate, with the class 1 error rate (error rate of predicting non-detects) very low, but error predicting detects still high (50/50)
#you don't typically extract anything from this model, rather, you use the trained model object to predict outcomes
# for example:
# predict(rf, newdata = new data you are assessing)
# to see which variables matter most:
importance(rf)
varImpPlot(rf)

#MeanDecreaseGini
# crassphage           4.051503
# Enterococcus         1.719782
# HF183                2.155645
# pmmov                5.973012
# Higher MeanDecreaseGini = more important predictor.

# More specifically:
#   
#   It measures how much a variable reduces node impurity (Gini impurity) when splitting trees.
# Variables that create better splits across many trees get higher values.
# The values are relative, not probabilities or percentages.

#MeanDecreaseGini ranks variables by importance, but it does not tell you whether the relationship is positive or negative, nor does it provide a usable equation.

#save random forest model:
#saveRDS(rf, "rf_model.rds")
