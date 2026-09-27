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
  filter(StationCode == c("UPP_ELY_REC")) %>% 
  select(-SampleDate, -StationCode)

#apply log transform, then add binary code
df2_num <- as.matrix(df_2)

df_log10 <- log10(df2_num) %>%
  as.data.frame() %>% 
  mutate(pathogen_binary = ifelse(pathogens_add <= 0.6, 0, 1)) %>% 
  rename(ecoli = "E. coli") 

#for predictive model, need vectors to be same length, so remove missing rows
df_2<- df_log10 %>% 
  filter(!is.na(ecoli)) %>% 
  filter(!is.na(HF183))

df_test_nolog <- df_2 %>% 
  mutate(pathogen_binary = ifelse(pathogens_add <= 3.4, 0, 1)) %>% 
  rename(ecoli = "E. coli") %>% 
  filter(!is.na(ecoli)) %>% 
  filter(!is.na(HF183))

# #multivariate predictors, binary outcome
# install.packages("mgcv")
# library(mgcv)
#  y <- df_log10$pathogen_binary
#  x1 <- df_log10$crassphage
#  x2 <- df_log10$`E. coli`
#  x3 <- df_log10$HF183
#  x4 <- df_log10$pmmov
#  
# #y is binary response (ND, Detect)
# #x1, x2, x3 = indicators (crassphage, e. coli, hf183, pmmov)
# model <- gam(y ~ s(x1) + s(x2) + s(x3) + s(x4),
#              family = binomial (link = "logit"),
#              data = df_long, method = "REML")
# 
# summary(model)


#classification detection modeling

#1) logistic regression:

#multiple variables
test1<- glm(pathogen_binary ~ crassphage + ecoli + HF183 + pmmov,
            data = df_2, family = binomial)
summary(test1)
print(test1)
#p-values: show which indicators help predict detection
#coefficients (odds): stregnth of association

#2) evaluate predictive ability (if coefficients exist, test performance)
library(pROC)
roc_obj <- roc(df_2$pathogen_binary, fitted(test1))
auc(roc_obj) #AUC = 0.7149
# AUC ~0.5 → no predictive ability
# AUC >0.7 → useful
# AUC >0.8 → strong

#compare to single variable
test2 <- glm(pathogen_binary ~ crassphage, family = binomial, data = df_2)
roc_obj_2 <- roc(df_2$pathogen_binary, fitted(test2))
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

library(randomForest)
rf <- randomForest(as.factor(pathogen_binary) ~ crassphage + ecoli + HF183 + pmmov, data = df_2, ntree = 500)
rf$confusion
tail(rf$err.rate)
#output from machine learning model
# [500,] 0.4545455 0.2068966 0.9333333
# 0.4545 (overall)= 45% error (weak model, near random)
# 0.2069 (class 0)= means the model predicts non-detects fairly well (% error)
# 0.9333 (class 1)= predicts detects VERY poorly (%error)
# model is biased towards non-detects, almost always misses actual detects, suggests indicators don't predict detection well


# Random forest tests many combinations of variables
# If one indicator (e.g., HF1) is truly predictive:
#   Trees will mostly split on that variable
# Other variables get used less or ignored
# Key idea- should I include all indicators in the model?
# 👉 You do not need to pre-choose
# 👉 RF automatically finds the strongest signal


#try when not log transformed data
df_test_nolog

rf2 <- randomForest(as.factor(pathogen_binary) ~ crassphage + ecoli + HF183 + pmmov, data = df_test_nolog, ntree = 500)
rf2$confusion
tail(rf2$err.rate)