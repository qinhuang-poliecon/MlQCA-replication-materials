##====================-- I Set up Environment ----
##----
## 1.1 load packages ----
rm(list=ls())

pkgs <- c( "tidyverse", "pals","ggpubr", "ggsci", "xgboost",
           "Ckmeans.1d.dp", "devtools", "R.utils", 
           "ggh4x", "SHAPforxgboost","scales", "uwot", "caret",
           "QCA", "SetMethods", "venn", "cna", "haven")

lapply(pkgs, require, character.only = TRUE)

## 1.2 load auxiliary script ----
source("Scripts/Utility/General_Utility.R")
source("Scripts/Utility/cluster_diag_hotfix.R")

## 1.3 load data ----
data <- readRDS("Data/data_for_ml_Aug03_2022.rds") %>%
  mutate_all(as.numeric) %>%
  select(-V502) %>%
  filter(!is.na(vote)) %>%
  mutate(vote = abs(vote-1))

dataNoNA <- data %>% drop_na()

##----
##====================-- II Build Xgboost model ----
##----
## 2.1 set up ----
set.seed(201)

## 2.2 train model ----
xgb <- buildXGboost(dataNoNA, "vote")
getAccXGboost(xgb, dataNoNA, "vote")
fDf <- getFeatureImp(xgb, dataNoNA, "vote")

## 2.4 write the feature ranking to disk ----
fDf <- fDf %>% arrange(desc(shap))
write.csv(fDf, "Results/feature_ranking_table_nonVoting.csv", row.names = F)

## 2.5 extract cutoff from xgboost models ----
cutDfFull <- getXGboostCut(xgb, full = T)
cutDf <- getXGboostCut(xgb, full = F)
write.csv(cutDfFull, "Results/cutoff_table_nonVoting.csv", row.names = F)

## 2.6 save data for reproducibility ----
xgbList <- list(xgb=xgb, fDf=fDf, data = data, dataNoNA = dataNoNA,
                cutDf = cutDf, cutDfFull = cutDfFull)
saveRDS(xgbList, "Results/xgb_dataList_nonVoting.rds")

##----
##====================-- III Crisp Calibration and 10 Choose 4 pathway ranking ----
##----
## 3.1 reload data ----
xgbList <- readRDS("Results/xgb_dataList_nonVoting.rds")
list2env(xgbList, .GlobalEnv)

## 3.2 re-calibrate using top xgboost threshold ----
subDataRC <- reCalibrateXGCut(dataNoNA, cutDf = cutDf, "vote", na.rm = F)

### 3.3 Conduct necessity test ----
qDf <- getQCAMetric(subDataRC, "vote")
rDf <- full_join(fDf, qDf, by = "feature") %>%
  arrange(desc(shap))
rDf <- rDf %>%mutate(fRank = c(1:nrow(rDf)))
write.csv(rDf, "Results/Crisp_FeatureTable_nonVoting.csv")

## 3.4 Do Sufficiency Test: Every interval of 10 choose random 4 ----
Nindex = 10
kterm = 4
stepSize = 10
Nfeature = ncol(data)-1
indList <- lapply(seq(1, (Nfeature-Nindex+1), stepSize), function(ii) c(ii:(ii+Nindex-1)))

QCAList <- lapply(indList, function(indVec){
  combnQCA(Nvec = indVec, k = kterm, dfInput = subDataRC, Lab = "vote", 
           fImpDf = fDf, incl.cut = 0.8)
})

## 3.5 clean up results ----
gpNames <- unlist(lapply(indList, function(ivec){
  paste0("Top ", ivec[1], "-", ivec[length(ivec)])
}))

dfx <- bind_rows(Map(function(gpx, qcaDfx){
  qcaDfx$group <- gpx
  qcaDfx
}, gpNames, QCAList))

dfx$group <- factor(dfx$group, levels = unique(dfx$group))

write.csv(dfx, "Results/iterativeQCA_Crisp_Every10Choose4_nonVoting.csv", row.names = F)


##----
##====================-- IV Fuzzy Calibration and 10 Choose 4 pathway ranking ----
##----
## 4.1 reload data ----
xgbList <- readRDS("Results/xgb_dataList_nonVoting.rds")
list2env(xgbList, .GlobalEnv)

## 4.2 re-calibrate using top xgboost threshold ----
subDataRC <- reCalibrateXGCutFuzzy(dfInput = dataNoNA, cutDf = cutDfFull, Lab = "vote", method = "direct")

## 4.3 Conduct necessity test ----
qDf <- getQCAMetric(subDataRC, "vote")
rDf <- full_join(fDf, qDf, by = "feature") %>%
  arrange(desc(shap))
rDf <- rDf %>%mutate(fRank = c(1:nrow(rDf)))
write.csv(rDf, "Results/Fuzzy_FeatureTable_nonVoting.csv")

## 4.4 Do Sufficiency Test: Every interval of 10 choose random 4 ----
Nindex = 10
kterm = 4
stepSize = 10
Nfeature = ncol(data)-1
indList <- lapply(seq(1, (Nfeature-Nindex+1), stepSize), function(ii) c(ii:(ii+Nindex-1)))

QCAList <- lapply(indList, function(indVec){
  combnQCA(Nvec = indVec, k = kterm, dfInput = subDataRC, Lab = "vote", 
           fImpDf = fDf, incl.cut = 0.8)
})

## 4.5 clean up results ----
gpNames <- unlist(lapply(indList, function(ivec){
  paste0("Top ", ivec[1], "-", ivec[length(ivec)])
}))

dfx <- bind_rows(Map(function(gpx, qcaDfx){
  qcaDfx$group <- gpx
  qcaDfx
}, gpNames, QCAList))

dfx$group <- factor(dfx$group, levels = unique(dfx$group))

write.csv(dfx, "Results/iterativeQCAFuzzy_Every10Choose4_nonVoting.csv", row.names = F)

