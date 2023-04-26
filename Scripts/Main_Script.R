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
  filter(!is.na(vote))

dataNoNA <- data %>% drop_na()

##----
##====================-- II Build XGboost Model ----
##----
## 2.1 set up ----
set.seed(201)

## 2.2 train model ----
xgb <- buildXGboost(dataNoNA, "vote")  # winning model #68 for Appendix 2
getAccXGboost(xgb, dataNoNA, "vote")
fDf <- getFeatureImp(xgb, dataNoNA, "vote")

## 2.4 write the feature ranking to disk ---- # Condition ranking for Table 1

fDf <- fDf %>% arrange(desc(shap))
write.csv(fDf, "Results/feature_ranking_table.csv", row.names = F)  
## 2.5 extract cutoff from xgboost models ---- # Cutoffs for Table 1
cutDfFull <- getXGboostCut(xgb, full = T)
cutDf <- getXGboostCut(xgb, full = F)
write.csv(cutDfFull, "Results/cutoff_table.csv", row.names = F)  

## 2.6 save data for reproducibility ----
xgbList <- list(xgb=xgb, fDf=fDf, data = data, dataNoNA = dataNoNA,
                cutDf = cutDf, cutDfFull = cutDfFull)
saveRDS(xgbList, "Results/xgb_dataList.rds")

##----
##====================-- III Crisp Calibration and 10 Choose 4 Model Ranking ----
##----
## 3.1 reload data ----
xgbList <- readRDS("Results/xgb_dataList.rds")
list2env(xgbList, .GlobalEnv)

## 3.2 re-calibrate using top xgboost threshold ----
subDataRC <- reCalibrateXGCut(dataNoNA, cutDf = cutDf, "vote", na.rm = F)

### 3.3 Conduct necessity test ---- # Necessity analysis for Appendix 1
qDf <- getQCAMetric(subDataRC, "vote")
rDf <- full_join(fDf, qDf, by = "feature") %>%
  arrange(desc(shap))
rDf <- rDf %>%mutate(fRank = c(1:nrow(rDf)))
write.csv(rDf, "Results/Crisp_FeatureTable.csv")  


## 3.4.1 Do Sufficiency Test: Every interval of 10 choose random 4 ----
Nindex = 10
kterm = 4
stepSize = 10
Nfeature = ncol(data)-1
indList <- lapply(seq(1, (Nfeature-Nindex+1), stepSize), function(ii) c(ii:(ii+Nindex-1)))

QCAList <- lapply(indList, function(indVec){
  combnQCA(Nvec = indVec, k = kterm, dfInput = subDataRC, Lab = "vote", 
           fImpDf = fDf, incl.cut = 0.8)
})

## 3.4.2 clean up results ----
gpNames <- unlist(lapply(indList, function(ivec){
  paste0("Top ", ivec[1], "-", ivec[length(ivec)])
}))

dfx <- bind_rows(Map(function(gpx, qcaDfx){
  qcaDfx$group <- gpx
  qcaDfx
}, gpNames, QCAList))

dfx$group <- factor(dfx$group, levels = unique(dfx$group))

write.csv(dfx, "Results/RadicalQCA.csv", row.names = F)    ## Results for Radical MlQCA in Table 2
saveRDS(dfx, "Results/Combine_iterative_QCA_TopNK_10chosse4.rds")

## 3.5 Plot results ---- #  Figure 2 Left Panel
dfx <- readRDS("Results/Combine_iterative_QCA_TopNK_10chosse4.rds")
### 3.5.1 coverage distribution
px <- ggplot(dfx, aes(x=group, y=convS, fill=group)) +
  geom_violin(trim=T) +
  geom_boxplot(width=0.1, fill="white") +
  theme_pubr() +
  labs(y="Coverage", fill="") +
  scale_fill_brewer(palette = "Blues", direction = -1) +
  theme(axis.text.x = element_text(angle = 60, hjust=1, size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=18),
        legend.text = element_text(size=15)) +
  guides(fill=guide_legend(override.aes = list(size=8)))
saveFig(px, prefix = "Results/Plot_Coverage_distribution_TopNK_10choose4", 
        wd = 2400, ht = 2400)              


### 3.5.2 distribution among top coverages candidates ----  # Figure 2 Right Panel

## by coverage rank
dfy <- dfx %>% 
  arrange(desc(convS)) %>%
  mutate(convS = 1:nrow(dfx)) %>%
  mutate(cGroup = case_when(
    convS > 0 & convS <= 10 ~ "Cov1-10",
    convS > 10 & convS <= 100 ~ "Cov11-100",
    convS > 100 & convS <= 500 ~ "Cov101-300",
    convS > 500 & convS <= 1000 ~ "Cov501-1000",
    convS >1000 ~"Cov>1000"
  )) %>%
  group_by(cGroup) %>%
  count(group) %>%
  mutate(frequency = n/sum(n)) %>%
  mutate(cGroup = factor(cGroup, 
                         levels = c("Cov1-10", "Cov11-100", 
                                    "Cov101-300", "Cov501-1000",
                                    "Cov>1000")))

p1b <- ggplot(dfy, aes(x=cGroup, y=frequency, fill=group)) +
  geom_bar(stat = "identity") +
  theme_pubr() +
  labs(y="Frequency", fill="") +
  scale_fill_brewer(palette = "Blues", direction = -1) +
  theme(axis.text.x = element_text(angle = 60, hjust=1, size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=18),
        legend.text = element_text(size=15)) +
  guides(fill=guide_legend(override.aes = list(size=8)))

saveFig(p1b, prefix = "Results/Plot_ShapGroup_distribution_Rank_byCoverage", 
        wd = 2400, ht = 2400)   


##----
##====================-- IV Fuzzy Calibration and 10 Choose 4 Model Ranking ----
##----
## 4.1 reload data ----
xgbList <- readRDS("Results/xgb_dataList.rds")
list2env(xgbList, .GlobalEnv)

## 4.2 re-calibrate using top xgboost threshold ----
subDataRC <- reCalibrateXGCutFuzzy(dfInput = dataNoNA, cutDf = cutDfFull, Lab = "vote", method = "direct")

## 4.3 Conduct necessity test ----
qDf <- getQCAMetric(subDataRC, "vote")
rDf <- full_join(fDf, qDf, by = "feature") %>%
  arrange(desc(shap))
rDf <- rDf %>%mutate(fRank = c(1:nrow(rDf)))
write.csv(rDf, "Results/Fuzzy_FeatureTable.csv")

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

write.csv(dfx, "Results/iterativeQCA_Fuzzy_Every10Choose4.csv", row.names = F)

##----
##====================-- V Conservative mlQCA ----  # Table 3
##----
## 5.1 reload data ----
xgbList <- readRDS("Results/xgb_dataList.rds")
list2env(xgbList, .GlobalEnv)

## 5.2 Recalibrate using Crisp methods  ----
crispRC <- reCalibrateXGCut(dataNoNA, cutDf = cutDf, "vote")


## 5.3 Conservative QCA----

### select variables
v1 <- "V141"
v2 <- c("V208", "P_V195", "P_V258")
vDf3 <- expand.grid(list(v1, v2), stringsAsFactors = F)
vDf3$Var3 <- "V648"
vDf3$Var4 <- "V380"
vDf <- vDf3

### conduct QCA
crispQCAList <- list()


crispQCAList <- lapply(1:nrow(vDf), function(ii){
  vDfRowi <- as.character(as.vector(vDf[ii,]))
  tryCatch({
    subsetSufficiency(dfInput=crispRC, features = vDfRowi, 
                      Lab="vote", fImpDf=fDf, incl.cut=0.8)
  }, error = function(e){
    return(NULL)
  })
})



### clean up results  
idxList <- lapply(1:nrow(vDf), function(ii) as.character(ii))
crispDf <- sumSuffList(idxList = idxList, resList = crispQCAList)
vDf$index <- as.character(1:nrow(vDf))
crispRes <- full_join(vDf, crispDf, by="index")
write.csv(crispRes, "Results/ConservativeQCA.csv", 
          quote = T, row.names = F)



##----
##====================-- VI Consistency Optimization Tests ---- # Appendix 5
##----
## 6.1 reload data ----
xgbList <- readRDS("Results/xgb_dataList.rds")
list2env(xgbList, .GlobalEnv)

## 6.2 re-calibrate using top xgboost threshold ----
subDataRC <- reCalibrateXGCut(dataNoNA, cutDf = cutDf, "vote", na.rm = F)

## 6.3 Do Consistency Distance Calculation: Every interval of 10 choose random 4 ----
Nindex = 10
kterm = 4
stepSize = 10
Nfeature = ncol(data)-1
indList <- lapply(seq(1, (Nfeature-Nindex+1), stepSize), function(ii) c(ii:(ii+Nindex-1)))

## iterate through
QCAList <- lapply(indList, function(indVec){
  combnQCAConsD(Nvec = indVec, k = kterm, dfInput = subDataRC, Lab = "vote", 
                fImpDf = fDf)
})

## 6.4 clean up results ----
###  names of groups
gpNames <- unlist(lapply(indList, function(ivec){
  paste0("Top ", ivec[1], "-", ivec[length(ivec)])
}))

### combine results
dfx <- bind_rows(Map(function(gpx, qcaDfx){
  qcaDfx$group <- gpx
  qcaDfx
}, gpNames, QCAList))

dfx$group <- factor(dfx$group, levels = unique(dfx$group))
saveRDS(dfx, "Results/iterativeQCA_Consistency_Distance_TopNK_10chosse4.rds")


## 6.5 Plot results ---- 
dfx <- readRDS("Results/iterativeQCA_Consistency_Distance_TopNK_10chosse4.rds")
### 3.5.1 coverage distribution
px <- ggplot(dfx, aes(x=group, y=meanConsDist, fill=group)) +
  geom_boxplot(width=0.5, color="#D96704") +
  theme_pubr() +
  labs(y="mean_Consis_Res", fill="") +
  scale_fill_brewer(palette = "Blues", direction = -1) +
  theme(axis.text.x = element_text(angle = 60, hjust=1, size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=18),
        legend.text = element_text(size=15)) +
  guides(fill=guide_legend(override.aes = list(size=8)))
saveFig(px, prefix = "Results/Plot_Consistency_distribution_TopNK_10choose4", 
        wd = 2400, ht = 2400)

##----
##====================-- VII Outcome-variable permutation test ---- # Appendix 6
##----
## 7.0 load previous models ----
xgbList <- readRDS("Results/xgb_dataList.rds")
list2env(xgbList, .GlobalEnv)
Top10Var <- head(fDf$feature,10)
fDf <- fDf %>% mutate(rank = 1:n())

## 7.1 utility function to retrain model ----
redoXGboost <- function(ShufDf){
  xgb <- buildXGboost(ShufDf, "vote")
  getAccXGboost(xgb, ShufDf, "vote")
  fDf <- getFeatureImp(xgb, ShufDf, "vote")
  fDf <- fDf %>% 
    arrange(desc(shap)) %>%
    mutate(rank = 1:n())
  return(fDf)
}


## 7.2 randomize one variable at a time ----
fDfCompareList <- list()
for (ii in 1:10){
  set.seed(ii)
  
  fDfCompare <- do.call("rbind", lapply(Top10Var, function(vx){
    dataPerm <- data %>%
      mutate_at(vx, sample)
    
    fDfPerm <- redoXGboost(dataPerm)
    
    fDfMerge <- fDf %>%
      full_join(fDfPerm, by="feature", suffix = c("_orig", '_perm')) %>%
      filter(feature == vx)
    
    fDfMerge
  }))
  
  fDfCompare <- fDfCompare %>% 
    mutate(Trial = ii)
  
  fDfCompareList[[ii]] <- fDfCompare
}
saveRDS(fDfCompareList, "Results/Permute_Trials_OneAtATime.rds")

## 7.3 combine all trials ----
fDfFinal <- do.call("rbind", fDfCompareList) %>%
  arrange(rank_orig, Trial)
write.csv(fDfFinal, "Results/Full_Results_Mutate_OneAtTime_10Trials.csv", row.names = F, quote = F)


## 7.4 plot results ----
dfx <- fDfFinal %>% 
  select(feature, rank_orig, rank_perm, Trial) %>%
  pivot_longer(cols = c(rank_orig, rank_perm), names_to = "group", values_to = "rank") %>%
  mutate(group = gsub("rank_", "",group)) %>%
  mutate(feature = factor(feature, levels = Top10Var)) %>%
  
  
  px <- ggplot(dfx, aes(x=group, y=rank, fill=group)) +
  geom_point() +
  geom_boxplot() +
  facet_wrap2(~feature, nrow = 2) +
  theme_bw() +
  labs(y="Rank") +
  theme(axis.text.x = element_text(size = 12, color="black"),
        axis.text.y = element_text(size=10),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=15),
        legend.text = element_text(size=12),
        legend.title = element_blank(),
        strip.background = element_rect(fill = "transparent", colour = "transparent"),
        strip.text = element_text(size=12, face = "bold"))
saveFig(px, prefix = "Results/Plot_Change_in_Rank_Mutate_OneAtTime_10Trials", wd = 2800, ht = 2400)


