### this script provides some custom helper functions, many of which were formalized into our mlQCA package later


##============================-- pure QCA utility  --=======================================##
## get QCA metrics
getQCAMetric <- function(dfx, labV){
  require(QCA)
  # labVec <- as.numeric(dfx[[labV]])
  # fDf <- dfx %>% select(-all_of(labV))
  # fDf <- as.data.frame(apply(fDf, 2, as.numeric))
  features <- setdiff(colnames(dfx), labV)
  metricDf <- do.call("rbind", lapply(features, function(fx) {
    subDfx <- dfx %>% select(all_of(c(fx, labV))) %>% drop_na()
    labVec <- as.numeric(subDfx[[labV]])
    valVec <- as.numeric(subDfx[[fx]])
    necDf <- as.data.frame( QCAfit(valVec, labVec, cond.lab = fx) )
    necDf$feature <- fx
    necDf
  }))
  return(metricDf)
}

## Annotate Label
AnnoLabel <- function(dfx, labV, posV, negV, labName=NULL){
  dfx$Label <- NA
  dfx$Label[dfx[[labV]] %in% posV] <- 0
  dfx$Label[dfx[[labV]] %in% negV] <- 1
  dfx[[labV]] <- NULL
  dfx <- dfx %>% filter(!is.na(Label))
  if(!is.null(labName)){
    colnames(dfx) <- gsub("Label", labName, colnames(dfx))
  }
  return(dfx)
}


## get QCA coverage for multiple conditions
getQCAIncl <- function(dfx, labV, condV){
  xx <- truthTable(data = dfx,
                      outcome=labV,
                      conditions = condV,
                      complete = TRUE) #shows all rows in the TT
  return(xx)
}

## generate index
getIndexSet <- function(inputVec, n){
  outList <- combn(inputVec, n, simplify = F)
  names(outList) <- paste0("index_", 1:length(outList))
  return(outList)
}

##============================-- xgboost + QCA utility  --=======================================##
## build XGBoost Model
buildXGboost <- function(dfx, Lab){
  dfx <- dfx %>% rename(Label = all_of(Lab))
  dfx$Label <- factor(dfx$Label)
  ## and all samples in every tree
  custGrid <- expand.grid(eta=c(0.2:0.4),
                          gamma=0,
                          max_depth=c(3:5),
                          colsample_bytree=c(0.6,0.8, 1),
                          min_child_weight=1,
                          subsample=c(0.75, 0.85,1),
                          nrounds=c(50,100,150))

  xgb = train(
    Label ~ .,
    data = dfx,
    method = "xgbTree",
    trControl = trainControl(method = "cv", number = 5),
    na.action = na.pass,
    tuneGrid = custGrid,
    verbose = 0
  )

  return(xgb)
}

## measure accuracy of model
getAccXGboost <- function(model, dfx, Lab){
  tst_lbl <- dfx[[Lab]]
  xgb_pred <- predict(model, newdata = dfx, na.action = na.pass)
  acc <- mean(tst_lbl==xgb_pred) ## 0.768 for Male; 0.760 for Female
  return(acc)
}

## extract feature importance of model
getFeatureImp <- function(model, dfx, Lab){
  bst <- model$finalModel
  impDf <- as.data.frame(xgb.importance(model = bst)) %>%
    rename(feature = Feature)

  ## 2.4 also add shap value
  cols <- setdiff(colnames(dfx), Lab)
  shapScore <- shap.values(bst, X=data.matrix(dfx[, cols]))$mean_shap_score
  shapDf <- data.frame(feature=names(shapScore), shap=shapScore)

  ## 2.5 combine with qca metric
  # qcaDf <- getQCAMetric(dfx, Lab)

  combDf <- impDf %>%
    full_join(shapDf, by="feature")

  return(combDf)
}

## extract threshold being used in model
getXGboostCut <- function(model, full=F, tIdx=NULL){
  bst <- model$finalModel
  Ntree <- bst$tuneValue$nrounds
  if (is.null(tIdx)){
    xx <- xgb.plot.tree(model=bst, trees = 0:(Ntree-1), render = F)
  } else {
    xx <- xgb.plot.tree(model=bst, trees = tIdx, render = F)
  }
  Edf <- xx$edges_df %>%
    filter(label!="") %>%
    select(from, label)
  ndf <- xx$nodes_df %>%
    select(id, data, label) %>%
    rename(metric=label)
  treeDf <- inner_join(Edf, ndf, by=c("from"="id")) %>%
    rename(feature=data) %>%
    select(-from) %>%
    mutate(cutoff = str_split(label, "<", simplify = T)[,2])
  cutDf <- treeDf %>%
    group_by(feature) %>%
    count(cutoff) %>%
    mutate(Prop = n/sum(n)) %>%
    ungroup() %>%
    arrange(feature, desc(Prop)) %>%
    mutate(cutoff = round(as.numeric(cutoff), 3))
  if (full) {
    return(cutDf)
  } else {
   finalDf <- cutDf %>%
     group_by(feature) %>%
     slice_head(n=1) %>%
     select(-n, -Prop) %>%
     ungroup() %>%
     mutate(cutoff = as.numeric(cutoff))
   return(finalDf)
  }
}

## do calibration with xgboost cut: crisp
reCalibrateXGCut <- function(dfInput, cutDf, Lab, na.rm=F){
  dfx <- dfInput %>%
    select(all_of(c(cutDf$feature, Lab)))

  dfx_Cali <- bind_cols(map2(cutDf$feature, cutDf$cutoff, function(fx, cx){
    vx <- dfx[[fx]]
    vx <- ifelse(vx>=cx, 1, 0)
    dfy <- data.frame(vx)
    colnames(dfy) <- fx
    dfy
  }))

  dfx_Cali[[Lab]] <- dfx[[Lab]]

  if (na.rm){
    dfx_Cali <- drop_na(dfx_Cali)
  }

  return(dfx_Cali)
}

## do calibration with xgboost cut: fuzzy
reCalibrateXGCutFuzzy <- function(dfInput, cutDf, Lab, method="direct", na.rm=F){
  ## check fuzziness
  if(!any(duplicated(cutDf$feature))){
    stop("not a full cutDf")
  }

  ## extract fuzzy and crisp features
  sumCut <- cutDf %>%
    count(feature) %>%
    arrange(desc(n))
  fCrisp <- sumCut %>% filter(n<=2) %>% pull(feature)
  fFuzzy <- sumCut %>% filter(n>2) %>% pull(feature)

  ## binarize crisp features
  dfxCrisp <- dfInput %>%
    select(all_of(fCrisp))
  cutDfCrisp <- cutDf %>%
    filter(feature %in% fCrisp) %>%
    group_by(feature) %>%
    arrange(desc(Prop)) %>%
    slice_head(n=1) %>%
    select(-n, -Prop) %>%
    ungroup() %>%
    mutate(cutoff = as.numeric(cutoff))

  ## calibrate crisp features
  dfxCrisp_Cali <- bind_cols(map2(cutDfCrisp$feature, cutDfCrisp$cutoff, function(fx, cx){
    vx <- dfxCrisp[[fx]]
    vx <- ifelse(vx>=cx, 1, 0)
    dfy <- data.frame(vx)
    colnames(dfy) <- fx
    dfy
  }))

  ## generate list of cutoffs for fuzzy features
  dfxFuzzy <- dfInput %>%
    select(all_of(fFuzzy))
  cutFuzzyDf <- cutDf %>% filter(feature %in% fFuzzy)
  cutFuzzyList <- lapply(split(cutFuzzyDf, f=cutFuzzyDf$feature), function(dd){
    dd %>%
      arrange(desc(Prop)) %>%
      slice_head(n=3) %>%
      mutate(cutoff = as.numeric(cutoff)) %>%
      pull(cutoff) %>%
      sort()
  })

  ## calibrate fuzzy features
  if (method=="direct" | method=="indirect"){
    dfxFuzzy_Cali <- bind_cols(lapply(names(cutFuzzyList), function(fx){
      vx <- dfxFuzzy[[fx]]
      cx <- cutFuzzyList[[fx]]
      vx <- calibrate(vx, thresholds = cx, type = "fuzzy", method = method)
      dfy <- data.frame(vx)
      colnames(dfy) <- fx
      dfy
    }))
  } else {
    stop("method must be either direct or indirect")
  }

  dfx_Cali <- bind_cols(dfxCrisp_Cali, dfxFuzzy_Cali)
  dfx_Cali[[Lab]] <- dfInput[[Lab]] ## label should have been binarized

  if (na.rm){
    dfx_Cali <- drop_na(dfx_Cali)
  }

  return(dfx_Cali)
}


## convert treedf to truth table
treeToTruthTable <- function(tDf){
  tDf <- tdfTreeFinal
  colnames(tDf) <- gsub("Level", "Feature", colnames(tDf))
  fN <- sum(grepl("Feature", names(tDf)))
  tDf$varList <- "NA"
  for (ii in 1:nrow(tDf)){
    ci <- tDf$cumCond[ii]
    if (grepl("NA", ci)){
      next()
    }
    vecI <- strsplit(ci, ";")[[1]]
    varI <- str_split(vecI, " ", simplify = T)[,1]
    signI <- str_split(vecI, " ", simplify = T)[,2]
    valueI <- plyr::mapvalues(signI, from = c(">", "<"), to = c(1, 0),
                              warn_missing = F)
    tDf <- as.data.frame(tDf)
    tDf$varList[ii] <- paste(varI, collapse = ";")
    for (jj in 1:fN){
      vx <- tDf[ii,jj]
      if(vx %in% varI){
        tDf[ii,jj] <- plyr::mapvalues(vx, from = varI, to = valueI, warn_missing = F)
      } else {
        tDf[ii,jj] <- NA
      }
    }
  }
  tDfOut <- tDf %>%
    filter(Feature_1 %in% c("0", "1")) %>%
    mutate_at(vars(starts_with("Feature")), as.numeric)
  return(tDfOut)
}

## sufficiency test based on a few selected variables
subsetSufficiency <- function(dfInput, rankVec, Lab, fImpDf=fDf, incl.cut=0.9, features=NULL){
  fDfsort <- fImpDf %>% arrange(desc(shap))
  if (is.null(features)){
    vSelect <- fDfsort$feature[rankVec]
  } else {
    vSelect <- features
  }
  subdata <- dfInput[,c(Lab, vSelect)] %>% drop_na()

  tTablePass1 <- truthTable(data = subdata,
                            outcome= Lab,
                            conditions = vSelect,
                            incl.cut = incl.cut, #consistency cut
                            pri.cut = 0.51, #PRI cut
                            n.cut = 10, #cases in a row before is classified as logical remainder (rows without enough empirical evidence)
                            sort.by = c("OUT", "incl"), #sort by output and consistency
                            complete = TRUE) #shows all rows in the TT


  tTableMin <- minimize(input =  tTablePass1,
                     use.tilde = TRUE,
                     details = TRUE)
  tTableList <- list(table=tTablePass1, path=tTableMin)

  return(tTableList)
}

## summarize sufficiency Results
sumSuffList <- function(idxList, resList){
  resDf <- bind_rows(Map(function(idx, resl){
    idset <- paste0(idx, collapse = ",")

    if (length(resl)==0){
      data.frame(index = idset, path = NA,
                 inclS = 0,
                 PRI = 0,
                 convS = 0)
    } else {
      resPath <- paste0(resl$path$essential, collapse = "+")
      resInc <- resl$path$IC

      if(length(resInc$sol.incl.cov)>0){
        incVec <- resInc$sol.incl.cov
      } else {
        incVec <- resInc$overall$sol.incl.cov
      }

      data.frame(index = idset, path = resPath,
                 inclS = incVec[1, 'inclS'],
                 PRI = incVec[1, 'PRI'],
                 convS = incVec[1, "covS"])
    }

  }, idxList, resList))

  resDfSorted <- resDf %>%
    arrange(desc(convS)) %>%
    mutate(rank = 1:n())

  return(resDfSorted)
}

## Combination run through
combnQCA <- function(Nvec, k, dfInput, Lab, fImpDf, incl.cut=0.9){
  idx_nk <- getIndexSet(Nvec, k)
  res_nk <- lapply(idx_nk, function(xv){
    tryCatch({
      subsetSufficiency(dfInput = dfInput, rankVec = xv,
                        Lab=Lab, incl.cut = incl.cut, fImpDf = fImpDf)
    }, error = function(e){
      return(NULL)
    })
  })
  df_nk <- sumSuffList(idx_nk, res_nk)
  return(df_nk)
}

## mean Consistency Distance test based on a few selected variables
subsetConsisDist <- function(dfInput, rankVec, Lab, fImpDf=fDf, features=NULL){
  fDfsort <- fImpDf %>% arrange(desc(shap))
  if (is.null(features)){
    vSelect <- fDfsort$feature[rankVec]
  } else {
    vSelect <- features
  }
  subdata <- dfInput[,c(Lab, vSelect)] %>% drop_na()

  tTablePass1 <- truthTable(data = subdata,
                            outcome= Lab,
                            conditions = vSelect,
                            incl.cut = 0, #consistency cut
                            pri.cut = 0, #PRI cut
                            n.cut = 5, #cases in a row before is classified as logical remainder (rows without enough empirical evidence)
                            sort.by = c("OUT", "incl"), #sort by output and consistency
                            complete = TRUE) #shows all rows in the TT


  minConsist <- tTablePass1$tt %>%
    select(-cases) %>%
    filter(OUT!="?") %>%
    mutate(incl = as.numeric(incl)) %>%
    mutate(minDistConsRaw = pmin(incl, 1-incl, na.rm = T)) %>%
    mutate(minDistCons = minDistConsRaw * n /sum(n)) %>%
    pull(minDistCons) %>%
    #mean(na.rm=T)
    sum(na.rm = T)

  return(minConsist)
}

## summarize sufficiency Results
sumConDList <- function(idxList, resList){
  resDf <- bind_rows(Map(function(idx, resl){
    idset <- paste0(idx, collapse = ",")

    if (length(resl)==0){
      data.frame(index = idset, meanConsDist = 1)
    } else {
      data.frame(index = idset, meanConsDist = resl)
    }
  }, idxList, resList))

  resDfSorted <- resDf %>%
    arrange(meanConsDist) %>%
    mutate(rank = 1:n())

  return(resDfSorted)
}

## Combination run through
combnQCAConsD <- function(Nvec, k, dfInput, Lab, fImpDf){
  idx_nk <- getIndexSet(Nvec, k)
  res_nk <- lapply(idx_nk, function(xv){
    tryCatch({
      subsetConsisDist(dfInput = dfInput, rankVec = xv,
                        Lab=Lab, fImpDf = fImpDf)
    }, error = function(e){
      return(NULL)
    })
  })
  df_nk <- sumConDList(idx_nk, res_nk)
  return(df_nk)
}


## functions to get all possible thresholds
getFthresh <- function(dataDf, fname, NcutMax=Inf){
  require(zoo)
  fVec <- sort(unique(dataDf[[fname]]))
  fVec <- fVec[!is.na(fVec)]
  N <- length(fVec)
  if (N == 2){
    return(mean(fVec))
  } else if (N > NcutMax+1) {
    fSub <- fVec[seq(1, N, length.out = NcutMax+1)]
    fCut <- rollapply(fSub, 2, mean)
    return(fCut)
  } else {
    fCut <- rollapply(fVec, 2, mean)
    return(fCut)
  }
}
crossThresh <- function(dataDf, fVec, NcutMax=Inf){
  threshList <- lapply(fVec, function(fx) getFthresh(dataDf = dataDf, fname = fx, NcutMax = NcutMax))
  names(threshList) <- fVec
  threshDf <- cross_df(threshList)
  return(threshDf)
}


##============================-- ggplot utility  --=======================================##
saveFig <- function(pobj, prefix, wd, ht){
  png(filename = paste0(prefix, ".png"), width = wd, height = ht, res = 300, units = "px")
  print(pobj)
  dev.off()

  ggsave(plot = pobj, filename = paste0(prefix, ".pdf"), width = wd, height = ht,
         units = "px")
}
