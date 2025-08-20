
#To plot a histgram of correlation of predicted & true metabolome
#read filename
#show and save histgram to Out 

histfunction <- function(filename, data) {
  x <- read.csv(paste0("../out/", filename), row.names = 1)
  subfilename <- substr(filename, 1, nchar(filename) - 4)
  correlation <- cor(x, data)
  correlation_diag <- diag(correlation)
  
  # Create a histogram plot
  hist(correlation_diag, xlim = c(-1, 1), ylim = c(0, 140),
       breaks = seq(-1, 1, length.out = 21), xlab = "Correlation",
       main = subfilename)
  
  # Save the plot as a PNG file
  png(filename = paste("../out/hist", subfilename, ".png", sep = ""),
      width = 640, height = 420, units = "px")
  
  dev.off()
}

getBLUP <- function(ZETAlist,y,Tv=NULL) {
  
  names<- names(ZETAlist)
  resultModelname <- names(ZETAlist)
  
  iterations=length(ZETAlist)
  cl <- makeCluster(4)
  registerDoSNOW(cl)
  pb <- txtProgressBar(max = iterations, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  xr<- foreach (zn = 1:iterations,.options.snow = opts,.export = c("EM3.cpp")) %dopar% {
    
    #xr <- foreach(zn = seq_along(ZETAlist)) %dopar% {
    
    zname <- names(ZETAlist)[zn]
    print(zname)
    
    var.pheno.all = c() #traitごと結果(要約)
    
    ZETA <-ZETAlist[[zn]]
    
    X0 <- NULL
    
    for (i in 2:ncol(y)){
      #for (i in 2){  
      trait.name <- colnames(y)[i]
      print(c(i,"/",ncol(y)))
      print(trait.name)
      # make the case for the result
      
      resultIndex <- c("Correlation", "R2", "RMSE","varU","varE","h2")
      # resultEachSeed <- matrix(NA, ncol = length(resultIndex), nrow = length(seedInd))
      # colnames(resultEachSeed) <- resultIndex
      # rownames(resultEachSeed) <- seedInd
      
      predictionDataRAINBOW <- rep(NA, nrow(y))
      
      for (t in 1:nrow(y)) { 
        #for (t in 1:3) {
        
        tf<- !(1:nrow(y) == t)
        if(!is.null(Tv)){
          ZETA[[1]][["K"]][tf,] <- Tv[tf,]
        }
        
        yNa <- y[,i]
        yNa[t] <- NA
        # RAINBOWR
        
        resEM3 <- EM3.cpp(y = yNa,
                          X0 = X0, n.core = 1, 
                          ZETA = ZETA)
        
        # input the predicted data
        predictionDataRAINBOW[t] <- resEM3$y.pred[t]
      }
      
      # resEM3$weights
      predictData <- (predictionDataRAINBOW)
      obsData <- (y[,i])
      
      # calculate the R2 and RMSE
      correlation <- cor(obsData, predictData)
      R2 <- 1 - sum((obsData - predictData) ^ 2) / sum((obsData - mean(obsData)) ^ 2)
      RMSE <- sqrt(sum((obsData - predictData) ^ 2) / length(obsData))
      VU <- resEM3$Vu
      VE <- resEM3$Ve
      h2 <-  VU/(VU+VE)
      
      resultsum<- c(correlation, R2, RMSE,VU,VE,h2)
      resultsum<- as.matrix(resultsum)
      colnames(resultsum) <-trait.name
      var.pheno.all <- cbind(var.pheno.all,resultsum)
    }
    
    rownames(var.pheno.all) <- resultIndex
    print(var.pheno.all)
    var.pheno.all
  }
  
  close(pb)
  stopCluster(cl) 
  
  resultModelname <- names(ZETAlist)
  resultEachModel <- matrix(NA, nrow = length(ZETAlist), ncol = ncol(y)-1)
  rownames(resultEachModel) <- names(ZETAlist)
  colnames(resultEachModel) <- colnames(y)[2:ncol(y)]
  
  #index p generator
  for (n in 1:length(names)){
    
    resultEachModel[n,]<- xr[[n]][1,]
    
    if (n == 1){
      p <- names[1]
    }
    else{
      p  <- paste(p,names[n],sep="")
    }
  }
  print(p)
  
  
  names(xr) <- names
  resultEachModel_all<- xr
  
  return(list(resultEachModel=resultEachModel,resultEachModel_all=resultEachModel_all,p=p))
  
}

getBLUP_simple <- function(ZETAlist, y, Tv = NULL) {
	
	names <- names(ZETAlist)
	resultModelname <- names(ZETAlist)
	
	iterations <- length(ZETAlist)
	
	xr <- list()  # 結果を格納するリスト
	
	for (zn in 1:iterations) {
		
		zname <- names(ZETAlist)[zn]
		print(zname)
		
		var.pheno.all <- c()  # 各traitごとの結果を格納
		
		ZETA <- ZETAlist[[zn]]
		
		X0 <- NULL
		
		for (i in 2:ncol(y)) {
			trait.name <- colnames(y)[i]
			print(c(i, "/", ncol(y)))
			print(trait.name)
			
			# 結果のインデックス
			resultIndex <- c("Correlation", "R2", "RMSE", "varU", "varE", "h2")
			
			predictionDataRAINBOW <- rep(NA, nrow(y))
			
			for (t in 1:nrow(y)) {
				print(c(t, "/", nrow(y)))
				
				tf <- !(1:nrow(y) == t)
				if (!is.null(Tv)) {
					ZETA[[1]][["K"]][tf, ] <- Tv[tf, ]
				}
				
				yNa <- y[, i]
				yNa[t] <- NA
				
				# RAINBOWR
				resEM3 <- EM3.cpp(y = yNa,
													X0 = X0, n.core = 1, 
													ZETA = ZETA)
				
				# 予測データの入力
				predictionDataRAINBOW[t] <- resEM3$y.pred[t]
			}
			
			predictData <- predictionDataRAINBOW
			obsData <- y[, i]
			
			# R2とRMSEを計算
			correlation <- cor(obsData, predictData)
			R2 <- 1 - sum((obsData - predictData) ^ 2) / sum((obsData - mean(obsData)) ^ 2)
			RMSE <- sqrt(sum((obsData - predictData) ^ 2) / length(obsData))
			VU <- resEM3$Vu
			VE <- resEM3$Ve
			h2 <- VU / (VU + VE)
			
			resultsum <- c(correlation, R2, RMSE, VU, VE, h2)
			resultsum <- as.matrix(resultsum)
			colnames(resultsum) <- trait.name
			var.pheno.all <- cbind(var.pheno.all, resultsum)
		}
		
		rownames(var.pheno.all) <- resultIndex
		print(var.pheno.all)
		xr[[zn]] <- var.pheno.all
	}
	
	resultModelname <- names(ZETAlist)
	resultEachModel <- matrix(NA, nrow = length(ZETAlist), ncol = ncol(y) - 1)
	rownames(resultEachModel) <- names(ZETAlist)
	colnames(resultEachModel) <- colnames(y)[2:ncol(y)]
	
	# 各モデルの結果をまとめる
	for (n in 1:length(names)) {
		resultEachModel[n, ] <- xr[[n]][1, ]
	}
	
	names(xr) <- names
	resultEachModel_all <- xr
	
	return(list(resultEachModel = resultEachModel, resultEachModel_all = resultEachModel_all))
}

getBLUP_simulation <- function(ZETAlist, y, Tv = NULL) {
  
  names <- names(ZETAlist)
  resultModelname <- names(ZETAlist)
  
  iterations <- length(ZETAlist)
  cl <- makeCluster(4)
  registerDoSNOW(cl)
  pb <- txtProgressBar(max = iterations, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  # foreach の結果をリストに保存
  results <- foreach(zn = 1:iterations, .options.snow = opts, .export = c("EM3.cpp")) %dopar% {
    zname <- names(ZETAlist)[zn]
    print(zname)
    
    var.pheno.all <- c() # trait ごとの結果(要約)
    ZETA <- ZETAlist[[zn]]
    
    
    X0 <- NULL
    predictDataEachModel <- matrix(NA, nrow = nrow(y), ncol = ncol(y)) # 保存用行列
    
    for (i in 1:ncol(y)) {
      trait.name <- colnames(y)[i]
      print(c(i, "/", ncol(y)))
      print(trait.name)
      
      resultIndex <- c("Correlation", "R2", "RMSE", "varU", "varE", "h2")
      predictionDataRAINBOW <- rep(NA, nrow(y))
      
      for (t in 1:nrow(y)) {
        tf <- !(1:nrow(y) == t)
        if (!is.null(Tv)) {
          ZETA[[1]][["K"]][tf, ] <- Tv[tf, ]
        }
        
        
        
        yNa <- y[, i]
        yNa[t] <- NA
        resEM3 <- EM3.cpp(y = yNa, X0 = X0, n.core = 1, ZETA = ZETA)
        
        # 予測データを入力
        predictionDataRAINBOW[t] <- resEM3$y.pred[t]
      }
      
      # 各 trait の結果を保存
      predictDataEachModel[, i] <- predictionDataRAINBOW
      obsData <- y[, i]
      correlation <- cor(obsData, predictionDataRAINBOW)
      R2 <- 1 - sum((obsData - predictionDataRAINBOW) ^ 2) / sum((obsData - mean(obsData)) ^ 2)
      RMSE <- sqrt(sum((obsData - predictionDataRAINBOW) ^ 2) / length(obsData))
      VU <- resEM3$Vu
      VE <- resEM3$Ve
      h2 <- VU / (VU + VE)
      
      resultsum <- c(correlation, R2, RMSE, VU, VE, h2)
      resultsum <- as.matrix(resultsum)
      colnames(resultsum) <- trait.name
      var.pheno.all <- cbind(var.pheno.all, resultsum)
    }
    
    rownames(var.pheno.all) <- resultIndex
    print(var.pheno.all)
    
    # 各モデルごとの予測データと結果をリストで返す
    list(varPhenoAll = var.pheno.all, predictData = predictDataEachModel)
  }
  
  close(pb)
  stopCluster(cl)
  
  # `results` から必要なデータを抽出
  varPhenoResults <- lapply(results, function(res) res$varPhenoAll) # 各結果
  predictDataList <- lapply(results, function(res) res$predictData) # 各予測データ
  
  # `predictDataList` を行列形式に統合
  predictDataMatrix <- do.call(rbind, predictDataList)
  
  # 結果モデル名と結果を保存
  resultModelname <- names(ZETAlist)
  resultEachModel <- matrix(NA, nrow = length(ZETAlist), ncol = ncol(y))
  rownames(resultEachModel) <- names(ZETAlist)
  colnames(resultEachModel) <- colnames(y)[2:ncol(y)]
  
  for (n in 1:length(names)) {
    resultEachModel[n, ] <- varPhenoResults[[n]][1, ]
    if (n == 1) {
      p <- names[1]
    } else {
      p <- paste(p, names[n], sep = "")
    }
  }
  
  print(p)
  
  # `resultEachModel_all` を作成
  resultEachModel_all <- varPhenoResults
  
  # 予測データを CSV に保存（オプション）
  #write.csv(predictDataMatrix, file = "predictData.csv", row.names = FALSE)
  
  return(list(
    resultEachModel = resultEachModel,
    resultEachModel_all = resultEachModel_all,
    p = p,
    predictDataList = predictDataList # 出力に追加
  ))
}

getBLUP_ab <- function(ZETAlist, y, B,Bpred,re=2) {
  
  Tv_list <- list() 
  
  
  for (s in 1:nrow(B)) {
    B[s, ] <- Bpred[s, ]  
    Tv <- bt(B) 
    Tv_list[[s]] <- Tv  
  }
  
  names <- names(ZETAlist)
  resultModelname <- names(ZETAlist)
  
  iterations <- length(ZETAlist)
  
  cl <- makeCluster(4)
  registerDoSNOW(cl)
  pb <- txtProgressBar(max = iterations, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  results <- foreach(zn = 1:iterations, .options.snow = opts, .export = c("EM3.cpp")) %dopar% {
    zname <- names(ZETAlist)[zn]
    print(zname)
    
    var.pheno.all <- c() 
    #ZETA <- ZETAlist[[zn]]
    ZETA <- ZETAlist[[zn]]
    
    X0 <- NULL
    predictDataEachModel <- matrix(NA, nrow = nrow(y), ncol = ncol(y)) 
    
    for (i in 1:ncol(y)) {
      trait.name <- colnames(y)[i]
      print(c(i, "/", ncol(y)))
      print(trait.name)
      
      resultIndex <- c("Correlation", "R2", "RMSE", "varU", "varE", "h2")
      predictionDataRAINBOW <- rep(NA, nrow(y))
      
      for (t in 1:nrow(y)) {
        tf <- !(1:nrow(y) == t)
        
        ZETA[[re]][["K"]] <- Tv_list[[t]]
        
        yNa <- y[, i]
        yNa[t] <- NA
        resEM3 <- EM3.cpp(y = yNa, X0 = X0, n.core = 1, ZETA = ZETA)
        
        # 予測データを入力
        predictionDataRAINBOW[t] <- resEM3$y.pred[t]
      }
      
      # 各 trait の結果を保存
      predictDataEachModel[, i] <- predictionDataRAINBOW
      obsData <- y[, i]
      correlation <- cor(obsData, predictionDataRAINBOW)
      R2 <- 1 - sum((obsData - predictionDataRAINBOW) ^ 2) / sum((obsData - mean(obsData)) ^ 2)
      RMSE <- sqrt(sum((obsData - predictionDataRAINBOW) ^ 2) / length(obsData))
      VU <- resEM3$Vu
      VE <- resEM3$Ve
      h2 <- VU / (VU + VE)
      
      resultsum <- c(correlation, R2, RMSE, VU, VE, h2)
      resultsum <- as.matrix(resultsum)
      colnames(resultsum) <- trait.name
      var.pheno.all <- cbind(var.pheno.all, resultsum)
    }
    
    rownames(var.pheno.all) <- resultIndex
    print(var.pheno.all)
    
    # 各モデルごとの予測データと結果をリストで返す
    list(varPhenoAll = var.pheno.all, predictData = predictDataEachModel)
  }
  
  close(pb)
  stopCluster(cl)
  
  # `results` から必要なデータを抽出
  varPhenoResults <- lapply(results, function(res) res$varPhenoAll) # 各結果
  predictDataList <- lapply(results, function(res) res$predictData) # 各予測データ
  
  # `predictDataList` を行列形式に統合
  predictDataMatrix <- do.call(rbind, predictDataList)
  
  # 結果モデル名と結果を保存
  resultModelname <- names(ZETAlist)
  resultEachModel <- matrix(NA, nrow = length(ZETAlist), ncol = ncol(y))
  rownames(resultEachModel) <- names(ZETAlist)
  colnames(resultEachModel) <- colnames(y)[1:ncol(y)]
  
  for (n in 1:length(names)) {
    resultEachModel[n, ] <- varPhenoResults[[n]][1, ]
    if (n == 1) {
      p <- names[1]
    } else {
      p <- paste(p, names[n], sep = "")
    }
  }
  
  print(p)
  
  # `resultEachModel_all` を作成
  resultEachModel_all <- varPhenoResults
  
  # 予測データを CSV に保存（オプション）
  #write.csv(predictDataMatrix, file = "predictData.csv", row.names = FALSE)
  
  return(list(
    resultEachModel = resultEachModel,
    resultEachModel_all = resultEachModel_all,
    p = p,
    predictDataList = predictDataList # 出力に追加
  ))
}


getRF <- function(RFlist,Y, dummy =NULL, Tv=NULL){

iterations=ncol(Y)
namesRF <- names(RFlist)

cl <- makeCluster(4)
registerDoSNOW(cl)
pb <- txtProgressBar(max = iterations, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

y.pred.list<- foreach (i = 1:iterations, .errorhandling = "pass",.options.snow = opts,.export = c("ranger","cvRF","cvRFb")) %dopar% {
  
  trait.name <- colnames(Y)[i]
  print(paste(i, "/", ncol(Y), trait.name))
  
  y <- Y[, i]
  
  if (is.null(dummy)){
    y.pred.each <- purrr::map(RFlist, ~ {
      X <- .x
      X<-do.call(cbind,X)
      y.pred <- cvRF(y, X)
      data.frame(pred = y.pred)
    })
    y.pred.all <- do.call(cbind, y.pred.each)
  } else {
    y.pred.each <- purrr::map(RFlist, ~ {
      X <- .x
      if (dummy %in% names(X)) {
        print("yea")
        y.pred <- cvRFb(y, X, dummy,Tv)
        data.frame(pred = y.pred)
      }
    })
    y.pred.all <- do.call(cbind, y.pred.each[!sapply(y.pred.each, is.null)])
    }
  
  y.pred.all <- cbind(y,y.pred.all)
  
  rownames(y.pred.all) <- rownames(Y)
  colnames(y.pred.all)[1] <- "y.obs"
  colnames(y.pred.all)[2:ncol(y.pred.all)] <- names(RFlist)
  return(y.pred.all)
}

close(pb)
stopCluster(cl) 

names(y.pred.list) <- colnames(Y)

rem <- matrix(NA,3,length(namesRF))
resultIndex <- c("Correlation", "R2", "RMSE")
rownames(rem) <- resultIndex
colnames(rem) <- namesRF


xr <- list()
t=1
for (yi in y.pred.list){
  
  predictData <- yi[,-1]
  obsData <- yi[,1]
  
  # calculate the R2 and RMSE
  correlation <- cor(obsData, predictData)
  
  if(length(namesRF) > 1){
    R2 <- (1 - apply((obsData - predictData) ^ 2,2,sum) / sum((obsData - mean(obsData)) ^ 2))
    RMSE <- sqrt(apply((obsData - predictData) ^ 2,2,sum) / length(obsData))
  } else {
    R2 <- (1 - sum((obsData - predictData) ^ 2)) / sum((obsData - mean(obsData)) ^ 2)
    RMSE <- sqrt(sum((obsData - predictData) ^ 2) / length(obsData))
  }
  
  rem[1,]<-correlation
  rem[2,]<-R2
  rem[3,]<-RMSE
  xr[[t]]<-rem
  t=t+1
}

resultModelname <- namesRF
resultEachModel <- matrix(NA, nrow = length(RFlist), ncol = ncol(Y))
rownames(resultEachModel) <- names(RFlist)
colnames(resultEachModel) <- colnames(Y)

#index p generator

for (n in 1:ncol(Y)) {
  resultEachModel[,n]<- xr[[n]][1,]
}
 
for (n in 1:length(namesRF)) {
   if (n == 1){
    p  <- namesRF[1]
  } else {
    p  <- paste(p,namesRF[n],sep="")
  }
}

print(p)
names(xr) <- colnames(Y)

return(list(y.pred.list=y.pred.list,resultEachModel=resultEachModel,p=p))

}


cvRF<- function(y, X) {
  dat <- data.frame(y, X)
  y.pred <- matrix(NA, nrow = length(y), ncol = 1)
  
  for(j in 1:nrow(y.pred)) {
    dat.train <- dat[-j, ]
    dat.test <- dat[j, ]
    
    model <- ranger(y ~ ., data = dat.train)
    model.pred <- predict(model, data = dat.test)
    y.pred[j, 1] <- model.pred$predictions 
  }
  y.pred
}




getRFab <- function(RFlist,Y, dummy =NULL, Tv=NULL){
  
  iterations=ncol(Y)
  namesRF <- names(RFlist)
  
  cl <- makeCluster(4)
  registerDoSNOW(cl)
  pb <- txtProgressBar(max = iterations, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  y.pred.list<- foreach (i = 1:iterations, .errorhandling = "pass",.options.snow = opts,.export = c("ranger","cvRF","cvRFb","cvRFab")) %dopar% {
    
    trait.name <- colnames(Y)[i]
    print(paste(i, "/", ncol(Y), trait.name))
    
    y <- Y[, i]
    
    if (is.null(dummy)){
      y.pred.each <- purrr::map(RFlist, ~ {
        X <- .x
        X<-do.call(cbind,X)
        y.pred <- cvRF(y, X)
        data.frame(pred = y.pred)
      })
      y.pred.all <- do.call(cbind, y.pred.each)
    } else {
      y.pred.each <- purrr::map(RFlist, ~ {
        X <- .x
        print("yea")
        y.pred <- cvRFab(y, X, dummy,Tv)
        data.frame(pred = y.pred)
        
      })
      y.pred.all <- do.call(cbind, y.pred.each[!sapply(y.pred.each, is.null)])
    }
    
    y.pred.all <- cbind(y,y.pred.all)
    
    rownames(y.pred.all) <- rownames(Y)
    colnames(y.pred.all)[1] <- "y.obs"
    colnames(y.pred.all)[2:ncol(y.pred.all)] <- names(RFlist)
    return(y.pred.all)
  }
  
  close(pb)
  stopCluster(cl) 
  
  names(y.pred.list) <- colnames(Y)
  
  rem <- matrix(NA,3,length(namesRF))
  resultIndex <- c("Correlation", "R2", "RMSE")
  rownames(rem) <- resultIndex
  colnames(rem) <- namesRF
  
  
  xr <- list()
  t=1
  for (yi in y.pred.list){
    
    predictData <- yi[,-1]
    obsData <- yi[,1]
    
    # calculate the R2 and RMSE
    correlation <- cor(obsData, predictData)
    
    if(length(namesRF) > 1){
      R2 <- (1 - apply((obsData - predictData) ^ 2,2,sum) / sum((obsData - mean(obsData)) ^ 2))
      RMSE <- sqrt(apply((obsData - predictData) ^ 2,2,sum) / length(obsData))
    } else {
      R2 <- (1 - sum((obsData - predictData) ^ 2)) / sum((obsData - mean(obsData)) ^ 2)
      RMSE <- sqrt(sum((obsData - predictData) ^ 2) / length(obsData))
    }
    
    rem[1,]<-correlation
    rem[2,]<-R2
    rem[3,]<-RMSE
    xr[[t]]<-rem
    t=t+1
  }
  
  resultModelname <- namesRF
  resultEachModel <- matrix(NA, nrow = length(RFlist), ncol = ncol(Y))
  rownames(resultEachModel) <- names(RFlist)
  colnames(resultEachModel) <- colnames(Y)
  
  #index p generator
  
  for (n in 1:ncol(Y)) {
    resultEachModel[,n]<- xr[[n]][1,]
  }
  
  for (n in 1:length(namesRF)) {
    if (n == 1){
      p  <- namesRF[1]
    } else {
      p  <- paste(p,namesRF[n],sep="")
    }
  }
  
  print(p)
  names(xr) <- colnames(Y)
  
  return(list(y.pred.list=y.pred.list,resultEachModel=resultEachModel,p=p))
  
}


cvRFb <- function(y, X, dummy, Tv) {
  
  Xcb <- do.call(cbind,X)
  dat <- data.frame(y, Xcb)
  
  if (dummy %in% names(X)){
    X[[dummy]] <- Tv
    Xcb<-do.call(cbind,X)
  }
  
  dat2<- data.frame(y,Xcb)
  
  y.pred <- matrix(NA, nrow = length(y), ncol = 1)
  
  for(j in 1:nrow(y.pred)) {
    dat.train <- dat2[-j, ]
    dat.test  <- dat[j,  ]
    model <- ranger(y ~ ., data = dat.train)
    model.pred <- predict(model, data = dat.test)
    y.pred[j, 1] <- model.pred$predictions 
  }
  y.pred
}


cvRFab <- function(y, X, dummy, Tv) {
  
  Xcb <- do.call(cbind,X)
  dat <- data.frame(y, Xcb)
  
  X[[dummy]] <- Tv
  Xcb2<-do.call(cbind,X)

  dat2<- data.frame(y,Xcb2)
  
  y.pred <- matrix(NA, nrow = length(y), ncol = 1)
  
  for(j in 1:nrow(y.pred)) {
    dat.train <- dat2[-j, ]
    dat.test  <- dat[j,  ]
    model <- ranger(y ~ ., data = dat.train)
    model.pred <- predict(model, data = dat.test)
    y.pred[j, 1] <- model.pred$predictions 
  }
  y.pred
}




scale_to_variance <- function(x, var_target) {
  scaled <- scale(x)
  scaled * sqrt(var_target)
}



unbiased_scale_2 <- function(x,m=NULL,s=NULL) {
  l <- nrow(x)
  if (is.null(m)){
    m <- apply(x, 2, mean)
  }
  if (is.null(s)){
    s <- apply(x, 2, sd) * sqrt((l - 1) / l)
  }
  xx <- sweep(x, 2, m, FUN = "-")
  xx <- sweep(xx, 2, s, FUN = "/")
  return(xx)
}


unbiased_scale <- function(x) {
  l <- nrow(x)
  m <- apply(x, 2, mean)
  s <- apply(x, 2, sd) * sqrt((l - 1) / l)
  xx <- sweep(x, 2, m, FUN = "-")
  xx <- sweep(xx, 2, s, FUN = "/")
  return(xx)
}

# Example sigmoid function
sigmoid <- function(x) {
 return(1 / (1 + exp(-x)))
}

# Example ReLU function
relu<- function(x) {
 return(matrix(pmax(0, x), nrow = nrow(x), ncol = ncol(x)))  # Ensures the output is a matrix
}

format_without_leading_zero <- function(x) {
  s <- formatC(x, format="f", digits=3)  # 3桁まで表示
  #s <- sub("^0\\.", ".", s)  # 0.を.に置き換え
  return(s)
}



binary <- function(matrix_input) {
      ifelse(matrix_input > 0, 1, 0)
  }

phl<- function(res,mo,cd,rot=45){
  p<-res$p
  pred.all<-res$resultEachModel
  pred.all<- pred.all[c(1,7,8,9,10,11,12),]
  model.name <- c("G","G+Micro","G+Micro+PredMetBLUPf19","G+Micro+PredMetBLUPf20","G+Micro+PredMetRF_f19","G+Micro+PredMetRF_f20","G+MIcro+TrueMet")
  trait.name <- colnames(pred.all) %>% sort() #sort
  pred.all<- data.frame(pred.all) %>% select(trait.name)
  #select(data.frame(pred.all),trait.name)
  df.all <- NULL
  for(i in 1:length(trait.name)) {
    pred <-  pred.all[,i]
    df <- data.frame(trait = trait.name[i],model = model.name, mean = pred)
    df.all <- rbind(df.all, df)
  }
  
  df.all[df.all < 0] <- 0
  df.all$model <- factor(df.all$model, levels = model.name)
  df.all$trait <- factor(df.all$trait)
  levels(df.all$trait) <- colnames(pred.all)
  g <- ggplot(df.all)
  g <- g + geom_bar(aes(x = trait, y = mean, fill = model), 
                    stat = "identity",  position = "dodge") +ylim(0,0.9)
  #g <- g + scale_x_discrete(limits=df.all$trait)
  g <- g + ggtitle(paste0("Yield Prediction ",mo," in 2019 ",cd))
  g <- g + ylab("COR") 
  g <- g + xlab("")
  g <- g +labs(fill = "Model")
  g <- g + theme(axis.text.x = element_text(angle =rot, hjust=0.8, vjust=0.9))
  #g<- g +theme(axis.text.x = element_text(vjust = 10))
  print(g)
  
}

phlg<- function(res,mo,cd,rot=45){
p<-res$p
pred.all<-res$resultEachModel

pred.all<- pred.all[c(1,7,8,9,10,11),]
model.name <- c("G","G+PredMetBlupf19","G+PredMetBlupf20","G+PredMetRFf19","G+PredMetRFf20","G+TrueMet")

trait.name <- colnames(pred.all) %>% sort() #sort
pred.all<- data.frame(pred.all) %>% select(trait.name)

df.all <- NULL

for(i in 1:length(trait.name)) {
  pred <-  pred.all[,colnames(pred.all)==trait.name[i]]
  df <- data.frame(trait = trait.name[i],model = model.name, mean = pred)
  df.all <- rbind(df.all, df)
}


df.all[df.all < 0] <- 0
df.all$model <- factor(df.all$model, levels = model.name)
df.all$trait <- factor(df.all$trait)
levels(df.all$trait) <- trait.name # trait
g <- ggplot(df.all)
g <- g + geom_bar(aes(x = trait, y = mean, fill = model), 
                  stat = "identity",  position = "dodge") +ylim(0,0.9) 
g <- g + ggtitle(paste0("Yield Prediction ",mo," in 2019 ",cd))
g <- g + ylab("COR")
g <- g + xlab("")
g <- g +labs(fill = "Model")
g <- g + theme(axis.text.x = element_text(angle =rot, hjust=0.8, vjust=0.9))

print(g)
}


phlb<- function(res_bpred,mo,cd,rot=45,b=TRUE){
p<-res_bpred$p
pred.all<-res_bpred$resultEachModel
model.name <- c("PredMetBLUPf19","PredMetBLUPf20","PredMetRF_f19","PredMetRF_f20")
trait.name <- colnames(pred.all) %>% sort() #sort
pred.all<- data.frame(pred.all) %>% select(trait.name)
df.all <- NULL
for(i in 1:length(trait.name)) {
  pred <-  pred.all[,colnames(pred.all)==trait.name[i]]
  df <- data.frame(trait = trait.name[i],model = model.name, mean = pred)
  df.all <- rbind(df.all, df)
}
df.all[df.all < 0] <- 0
df.all$model <- factor(df.all$model, levels = model.name)
df.all$trait <- factor(df.all$trait)
levels(df.all$trait) <- trait.name # trait
g <- ggplot(df.all)
g <- g + geom_bar(aes(x = trait, y = mean, fill = model), 
                  stat = "identity",  position = "dodge") +ylim(0,0.9) 
if (b==TRUE){
  g <- g + ggtitle(paste0("Yield Prediction ",mo," in 2019 ",cd))
} else{
  g <- g + ggtitle(paste0("Yield Prediction ",mo," in 2019 ",cd,"(pred by True Met)"))
}

g <- g +labs(fill = "Model")
g <- g + ylab("COR")
g <- g + xlab("")
g <- g + theme(axis.text.x = element_text(angle =rot, hjust=0.8, vjust=0.9))
print(g)
}


phlrf<- function(res,mo,cd,rot=45){
 p<-res$p
 pred.all<-res$resultEachModel
 pred.all<- pred.all[c(7,10,11),]
 model.name <- c("G+Micro","G+Micro+PredMetRF_f19","G+Micro+PredMetRF_f20")
 trait.name <- colnames(pred.all) %>% sort() #sort
 pred.all<- data.frame(pred.all) %>% select(trait.name)
 #select(data.frame(pred.all),trait.name)
 df.all <- NULL
 for(i in 1:length(trait.name)) {
  pred <-  pred.all[,i]
  df <- data.frame(trait = trait.name[i],model = model.name, mean = pred)
  df.all <- rbind(df.all, df)
 }
 
 df.all[df.all < 0] <- 0
 df.all$model <- factor(df.all$model, levels = model.name)
 df.all$trait <- factor(df.all$trait)
 levels(df.all$trait) <- colnames(pred.all)
 g <- ggplot(df.all)
 g <- g + geom_bar(aes(x = trait, y = mean, fill = model), 
                   stat = "identity",  position = "dodge") +ylim(0,0.9)
 #g <- g + scale_x_discrete(limits=df.all$trait)
 g <- g + ggtitle(paste0("G + Micro + Met model : Yield Prediction ",mo," in 2019 ",cd))
 g <- g + ylab("COR") 
 g <- g + xlab("")
 g <- g +labs(fill = "Model")
 #g <- g + theme(axis.text.x = element_text(angle =rot, hjust=0.8, vjust=0.9))
 g <- g +  theme(
 	axis.text.x = element_text(angle = rot, hjust = 0.8, vjust = 0.9),
 	legend.position = c(0.5, -0.40),  # Custom position
 	legend.box = "horizontal",        # Horizontal arrangement
 	legend.direction = "horizontal",   # Ensures items are horizontal
 	legend.key.size = unit(0.5, "cm"),
 	legend.text = element_text(size = 10),
 	legend.justification = "center",   # Center the legend at the specified position
 	legend.margin = margin(t = -5),
 	plot.title = element_text(hjust = 0.5) 
 )
 #g<- g +theme(axis.text.x = element_text(vjust = 10))
 gentle_colors <- c(
   
   rgb(198 / 255, 79 / 255, 83 / 255),  # Red
   rgb(249 / 255, 177 / 255, 25 / 255), # Yellow
   rgb(95 / 255, 158 / 255, 160 / 255)  # Soft cyan for the 
                    
                    )  # Light Blue
 

 
 #gentle_colors <- c("#D95F02", "#1B9E77", "#7570B3")  # Softer tones of red, green, blue
 
 g <- g + scale_fill_manual(values = gentle_colors)
}


phlrf4 <- function(res, mo, cd, rot = 45) {
  # Extract data from the result
  p <- res$p
  pred.all <- res$resultEachModel
  pred.all <- pred.all[c(1, 7, 10, 11), ] # Select the four relevant models
  model.name <- c("G", "G+Micro", "G+Micro+PredMetRF_f19", "G+Micro+PredMetRF_f20")
  trait.name <- colnames(pred.all) %>% sort() # Sort trait names
  
  # Ensure pred.all has columns ordered by trait.name
  pred.all <- data.frame(pred.all) %>% select(all_of(trait.name))
  
  # Build a combined data frame for plotting
  df.all <- NULL
  for (i in seq_along(trait.name)) {
    pred <- pred.all[, i]
    df <- data.frame(trait = trait.name[i], model = model.name, mean = pred)
    df.all <- rbind(df.all, df)
  }
  
  # Replace negative values with 0 and ensure factors are ordered
  df.all$mean[df.all$mean < 0] <- 0
  df.all$model <- factor(df.all$model, levels = model.name)
  df.all$trait <- factor(df.all$trait, levels = trait.name)
  
  # Create the ggplot object
  g <- ggplot(df.all) +
    geom_bar(
      aes(x = trait, y = mean, fill = model),
      stat = "identity",
      position = "dodge"
    ) +
    ylim(0, 0.9) + # Adjust y-axis limit
    ggtitle(paste0("G + Micro + Met model : Yield Prediction ", mo, " in 2019 ", cd)) +
    ylab("COR") +
    xlab("") +
    labs(fill = "Model") +
    theme(
      axis.text.x = element_text(angle = rot, hjust = 0.8, vjust = 0.9),
      legend.position = c(0.5, -0.4),
      legend.box = "horizontal",
      legend.direction = "horizontal",
      legend.key.size = unit(0.5, "cm"),
      legend.text = element_text(size = 10),
      legend.justification = "center",
      legend.margin = margin(t = -5),
      plot.title = element_text(hjust = 0.5)
    )
  
  # Define gentle color palette
  gentle_colors <- c(
    rgb(198 / 255, 79 / 255, 83 / 255),  # Red
    rgb(82 / 255, 182 / 255, 139 / 255), # Green
    rgb(249 / 255, 177 / 255, 25 / 255), # Yellow
    rgb(95 / 255, 158 / 255, 160 / 255)  # Soft cyan for the 4th element
  )
  
  # Apply the color palette to the plot
  g <- g + scale_fill_manual(values = gentle_colors)
  
  return(g)
}

phlrf4_2 <- function(res, mo, cd, rot = 45) {
  # Extract data from the result
  p <- res$p
  pred.all <- res$resultEachModel
  pred.all <- pred.all[c(1, 7, 10, 11), ] # Select the four relevant models
  model.name <- c("G", "G+Micro", "G+Micro+PredMetRF_f19", "G+Micro+PredMetRF_f20")
  trait.name <- colnames(pred.all) %>% sort() # Sort trait names
  
  # Ensure pred.all has columns ordered by trait.name
  pred.all <- data.frame(pred.all) %>% select(all_of(trait.name))
  
  # Build a combined data frame for plotting
  df.all <- NULL
  for (i in seq_along(trait.name)) {
    pred <- pred.all[, i]
    df <- data.frame(trait = trait.name[i], model = model.name, mean = pred)
    df.all <- rbind(df.all, df)
  }
  
  # Replace negative values with 0 and ensure factors are ordered
  df.all$mean[df.all$mean < 0] <- 0
  df.all$model <- factor(df.all$model, levels = model.name)
  df.all$trait <- factor(df.all$trait, levels = trait.name)
  
  # Create the ggplot object
  g <- ggplot(df.all) +
    geom_bar(
      aes(x = trait, y = mean, fill = model),
      stat = "identity",
      position = "dodge"
    ) +
    ylim(0, 0.9) + # Adjust y-axis limit
    ggtitle(paste0("Yield Prediction : G-Micro-Proxy model : ", cd)) +
    ylab("COR") +
    xlab("") +
    labs(fill = "Model") +
    theme(
      axis.text.x = element_text(angle = rot, hjust = 0.8, vjust = 0.9),
      legend.position = c(0.5, -0.35),
      legend.box = "horizontal",
      legend.direction = "horizontal",
      legend.key.size = unit(0.5, "cm"),
      legend.text = element_text(size = 10),
      legend.justification = "center",
      legend.margin = margin(t = -5),
      plot.title = element_text(hjust = 0.5)
    )
  
  # Define gentle color palette
  gentle_colors <- c(
    rgb(198 / 255, 79 / 255, 83 / 255),  # Red
    rgb(82 / 255, 182 / 255, 139 / 255), # Green
    rgb(249 / 255, 177 / 255, 25 / 255), # Yellow
    rgb(95 / 255, 158 / 255, 160 / 255)  # Soft cyan for the 4th element
  )
  
  # Apply the color palette to the plot
  g <- g + scale_fill_manual(values = gentle_colors)
  
  return(g)
}

phlgrf<- function(res,mo,cd,rot=45){
 p<-res$p
 pred.all<-res$resultEachModel
 
 pred.all<- pred.all[c(1,9,10),]
 model.name <- c("G","G+PredMetRFf19","G+PredMetRFf20")
 
 trait.name <- colnames(pred.all) %>% sort() #sort
 pred.all<- data.frame(pred.all) %>% select(trait.name)
 
 df.all <- NULL
 
 for(i in 1:length(trait.name)) {
  pred <-  pred.all[,colnames(pred.all)==trait.name[i]]
  df <- data.frame(trait = trait.name[i],model = model.name, mean = pred)
  df.all <- rbind(df.all, df)
 }
 
 
 df.all[df.all < 0] <- 0
 df.all$model <- factor(df.all$model, levels = model.name)
 df.all$trait <- factor(df.all$trait)
 levels(df.all$trait) <- trait.name # trait
 g <- ggplot(df.all)
 g <- g + geom_bar(aes(x = trait, y = mean, fill = model), 
                   stat = "identity",  position = "dodge") +ylim(0,0.9) 
 g <- g + ggtitle(paste0("G-Proxy model : Yield Prediction ",mo," in 2019 ",cd))
 g <- g + ylab("COR")
 g <- g + xlab("")
 g <- g +labs(fill = "Model")
 #g <- g + theme(axis.text.x = element_text(angle =rot, hjust=0.8, vjust=0.9))
 g <- g +  theme(
 	axis.text.x = element_text(angle = rot, hjust = 0.8, vjust = 0.9),
 	legend.position = c(0.5, -0.32),  # Custom position
 	legend.box = "horizontal",        # Horizontal arrangement
 	legend.direction = "horizontal",   # Ensures items are horizontal
 	legend.key.size = unit(0.5, "cm"),
 	legend.text = element_text(size = 10),
 	legend.justification = "center",   # Center the legend at the specified position
 	legend.margin = margin(t = -5),
 	plot.title = element_text(hjust = 0.5) 
 )
 # gentle_colors <- c(rgb(198/255, 79/255, 83/255),  # Light Red
 #                    rgb(82/255, 182/255, 139/255),  # Light Green
 #                    rgb(249/255, 177/255, 25/255))  # Light Blue
 # 
 gentle_colors <- c(
   rgb(198 / 255, 79 / 255, 83 / 255),  # Red
   rgb(249 / 255, 177 / 255, 25 / 255), # Yellow
   rgb(95 / 255, 158 / 255, 160 / 255)  # Soft cyan for the 4th element
 )
 
 #gentle_colors <- c("#D95F02", "#1B9E77", "#7570B3")  # Softer tones of red, green, blue
 
 g <- g + scale_fill_manual(values = gentle_colors)
}

fi <- function(x){
 x/x[1] *100
}

bt <-function(M){
  K<-tcrossprod(M) / ncol(M)
  return(K)
}


phlrf100<- function(res,mo,cd,rot=45){
 p<-res$p
 pred.all<-res$resultEachModel
 
 pred.all<- pred.all[c(7,10,11),]
 model.name <- c("G+Micro","G+Micro+PredMetRFf19","G+Micro+PredMetRFf20")
 
 pred.all<- apply(pred.all,2,fi)
 
 trait.name <- colnames(pred.all) %>% sort() #sort
 pred.all<- data.frame(pred.all) %>% select(trait.name)
 
 df.all <- NULL
 
 for(i in 1:length(trait.name)) {
  pred <-  pred.all[,colnames(pred.all)==trait.name[i]]
  df <- data.frame(trait = trait.name[i],model = model.name, mean = pred)
  df.all <- rbind(df.all, df)
 }
 
 
 df.all[df.all < 0] <- 0
 df.all$model <- factor(df.all$model, levels = model.name)
 df.all$trait <- factor(df.all$trait)
 levels(df.all$trait) <- trait.name # trait
 g <- ggplot(df.all)
 g <- g + geom_bar(aes(x = trait, y = mean, fill = model), 
                   stat = "identity",  position = "dodge") 
 #+ylim(0,0.9) 
 g <- g + ggtitle(paste0("G + Micro + Met model : Yield Prediction in 2019 ",cd))
 g <- g + ylab("COR (%)")
 g <- g + xlab("")
 g <- g +labs(fill = "Model")
 g <- g +  theme(
  axis.text.x = element_text(angle = rot, hjust = 0.8, vjust = 0.9),
  legend.position = c(0.5, -0.40),  # Custom position
  legend.box = "horizontal",        # Horizontal arrangement
  legend.direction = "horizontal",   # Ensures items are horizontal
  legend.key.size = unit(0.5, "cm"),
  legend.text = element_text(size = 10),
  legend.justification = "center",   # Center the legend at the specified position
  legend.margin = margin(t = -5),
  plot.title = element_text(hjust = 0.5) 
 )
 
 g <- g + geom_hline(yintercept = 100, color = "black")
 #gentle_colors <- c("#D95F02", "#1B9E77", "#7570B3")  # Softer tones of red, green, blue
 gentle_colors <- c(rgb(198/255, 79/255, 83/255),  # Light Red
                    rgb(82/255, 182/255, 139/255),  # Light Green
                    rgb(249/255, 177/255, 25/255))  # Light Blue
 
 g <- g + scale_fill_manual(values = gentle_colors)
 print(g)
}


phlgrf100<- function(res,mo,cd,rot=45){
 p<-res$p
 pred.all<-res$resultEachModel
 
 pred.all<- pred.all[c(1,9,10),]
 model.name <- c("G","G+PredMetRFf19","G+PredMetRFf20")
 

 pred.all<- apply(pred.all,2,fi)

 trait.name <- colnames(pred.all) %>% sort() #sort
 pred.all<- data.frame(pred.all) %>% select(trait.name)
 
 df.all <- NULL
 
 for(i in 1:length(trait.name)) {
  pred <-  pred.all[,colnames(pred.all)==trait.name[i]]
  df <- data.frame(trait = trait.name[i],model = model.name, mean = pred)
  df.all <- rbind(df.all, df)
 }
 
 
 df.all[df.all < 0] <- 0
 df.all$model <- factor(df.all$model, levels = model.name)
 df.all$trait <- factor(df.all$trait)
 levels(df.all$trait) <- trait.name # trait
 g <- ggplot(df.all)
 g <- g + geom_bar(aes(x = trait, y = mean, fill = model), 
                   stat = "identity",  position = "dodge") 
 #+ylim(0,0.9) 
 g <- g + ggtitle(paste0("G + Met model : Yield Prediction in 2019 ",cd))
 g <- g + ylab("COR (%)")
 g <- g + xlab("")
 g <- g +labs(fill = "Model")
 g <- g +  theme(
  axis.text.x = element_text(angle = rot, hjust = 0.8, vjust = 0.9),
  legend.position = c(0.5, -0.40),  # Custom position
  legend.box = "horizontal",        # Horizontal arrangement
  legend.direction = "horizontal",   # Ensures items are horizontal
  legend.key.size = unit(0.5, "cm"),
  legend.text = element_text(size = 10),
  legend.justification = "center",   # Center the legend at the specified position
  legend.margin = margin(t = -5),
  plot.title = element_text(hjust = 0.5) 
 )
 g <- g + geom_hline(yintercept = 100, color = "black")
 gentle_colors <- c(rgb(198/255, 79/255, 83/255),  # Light Red
                    rgb(82/255, 182/255, 139/255),  # Light Green
                    rgb(249/255, 177/255, 25/255))  # Light Blue
 
 #gentle_colors <- c("#D95F02", "#1B9E77", "#7570B3")  # Softer tones of red, green, blue

 g <- g + scale_fill_manual(values = gentle_colors)
}



#' @title Estimate an appropriate value for the ridge penalty (lambda).
#' @description Estimate an appropriate value for the ridge penalty (lambda).
#' @param Y A numeric matrix of response variables.
#' @param X A numeric matrix of explanatory variables.
#' @return A numeric vector of the range of lambda values.
glam <- function(Y, X) {
 y_rate<- 1
 Y<-as.matrix(Y)
 X<-as.matrix(X)
 if(ncol(X)>1000){
  X<-X[,sample(ncol(X), 1000, replace = FALSE)]
 }
 if(ncol(Y)>1000){
  y_rate<- ncol(Y)/1000
  Y<-Y[,sample(ncol(Y), 1000, replace = FALSE)]
 }
 usx<-unb(X)
 usy<-unb(Y)
 
 tmp_sum<- colSums((crossprod(usy,usx))^2) * y_rate
 result <- max(sqrt(tmp_sum)) / (nrow(Y) * 10^(-3))
 e<- 10^seq(-4,0, length.out = 50)
 lambda_def<- result * e
 lambda_def<- signif(lambda_def,4)
 return(lambda_def)
}


#' @title Scale a matrix using unbiased estimators for the mean and standard deviation.
#' @description Scale a matrix using unbiased estimators for the mean and standard deviation.
#' @param x A numeric matrix to be scaled.
#' @importFrom stats sd
#' @return A scaled numeric matrix.
unb <- function(x) {
 l <- nrow(x)
 xx<-scale(x) / sqrt((l-1)/l)
 xx[is.na(xx)] <- 0
 return(xx)
}



glam2 <- function(Y, X) {
 y_rate<- 1
 Y<-as.matrix(Y)
 X<-as.matrix(X)
 
 if(ncol(X)>1000){
  X<-X[,sample(ncol(X), 1000, replace = FALSE)]
 }
 if(ncol(Y)>1000){
  y_rate<- ncol(Y)/1000
  Y<-Y[,sample(ncol(Y), 1000, replace = FALSE)]
 }
 usx<-unb(X)
 usy<-Y
 
 if(ncol(Y)==1){
  result <- max(abs(colSums(usx*Y[,1])))/nrow(Y)
 } else{
  tmp_sum<- colSums((crossprod(usy,usx))^2) * y_rate
  result <- max(sqrt(tmp_sum)) / (nrow(Y) * 10^(-3))
 }

 e<- 10^seq(-4,0, length.out = 50)
 lambda_def<- result * e
 lambda_def<- signif(lambda_def,4)
 return(lambda_def)
}

glam3 <- function(Y, X, scale.X = FALSE) {
 Y<-as.matrix(Y)
 X<-as.matrix(X)
 y_rate<- 1
 if(ncol(X)>1000){
  X<-X[,sample(ncol(X), 1000, replace = FALSE)]
 }
 if(ncol(Y)>1000){
  y_rate<- ncol(Y)/1000
  Y<-Y[,sample(ncol(Y), 1000, replace = FALSE)]
 }
 
 if (scale.X){
  X<-unb(X)
 }
 
 if(ncol(Y)==1){
  result <- max(abs(colSums(X*Y[,1])))/nrow(Y)
 } else{
  tmp_sum<- colSums((crossprod(Y,X))^2) * y_rate
  result <- max(sqrt(tmp_sum)) / (nrow(Y) * 10^(-3))
 }
 
 e<- 10^seq(-4,0, length.out = 50)
 lambda_def<- result * e
 lambda_def<- signif(lambda_def,4)
 return(lambda_def)
}

rsbc_stars <- function(Y,X,num,maxrank,save,nfold=5){
 
 # rrr_result <- rrr2(Y, X, penaltySVD = c("ann"))
 # Lambda<-rrr_result$lambda
 lm<-readRDS(file.path(base_path, "aic_lambda.RDS"))
 results <- list()
 
 pb <- txtProgressBar(min = 0, max = num, style = 3)
 for(i in 1:num){
  Lambda<-lm[,i]
  rn<-sample(c(1:nrow(X)),sample,rep=F)
  
  Y_train<-Y[-rn,]
  Y_test<-Y[rn,]
  X_train<-X[-rn,]
  X_test<-X[rn,]
  
  sink(nullfile())
  plan(multisession)
  tic()
  cv_result<- rrda.cv(Y = Y_train, X = X_train, lambda=Lambda,maxrank =maxrank,nfold = nfold,scale.X = F,scale.Y = F)
  a<-toc()
  plan(sequential)
  sink()
  K<-cv_result$opt_min$rank
  L<-cv_result$opt_min$lambda
  K1<-cv_result$opt_rank.1se$rank
  L1<-cv_result$opt_lambda.1se$lambda
  
  et<-a$toc-a$tic
  cv_result$et<-et
  results[[i]] <- list(
   rn = rn,
   cv_result = cv_result
  )
  
  setTxtProgressBar(pb, i)
 }
 close(pb)
 
 file_path <- file.path(save, "cv_stars.RDS")
 saveRDS(results,file_path)
 
 results <-readRDS(file_path)
 b_l <- c()
 b_l.1se <-c()
 b_r <- c()
 b_r.1se <- c()
 time_cv <- c()
 
 for (i in 1:num) {
  b_l <- c(b_l, results[[i]]$cv_result$opt_min$lambda)
  b_l.1se <- c(b_l.1se, results[[i]]$cv_result$opt_lambda.1se$lambda)
  b_r <- c(b_r, results[[i]]$cv_result$opt_min$rank)
  b_r.1se <- c(b_r.1se, results[[i]]$cv_result$opt_rank.1se$rank)
  time_cv <- c(time_cv, results[[i]]$cv_result$et)
 }
 
 print(paste("estimated maxrank",max(b_r)))
 
 possible_ranks <- 1:max(b_r)
 rank_table <- table(factor(b_r, levels = possible_ranks))
 barplot(rank_table, main = paste0("CV min Estimated rank"), xlab = "Rank", ylab = "Frequency")
 
 possible_ranks <- 1:max(b_r)
 rank_table <- table(factor(b_r.1se, levels = possible_ranks))
 barplot(rank_table, main = paste0("CV 1se Estimated rank"), xlab = "Rank", ylab = "Frequency")
 
 ## mspe
 #```{r,eval=F}
 set.seed(123)
 results <-readRDS(file_path)
 
 # n_lam <- 50
 # Lambda <- 10^seq(0, 6, length.out = n_lam)
 # 
 results_mspe <- list()
 
 
 Y<-dna[,dch==chr]
 X<-rna[,rch==chr]
 
 pb <- txtProgressBar(min = 0, max = num, style = 3)
 for(i in 1:num){
  
  rn<-results[[i]]$rn
  
  cv_result<-results[[i]]$cv_result
  #rrda.plot(cv_result = cv_result,show_error_bar = FALSE)
  
  K<-cv_result$opt_min$rank
  L<-cv_result$opt_min$lambda
  K1<-cv_result$opt_rank.1se$rank
  L1<-cv_result$opt_lambda.1se$lambda
  K2<- min(dim(Y_train),dim(X_train))

  Y_train<-Y[-rn,]
  Y_test<-Y[rn,,drop=F]
  X_train<-X[-rn,]
  X_test<-X[rn,,drop=F]
  
  mspe <- c(NA,NA,NA,NA,NA,NA)
  cor <- c(NA,NA,NA,NA,NA,NA)
  rms<- c(K,K1,K,K2,K,K2)
  lms<- c(L,L,L1,L,0,0)
  
  
  for (ms in 1: length(mspe)){
   
   best_Bhat <- rrda.fit(Y = Y_train, X = X_train, nrank = rms[ms],lambda = lms[ms],scale.X = F,scale.Y = F)
   best_Pred <- rrda.predict(Bhat = best_Bhat, X = X_test)[[1]][[1]][[1]]
   
   mspe[ms] <- 100 * sum((Y_test - best_Pred)^2) / (ncol(Y_test) * nrow(Y_test))
   
   cor[ms] <- mean(diag(cor(Y_test, best_Pred)))
  }
  
  
  results_mspe[[i]] <- list(
   rn = rn,
   mspe = mspe,
   cor = cor,
   L = L,
   K = K,
   L1 = L1,
   K1 = K1,
   K2 = K2
  )
  
  setTxtProgressBar(pb, i)
 }
 close(pb)
 
 file_path <- file.path(save, "test_stars.RDS")
 saveRDS(results_mspe,file_path)
 
}

# 
# 
# simdata<-function(n,a,b,c,d,ws,w,model=c("p2","abs","sigmoid","relu","linear"),noise=0.1){
# 
#  theta.ab <- matrix(runif(a*b,ws,w),a,b)
#  theta.ac <- matrix(runif(a*c,ws,w),a,c)
#  
#  theta.bc <- matrix(runif(b*c, ws,w),b,c)
#  
#  theta.db <- matrix(runif(d*b,ws,w),d,b)
#  theta.dc <- matrix(runif(d*c,ws,w),d,c)
#  
#  noise.B <- matrix(rnorm(n = n*b, mean = 0, sd = noise),n,b)
#  noise.C <- matrix(rnorm(n = n*c, mean = 0, sd = noise),n,c)
#  
#  A <- matrix(rnorm(n = n*a, mean = 0, sd = 1 / (a**(1/2))),n,a)
#  D <- matrix(rnorm(n = n*d, mean = 0, sd = 1 / (d**(1/2))),n,d)
#  
#  # Choose the model to generate B
#  if (model == "sigmoid") {
#   B <- scale(scale(sigmoid(A) %*% theta.ab + D %*% theta.db) + noise.B)
#  } else if (model == "relu") {
#   B <- scale(scale(relu(A) %*% theta.ab + D %*% theta.db) + noise.B)
#  } else if (model == "abs") {
#   B <- scale(scale(abs(A) %*% theta.ab + D %*% theta.db) + noise.B)
#  } else if (model == "p2") {
#   B <- scale(scale((A)^2 %*% theta.ab + D %*% theta.db) + noise.B)
#  } else if (model == "linear") {
#   B <- scale(scale(A %*% theta.ab + D %*% theta.db) + noise.B)
#  } else {
#   # Stop the function for any unrecognized model
#   stop("Error: Unrecognized model. Please use one of the following: 'p2', 'abs', 'sigmoid', or 'relu'.")
#  }
#  
#  B <- B / (b**(1/2))
#  C <- scale(A %*% theta.ac + B %*% theta.bc + D %*% theta.dc) + noise.C
#  #C <- A %*% theta.ac + B %*% theta.bc + D %*% theta.dc + noise.C
#  return(list(A = A, B = B, C = C, D = D))
# }



# 
# simdata2<-function(n,a,b,c,d,ws,w,model=c("p2","abs","sigmoid","relu","linear"),noise=0.1){
#  
#  theta.ab <- matrix(runif(a*b,ws,w),a,b)
#  theta.ac <- matrix(runif(a*c,ws,w),a,c)
#  
#  theta.bc <- matrix(runif(b*c, ws,w),b,c)
#  
#  theta.db <- matrix(runif(d*b,ws,w),d,b)
#  theta.dc <- matrix(runif(d*c,ws,w),d,c)
#  
#  noise.B <- matrix(rnorm(n = n*b, mean = 0, sd = noise),n,b)
#  noise.C <- matrix(rnorm(n = n*c, mean = 0, sd = noise),n,c)
#  
#  A <- matrix(rnorm(n = n*a, mean = 0, sd = 1 / (a**(1/2))),n,a)
#  D <- matrix(rnorm(n = n*d, mean = 0, sd = 1 / (d**(1/2))),n,d)
#  
#  # Choose the model to generate B
#  if (model == "sigmoid") {
#   B <- ((sigmoid(A) %*% theta.ab + D %*% theta.db) + noise.B)
#  } else if (model == "relu") {
#   B <- ((relu(A) %*% theta.ab + D %*% theta.db) + noise.B)
#  } else if (model == "abs") {
#   B <- ((abs(A) %*% theta.ab + D %*% theta.db) + noise.B)
#  } else if (model == "p2") {
#   B <- (((A)^2 %*% theta.ab + D %*% theta.db) + noise.B)
#  } else if (model == "linear") {
#   B <- ((A %*% theta.ab + D %*% theta.db) + noise.B)
#  } else {
#   # Stop the function for any unrecognized model
#   stop("Error: Unrecognized model. Please use one of the following: 'p2', 'abs', 'sigmoid', or 'relu'.")
#  }
#  
#  B <- B / (b**(1/2))
#  C <- (A %*% theta.ac + B %*% theta.bc + D %*% theta.dc) + noise.C
#  
#  return(list(A = A, B = B, C = C, D = D))
# }
# 


# simdata3<-function(n,a,b,c,d,ws,w,model=c("p2","abs","sigmoid","relu","linear"),noise=0.1){
#  
#  theta.ab <- matrix(rnorm(a * b, mean = 0, sd = 1), nrow = a, ncol = b)
#  theta.ac <- matrix(rnorm(a * c, mean = 0, sd = 1), nrow = a, ncol = c)
#  theta.bc <- matrix(rnorm(b * c, mean = 0, sd = 1), nrow = b, ncol = c)
#  theta.db <- matrix(rnorm(d * b, mean = 0, sd = 1), nrow = d, ncol = b)
#  theta.dc <- matrix(rnorm(d * c, mean = 0, sd = 1), nrow = d, ncol = c)
#  
#  noise.B <- matrix(rnorm(n = n*b, mean = 0, sd = noise),n,b)
#  noise.C <- matrix(rnorm(n = n*c, mean = 0, sd = noise),n,c)
#  
#  A <- matrix(rnorm(n = n*a, mean = 0, sd = 1 / (a**(1/2))),n,a)
#  D <- matrix(rnorm(n = n*d, mean = 0, sd = 1 / (d**(1/2))),n,d)
#  
#  # Choose the model to generate B
#  if (model == "sigmoid") {
#   B <- scale(scale(sigmoid(A) %*% theta.ab + D %*% theta.db) + noise.B)
#  } else if (model == "relu") {
#   B <- scale(scale(relu(A) %*% theta.ab + D %*% theta.db) + noise.B)
#  } else if (model == "abs") {
#   B <- scale(scale(abs(A) %*% theta.ab + D %*% theta.db) + noise.B)
#  } else if (model == "p2") {
#   B <- scale(scale((A)^2 %*% theta.ab + D %*% theta.db) + noise.B)
#  } else if (model == "linear") {
#   B <- scale(scale(A %*% theta.ab + D %*% theta.db) + noise.B)
#  } else {
#   # Stop the function for any unrecognized model
#   stop("Error: Unrecognized model. Please use one of the following: 'p2', 'abs', 'sigmoid', or 'relu'.")
#  }
#  
#  B <- B / (b**(1/2))
#  C <- scale(A %*% theta.ac + B %*% theta.bc + D %*% theta.dc) + noise.C
#  #C <- A %*% theta.ac + B %*% theta.bc + D %*% theta.dc + noise.C
#  return(list(A = A, B = B, C = C, D = D))
# }


simdata4<-function(n,a,b,c,d,ws,w,model=c("p2","p4","abs","sigmoid","relu","linear","exp","tanh","binary"),noise=0.1){
 
 theta.ab <- matrix(rnorm(a * b, mean = 0, sd = 1), nrow = a, ncol = b)
 theta.ac <- matrix(rnorm(a * c, mean = 0, sd = 1), nrow = a, ncol = c)
 theta.bc <- matrix(rnorm(b * c, mean = 0, sd = 1), nrow = b, ncol = c)
 theta.db <- matrix(rnorm(d * b, mean = 0, sd = 1), nrow = d, ncol = b)
 theta.dc <- matrix(rnorm(d * c, mean = 0, sd = 1), nrow = d, ncol = c)
 
 noise.B <- matrix(rnorm(n = n*b, mean = 0, sd = noise),n,b)
 noise.C <- matrix(rnorm(n = n*c, mean = 0, sd = noise),n,c)
 
 A <- matrix(rnorm(n = n*a, mean = 0, sd = 1 / (a**(1/2))),n,a)
 D <- matrix(rnorm(n = n*d, mean = 0, sd = 1 / (d**(1/2))),n,d)
 
 if (model == "sigmoid") {
  An <- sigmoid(A)
 } else if (model == "relu") {
  An <- relu(A)
 } else if (model == "abs") {
  An <- abs(A)
 } else if (model == "p2") {
  An <- (A)^2
 } else if (model == "p4") {
  An <- (A)^4
 } else if (model == "exp") {
  An <- exp(A)
 } else if (model == "tanh") {
  An <- tanh(A)
 } else if (model == "binary") {
  An <- binary(A)
 } else if (model == "linear") {
  An<-A
 } else {
  # Stop the function for any unrecognized model
  stop("Error: Unrecognized model.")
 }
 
 # Calculate B after transforming A
 B <- (An %*% theta.ab + D %*% theta.db + noise.B)
 
 B <- B / (b**(1/2))
 C <- (A %*% theta.ac + B %*% theta.bc + D %*% theta.dc) + noise.C
 #C <- A %*% theta.ac + B %*% theta.bc + D %*% theta.dc + noise.C
 return(list(A = A, B = B, C = C, D = D))
}


simdata6<-function(n,a,b,c,d,ws,w,model=c("p2","p4","abs","sigmoid","relu","linear","exp","tanh","binary"),noise=0.1){
  
  theta.ab <- matrix(rnorm(a * b, mean = 0, sd = 1), nrow = a, ncol = b)
  theta.ac <- matrix(rnorm(a * c, mean = 0, sd = 1), nrow = a, ncol = c)
  theta.bc <- matrix(rnorm(b * c, mean = 0, sd = 1), nrow = b, ncol = c)
  theta.db <- matrix(rnorm(d * b, mean = 0, sd = 1), nrow = d, ncol = b)
  theta.dc <- matrix(rnorm(d * c, mean = 0, sd = 1), nrow = d, ncol = c)
  
  noise.B <- matrix(rnorm(n = n*b, mean = 0, sd = noise),n,b)
  noise.C <- matrix(rnorm(n = n*c, mean = 0, sd = noise),n,c)
  
  A <- matrix(rnorm(n = n*a, mean = 0, sd = 1 / (a**(1/2))),n,a)
  D <- matrix(rnorm(n = n*d, mean = 0, sd = 1 / (d**(1/2))),n,d)
  
  if (model == "sigmoid") {
    An <- sigmoid(A)
  } else if (model == "relu") {
    An <- relu(A)
  } else if (model == "abs") {
    An <- abs(A)
  } else if (model == "p2") {
    An <- (A)^2
  } else if (model == "p4") {
    An <- (A)^4
  } else if (model == "exp") {
    An <- exp(A)
  } else if (model == "tanh") {
    An <- tanh(A)
  } else if (model == "binary") {
    An <- binary(A)
  } else if (model == "linear") {
    An<-A
  } else {
    # Stop the function for any unrecognized model
    stop("Error: Unrecognized model.")
  }
  
  # Calculate B after transforming A
  B <- (An %*% theta.ab + D %*% theta.db + noise.B)
  
  C <- (A %*% theta.ac + B %*% theta.bc + D %*% theta.dc) + noise.C
  #C <- A %*% theta.ac + B %*% theta.bc + D %*% theta.dc + noise.C
  return(list(A = A, B = B, C = C, D = D))
}


#
simdata5<-function(n,a,b,c,d,ws,w,model=c("p2","p4","abs","sigmoid","relu","linear","exp","tanh","binary"),noise=0.1){

 theta.ab <- matrix(rnorm(a * b, mean = 0, sd = 1), nrow = a, ncol = b)
 theta.ac <- matrix(rnorm(a * c, mean = 0, sd = 1), nrow = a, ncol = c)
 theta.bc <- matrix(rnorm(b * c, mean = 0, sd = 1), nrow = b, ncol = c)
 theta.db <- matrix(rnorm(d * b, mean = 0, sd = 1), nrow = d, ncol = b)
 theta.dc <- matrix(rnorm(d * c, mean = 0, sd = 1), nrow = d, ncol = c)

 noise.B <- matrix(rnorm(n = n*b, mean = 0, sd = noise),n,b)
 noise.C <- matrix(rnorm(n = n*c, mean = 0, sd = noise),n,c)

 A <- matrix(rnorm(n = n*a, mean = 0, sd = 1 / (a**(1/2))),n,a)
 D <- matrix(rnorm(n = n*d, mean = 0, sd = 1 / (d**(1/2))),n,d)

 # Apply the nonlinear transformation to A
 if (model == "sigmoid") {
  An <- sigmoid(A)
 } else if (model == "relu") {
  An <- relu(A)
 } else if (model == "abs") {
  An <- abs(A)
 } else if (model == "p2") {
  An <- (A)^2
 } else if (model == "p4") {
  An <- (A)^4
 } else if (model == "exp") {
  An <- exp(A)
 } else if (model == "tanh") {
  An <- tanh(A)
 } else if (model == "binary") {
  An <- binary(A)
 } else if (model == "linear") {
  An<-A
 } else {
  # Stop the function for any unrecognized model
  stop("Error: Unrecognized model.")
 }

 # Calculate B after transforming A
 B <- (An %*% theta.ab + D %*% theta.db + noise.B)

 B <- B / (b**(1/2))
 C <- (An %*% theta.ac + B %*% theta.bc + D %*% theta.dc) + noise.C
 #C <- A %*% theta.ac + B %*% theta.bc + D %*% theta.dc + noise.C
 return(list(A = A, B = B, C = C, D = D))
}



#. 

#stars is made ti be functions
# lambda option and rn options



