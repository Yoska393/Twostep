

getBLUP <- function(ZETAlist,y,Tv=NULL) {
	
	names<- names(ZETAlist)
	resultModelname <- names(ZETAlist)
	
	iterations=length(ZETAlist)
	cl <- makeCluster(10)
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

getRF <- function(RFlist,Y, dummy =NULL, Tv=NULL){
	
	iterations=ncol(Y)
	namesRF <- names(RFlist)
	
	cl <- makeCluster(10)
	registerDoSNOW(cl)
	pb <- txtProgressBar(max = iterations, style = 3)
	progress <- function(n) setTxtProgressBar(pb, n)
	opts <- list(progress = progress)
	
	y.pred.list<- foreach (i = 1:iterations,.options.snow = opts,.export = c("ranger","cvRF","cvRFb")) %dopar% {
		
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

# Example ReLU function
relu<- function(x) {
	return(matrix(pmax(0, x), nrow = nrow(x), ncol = ncol(x)))  # Ensures the output is a matrix
}

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
	#C <- (An %*% theta.ac + B %*% theta.bc + D %*% theta.dc) + noise.C
	C <- A %*% theta.ac + B %*% theta.bc + D %*% theta.dc + noise.C
	return(list(A = A, B = B, C = C, D = D))
}
