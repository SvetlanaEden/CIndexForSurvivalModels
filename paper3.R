#################### Third paper proposal:
### Read Uno & Pencina
### Compared to their statistic, ours:
### 1. Does not depend on the censoring distribution (becase it does not compute its distribution)
### 2. Does not require semiparametric model with beta coefficients in order to compute C-statistics
####################
### MY PLAN:
### Simulate the data using a finite tau and compare the two methods
### Simulate the data where censoring distribution depends on covariates and compare two models.
####################
### QUESTION:
### When we introduce tau, what is the population parameter of C|tau, when we have type-I censoring?
### How to evaluate the above population parameter? Run a model on full uncensored data and then compute C-index only on time < tau?
### Is it possible to weigh the correlation of the probability scale residuals by censoring distribution?
####################
### 20200615: QUESTION/ANSWER via sumulations:
### I could not fall asleep the other day, so I had the following thought. Suppose time to outcome (exponentially distributed) depends only on one covariate; and we estimate this covariate's beta using exponetial distribution.
### Suppose we also have censoring, which has the same marginal distribution (same covariate, same beta)
### Suppose we also estimating survival function using exponential model
### As long as beta is the same sing as the real beta, censoring should not be a problem, and the C-index should be unbiased.
### Rationale of this thought is that the predicted values will be concordant with beta*Z because we predict survival probability using exp model.
### This might not be true if we use non-parametric estimator of the c-index
### ANSWER: The estimator of Uno & Pencina is quite good under correctly specified model. In rare cases it does not compute, but is very precise event under dependent (conditionally independent) censoring.

####################
### 20200616: QUESTIONs/ACTION ITMES:
### Construct the parametric c-index
### Why Uno & Pencina did not adjust censoring distribution for covariantes?
### Compute from scratch c-index of Uno & Pencina and ajdust the censoring distr. for cov.

simSurvDataWithCensDepOnZ = function(subjNum, beta1, beta2, betaForCensoring, theta, family){
  ### it is not a causal association, but this is temporary
  designMatr1 = matrix(1, ncol=length(beta1), nrow = subjNum)
  designMatr1[,1] = rnorm(subjNum)
  designMatr1[,2] = sample(c(0, 1), size = subjNum, replace = TRUE)
  designMatr2 = designMatr1
  mean1 = exp(-designMatr1 %*% matrix(beta1, ncol=1))
  mean2 = exp(-designMatr2 %*% matrix(beta2, ncol=1))
  mean3 = exp(-designMatr2 %*% matrix(betaForCensoring, ncol=1))
  #### simulate the outcome and the covariate
  res12 = rclayton(nSim = subjNum, thetaPar = theta, lambda1 = mean1, lambda2 = mean2, censoringProb1 = 0, censoringProb2 = 0, independentCensoring = TRUE, family = family)
  #### simulate censoring variable res34[, "t1"]
  res34 = rclayton(nSim = subjNum, thetaPar = theta, lambda1 = mean1, lambda2 = mean3, censoringProb1 = 0, censoringProb2 = 0, independentCensoring = TRUE, family = family)
  #### make censored data:
  data = data.frame(time = res12[, "t1"], event = res12[, "delta1"], V1 = res12[, "t2"])
  if(any(betaForCensoring!= 0)){
    data$event = as.numeric(data$time <= res34[, "t1"])
    data$time[data$event == 0] = res34[data$event == 0, "t1"]
  }
  data
}

simSurvDataWithCensDepOnZ_take_two = function(subjNum, betaEvent, betaCens = NULL, indepCens = FALSE, theta, family){
  ### it is not a causal association, but this is temporary
  ### betaEvent = c(2, 3) and betaCens = c(3, 4) gives about 38% censored
  designMatr1 = matrix(1, ncol=length(betaEvent), nrow = subjNum)
  designMatr1[,1] = rnorm(subjNum)
  designMatr1[,2] = sample(c(0, 1), size = subjNum, replace = TRUE)
  
  meanEvent = exp(-designMatr1 %*% matrix(betaEvent, ncol=1))
  time = rexp(n = subjNum, rate = 1/meanEvent)
  data = data.frame(timeNoCens = time, time  = time, event = 1, V1 = designMatr1[,1], V2 = designMatr1[,2])

  # data1 = data.frame(time, event = 1, V1 = designMatr1[,1], V2 = designMatr1[,2])
  # survObj = Surv(data1$time, data1$event)
  # m1 = cph(survObj ~ V1 + V2, data = data1, x = TRUE, y = TRUE)
  
  if(!is.null(betaCens)){
    if(indepCens){
      ### define a new matrix for censoring if censoring is independent
      designMatr1 = matrix(1, ncol=length(betaEvent), nrow = subjNum)
      designMatr1[,1] = rnorm(subjNum)
      designMatr1[,2] = sample(c(0, 1), size = subjNum, replace = TRUE)
    }
    meanCens = exp(-designMatr1 %*% matrix(betaCens, ncol=1))
    censTime = rexp(n = subjNum, rate = 1/meanCens)
    data$censTime = censTime
    data$event = as.numeric(data$time <= censTime)
#    data$x = data$time
    data$time[data$event == 0] = censTime[data$event == 0]
  }
  data
}

hazardStretcher_take_two = function(data, time){
  ### data - data with hazard, time  what is supplied by cox model
  ### time - is the original time variable (with duplicates and everything...)
  
  ### data = data.frame(hazard = c(0, 0, .11, .11, .11, .22, .33, .33), time = c(1.1, 1.2, 2.2, 2.22, 2.3, 3, 4, 5.5)); time = c(1.1, 1.2, 2.2, 2.22, 2.22, 2.3, 3, 4, 4, 4, 5.5)
  ### data = data.frame(hazard = c(0, 0, .11, .11, .11, .22, .33, .44), delta = c(0, 0, 1, 0, 0, 1, 1, 1), time = c(1.1, 1.2, 2, 2.2, 2.3, 3, 4, 5.5))
  ### data = data.frame(hazard = c(.11, .11, .22, .33, .44), delta = c(1, 0, 1, 1, 1), time = c(2, 2.2, 3, 4, 5.5))
  ### data = data.frame(hazard = c(.11, .11, .22, .33, .44), delta = c(1, 0, 1, 1, 1), time = c(2, 2.2, 3, 4, 5.5))
  
  ######################### important note:
  ### when this funciton is applied to coxph baseline hazard (model1), 
  ### it is important that time comes from object$y[,1] because of issue with
  ### number representation (see argument timefix in coxph.control and Terri Therneau's email)

  if(!all(c("hazard", "time") %in% names(data))){
    stop("The data has to have the following names: hazard", "time")
  }
  
  orderedTime = time[order(time)]
  data = data[order(data$time),]
  uniqueHazard = unique(data[, c("hazard", "time")])
  #uniqueHazard = rbind(c(0, 1), uniqueHazard)
  uniqueHazard = rbind(c(0, 0), uniqueHazard)
  data$Hminus = NA
  
  ##################### stretch the data:
  ##################### define hazard for every point of the original time:
  orig_i = 1
  stretchedData = data.frame(time = orderedTime, hazard = NA)

  for(i in 1:nrow(uniqueHazard)){
    while((uniqueHazard$time[i] == orderedTime[orig_i]) & (orig_i <= length(orderedTime))){
      stretchedData$hazard[orig_i] = uniqueHazard$hazard[i]
      orig_i = orig_i + 1
    }
  }
  
  i_minus = 1
  i_unique = 2
  while(stretchedData$time[i_minus] <= uniqueHazard$time[nrow(uniqueHazard)] & i_minus <= nrow(stretchedData) & i_unique < nrow(uniqueHazard) ){
    ### fill in censored observations 
    if(stretchedData$time[i_minus] < uniqueHazard$time[i_unique+1]){
      stretchedData$Hminus[i_minus] = uniqueHazard$hazard[i_unique-1]
      i_minus = i_minus + 1
    }else{
      stretchedData$Hminus[i_minus] = uniqueHazard$hazard[i_unique]
      i_minus = i_minus + 1
      i_unique = i_unique + 1
    }
    #cat(i_minus, i_unique, "\n")
  }
  
  ### if there are any observations that are censored after the last failure point.
  while(i_minus <= nrow(stretchedData)){
    stretchedData$Hminus[i_minus] = uniqueHazard$hazard[nrow(uniqueHazard) - 1]
    i_minus = i_minus + 1
  }
  
  stretchedData
}

myOwnKM = function(time, delta, returnToOriginalOrder = TRUE){
  ### returnToOriginalOrder = TRUE - the order of time values is the same as in original data
  ###            it also returns the delta (event indicator)
  ### returnToOriginalOrder = FALSE - the time is unique and ordered and the resulting data
  ###            frame does not contain delta (event indicator).
  uniqueAndOrderedTime = unique(time)[order(unique(time))]
  if(TRUE){ 
    ###-----------------------------------------------
    ### computes KM (the one from survival package was hard to use)
    ### the one from survival package returned rounded times values,
    ###    which would not allow to perform precise PSR calculations
    #dataKM = data.frame(time=c(1, 2, 3, 5, 2, 6, 6, 6), delta=c(1, 0, 0, 1, 1, 0, 1, 1))
    # dataKM = data.frame(time=time, delta=delta)
    # dataKM = dataKM[order(dataKM$time),]
    # #KM <- survfit(Surv(dataKM$time, dataKM$delta) ~ 1, type="kaplan-meier", conf.type="log")
    # nEvents = tapply(dataKM$delta, dataKM$time, sum)
    # nDrop = tapply(dataKM$delta, dataKM$time, length)
    # atRisk = c(length(dataKM$time), length(dataKM$time) - cumsum(nDrop))[1:length(nDrop)]
    # probForEachTime = (1-nEvents/atRisk)
    nEvents = tapply(delta, time, sum)
    nDrop = tapply(delta, time, length)
    atRisk = c(length(time), length(time) - cumsum(nDrop))[1:length(nDrop)]
    probForEachTime = (1-nEvents/atRisk)
    dataKM = data.frame(time=uniqueAndOrderedTime, nEvents = nEvents, atRisk = atRisk, KM = cumprod(probForEachTime))
    fit <- survfit(Surv(time, delta) ~ 1)  ### after some time this line should be removed
    ### together with the next line (temporary check for my code against R's)
    if(length(dataKM$KM) == length(fit$surv)){  ### sometimes myKM is longer than survfit
      if( max(abs(dataKM$KM - fit$surv))>0.0001) {
        problemData = list(uniqueAndOrderedTime = uniqueAndOrderedTime, df = data.frame(time = time, delta = delta), fit = fit, dataKM = dataKM)
        fileName = gsub("[ :]", "_", paste("./RESULTS/problemData", Sys.time(),".rda"))
        save(problemData, file = fileName, compress = TRUE)
      }
    }
  }else{
    ### not using this for now b/c it gave me trouble on ACCRE
    ### it rounded up some of the time values, which caused problems when "stretching" KM
    ### over time...
    fit <- survfit(Surv(time, delta) ~ 1)
    if(length(uniqueAndOrderedTime) != length(fit$n.event)){
      problemData = list(uniqueAndOrderedTime = uniqueAndOrderedTime, fitN = fit$n.event, df = data.frame(time = time, delta = delta), survFit = fit)
      save(problemData, file = "./RESULTS/problemData.rda", compress = TRUE)
      stop("Saved problematic data\n")
    }
    dataKM = data.frame(time=uniqueAndOrderedTime, nEvents = fit$n.event, atRisk = fit$n.risk, KM = fit$surv)
  }
  dataKM$CDF = 1 - dataKM$KM
  # cat(nrow(dataKM), "\n")
  # cat(paste(dataKM$CDF, collapse = ", "), "\n")
  dataKM$CDF_M = c(0, dataKM$CDF[1:(nrow(dataKM)-1)])
  if(returnToOriginalOrder){
    rownames(dataKM) = uniqueAndOrderedTime
    ### let's order the output according to the original order of time
    dataKM = dataKM[as.character(time),]
    dataKM$delta = delta
  }
  rownames(dataKM) = 1:nrow(dataKM)
  dataKM
}

survalProbabilityCond <- function(object, newdata){
  ### Check: the newdata has to be the same as data b/c
  
  ### Restricted cubic splines are not implemented
  
  ### ONLY CPH OBJECT IS IMPLEMENTED
  ### object - object returned by a surival model
  ### designMatrix - which patients I want survival functions to be computed.
  ###           it should contain the same variables as the fitted model's desing matrix
  ###           if NULL, the original data will be used
  ### newdata is the data that survival curves are computed for
  
  if(FALSE){
    data = data.frame(time = c(2,2, 4,4, 5,5, 1, 6), delta = c(1,0, 0,0, 0,1, 1, 0), z = 1:8, w = c(0.31, 0.24, 0.42, 0.34, 0.13, 0.9, 0.49, 0.86))
#    object = coxph(Surv(data$time, data$delta) ~ z + w, data = data, x = TRUE, y = TRUE)
    object = cph(Surv(data$time, data$delta) ~ z + w, data = data, x = TRUE, y = TRUE)
    newdata = data[1:4,]
    newdata = data
  }
  
  if(is.null(object$y))stop("Please fit the regression with x = TRUE y = TRUE")

  if(class(object)[1] %in% c("cph", "coxph")){
    ### compute the hazard
    # if(class(object$y) == "Surv"){
    #   time = object$y[, "time"]
    # }else{
    #   time = exp(object$y[, "time"])
    # }
    # orderedTime = time[order(time)]
    # delta = object$y[, "status"]
    # data = data.frame(time = time, delta = delta)
    
    H0 = basehaz(object, centered=FALSE)
    toStretchHaz = FALSE
    ### I do not think the hazard has to be stretcheds
    if(toStretchHaz){
      H0 = hazardStretcher_take_two(H0, orderedTime)
    }else{
      H0 = data.frame(time = H0$time, hazard = H0$hazard, Hminus = c(0, H0$hazard[1:(nrow(H0)-1)]))
    }
    dH0 = H0$hazard - H0$Hminus

    ### compute linear predictors
    varsToInclude = names(object$coefficients)
    designMatrix = as.matrix(newdata[, varsToInclude], nrow = nrow(newdata), ncol = length(varsToInclude))
    linPred = as.vector(designMatrix %*% object$coefficients)
    ### another way to obtain linear predictors
    ### BEWARE:  these linear predictors differ from the ones above by 
    ###          a constant value
    # if(FALSE){
    #   dd = datadist(data)
    #   options(datadist = "dd")
    #   linPred = predict(object, newdata = newdata, type = "lp")
    #   linPred =
    # }
    newdata = cbind(newdata, linPred)
    
    ### matrix of survival probabilities:
    ###   row is time
    ###   col is subject
    survProbMatrix = matrix(NA, nrow = nrow(H0), ncol = nrow(newdata))
    rownames(survProbMatrix) = H0$time
    for(i in 1:nrow(newdata)){
      adjHaz = dH0 * exp(linPred[i])
      adjCumHaz = cumsum(adjHaz)
      survProbMatrix[, i] = exp(-adjCumHaz)
    }
    # survProbMatrix = rbind(rep(1, ncol(survProbMatrix)), survProbMatrix)
    # rownames(survProbMatrix)[1] = "0"
    colnames(survProbMatrix) = as.character(newdata$linPred)
    
    # plot(0, 0, xlim = range(c(0,time)), ylim = c(0.3, 1))
    # for(i in 1:nrow(data)){
    #   lines(rownames(survProbMatrix), survProbMatrix[, i], type = "s")
    # }
    
    survProbMatrix
  }
}

### Not working... it used to work, but then .... just stoppe working:
###    survfit function stopped returning a matrix for each subject
survalProbabilityTherneau = function(object, newdata = NULL){
  ### object - object returned by a surival model
  ### newdata is the data that survival curves are computed for
  
  if(class(object)[1] %in% c("cph", "coxph")){
    ### compute the hazard
    if(is.null(newdata)){
      survFitCurves <- survfit(object, conf.type="plain", conf.int=0.95)
    }else{
      survFitCurves <- survfit(object, newdata=newdata, conf.type="plain", conf.int=0.95)
    }
  }
  res = list(est = survFitCurves$surv, lower = survFitCurves$lower, upper = survFitCurves$upper)
  res
}

### not in use ->  NOW USING paramBivarSurvProbWithThetasAndGammas()
paramBivarSurvProb <- function(object, newdata){
  if(FALSE){
    data = data.frame(time = c(2,2, 4,4, 5,5, 1, 6), delta = c(1,0, 0,0, 0,1, 1, 0), z = 1:8, w = c(0.31, 0.24, 0.42, 0.34, 0.13, 0.9, 0.49, 0.86))
    object = coxph(Surv(data$time, data$delta) ~ z + w, data = data, x = TRUE, y = TRUE)
    newdata = data[1:4,]
  }
  
  ### get conditional surv prob
  survSurfCondOnCovar =  survalProbabilityCond(object, newdata)
  n = nrow(survSurfCondOnCovar)

  ### sort the matrix of cond probability with in the order of increasing
  ###       value of linear predictor
  ### survSurfCondOnCovar is NOT a bivariate survival function, therefore,
  ###       survival prob. does NOT have to decrease with increasing linear predictor
  linPred = as.numeric(colnames(survSurfCondOnCovar))
  linPred = -linPred
  colnames(survSurfCondOnCovar) = as.character(linPred)
  survSurfCondOnCovar = survSurfCondOnCovar[, order(linPred)]
  
  cdfOfLinPred = myOwnKM(time = linPred, delta = rep(1, length(linPred)), returnToOriginalOrder = FALSE)$CDF
  dCDF = diff(c(0, cdfOfLinPred))
  dCDFmatrix = matrix(dCDF, ncol = ncol(survSurfCondOnCovar), nrow = n, byrow = TRUE)

  S_ti_dw1 = -survSurfCondOnCovar[,1] * cdfOfLinPred[1]

  ### marginal S(t):
  margSt = apply(survSurfCondOnCovar*dCDFmatrix, 1, sum)
  
  ### Compute the joint survival probability
  ### the rows are times, the columns are linear predictors
  jointSurvProb = matrix(NA, nrow = nrow(survSurfCondOnCovar), ncol = ncol(survSurfCondOnCovar))

  jointSurvProb[ , 1] = S_ti_dw1 + margSt
  
  # ### the last column is zero because the prob of being greater than the last
  # ###   always uncensored linear predictor is zero (because it is the maximum)
  # colNum = ncol(jointSurvProb)
  # jointSurvProb[, colNum] = 0
  # jointSurvProb[, colNum - 1] = survSurfCondOnCovar[, (colNum - 1)]*dKM_matrix[, (colNum - 1)]
  # if(colNum > 2){
    for(j in 2:ncol(jointSurvProb)){
      S_t_dwj = -survSurfCondOnCovar[,j] * dCDFmatrix[, j]
      jointSurvProb[, j] = S_t_dwj  + jointSurvProb[, j-1]
      # jointSurvProb[, j] = apply(survSurfCondOnCovar[, j:(colNum - 1)]*dKM_matrix[, j:(colNum - 1)], 1, sum)
    }
  # }
  
  ### add marginal survival probabilities:
  jointSurvProb = cbind(margSt, jointSurvProb)
  jointSurvProb = rbind(c(1, 1 - cdfOfLinPred), jointSurvProb)
  colnames(jointSurvProb) = c("-Inf", colnames(survSurfCondOnCovar))
  rownames(jointSurvProb) = c("0", rownames(survSurfCondOnCovar))
  jointSurvProb
}

paramBivarSurvProbWithThetasAndGammas <- function(object, newdata){
  if(FALSE){
    data = data.frame(time = c(2,2, 4,4, 5,5, 1, 6), delta = c(1,0, 0,0, 0,1, 1, 0), z = 1:8, w = c(0.31, 0.24, 0.42, 0.34, 0.13, 0.9, 0.49, 0.86))
    object = coxph(Surv(data$time, data$delta) ~ z + w, data = data, x = TRUE, y = TRUE)
    newdata = data[1:4,]
  }
  
  ### get conditional surv prob
  survSurfCondOnCovar =  survalProbabilityCond(object, newdata)
  n = nrow(survSurfCondOnCovar)

  ### sort the matrix of cond probability with in the order of increasing
  ###       value of linear predictor
  ### survSurfCondOnCovar is NOT a bivariate survival function, therefore,
  ###       survival prob. does NOT have to decrease with increasing linear predictor
  linPred = as.numeric(colnames(survSurfCondOnCovar))
  linPred = -linPred
  colnames(survSurfCondOnCovar) = as.character(linPred)
  ###       Order the conditional survival probability by the linear predictor
  survSurfCondOnCovar = survSurfCondOnCovar[, order(linPred)]
  ### get rid of columns generated by duplicated values of linear predictors:
  survCondUnique = t(unique(t(survSurfCondOnCovar)))
    
  ###       Compute the ordered marginal survival probability of the linear predictor
  cdfOfLinPred = myOwnKM(time = linPred, delta = rep(1, length(linPred)), returnToOriginalOrder = FALSE)$CDF
  dCDF = diff(c(0, cdfOfLinPred))
  dCDFmatrix = matrix(dCDF, ncol = ncol(survCondUnique), nrow = n, byrow = TRUE)

  ########## The following quantities are unordered:
  ########## Define gamma_{ij} and theta_{j} so you could do M-estimation and compute
  ########## the C-index without building a survival surface
  thetaj = dCDFmatrix
  gammaij = survCondUnique
  gamma0ij = rbind(rep(1, ncol(gammaij)), gammaij)
  gamma0ij = gamma0ij[1:(nrow(gamma0ij)-1), ]
  
  ########## St_dw = gamma_ij * thetaj
  St_dw = -gammaij * thetaj

  ########## Sdt_dw = (gamma_{i-1,j} - gamma_{ij}) * thetaj
  Sdt_dw = (gamma0ij - gammaij) * thetaj

  ########## Stw = sum_1^n(gamma_{i,j}*thetaj) - sum_1^j(gamma_{i,j}*thetaj)  
  #              =  sum_{j+1}^n(gamma_{i,j}*thetaj)
  Stw = thetaj * 0
  for(i in 1:nrow(Stw)){
    Stw[i, ] = sum(gammaij[i, ]*thetaj[i, ])  - cumsum(gammaij[i, ]*thetaj[i, ])
  }
  
  ########## I've checked that for continuous time Sdtw it is good.
  Sdt_w = thetaj * 0
  for(i in 1:nrow(Sdt_w)){
    Sdt_w[i, ] = sum(gammaij[i, ]*thetaj[i, ])  - cumsum(gammaij[i, ]*thetaj[i, ])   -    (sum(gamma0ij[i, ]*thetaj[i, ])  - cumsum(gamma0ij[i, ]*thetaj[i, ]))
  }  
  
  # jointSurvProb = paramBivarSurvProb(object, newdata)
  # r = cIndexProbOfConcAndDisc(newdata, jointSurvProb)
  
  # ### We need to check these values:
  # Sdxdy = r$someOther$Sdxdy
  # Sxy = r$someOther$Sxy
  # Sdx_y = -r$someOther$Sdx_y
  # Sx_dy = -r$someOther$Sx_dy
  
  probOfConc = sum(Sdt_dw * Stw)
  probOfDisc = sum(Sdt_w * St_dw)
  probOfConcByT = apply(Sdt_dw * Stw, 1, sum)
  probOfDiscByT = apply(Sdt_w * St_dw, 1, sum)
  
  cInd = probOfConc / (probOfConc + probOfDisc)
  cIndByT = data.frame(time = as.numeric(rownames(survSurfCondOnCovar)), cInd = probOfConcByT / (probOfConcByT + probOfDiscByT))
  
  ### marginal S(t):
  #margSt = apply(survSurfCondOnCovar*dCDFmatrix, 1, sum)
  margSt = apply(-St_dw, 1, sum)
  margStMatrix = matrix(margSt, nrow = nrow(Stw), ncol = ncol(Stw), byrow = FALSE)
  margSdtMatrix = matrix(diff(c(1, margSt)), nrow = nrow(Stw), ncol = ncol(Stw), byrow = FALSE)
  
  ### compute incident sensitivity:
  ### S(dt,w)/(-margSt(dt))
  incSens = Sdt_w/margSdtMatrix
  colnames(incSens) = c(colnames(survCondUnique))
  rownames(incSens) = c(rownames(survCondUnique))
  
  # ### making Survival surface
  # S_ti_dw1 = -survSurfCondOnCovar[,1] * cdfOfLinPred[1]
  #
  # ### Compute the joint survival probability
  # ### the rows are times, the columns are linear predictors
  # jointSurvProb = matrix(NA, nrow = nrow(survSurfCondOnCovar), ncol = ncol(survSurfCondOnCovar))
  #
  # jointSurvProb[ , 1] = S_ti_dw1 + margSt
  # # ### the last column is zero because the prob of being greater than the last
  # # ###   always uncensored linear predictor is zero (because it is the maximum)
  #   for(j in 2:ncol(jointSurvProb)){
  #     S_t_dwj = -survSurfCondOnCovar[,j] * dCDFmatrix[, j]
  #     jointSurvProb[, j] = S_t_dwj  + jointSurvProb[, j-1]
  #   }
  #
  # ### add marginal survival probabilities:
  # jointSurvProb = cbind(margSt, jointSurvProb)
  # jointSurvProb = rbind(c(1, 1 - cdfOfLinPred), jointSurvProb)
  # colnames(jointSurvProb) = c("-Inf", colnames(survSurfCondOnCovar))
  # rownames(jointSurvProb) = c("0", rownames(survSurfCondOnCovar))
  # jointSurvProb
  
  Stw_WithMarg = cbind(margSt, Stw)
  Stw_WithMarg = rbind(c(1, 1 - cdfOfLinPred), Stw_WithMarg)
  colnames(Stw_WithMarg) = c("-Inf", colnames(survCondUnique))
  rownames(Stw_WithMarg) = c("0", rownames(survCondUnique))

  list(incSens = incSens, cInd = cInd, cIndByT = cIndByT, Stw = Stw_WithMarg)  
}

checkDerivationsForMEstimation <- function(){
  ### Discrete time case:
  data = data.frame(time = c(2,2, 4,4, 5,5, 1, 6), delta = c(1,0, 0,0, 0,1, 1, 0), z = 1:8, w = c(0.31, 0.24, 0.42, 0.34, 0.13, 0.9, 0.49, 0.86))
  ### Discrete time and linear predictor case:
  # data = data.frame(time = c(2,2, 4,4, 5,5, 1, 6), delta = c(1,0, 0,0, 0,1, 1, 0), z = 1:8, w = c(0.31, 0.24, 0.42, 0.34, 0.13, 0.9, 0.49, 0.86))
  ### Continuous time example:
  # data = data.frame(time = c(2,2.1, 4,4.1, 5,5.1, 1, 6), delta = c(1,0, 0,0, 0,1, 1, 0), z = 1:8, w = c(0.31, 0.24, 0.42, 0.34, 0.13, 0.9, 0.49, 0.86))
#    object = coxph(Surv(data$time, data$delta) ~ z + w, data = data, x = TRUE, y = TRUE)
  object = cph(Surv(data$time, data$delta) ~ z + w, data = data, x = TRUE, y = TRUE)
  # newdata = data[1:4,]
  newdata = data
  
  gamma_ij =  survalProbabilityCond(object, newdata)
  
  KM_W = myOwnKM(time = data$w, delta = rep(1, nrow(data)), returnToOriginalOrder = FALSE)$KM
  theta_j = matrix(-diff(c(1, KM_W)), nrow = nrow(gamma_ij), ncol = ncol(gamma_ij), byrow = TRUE)
  
  jointSurvProb = paramBivarSurvProb(object, newdata)
  
  r = cIndexProbOfConcAndDisc(newdata, jointSurvProb)$someOther
  
  ### We need to check these values:
  
  Sdtdw = r$someOther$Sdxdy
  Stw = r$someOther$Sxy
  Sdt_w = r$someOther$Sdx_y
  St_dw = r$someOther$Sx_dy
  
  
  Sti_w0 = apply(gamma_ij/8, 1, sum)
  #Sdti_w0 = apply(gamma_ij/8 - gamma_ijWithZero[1:(nrow(gamma_ijWithZero)-1), ]/8, 1, sum)
  
  ### I've checked that for continuous time St_dw = gamma_ij / n
  # > range(St_dw - gamma_ij/8)
  # [1] -8.326673e-17  8.326673e-17
  
  ### Discrete
  range(St_dw + gamma_ij * theta_j)
  
  ### Compute Kaplan-Meier for 
  
  ### I've checked that for continuous time Sdtdw = (gamma_{i-1,j} - gamma_{ij}) / n
  # gamma_ijWithZero = rbind(rep(1, ncol(gamma_ij)), gamma_ij)
  # gamma_ijWithZero = gamma_ijWithZero[1:(nrow(gamma_ijWithZero)-1), ]
  # >   tmp = gamma_ijWithZero - gamma_ij
  # >   range(Sdtdw - tmp/8)
  # [1] -6.938894e-17  6.245005e-17
  
  ### I've checked that for continuous time Stw = sum_{j+1}^n(gamma_{i,j}*theta_j) / n
  # tmp1 = gamma_ij*0
  # for(i in 1:nrow(tmp1)){
  #   tmp1[i, ] = sum(gamma_ij[i, ]/8)  - cumsum(gamma_ij[i, ]/8)
  # }
  # range(Stw - tmp1)

  ### I've checked that for continuous time Sdtw it is good.
  # tmp = gamma_ij*0
  # for(i in 1:nrow(tmp)){
  #   tmp[i, ] = sum(gamma_ij[i, ]/8)  - cumsum(gamma_ij[i, ]/8)   -    (sum(gamma_ijWithZero[i, ]/8)  - cumsum(gamma_ijWithZero[i, ]/8))
  # }
  # range(Sdt_w[1:nrow(tmp)] + tmp[1:nrow(tmp)])

  ### NEXT: check for discrete time, and discrete linear predictor
  
}

paramCIndex = function(object, data){
  if(FALSE){
    data = data.frame(time = c(2,2, 4,4, 5,5, 1, 6), delta = c(1,0, 0,0, 0,1, 1, 0), z = 1:8, w = c(0.31, 0.24, 0.42, 0.34, 0.13, 0.9, 0.49, 0.86))
    object = coxph(Surv(data$time, data$delta) ~ z + w, data = data, x = TRUE, y = TRUE)
  }
  
  # linPred = predict(object, data, "lp")
  # bivarData = data.frame(X = as.matrix(object$y)[, "time"], delta = as.matrix(object$y)[, "status"], Y = exp(-linPred), epsilon = rep(1, length(linPred)))
  # dabrSurf = SurvSurfaceOfDabrowska_RecursiveEvenFaster(bivarData)$DabrowskaEst
  # tmp = cIndexProbOfConcAndDisc(data, dabrSurf)$res
  
  survSurf = paramBivarSurvProb(object, newdata = data)
  survSurf = rbind(survSurf, rep(0, ncol(survSurf)))
  rownames(survSurf)[nrow(survSurf)] = "Inf"
  
  cIndexProbOfConcAndDisc(data, survSurf)$res
}

### likely correct:
censDepOnOneZ = function(subjNum, betaEvent, betaCens = NULL, indepCens = FALSE){
  ### it is not a causal association, but this is temporary
  ### betaEvent = c(2, 3) and betaCens = c(3, 4) gives about 38% censored
  
  if(!indepCens){
    if(length(betaEvent) != length(betaCens)){
      stop("For dependent censoring (indepCens = FALSE) arguments betaEvent and betaCens should be the same length")
    }
  }
  
  designMatr1 = matrix(1, ncol=length(betaEvent), nrow = subjNum)
  colnames(designMatr1) = paste0("z", 1:length(betaEvent))
  for(i in 1:length(betaEvent)){
    designMatr1[, i] = rnorm(subjNum)
  }
  
  meanEvent = exp(-designMatr1 %*% matrix(betaEvent, ncol=1))
  time = rexp(n = subjNum, rate = 1/meanEvent)
  # data = data.frame(time, censTime = Inf, x = time, event = 1, Z = designMatr1[,1])
  data = data.frame(time, censTime = Inf, x = time, event = 1)
  data = cbind(data, designMatr1)

  if(!is.null(betaCens)){
    if(indepCens){
      designMatr2 = matrix(1, ncol=length(betaCens), nrow = subjNum)
      for(i in 1:length(betaCens)){
        designMatr2[,i] = rnorm(subjNum)
      }
      meanCens = exp(-designMatr2 %*% matrix(betaCens, ncol=1))
    }else{
      meanCens = exp(-designMatr1 %*% matrix(betaCens, ncol=1))
    }
    data$censTime = rexp(n = subjNum, rate = 1/meanCens)
    data$event = as.numeric(data$time < data$censTime)
    data$x[data$event == 0] = data$censTime[data$event == 0]
  }
  data
}

answer20200615 = function(simNum = 100, subjNum = 80, indepCens = FALSE){
  #simNum = 100; subjNum = 80
  #set.seed(20190708)
  seeds = sample(1:simNum*1000, simNum) 
  labelsForStuff = c("Harrell", "Uno", "Our and Exp", "Our and Cox", "Out semipar")
  resMatr = matrix(NA, ncol = length(labelsForStuff) + 1, nrow = simNum)
  colnames(resMatr) = c("noCensValue", labelsForStuff)
  betaEvent = betaCens = 2
  for(i in 1:simNum){
    set.seed(seeds[i])
    # betaEvent = runif(2, min = -4, max = 4)
    # betaCens = runif(2, min = -4, max = 4)
    # betaEvent = runif(2, min = -4, max = -2)
    # betaCens = runif(2, min = 2, max = 4)
    betaEvent = runif(2, min = 2, max = 4)
    betaCens = runif(2, min = -4, max = -2)
    # betaEvent = runif(2, min = c(-2, 0), max = c(0, 2))
    # betaCens = runif(2, min = c(0, -2), max = c(2, 0))
    data = censDepOnOneZ(subjNum = subjNum, betaEvent = betaEvent, betaCens = betaCens, indepCens = TRUE)
    survObj0 = Surv(data$time, data$event)
    survObj1 = Surv(data$x, data$event)
    
    dd = datadist(data)
    options(datadist = "dd")
    # covars = c("z1", "z2")
    covars = c("z1")
    m0 = cph(as.formula(paste("survObj0 ~", paste(covars, collapse = " + "))), data = data, x = TRUE, y = TRUE)
    m1 = cph(as.formula(paste("survObj1 ~", paste(covars, collapse = " + "))), data = data, x = TRUE, y = TRUE)
    m2 = survreg(as.formula(paste("survObj1 ~", paste(covars, collapse = " + "))), data = data, dist = "exponential")
    # m0 = cph(survObj0 ~ z1, data = data, x = TRUE, y = TRUE)
    # m1 = cph(survObj1 ~ z1, data = data, x = TRUE, y = TRUE)
    # m2 = survreg(survObj1 ~ z1, data = data, dist = "exponential")

    ###########################################################
    ### Harrell, no censoring
    resMatr[i, "noCensValue"] = survConcordance(survObj0 ~ predict(m0))[["concordance"]]
  
    ###########################################################
    ### Harrell, with censoring
    resMatr[i, "Harrell"] = survConcordance(survObj1 ~ predict(m1))[["concordance"]]
  
    ###########################################################
    ### Uno & Pencina, ...
    try({
      unotmp = Inf.Cval(data[, c("time", "event", covars)],  tau = 1000, itr = 10);
      resMatr[i, "Uno"] = unotmp$Dhat
    }, silent = TRUE)
    
    ###########################################################
    ### Eden using Dabrowska
    desMatr = predict(m1, newdata = data, type = "x")
    cumHaz = basehaz(m1, centered = FALSE)
    regHaz = diff(c(0, cumHaz$hazard))
    names(regHaz) = cumHaz$time
    linPred = as.vector(desMatr %*% m1$coefficients)
    names(linPred) = data$time
    adjHaz = regHaz[as.character(data$time)] * exp(linPred[as.character(data$time)])
    adjCumHaz = cumsum(adjHaz)
    survProb = exp(-adjCumHaz)
  
    bivarData = data.frame(X = data$x, delta = data$event, Y = exp(-linPred), epsilon = rep(1, subjNum))
    dabrSurf = SurvSurfaceOfDabrowska_RecursiveEvenFaster(bivarData)$DabrowskaEst
    tmp = cIndexProbOfConcAndDisc(data, dabrSurf)$res
    # concord = tmp["ConcordProb"]
    # discord = tmp["DiscordProb"]
    myCIndex = tmp["CIndex"]
    resMatr[i, "Our and Cox"] = c(myCIndex)
    
    
    ###########################################################
    ### Eden using Cox model
    resMatr[i, "Out semipar"] = paramCIndex(m1, data)["CIndex"]
 
    # ###########################################################
    # ### Eden using Exponential model
    # tmpExp = predict(m2, type = "response")
    # bivarData = data.frame(X = data$x, delta = data$event, Y = tmpExp, epsilon = rep(1, subjNum))
    # dabrSurf = SurvSurfaceOfDabrowska_RecursiveEvenFaster(bivarData)$DabrowskaEst
    # tmp = cIndexProbOfConcAndDisc(data, dabrSurf)$res
    # concord = tmp["ConcordProb"]
    # discord = tmp["DiscordProb"]
    # myCIndex = tmp["CIndex"]
    #
    # # resMatr[i, ] = c(cIndex2, myCIndex, concord)
    # resMatr[i, "Our and Exp"] = c(myCIndex)
  }
  
  colToPrint = c(1:4,5)
  whatToPrint = c("Harrell", "Uno", "Our and Cox", "Out semipar")
  par(mfrow = c(length(whatToPrint), 1), mar = c(3, 0, 0, 0))
#  xlim = range(resMatr, na.rm = TRUE)
  # xlim = c(0.5, 1)
  xlim = c(0.2, 1)
  #for(i in colToPrint){
  for(n_i in whatToPrint){
    #plot(0, 0, type = "n", xlim = xlim)
    hist(resMatr[, n_i], xlim = xlim, breaks = 25, col = "gray", border = "white", xlab = n_i)
    abline(v = c( mean(resMatr[, n_i], na.rm = TRUE),  mean(resMatr[, "noCensValue"], na.rm = TRUE)), col = c("black", "red" ))
    #lines(x = mean(resMatr[, i]) + c(-3, 3)*sqrt(var(resMatr[, i])), y = c(0, 0), col = "green")
    text(mean(resMatr[, "noCensValue"], na.rm = TRUE), 0, paste0(n_i, ", ", round(sd(resMatr[, n_i], na.rm = TRUE), 4)), pos = 4, srt = 45)
  }
  
  par(mfrow = c(1, length(whatToPrint)), mar = c(3, 0, 0, 0))
#  xlim = range(resMatr, na.rm = TRUE)
  xlim = c(0.5, 1)
  #for(n_i in c("Uno", "Our and Cox")){
  for(n_i in whatToPrint){
    plot(0, 0, xlim = c(0, 1), ylim = c(0, 1), ylab = n_i, xlab = "noCensValue")
    points(resMatr[, "noCensValue"], resMatr[, n_i])
    abline(a = 0, b = 1)
    # text(mean(resMatr[, "noCensValue"], na.rm = TRUE), 0, paste0(colnames(resMatr)[i], ", ", round(sd(resMatr[, i], na.rm = TRUE), 4)), pos = 4, srt = 45)
  }
  
  plot(resMatr[, "Uno"], resMatr[, "Our and Cox"], xlim = c(0, 1), ylim = c(0, 1))
  abline(a = 0, b = 1)
  
}

plotResults = function(resMatr, trueCIndex, pch, col, cex, ...){
  plot(range(dim(resMatr)), range(resMatr), type = "n", ...)
  for(j in colnames(resMatr)){
    points(1:nrow(resMatr), resMatr[, j], pch = pch[j], col = col[j], cex = cex[j])
    points(x = 1:nrow(resMatr), y = rep(mean(resMatr[, j]), nrow(resMatr)), col = col[j], pch = pch[j], cex = cex[j])
  }
  abline(h = trueCIndex, col = "black", lty = 5)
  legend(x = "bottomleft", pch = pch, bg = "white", pt.cex = cex, col = col, legend = colnames(resMatr))
  legend(x = "bottomright", lty = 2, bg = "white", pt.cex = cex, col = "black", legend = "True value")
  list(mean = apply(resMatr, 2, mean), var = apply(resMatr, 2, var))
}

reproduceHeagertyAndZheng <- function(object){
  data(ovarian)
  simdata = ovarian[1:10, ]
  simdata$fustat = 1
  object = coxph(Surv(futime, fustat) ~ resid.ds + rx + ecog.ps, simdata)
  res  = paramBivarSurvProbWithThetasAndGammas(object, simdata)
  plot(res$cIndByT$time/365.25, res$cIndByT$cInd)
  abline(h = res$cInd, col = "red")
  
  ### Their incidence sensitivity:
  newdata = simdata
  varsToInclude = names(object$coefficients)
  designMatrix = as.matrix(newdata[, varsToInclude], nrow = nrow(newdata), ncol = length(varsToInclude))
  linPred = as.vector(designMatrix %*% object$coefficients)
  linPredMatr = matrix(linPred, nrow = length(time), ncol = length(linPred), byrow = TRUE)
  uniqueLinPred = unique(linPred)
  
  ### time to event or censoring:
  time = object$y[,1]
  ### Compute W(t) = sum(  1(xi>=t) exp(beta zi)  )
  
  Rjt = matrix(0, nrow = length(time), ncol = length(time))
  for(j in 1:length(linPred)){
    Rjt[, j] = as.numeric(time[j] >= time)
  }
  
  Wt = apply(Rjt * exp(linPredMatr), 1, sum)
  WtMatr = matrix(Wt, nrow = length(Wt), ncol = length(linPred), byrow = FALSE)
  PIjt = Rjt * exp(linPredMatr)/WtMatr
  
  indMkGreaterThanC = matrix(NA, nrow = length(Wt), ncol = length(linPred))
  
  
  incSenHZ = matrix(NA, nrow = length(Wt), ncol = length(linPred))
  vecOfC = c(linPred, Inf)
  for(i in 1:length(Wt)){
    for(j in 1:length(linPred)){
      indMxGreaterThanC = as.numeric(linPred > vecOfC[j+1])
      incSenHZ[i, j] = sum(indMxGreaterThanC   *  PIjt[i, ])
    }
  }
  incSenHZ = incSenHZ[order(time), ]
  incSenHZ = incSenHZ[, order(-linPred)]
  
}

#################### Uno, Pencina ####################
setwd("/Users/svetlanaeden/stuff/StatStuff/BRYAN/Simulations/PAPER3")
source("../tools.R")
source("../simFunctions.R")
library(survC1)   ### for Uno's C-index
library(rms)

#load("../data/coxdata.rda")
# missing = apply(coxdata[, c("timetoevent", "coxevent", "age", "bmi", "cd4", "viralload")], 1, function(x){any(is.na(x))})
# coxdata = coxdata[!missing,]

# survObj = Surv(res34$X, res34$delta)
# mod = coxph(survObj ~ Y, data = res34)
# risk = plogis(predict(mod))
#
# ### that huno guy (and pencina and d'agostino)
# conc(X = res34$X, D = res34$delta, W=rep(1,nrow(res34)), R=risk)
# #Inf.Cval(coxdata[, c("timetoevent", "coxevent", "age", "bmi", "cd4", "viralload")],  tau = max(coxdata$timetoevent), itr = 10)
# unotmp = Inf.Cval(res34[, c("X", "delta", "Y")],  tau = 1000, itr = 10)

labelsForStuff = c("Harrell", "Uno", "Our estimator")
# col = c("#FF000044", "#00FF0044", "#0000FF44"); names(col) = c("Uno", "Harrell", "Eden")
col = c("blue", "red", "green"); names(col) = labelsForStuff
# pch = c(23, 1, 4); names(pch) = labelsForStuff
pch = c(2, 1, 4); names(pch) = labelsForStuff
lwd = c(12, 7, 1); names(lwd) = labelsForStuff
cex = c(0.6, 1.2, 0.8); names(cex) = labelsForStuff

################## no tau, random censoring
################## no tau, random censoring
################## no tau, random censoring
simNum = 100
#set.seed(20190708)
seeds = sample(1:1000, simNum)
sampleS = 80
resMatr = matrix(NA, ncol = 4, nrow = simNum)
colnames(resMatr) = c("noCensValue", labelsForStuff)
for(i in 1:simNum){
  set.seed(seeds[i])
  res12 = simulateRealisticData(nSubj = sampleS, thetaPar = 4.426265, family = "frank", censoringProb1 = 0, censoringProb2 = 0, independentCensoring = TRUE, restrictedTimeX = Inf, restrictedTimeY = Inf)
  set.seed(seeds[i])
  res34 = simulateRealisticData(nSubj = sampleS, thetaPar = 4.426265, family = "frank", censoringProb1 = 0.5, censoringProb2 = 0, independentCensoring = TRUE, restrictedTimeX = Inf, restrictedTimeY = Inf)
  names(res12) = gsub("delta$", "deltaX", names(res12)); names(res12) = gsub("epsilon$", "deltaY", names(res12))
  names(res34) = gsub("delta$", "deltaX", names(res34)); names(res34) = gsub("epsilon$", "deltaY", names(res34))
  data0 = data.frame(time = res12$X, event = res12$deltaX, V1 = res12$Y)
  data1 = data.frame(time = res34$X, event = res34$deltaX, V1 = res34$Y)
  survObj0 = Surv(data0$time, data0$event)
  survObj1 = Surv(data1$time, data1$event)

  ###########################################################
  ### Uno & Pencina, ...
  mod = coxph(survObj1 ~ V1, data = data1)
  risk = plogis(predict(mod))
  unotmp = Inf.Cval(data1[, c("time", "event", "V1")],  tau = 1000, itr = 10)
  resMatr[i, "Uno"] = unotmp$Dhat

  ###########################################################
  ### Harrell ...
  dd = datadist(data0)
  options(datadist = "dd")
  m0 = cph(survObj0 ~ V1, data = data0, x = TRUE, y = TRUE)
  sum.surv <- summary(m0)
  resMatr[i, "noCensValue"] = survConcordance(survObj0 ~ predict(m0))[["concordance"]]
  
  dd = datadist(data1)
  options(datadist = "dd")
  m1 = cph(survObj1 ~ V1, data = data1, x = TRUE, y = TRUE)
  sum.surv <- summary(m1)
  resMatr[i, "Harrell"] = survConcordance(survObj1 ~ predict(m1))[["concordance"]]
  
  ###########################################################
  ### Eden ...
  desMatr = predict(m1, newdata = data1, type = "x")
  cumHaz = basehaz(m1, centered = FALSE)
  regHaz = diff(c(0, cumHaz$hazard))
  names(regHaz) = cumHaz$time
  linPred = as.vector(desMatr %*% m1$coefficients)
  names(linPred) = data1$time
  adjHaz = regHaz[as.character(data1$time)] * exp(linPred[as.character(data1$time)])
  adjCumHaz = cumsum(adjHaz)
  survProb = exp(-adjCumHaz)
  
  bivarData = data.frame(X = data1$time, delta = data1$event, Y = exp(-linPred), epsilon = rep(1, sampleS))
  dabrSurf = SurvSurfaceOfDabrowska_RecursiveEvenFaster(bivarData)$DabrowskaEst
  tmp = cIndexProbOfConcAndDisc(data1, dabrSurf)$res
  concord = tmp["ConcordProb"]
  discord = tmp["DiscordProb"]
  myCIndex = tmp["CIndex"]

  # resMatr[i, ] = c(cIndex2, myCIndex, concord)
  resMatr[i, "Our estimator"] = c(myCIndex)
}

resMatr[, "Our estimator"][resMatr[, 3]>1] = 1

trueCIndex = 0.7089829
plotResults(resMatr[,2:4], trueCIndex = trueCIndex, pch = pch, col = col, cex = cex)
plot(range(resMatr), range(resMatr), type = "n")

plot(range(resMatr), range(resMatr), type = "n")
points(resMatr[,"noCensValue"], resMatr[, "Uno"], col = "blue", cex  = 0.7)
points(resMatr[,"noCensValue"], resMatr[, "Harrell"], col = "red", cex  = 0.7)
points(resMatr[,"noCensValue"], resMatr[, "Our estimator"], col = "green", cex  = 0.7)
abline(a = 0, b = 1)

par(mfrow = c(4, 1))
hist(resMatr[,"noCensValue"], xlim = c(range(resMatr)),  col = "gray", breaks = 25)
hist(resMatr[,"Uno"], xlim = c(range(resMatr)),  col = "blue", breaks = 25)
hist(resMatr[,"Harrell"], xlim = c(range(resMatr)),  col = "red", breaks = 25)
hist(resMatr[,"Our estimator"], xlim = c(range(resMatr)),  col = "green", breaks = 25)



##################### with tau, type-I censoring
##################### with tau, type-I censoring
##################### with tau, type-I censoring
simNum = 100
#set.seed(20190708)
sampleS = 80
resMatr = matrix(NA, ncol = 3, nrow = simNum)
colnames(resMatr) = labelsForStuff
tau075 = qexp(0.5)
trueCIndexTau075 = 0.7413919 ### for 0.5
#################### main loop
for(i in 1:simNum){
  res34 = simulateRealisticData(nSubj = sampleS, thetaPar = 4.426265, family = "frank", censoringProb1 = 0.3, censoringProb2 = 0, restrictedTimeX = tau075, restrictedTimeY = Inf)
  names(res34) = gsub("delta$", "deltaX", names(res34)); names(res34) = gsub("epsilon$", "deltaY", names(res34))
  data1 = data.frame(time = res34$X, event = res34$deltaX, V1 = res34$Y)
  onlyLessThanTau = data1$time < tau075

  ###########################################################
  ### Uno & Pencina, ...
  survObj = Surv(data1$time, data1$event)
  mod = coxph(survObj ~ V1, data = data1)
  risk = plogis(predict(mod))
  unotmp = Inf.Cval(data1[, c("time", "event", "V1")],  tau = tau075, itr = 10)
  resMatr[i, "Uno"] = unotmp$Dhat

  ###########################################################
  ### Harrell ...
  dd = datadist(data1)
  options(datadist = "dd")
  survObj = Surv(data1$time, data1$event)
  m1 = cph(survObj ~ V1, data = data1, x = TRUE, y = TRUE)
  sum.surv <- summary(m1)
  resMatr[i, "Harrell"] = survConcordance(survObj[onlyLessThanTau] ~ predict(m1)[onlyLessThanTau])[["concordance"]]
  
  ###########################################################
  ### Eden ...
  desMatr = predict(m1, newdata = data1, type = "x")
  cumHaz = basehaz(m1, centered = FALSE)
  regHaz = diff(c(0, cumHaz$hazard))
  names(regHaz) = cumHaz$time
  linPred = as.vector(desMatr %*% m1$coefficients)
  names(linPred) = data1$time
  adjHaz = regHaz[as.character(data1$time)] * exp(linPred[as.character(data1$time)])
  adjCumHaz = cumsum(adjHaz)
  survProb = exp(-adjCumHaz)
  
  medSurvTime = survfit(m1, newdata= data1)$median
  bivarData = data.frame(X = data1$time, delta = data1$event, Y = exp(-linPred), epsilon = rep(1, sampleS))
  dabrSurf = SurvSurfaceOfDabrowska_RecursiveEvenFaster(bivarData)$DabrowskaEst
  tmp = cIndexProbOfConcAndDisc(bivarData, dabrSurf)$res
  concord = tmp["ConcordProb"]
  discord = tmp["DiscordProb"]
  myCIndex = tmp["CIndex"]

  # resMatr[i, ] = c(cIndex2, myCIndex, concord)
  resMatr[i, "Our estimator"] = c(myCIndex)
}

resMatr[, "Our estimator"][resMatr[, 3]>1] = 1
plotResults(resMatr, trueCIndex = trueCIndexTau075, pch = pch, col = col, cex = cex)

##################### censoring distribution depends on Z
##################### censoring distribution depends on Z
##################### censoring distribution depends on Z
simNum = 100
#set.seed(20190708)
sampleS = 160
resMatr = matrix(NA, ncol = 3, nrow = simNum)
colnames(resMatr) = labelsForStuff
tau075 = qexp(0.5)
trueCIndexTau075 = 0.8586082 ### for 0.5
#################### main loop
for(i in 1:simNum){
  # data1 = simSurvDataWithCensDepOnZ_take_two(subjNum = sampleS, betaEvent = c(2, 3), betaCens = c(1, 4), theta = 4.426265, family = "frank")
  data1 = simSurvDataWithCensDepOnZ_take_two(subjNum = sampleS, betaEvent = c(2, 3), betaCens = c(3, 4), theta = 4.426265, family = "frank")
  # data1 = simSurvDataWithCensDepOnZ_take_two(subjNum = sampleS, betaEvent = c(2, 3), betaCens = c(3, -1), theta = 4.426265, family = "frank")
  # data1 = simSurvDataWithCensDepOnZ_take_two(subjNum = sampleS, betaEvent = c(2, 3), betaCens = c(-1, -1), theta = 4.426265, family = "frank")
  # data1 = simSurvDataWithCensDepOnZ_take_two(subjNum = sampleS, betaEvent = c(2, 3), betaCens = NULL, theta = 4.426265, family = "frank")

  data1$censEvent = 1 - data1$event
  ###########################################################
  ### Uno & Pencina, ...
  survObj = Surv(data1$time, data1$event)
  mod = coxph(survObj ~ V1 + V2, data = data1)
  risk = plogis(predict(mod))
  # resMatr[i, "Uno"] = conc(X = data1$time, D = data1$event, W=rep(1,nrow(data1)),R=risk)
  unotmp = Inf.Cval(data1[, c("time", "event", "V1", "V2")],  tau = 1000, itr = 10)
  resMatr[i, "Uno"] = unotmp$Dhat

  ###########################################################
  ### Harrell ...
  dd = datadist(data1)
  options(datadist = "dd")
  survObj = Surv(data1$time, data1$event)
  m1 = cph(survObj ~ V1 + V2, data = data1, x = TRUE, y = TRUE)
  sum.surv <- summary(m1)
  resMatr[i, "Harrell"] = survConcordance(survObj ~ predict(m1))[["concordance"]]
  
  ###########################################################
  ### Eden ...
  if(FALSE){
    modForCens = cph(Surv(data1$time, data1$censEvent) ~ V1 + V2, data = data1, x = TRUE, y = TRUE)
    desMatr = predict(modForCens, newdata = data1, type = "x")
    cumHaz = basehaz(modForCens, centered = FALSE)
    regHaz = diff(c(0, cumHaz$hazard))
    names(regHaz) = cumHaz$time  
    linPred = as.vector(desMatr %*% modForCens$coefficients)
    names(linPred) = data1$time
    adjHaz = regHaz[as.character(data1$time)] * exp(linPred[as.character(data1$time)])
    adjCumHaz = cumsum(adjHaz)
    survProbForCens = exp(-adjCumHaz)
    m3 = cph(survObj ~ V1 + V2, data = data1, weights = 1/survProbForCens, x = TRUE, y = TRUE)
  }else{
    m3 = cph(survObj ~ V1 + V2, data = data1, x = TRUE, y = TRUE)
  }
  desMatr = predict(m3, newdata = data1, type = "x")
  cumHaz = basehaz(m3, centered = FALSE)
  regHaz = diff(c(0, cumHaz$hazard))
  names(regHaz) = cumHaz$time
  linPred = as.vector(desMatr %*% m3$coefficients)
  #linPred = as.vector(desMatr %*% c(2,3))
  names(linPred) = data1$time
  adjHaz = regHaz[as.character(data1$time)] * exp(linPred[as.character(data1$time)])
  adjCumHaz = cumsum(adjHaz)
  survProb = exp(-adjCumHaz)
  # survProb = exp(-linPred[as.character(data1$time)])
  
  medSurvTime = survfit(m1, newdata= data1)$median
  bivarData = data.frame(X = data1$time, delta = data1$event, Y = exp(-linPred), epsilon = rep(1, sampleS))
  dabrSurf = SurvSurfaceOfDabrowska_RecursiveEvenFaster(bivarData)$DabrowskaEst
  tmp = cIndexProbOfConcAndDisc(data1, dabrSurf)$res
  concord = tmp["ConcordProb"]
  discord = tmp["DiscordProb"]
  myCIndex = tmp["CIndex"]

  # resMatr[i, ] = c(cIndex2, myCIndex, concord)
  resMatr[i, "Our estimator"] = c(myCIndex)
}

resMatr[, "Our estimator"][resMatr[, 3]>1] = 1
plotResults(resMatr, trueCIndex = trueCIndexTau075, pch = pch, col = col, cex = cex, ylim = c(0.5, 1))


##################### censoring distribution depends on random beta * Z
##################### censoring distribution depends on random beta * Z
##################### censoring distribution depends on random beta * Z
simNum = 100
# set.seed(20190708)
seeds = sample(1:1000, simNum)
sampleS = 160
resMatr = matrix(NA, ncol = 5, nrow = simNum)
colnames(resMatr) = c("noCensValue", labelsForStuff, "EdenPlusSurvProb")
auxMatr = matrix(NA, ncol = 1, nrow = simNum)
colnames(auxMatr) = c("percCens")
tau075 = qexp(0.5)
trueCIndexTau075 = 0.8586082 ### for 0.5
#################### main loop
for(i in 1:simNum){
  betaEvent = runif(2, min(-2, 2), c(-2, 2))
  betaCens = runif(2, min(-2, 2), c(-2, 2))
  set.seed(seeds[i])
  data1_no_missing = simSurvDataWithCensDepOnZ_take_two(subjNum = sampleS, betaEvent = betaEvent, betaCens = NULL, theta = 4.426265, family = "frank")
  set.seed(seeds[i])
  data1 = simSurvDataWithCensDepOnZ_take_two(subjNum = sampleS, betaEvent = betaEvent, betaCens = betaCens, theta = 4.426265, family = "frank")
  data1 = simSurvDataWithCensDepOnZ_take_two(subjNum = sampleS, betaEvent = betaEvent, betaCens = betaCens, indepCens = TRUE, theta = 4.426265, family = "frank")
  
  auxMatr[i, "percCens"] = mean(data1$event == 0)

  data1$censEvent = 1 - data1$event
  ###########################################################
  ### Uno & Pencina, ...
  mod = coxph(Surv(data1$time, data1$event) ~ V1 + V2, data = data1)
  risk = plogis(predict(mod))
  unotmp = Inf.Cval(data1[, c("time", "event", "V1", "V2")],  tau = 1000, itr = 10)
  resMatr[i, "Uno"] = unotmp$Dhat

  ###########################################################
  ### Harrell ...
  dd = datadist(data1)
  options(datadist = "dd")
  survObj0 = Surv(data1_no_missing$time, data1_no_missing$event)
  survObj1 = Surv(data1$time, data1$event)
  m0 = cph(survObj0 ~ V1 + V2, data = data1_no_missing, x = TRUE, y = TRUE)
  m1 = cph(survObj1 ~ V1 + V2, data = data1, x = TRUE, y = TRUE)
  sum.surv <- summary(m1)
  resMatr[i, "noCensValue"] = survConcordance(survObj0 ~ predict(m0))[["concordance"]]
  resMatr[i, "Harrell"] = survConcordance(survObj1 ~ predict(m1))[["concordance"]]
  
  ###########################################################
  ### Eden ...
  m3 = cph(Surv(data1$time, data1$event) ~ V1 + V2, data = data1, x = TRUE, y = TRUE)
  desMatr = predict(m3, newdata = data1, type = "x")
  cumHaz = basehaz(m3, centered = FALSE)
  regHaz = diff(c(0, cumHaz$hazard))
  names(regHaz) = cumHaz$time
  linPred = as.vector(desMatr %*% m3$coefficients)
  names(linPred) = data1$time
  adjHaz = regHaz[as.character(data1$time)] * exp(linPred[as.character(data1$time)])
  adjCumHaz = cumsum(adjHaz)
  survProb = exp(-adjCumHaz)
  
  medSurvTime = survfit(m1, newdata= data1)$median
  bivarData = data.frame(X = data1$time, delta = data1$event, Y = exp(-linPred), epsilon = rep(1, sampleS))
  dabrSurf = SurvSurfaceOfDabrowska_RecursiveEvenFaster(bivarData)$DabrowskaEst
  tmp = cIndexProbOfConcAndDisc(data1, dabrSurf)$res
  concord = tmp["ConcordProb"]
  discord = tmp["DiscordProb"]
  myCIndex = tmp["CIndex"]
  resMatr[i, "Our estimator"] = c(myCIndex)

  bivarData = data.frame(X = data1$time, delta = data1$event, Y = survProb, epsilon = rep(1, sampleS))
  dabrSurf = SurvSurfaceOfDabrowska_RecursiveEvenFaster(bivarData)$DabrowskaEst
  tmp = cIndexProbOfConcAndDisc(data1, dabrSurf)$res
  resMatr[i, "EdenPlusSurvProb"] = tmp["CIndex"]
}

#resMatr[, "Our estimator"][resMatr[, 3]>1] = 1
# plotResults(resMatr, trueCIndex = trueCIndexTau075, pch = pch, col = col, cex = cex, ylim = c(0.5, 1))

pch = c(noCensValue = 5, Harrell = 2, Uno = 1, "Our estimator" = 4)
cex = c(noCensValue = 1, Harrell = 0.6, Uno = 1.2, "Our estimator" = 0.8)
col = c(noCensValue = "gray", Harrell = "blue", Uno = "red", "Our estimator" = "green")
plotResults(resMatr[, setdiff(colnames(resMatr), "EdenPlusSurvProb")], trueCIndex = mean(resMatr[,"noCensValue"]), pch = pch, col = col, cex = cex)



### the other guys (gonen, heller):
library(clinfun)
coxphCPE(mod)


################## C-INDEX stuff
################## C-INDEX stuff
################## C-INDEX stuff
# library(prodlim)
library(pec)
library(riskRegression)
library(cmprsk)
library(rms)

### simulate data:
set.seed(1)
sampleS = 100
X1 = sample(c(0, 1), sampleS, prob = c(.2, .8), replace = TRUE)
X2 = rnorm(sampleS)
time = sapply(exp(-2*X1 + 1.5*X2), function(x){rexp(1, rate = x)})
#outcome[outcome > 3000] = rexp(1, rate = (-2*X1 + 1.5*X2)[1])
event = rep(1, sampleS)
event[1] = 0

data1 = data.frame(time = time, event = event, X1 = X1, X2 = X2)
dd = datadist(data1)
options(datadist = "dd")
m1 = coxph(Surv(time, event) ~ X1 + X2, data = data1)
predM1 = predict(m1)

# pec::cindex(object = list(m1 = m1), formula = Surv(outcome, event) ~ X1 + X2, data = data.frame(outcome = outcome, X1 = X1, X2 = X2), eval.times = outcome, pred.times = predM1, cause, lyl = FALSE, cens.model = "marginal")

data()
library(dynpred)
cindex(Surv(time, event) ~ X1 + X2, data = data1)

library(survival)
library(rms)

##### Approach 1:
surv <- Surv(time, event)
rcorrcens(surv ~ X1, data = data1)

##### Approach 2:
surv <- Surv(time, event)
sum.surv <- summary(coxph(surv ~ X1))
c_index <- sum.surv$concordance

##### Approach 3:
surv <- Surv(time, event)
fit <- coxph(surv ~ X1)
survConcordance(surv ~ predict(fit))

##### Approach 4:
source("http://bioconductor.org/biocLite.R")
biocLite("survcomp")
library(survcomp)
surv <- Surv(surv, censor) 
fit <- coxph(surv ~ X1, data= data1)
coxPredict <- predict(fit, data=sample.data, type="risk")  
concordance.index(x=coxPredict, surv.time=time, surv.event=event, method="noether")

##### Approach 5:
library(rms)
set.seed(1)
dd = datadist(data1)
options(datadist = "dd")
fit.cph <- cph(surv ~ X1, data= data1, x=TRUE, y=TRUE, surv=TRUE)
  
# Get the Dxy
v <- validate.cph(fit.cph, dxy=TRUE, B=1000)
Dxy = v[rownames(v)=="Dxy", colnames(v)=="index.corrected"]
orig_Dxy = v[rownames(v)=="Dxy", colnames(v)=="index.orig"]
  
# The c-statistic according to Dxy=2(c-0.5)
bias_corrected_c_index  <- abs(Dxy)/2+0.5
orig_c_index <- abs(Orig.Dxy)/2+0.5

##### Yet another approach:
surv <- Surv(survival, censor)
c_index <- function(group, ties=TRUE){
  fit <- coxph(surv ~ group, data=sample.data)
  coxPredict <- predict(fit, data=sample.data, type="risk")  
  
  # Approaches 4/5
  if (ties==F) {
  concordance.index(x=coxPredict, surv.time=survival, surv.event=censor, method="noether")
  }
  # Approaches 2/3
  else if (ties==T) {
  survConcordance(surv ~ coxPredict, data=sample.data)
  }
}
c_index_ties1 <- c_index(group=group1, ties=TRUE)
c_index_ties2 <- c_index(group=group2, ties=TRUE)

c_index_no_ties1 <- c_index_ties(group=group1, ties=F)
c_index_no_ties2 <- c_index_ties(group=group2, ties=F)

# p-value of testing two c-indices ignoring ties
round(cindex.comp(c_index_no_ties1, c_index_no_ties2)$p.value,3)

# Function for p-value of testing two c-indices accounting for ties
# t-test for dependent variables is used for significance 
# Input variables are objects obtained from the first function

cindex.p.ties <- function(c_index_ties1, c_index_ties2, c_index_no_ties1, c_index_no_ties2) {
    eps <- 1E-15
    n <- c_index_no_ties1$n
    r <- cor(c_index_no_ties1$data$x, c_index_no_ties2$data$x, use="complete.obs", method="spearman")
    if ((1 - abs(r)) > eps) {
      t.stat <- (c_index_ties1$concordance - c_index_ties2$concordance) / sqrt(c_index_ties1$std.err^2 + c_index_ties2$std.err^2 - 2 * r * c_index_ties1$std.err * c_index_ties2$std.err)
      diff.ci.p <- pt(q=t.stat, df=n - 1, lower.tail=FALSE)
    } else { diff.ci.p <- 1 }
    return(list("p.value"=diff.ci.p))
  }
  
cindex.p.ties(c_index_ties1=c_index_ties1, c_index_ties2=c_index_ties2, c_index_no_ties1=c_index_no_ties1, c_index_no_ties2=c_index_no_ties2)


sampleS = 10000
X1 = rnorm(sampleS)
X2 = rnorm(sampleS)
time = sapply(exp(-2*X1 + 1.5*X2), function(x){rexp(1, rate = x)})
event = rep(1, sampleS)
data1 = data.frame(time = time, event = event, X1 = X1, X2 = X2)
dd = datadist(data1)
options(datadist = "dd")
survObj = Surv(time, event)
m1 = cph(survObj ~ X1 + X2, data = data1, x = TRUE, y = TRUE)
sum.surv <- summary(m1)
# cIndex1 <- sum.surv$concordance[["C"]]
cIndex2 = survConcordance(survObj ~ predict(m1))[["concordance"]]


simNum = 200
set.seed(20190708)
sampleS = 100
resMatr = matrix(NA, ncol = 3, nrow = simNum)
for(i in 1:simNum){
  X1 = rnorm(sampleS)
  # X1 = sample(c(0, 1), sampleS, prob = c(.2, .8), replace = TRUE)
  X2 = rnorm(sampleS)
  time = sapply(exp(-2*X1 + 1.5*X2), function(x){rexp(1, rate = x)})
  event = rep(1, sampleS)
  ### sim censoring
  if(TRUE){
    cenTime = rexp(sampleS, rate = mean(exp(-2*X1 + 1.5*X2))/5)
    event[time>cenTime] = 0
    time[time>cenTime] = cenTime[time>cenTime]
  }
  data1 = data.frame(time = time, event = event, X1 = X1, X2 = X2)
  dd = datadist(data1)
  options(datadist = "dd")
  survObj = Surv(time, event)
  m1 = cph(survObj ~ X1 + X2, data = data1, x = TRUE, y = TRUE)
  # predM1 = predict(m1, newdata = data1, type = "data.frame")
  
  sum.surv <- summary(m1)
  # cIndex1 <- sum.surv$concordance[["C"]]
  cIndex2 = survConcordance(survObj ~ predict(m1))[["concordance"]]
  
  desMatr = predict(m1, newdata = data1, type = "x")
  cumHaz = basehaz(m1, centered = FALSE)
  regHaz = diff(c(0, cumHaz$hazard))
  names(regHaz) = cumHaz$time
  linPred = as.vector(desMatr %*% m1$coefficients)
  names(linPred) = data1$time
  adjHaz = regHaz[as.character(data1$time)] * exp(linPred[as.character(data1$time)])
  adjCumHaz = cumsum(adjHaz)
  survProb = exp(-adjCumHaz)
  
  medSurvTime = survfit(m1,newdata= data1)$median
  bivarData = data.frame(X = time, delta = event, Y = exp(-linPred), epsilon = rep(1, sampleS))
  dabrSurf = SurvSurfaceOfDabrowska_RecursiveEvenFaster(bivarData)$DabrowskaEst
  
  tmp = cIndexProbOfConcAndDisc(bivarData, dabrSurf)$res
  concord = tmp["ConcordProb"]
  discord = tmp["DiscordProb"]
  myCIndex = tmp["CIndex"]

  # resMatr[i, ] = c(cIndex1, cIndex2, concord)
  resMatr[i, ] = c(cIndex2, myCIndex, concord)
}
resMatr[, 3][resMatr[, 3]>1] = 1

trueCIndex = 0.857143
col = c("orange", "blue", "green")
plot(range(dim(resMatr)), range(resMatr[,1:2]), type = "n")
points(1:nrow(resMatr), resMatr[, 1], pch = 23, col = col[1], cex = 1)
points(1:nrow(resMatr), resMatr[, 2], pch = 19, col = col[2], cex = .5)
abline(h = apply(resMatr[,1:2], 2, mean), col = col[1:2], lwd = c(3, 1))
abline(h = trueCIndex, col = "black", lty = 3)
legend(x = "bottomleft", pch = c(23, 19, 21), col = col, legend = c("reg CInd", "Your concord", "Nothing for now"))
bias = trueCIndex - apply(resMatr[,1:2], 2, mean)
50  -0.03879750  -0.01793876
100 -0.03513393  -0.009500934
    -0.03558022  -0.01384252
200 -0.03727289  -0.01681142
300 -0.036113758 -0.009874169
400 -0.03524658  -0.01070548
# restrictedRegionSpearmanRhoAndTauUsingSurvSurf_close_to_final(bivarData, dabrSurf)$res
