
#' @import stats survival geepack
#' @importFrom frailtypack frailtyPenal
#' @importFrom stats runif formula model.matrix predict terms
#' @importFrom graphics plot points segments
#' @importFrom data.table data.table setkey rbindlist


lagby1.1person <- function(x){
	if(length(x)>1){ return(c(NA,x[-length(x)]))}
	else return(NA)
}


lagby1.1var <- function(x,id,time,lagfirst=NA){
  if(length(x)>1){
    lagx <- c(lagfirst,x[-length(x)])
    rows <- 1:length(x)

    data <- cbind(rows,x,id,time,lagx)
    data <- as.data.frame(data)
    names(data) <- c("rows","x","id","time","lagx")
    DT <- data.table(data)
    setkey(DT, id)
    dt1 <- data.frame(DT[data.table(unique(id)), mult = "first"])
    dt1$lagx <- lagfirst

    lagx[dt1$rows] <- lagfirst
    return(lagx)
  }
  else return(NA)
}



# data must be ordered by id and time.
# id variable should be called id, time variable should be called time

#' Create lagged versions the variables in data
#'
#' @param data The data to be lagged
#' @param lagvars The names of the columns in the data to be lagged
#' @param id A character indicating which column of the data contains subject identifiers. ids are assumed to be consecutive integers, with the first subject having id 1
#' @param time A character indicating which column of the data contains the times at which each of the observations in data was made
#' @param lagfirst A vector giving the value of each lagged variable for the first time within each subject. This is helpful if, for example, time is the variable to be lagged and you know that all subjects entered the study at time zero
#' @return The original data frame with lagged variables added on as columns. For example, if the data frame contains a variable named x giving the value of x for each subject i at each visit j, the returned data frame will contain a column named x.lag containing the value of x for subject i at visit j-1. If j is the first visit for subject i, the lagged value is set to NA
#' @examples
#' library(nlme)
#' data(Phenobarb)
#' head(Phenobarb)
#'
#' data <- lagfn(Phenobarb,"time","Subject","time")
#' head(data)
#' @export


lagfn <- function(data,lagvars,id,time,lagfirst=NA){

#	lagged <- sapply(data[,columns],lagby1.1var,data[,names(data)%in%id],data[,names(data)%in%time])
  columns <- (1:ncol(data))[is.finite(match(names(data), lagvars))]
  for(col in 1:length(columns)){
  	if(col==1) lagged <- lagby1.1var(x=data[,columns[col]],id=data[,names(data)%in%id],time=data[,names(data)%in%time],lagfirst=lagfirst[col])
  	if(col>1) lagged <- cbind(lagged,lagby1.1var(x=data[,columns[col]],id=data[,names(data)%in%id],time=data[,names(data)%in%time],lagfirst=lagfirst[col]))
  }
  	lagged <- as.data.frame(lagged)
	for(col in 1:length(columns)) names(lagged)[col] <- paste(names(data)[columns][col],".lag",sep="")
	data <- cbind(data,lagged)
	return(data)
}

#' Add rows corresponding to censoring times to a longitudinal dataset
#'
#' @param data The dataset to which rows are to be added. The data should have one row per observation
#' @param maxfu The maximum follow-up time per subject. If all subjects have the same follow-up time, this can be supplied as a single number. Otherwise, maxfu should be a dataframe with the first column specifying subject identifiers and the second giving the follow-up time for each subject.
#' @param tinvarcols A vector of column numbers corresponding to variables in data that are time-invariant.
#' @param id character string indicating which column of the data identifies subjects
#' @param time character string indicating which column of the data contains the time at which the visit occurred
#' @param event character string indicating which column of the data indicates whether or not a visit occurred. If every row corresponds to a visit, then this column will consist entirely of ones
#' @examples
#' x <- c(1:3,1:2,1:5)
#' x0 <- c(rep(2,3),rep(0,2),rep(1,5))
#' id <- c(rep(1,3),rep(2,2),rep(3,5))
#' time <- c(0,4,6,2,3,1,3,5,6,7)
#' event <- c(1,1,1,0,1,0,1,1,1,1)
#' data <- as.data.frame(cbind(x,id,time,event,x0))
#' addcensoredrows(data,maxfu=8,id="id",time="time",tinvarcols=5,event="event")
#'
#'
#' x <- c(1:3,1:2,1:5)
#' x0 <- c(rep(2,3),rep(0,2),rep(1,5))
#' id <- c(rep(1,3),rep(2,2),rep(3,5))
#' time <- c(0,4,6,2,3,1,3,5,6,7)
#' event <- c(1,1,1,0,1,0,1,1,1,1)
#' data <- as.data.frame(cbind(x,id,time,event,x0))
#' maxfu.id <- 1:3
#' maxfu.time <- c(6,5,8)
#' maxfu <- cbind(maxfu.id,maxfu.time)
#' maxfu <- as.data.frame(maxfu)
#' addcensoredrows(data,maxfu=maxfu,id="id",time="time",tinvarcols=5,event="event")
#' @return The original dataset with extra rows corresponding to censoring times
#' @export

addcensoredrows <- function(data,maxfu,tinvarcols,id,time,event){
  if(length(maxfu)==1){ maxfu <- rep(maxfu,length(table(data[,names(data)%in%id]))) } else maxfu <- maxfu[,2]

  data <- data[order(data[,names(data)%in%id],data[,names(data)%in%time]),]

  extrarows <- array(dim=c(length(table(data[,names(data)%in%id])),ncol(data)))
  extrarows <- as.data.frame(extrarows)
  names(extrarows) <- names(data)

  DT <- data.table(data)
  setkey(DT, id)

  dt1 <- data.frame(DT[data.table(unique(id)), mult = "first"])
  extrarows[,tinvarcols] <- dt1[,tinvarcols]
  extrarows[,names(data)%in%id] <- dt1[,names(dt1)%in%id]
  extrarows[,names(data)%in%time] <- maxfu
  extrarows[,names(data)%in%event] <- 0

  maxobstime <- tapply(data[,names(data)%in%time],data[,names(data)%in%id],max)
  alreadythere <- sapply(maxobstime-maxfu,identical,0)
  extrarows <- extrarows[!alreadythere,]

  data <- rbind(data.matrix(data),data.matrix(extrarows))
  data <- as.data.frame(data)
  data <- data[order(data[,names(data)%in%id],data[,names(data)%in%time]),]
  return(data)
}

# assumes no missing data
phfn <- function(datacox,regcols,data){
	X <- datacox[,regcols]
	datacox <- datacox[is.finite(datacox$time.lag),]
	response <- names(datacox[,regcols])[1]
	if(length(regcols)>1) for(p in 2:length(regcols)) response <- paste(response,"+",names(datacox[,regcols])[p])
	m <- coxph(Surv(time.lag,time,event) ~ response,data=datacox)
	print(summary(m))
	hr <- predict(m,newdata=data,type="lp")
	iiw.weight <- 1/hr
	return(iiw.weight)
}



#' Given a proportional hazards model for visit intensities, compute inverse-intensity weights.
#'
#' For a longitudinal dataset subject to irregular observation, use a Cox proportional hazards model for visit intensities to compute inverse intensity weights
#' @family iiw
#' @param phfit coxph object for the visit process
#' @param data The dataset featuring longitudinal data subject to irregular observation for which inverse-intensity weights are desired
#' @param time character string indicating which column of the data contains the time at which the visit occurred
#' @param id character string indicating which column of the data identifies subjects
#' @param first logical variable. If TRUE, the first observation for each individual is assigned an intensity of 1. This is appropriate if the first visit is a baseline visit at which recruitment to the study occurred; in this case the baseline visit is observed with probability 1.
#' @return A vector of inverse-intensity weights for each row of the dataset. The first observation for each subject is assumed to have an intensity of 1.
#' @examples
#' library(nlme)
#' data(Phenobarb)
#' library(survival)
#' library(geepack)
#' data(Phenobarb)
#' Phenobarb$event <- 1-as.numeric(is.na(Phenobarb$conc))
#' data <- lagfn(Phenobarb, lagvars="dose", id="Subject", time="time", lagfirst = 0)
#' data <- lagfn(data, lagvars="dose.lag", id="Subject", time="time", lagfirst = 0)
#' data <- lagfn(data, lagvars="dose.lag.lag", id="Subject", time="time", lagfirst = 0)
#' data$dose.lag[is.na(data$dose.lag)] <- data$dose.lag.lag[is.na(data$dose.lag)]
#' data$dose.lag[is.na(data$dose.lag)] <- data$dose.lag.lag.lag[is.na(data$dose.lag)]
#' data <- data[data$event==1,]
#' data$id <- as.numeric(data$Subject)
#' data <- data[data$time<16*24,]
#' data <- lagfn(data, lagvars=c("time","conc"), id="Subject", time="time", lagfirst = NA)
#' head(data)
#'
#' mph <- coxph(Surv(time.lag,time,event)~I(conc.lag>0) + conc.lag + cluster(id),data=data)
#' summary(mph)
#' data$weight <- iiw(mph,data,"id","time",TRUE)
#' head(data)
#' @export

iiw <- function(phfit,data,id,time,first){
	data$iiw.weight <- exp(-predict(phfit,type="lp",newdata=data))
	tmin <- tapply(data[,names(data)%in%time],data[,names(data)%in%id],min)
	firstvisit <- as.numeric(data[,names(data)%in%time]==tmin[data[,names(data)%in%id]])
	if(first) data$iiw.weight[firstvisit==1] <- 1
	return(data$iiw.weight)
}

#' Fit an inverse-intensity weighted GEE.
#'
#' Implements inverse-intensity weighted GEEs as first described by Lin, Scharfstein and Rosenheck (2004). A Cox proportional hazards model is applied to the visit intensities, and the hazard multipliers are used to compute inverse-intensity weights. Using the approach described by Buzkova and Lumley (2007) avoids the need to compute the baseline hazard.
#'
#' @details
#' Let the outcome of interest be \eqn{Y} and suppose that subject i has \eqn{j^{th}} observation at \eqn{T_{ij}}. Let \eqn{N_i(t)} be a counting process for the number of observations for subject i up to and including time t. Suppose that \eqn{N_i} has intensity \eqn{\lambda} given by \deqn{\lambda_i(t)=\lambda0(t)exp(Z_i(t)\gamma).} Then the inverse-intensity weights are \deqn{exp(-Z_i(t)\gamma).} If \eqn{Y_i} is the vector of observations for subject \eqn{i}, to be regressed onto \eqn{X_i} (i.e. \eqn{E(Y_i|X_i)=\mu(X_i;\beta)} with \eqn{g(\mu(X_i;beta)=X_i\beta}, then the inverse-intensity weighted GEE equations are \deqn{\sum_i \frac{\partial\mu_i}{\partial\beta}V_i^{-1}\Delta_i(Y_i X_i\beta)=0}, where \eqn{\Delta_i} is a diagonal matrix with \eqn{j^{th}} entry equal to \eqn{\exp(-Z_i(T_{ij})\gamma)} and $V_i$ is the working variance matrix.
#' Warning: Due to the way some gee functions incorporate weights, if using inverse-intensity weighting you should use working independence.
#' @param formulagee the formula for the GEE model to be fit. The syntax used is the same as in geeglm
#' @param formulaph the formula for the proportional hazards model for the visit intensity that will be used to derive inverse-intensity weights. The formula should usually use the counting process format (i.e. Surv(start,stop,event))
#' @family iiw
#' @param formulanull if stabilised weights are to be used, the formula for the null model used to stabilise the weights
#' @param data data frame containing the variables in the model
#' @param id character string indicating which column of the data identifies subjects
#' @param time character string indicating which column of the data contains the time at which the visit occurred
#' @param event character string indicating which column of the data indicates whether or not a visit occurred. If every row corresponds to a visit, then this column will consist entirely of ones
#' @param family family to be used in the GEE fit. See geeglm for documentation
#' @param lagvars a vector of variable names corresponding to variables which need to be lagged by one visit to fit the visit intensity model. Typically time will be one of these variables. The function will internally add columns to the data containing the values of the lagged variables from the previous visit. Values of lagged variables for a subject's first visit will be set to NA. To access these variables in specifying the proportional hazards formulae, add ".lag" to the variable you wish to lag. For example, if time is the variable for time, time.lag is the time of the previous visit
#' @param invariant a vector of variable names corresponding to variables in data that are time-invariant. It is not necessary to list every such variable, just those that are invariant and also included in the proportional hazards model
#' @param maxfu the maximum follow-up time(s). If everyone is followed for the same length of time, this can be given as a single value. If individuals have different follow-up times, maxfu should have the same number of elements as there are rows of data
#' @param lagfirst A vector giving the value of each lagged variable for the first time within each subject. This is helpful if, for example, time is the variable to be lagged and you know that all subjects entered the study at time zero
#' @param first logical variable. If TRUE, the first observation for each individual is assigned an intensity of 1. This is appropriate if the first visit is a baseline visit at which recruitment to the study occurred; in this case the baseline visit is observed with probability 1.
#' @return a list, with the following elements:
#' \item{geefit}{the fitted GEE, see documentation for geeglm for details}
#' \item{phfit}{the fitted proportional hazards model, see documentation for coxph for details}
#' @references
#' \itemize{
#' \item Lin H, Scharfstein DO, Rosenheck RA. Analysis of Longitudinal data with Irregular, Informative Follow-up. Journal of the Royal Statistical Society, Series B (2004), 66:791-813
#' \item Buzkova P, Lumley T. Longitudinal data analysis for generalized linear models with follow-up dependent on outcome-related variables. The Canadian Journal of Statistics 2007; 35:485-500.}
#' @examples
#' library(nlme)
#' data(Phenobarb)
#' library(survival)
#' library(geepack)
#' Phenobarb$event <- 1-as.numeric(is.na(Phenobarb$conc))
#' data <- lagfn(Phenobarb, lagvars="dose", id="Subject", time="time", lagfirst = 0)
#' data <- lagfn(data, lagvars="dose.lag", id="Subject", time="time", lagfirst = 0)
#' data <- lagfn(data, lagvars="dose.lag.lag", id="Subject", time="time", lagfirst = 0)
#' data$dose.lag[is.na(data$dose.lag)] <- data$dose.lag.lag[is.na(data$dose.lag)]
#' data$dose.lag[is.na(data$dose.lag)] <- data$dose.lag.lag.lag[is.na(data$dose.lag)]
#' data <- data[data$event==1,]
#' data$id <- as.numeric(data$Subject)
#' data <- data[data$time<16*24,]
#' miiwgee <- iiwgee(conc ~ I(time^3) + log(time),Surv(time.lag,time,event)~Wt *dose.lag +
#' dose.lag*(I(conc.lag>0) + conc.lag) +cluster(Subject),
#' formulanull=NULL,id="id",time="time",event="event",data=data,
#' invariant="id",lagvars=c("time","conc"),maxfu=16*24,lagfirst=0,first=TRUE)
#' summary(miiwgee$geefit)
#' summary(miiwgee$phfit)
#'
#' # compare to results without weighting
#' m <- geeglm(conc ~ I(time^3) + log(time) , id=Subject, data=data); print(summary(m))
#' time <- (1:200)
#' unweighted <- cbind(rep(1,200),time^3,log(time))%*%m$coefficients
#' weighted <- cbind(rep(1,200),time^3,log(time))%*%miiwgee$geefit$coefficients
#' plot(data$time,data$conc,xlim=c(0,200),pch=16)
#' lines(time,unweighted,type="l")
#' lines(time,weighted,col=2)
#' legend (0,60,legend=c("Unweighted","Inverse-intensity weighted"),col=1:2,bty="n",lty=1)
#' @export

iiwgee <- function(formulagee,formulaph,formulanull=NULL,data,id,time,event,family=gaussian,lagvars,invariant=NULL,maxfu,lagfirst,first){
# id is the id variable
# lagvars are the variables to be lagged
# invariant are the variables that are invariant. Only need to be entered if they are in formulaph
# data should include a variable named event, which should be 1 if the event occurred ################ need to create an event argument for the formula
# maxfu is either the administrative censoring time (same for everyone) or else is a dataframe one row for each subject, columns for ids and the censoring times. The name of the id variable should be the same as in data

	# sort the data on id then time
	data <- data[order(data[,names(data)%in%id],data[,names(data)%in%time]),]
	weights <- iiw.weights(formulaph, formulanull,data=data,id=id,time=time,event=event,lagvars=lagvars,invariant=invariant,maxfu=maxfu,first=first,lagfirst=lagfirst,frailty=FALSE)
	data$useweight <- weights$iiw.weight
	m <- weights$m

	# fit IIW-GEE
	data$iddup <- data[,names(data)%in%id]
	useweight <- data$useweight
	iddup <- data$iddup

	mgee <- geeglm(formulagee,id=iddup,data=data,corstr="independence",weights=useweight,family=family)
	return(list(geefit=mgee,phfit=m))
}

#' Compute inverse-intensity weights.
#'
#' Given longitudinal data with irregular visit times, fit a Cox proportional hazards model for the visit intensity, then use it to compute inverse-intensity weights
#' @param formulaph the formula for the proportional hazards model for the visit intensity that will be used to derive inverse-intensity weights. The formula should usually use the counting process format (i.e. Surv(start,stop,event)). If a frailty model is used, the cluster(id) term should appear before other covariates
#' @family iiw
#' @param formulanull if stabilised weights are to be used, the formula for the null model used to stabilise the weights
#' @param data data frame containing the variables in the model
#' @param id character string indicating which column of the data identifies subjects
#' @param time character string indicating which column of the data contains the time at which the visit occurred
#' @param event character string indicating which column of the data indicates whether or not a visit occurred. If every row corresponds to a visit, then this column will consist entirely of ones
#' @param lagvars a vector of variable names corresponding to variables which need to be lagged by one visit to fit the visit intensity model. Typically time will be one of these variables. The function will internally add columns to the data containing the values of the lagged variables from the previous visit. Values of lagged variables for a subject's first visit will be set to NA. To access these variables in specifying the proportional hazards formulae, add ".lag" to the variable you wish to lag. For example, if time is the variable for time, time.lag is the time of the previous visit
#' @param invariant a vector of variable names corresponding to variables in data that are time-invariant. It is not necessary to list every such variable, just those that are invariant and also included in the proportional hazards model
#' @param maxfu the maximum follow-up time(s). If everyone is followed for the same length of time, this can be given as a single value. If individuals have different follow-up times, maxfu should have the same number of elements as there are rows of data
#' @param lagfirst A vector giving the value of each lagged variable for the first time within each subject. This is helpful if, for example, time is the variable to be lagged and you know that all subjects entered the study at time zero
#' @param first logical variable. If TRUE, the first observation for each individual is assigned an intensity of 1. This is appropriate if the first visit is a baseline visit at which recruitment to the study occurred; in this case the baseline visit is observed with probability 1.
#' @param frailty logical variable. If TRUE, a frailty model is fit to calculate the inverse intensity weights. If FALSE, a marginal semi-parametric model is fit. Frailty models are helpful when fitting semi-parametric joint models.
#' @return a vector of inverse-intensity weights, ordered on id then time
#' @description Since the vector of weights is ordered on id and time, if you intend to merge these weights onto your original dataset it is highly recommended that you sort the data before running iiw.weights
#' @references
#' \itemize{
#' \item Lin H, Scharfstein DO, Rosenheck RA. Analysis of Longitudinal data with Irregular, Informative Follow-up. Journal of the Royal Statistical Society, Series B (2004), 66:791-813
#' \item Buzkova P, Lumley T. Longitudinal data analysis for generalized linear models with follow-up dependent on outcome-related variables. The Canadian Journal of Statistics 2007; 35:485-500.}
#' @examples
#' library(nlme)
#' data(Phenobarb)
#' library(survival)
#' library(geepack)
#' Phenobarb$event <- 1-as.numeric(is.na(Phenobarb$conc))
#' data <- lagfn(Phenobarb, lagvars="dose", id="Subject", time="time", lagfirst = 0)
#' data <- lagfn(data, lagvars="dose.lag", id="Subject", time="time", lagfirst = 0)
#' data <- lagfn(data, lagvars="dose.lag.lag", id="Subject", time="time", lagfirst = 0)
#' data$dose.lag[is.na(data$dose.lag)] <- data$dose.lag.lag[is.na(data$dose.lag)]
#' data$dose.lag[is.na(data$dose.lag)] <- data$dose.lag.lag.lag[is.na(data$dose.lag)]
#' data <- data[data$event==1,]
#' data$id <- as.numeric(data$Subject)
#' data <- data[data$time<16*24,]
#' i <- iiw.weights(Surv(time.lag,time,event)~Wt *dose.lag + dose.lag*(I(conc.lag>0) + conc.lag) +
#' cluster(Subject),id="id",time="time",event="event",data=data,
#' invariant="Subject",lagvars=c("time","conc"),maxfu=16*24,lagfirst=0,first=TRUE)
#' data$weight <- i$iiw.weight
#' summary(i$m)
#' # can use to fit a weighted GEE
#' mw <- geeglm(conc ~ I(time^3) + log(time) , id=Subject, data=data, weights=weight)
#' summary(mw)
#' # agrees with results through the single command iiwgee
#' miiwgee <- iiwgee(conc ~ I(time^3) + log(time),Surv(time.lag,time,event)~Wt *dose.lag +
#' dose.lag*(I(conc.lag>0) + conc.lag) +cluster(Subject),
#' formulanull=NULL,id="id",time="time",event="event",data=data,
#' invariant="id",lagvars=c("time","conc"),maxfu=16*24,lagfirst=0,first=TRUE)
#' summary(miiwgee$geefit)
#' @family iiw
#' @export

iiw.weights <- function(formulaph,formulanull=NULL,data,id,time,event,lagvars,invariant=NULL,maxfu,lagfirst=lagfirst,first,frailty=FALSE){
	if(is.null(maxfu) & frailty){ maxtable <- tapply(data[,names(data)%in%time],data[,names(data)%in%id],max); maxfu <- cbind(1:length(maxtable),maxtable + max(maxtable)*0.001)}
	data <- data[order(data[,names(data)%in%id],data[,names(data)%in%time]),]

  if(is.null(invariant) & !is.null(maxfu)) print("invariant must be specified if maxfu is non-null")

	#	lagcols <- (1:ncol(data))[is.finite(match(names(data), lagvars))]
	invarcols <- (1:ncol(data))[is.finite(match(names(data), invariant))]

	stabilize <- !is.null(formulanull)

	# redo id variable so ids are numbered using consecutive integers
	ids <- names(table(data[,names(data)%in%id]))
	idnum <- array(dim=nrow(data))
	for(i in 1:nrow(data)) idnum[i] <- (1:length(ids))[data[i,names(data)%in%id]==ids]

	if(!is.null(maxfu)){if(length(maxfu)>1 | sum(is.na(maxfu))==0){ if(is.data.frame(maxfu)){ for(i in 1:nrow(maxfu)){ maxfu[i,names(maxfu)%in%id] <- (1:length(ids))[maxfu[i,names(maxfu)%in%id]==ids]}}}}
	data[,names(data)%in%id] <- idnum


	# create weights
	if(!is.null(maxfu)|frailty){if(length(maxfu)==1) datacox <- data else datacox <- addcensoredrows(data=data,maxfu=maxfu,tinvarcols=invarcols,id=id,time=time,event=event)}
	if(is.null(maxfu)) datacox <- data
	datacox <- lagfn(datacox,lagvars,id,time,lagfirst)
	if(!frailty){
	m <- coxph(formula=formulaph,data=datacox)
	data$iiw.weight <- iiw(m,lagfn(data,lagvars,id,time,lagfirst),id,time,first)
	if(stabilize){
		m0 <- coxph(formula=formulanull,data=datacox)
		data$nullweight <- iiw(m0,data,id,time,first)
		data$useweight <- data$iiw.weight/data$nullweight
	} else{ data$useweight <- data$iiw.weight}
	}
	if(frailty){
    row.names(datacox)[datacox$event==1] <- row.names(data)
	  use <- (1:nrow(datacox))[!is.na(datacox[,names(datacox)%in%paste(time,".lag",sep="")])]
	  		m <- frailtyPenal(formula=formulaph,data=datacox[use,], recurrentAG=TRUE, n.knots=6, kappa=10000,cross.validation=TRUE)

		omit <- c("(Intercept)",paste("cluster(",id,")",sep=""))
		Xmat <- model.matrix(terms(formula(m)),data=datacox)
		columns <- colnames(Xmat)
		Xmat.use <- Xmat[row.names(Xmat)%in%row.names(data),!columns%in%omit]
		Xmat.use <- array(Xmat.use,dim=c(nrow(Xmat[row.names(Xmat)%in%row.names(data),]),length(columns[!columns%in%omit])))
#		data$iiw.weight[row.names(data)%in%row.names(Xmat)] <- exp(-Xmat[row.names(Xmat)%in%row.names(data),!columns%in%omit]%*%m$coef)
		data$iiw.weight[row.names(data)%in%row.names(Xmat)] <- exp(-Xmat.use%*%m$coef)

	if(stabilize){
		m0 <- frailtyPenal(formula=formulanull,data=datacox[use,],recurrentAG=TRUE,n.knots=6,kappa=10000,cross.validation=TRUE)
		omit <- c("(Intercept)",paste("cluster(",id,")",sep=""))

		Xmat <- model.matrix(terms(formula(m0)),data=datacox)
		columns <- colnames(Xmat)
		Xmat.use <- Xmat[row.names(Xmat)%in%row.names(data),!columns%in%omit]
		Xmat.use <- array(Xmat.use,dim=c(nrow(Xmat[row.names(Xmat)%in%row.names(data),]),length(columns[!columns%in%omit])))

		data$null.weight[row.names(data)%in%row.names(Xmat)] <- exp(-Xmat.use%*%m0$coef)
				data$useweight <- data$iiw.weight/data$null.weight
	} else{ data$useweight <- data$iiw.weight; m0 <- NULL}

	}
	if(stabilize==FALSE) m0 <- NULL
	tmin <- tapply(data[,names(data)%in%time],data[,names(data)%in%id],min)
	firstvisit <- as.numeric(data[,names(data)%in%time]==tmin[data[,names(data)%in%id]])
	if(first){data$useweight[firstvisit==1] <- 1}

	return(list(iiw.weight=data$useweight,m=m,m0=m0,datacox=datacox))
}



#' Create an outputted dataset for use with multiple outputation.
#'
#' Multiple outputation is a procedure whereby excess observations are repeatedly randomly sampled and discarded. The method was originally developed to handle clustered data where cluster size is informative, for example when studying pups in a litter. In this case, analysis that ignores cluster size results in larger litters being over-represented in a marginal analysis. Multiple outputation circumvents this problem by randomly selecting one observation per cluster. Multiple outputation has been further adapted to handle longitudinal data subject to irregular observation; here the probability of being retained on any given outputation is inversely proportional to the visit intensity. This function creates a single outputted dataset.

#' @param data the original dataset on which multiple outputation is to be performed
#' @param weights the weights to be used in the outputation, i.e. the inverse of the probability that a given observation will be selected in creating an outputted dataset. Ignored if singleobs=TRUE
#' @param singleobs logical variable indicating whether a single observation should be retained for each subject
#' @param id character string indicating which column of the data identifies subjects
#' @param time character string indicating which column of the data contains the time at which the visit occurred
#' @param keep.first logical variable indicating whether the first observation should be retained with probability 1. This is useful if the data consists of an observation at baseline followed by follow-up at stochastic time points.
#' @return the outputted dataset.
#' @references
#' \itemize{
#' \item Hoffman E, Sen P, Weinberg C. Within-cluster resampling. Biometrika 2001; 88:1121-1134
#' \item Follmann D, Proschan M, Leifer E. Multiple outputation: inference for complex clustered data by averaging analyses from independent data. Biometrics 2003; 59:420-429
#' \item Pullenayegum EM. Multiple outputation for the analysis of longitudinal data subject to irregular observation. Statistics in Medicine (in press).}
#' @family mo
#' @examples
#' library(nlme)
#' data(Phenobarb)
#' library(survival)
#' library(geepack)
#' Phenobarb$event <- 1-as.numeric(is.na(Phenobarb$conc))
#' data <- lagfn(Phenobarb, lagvars="dose", id="Subject", time="time", lagfirst = 0)
#' data <- lagfn(data, lagvars="dose.lag", id="Subject", time="time", lagfirst = 0)
#' data <- lagfn(data, lagvars="dose.lag.lag", id="Subject", time="time", lagfirst = 0)
#' data$dose.lag[is.na(data$dose.lag)] <- data$dose.lag.lag[is.na(data$dose.lag)]
#' data$dose.lag[is.na(data$dose.lag)] <- data$dose.lag.lag.lag[is.na(data$dose.lag)]
#' data <- data[data$event==1,]
#' data$id <- as.numeric(data$Subject)
#' data <- data[data$time<16*24,]
#' i <- iiw.weights(Surv(time.lag,time,event)~Wt *dose.lag + dose.lag*(I(conc.lag>0) + conc.lag) +
#' cluster(Subject),id="id",time="time",event="event",data=data,
#' invariant="Subject",lagvars=c("time","conc"),maxfu=16*24,lagfirst=0,first=TRUE)
#' data$weight <- i$iiw.weight
#' head(data)
#' data.output1 <-   outputation(data,data$weight,singleobs=FALSE,
#' id="id",time="time",keep.first=FALSE)
#' head(data.output1)
#' data.output2 <-   outputation(data,data$weight,singleobs=FALSE,
#' id="id",time="time",keep.first=FALSE)
#' head(data.output2)
#' data.output3 <-   outputation(data,data$weight,singleobs=FALSE,
#' id="id",time="time",keep.first=FALSE)
#' head(data.output3)
#' # Note that the outputted dataset varies with each command run; outputation is done at random
#' @export


outputation <- function(data,weights,singleobs,id,time,keep.first){
	if(singleobs==TRUE){
		ids <- names(table(data[,names(data)%in%id]))
		idnum <- array(dim=nrow(data))
		for(i in 1:nrow(data)) idnum[i] <- (1:length(ids))[data[i,names(data)%in%id]==ids]
#		if(is.data.frame(maxfu)){ for(i in 1:nrow(maxfu)){ maxfu[i,names(maxfu)%in%id] <- (1:length(ids))[maxfu[i,names(maxfu)%in%id]==ids]}}
		data[,names(data)%in%id] <- idnum

		obsnum <- array(1,nrow(data))
		for(row in 2:nrow(data)){
			if(data[row,names(data)%in%id]==data[row-1,names(data)%in%id]) obsnum[row] <- obsnum[row-1] + 1
		}

		nobs <- tapply(data[,names(data)%in%id],data[,names(data)%in%id],length)
		samples <- sapply(nobs,sample.int,size=1)
		include <- obsnum==samples[data[,names(data)%in%id]]
		datamo <- data[include,]
	}
	if(singleobs==FALSE){
		unif <- runif(nrow(data))
		choosevec <- as.numeric(unif<(weights)/max(weights))


		if(keep.first==1){
		ids <- names(table(data[,names(data)%in%id]))
		idnum <- array(dim=nrow(data))
		for(i in 1:nrow(data)) idnum[i] <- (1:length(ids))[data[i,names(data)%in%id]==ids]
#		if(is.data.frame(maxfu)){ for(i in 1:nrow(maxfu)){ maxfu[i,names(maxfu)%in%id] <- (1:length(ids))[maxfu[i,names(maxfu)%in%id]==ids]}}
		datanew <- data
		datanew[,names(data)%in%id] <- idnum
		datanew$first <- 0
		mintime <- tapply(datanew[,names(datanew)%in%time],datanew[,names(datanew)%in%id],min)
		datanew$mintime <- mintime[datanew[,names(datanew)%in%id]]
		datanew$first <- as.numeric(datanew[,names(datanew)%in%time]==datanew$mintime)
		choosevec[datanew$first==1] <- 1
		}
		datamo <- data[choosevec==1,]

	}
	return(datamo)
}

outputanalfn <- function(fn,data,weights,singleobs,id,time,keep.first,...){
	datamo <- outputation(data,weights,singleobs,id,time,keep.first)
	ans <- fn(datamo,...)
	return(ans)
}


#' Multiple outputation for longitudinal data subject to irregular observation.
#'
#' Multiple outputation is a procedure whereby excess observations are repeatedly randomly sampled and discarded. The method was originally developed to handle clustered data where cluster size is informative, for example when studying pups in a litter. In this case, analysis that ignores cluster size results in larger litters being over-represented in a marginal analysis. Multiple outputation circumvents this problem by randomly selecting one observation per cluster. Multiple outputation has been further adapted to handle longitudinal data subject to irregular observation; here the probability of being retained on any given outputation is inversely proportional to the visit intensity. This function creates multiply outputted datasets, analyses each separately, and combines the results to produce a single estimate.
#' @param noutput the number of outputations to be used
#' @param fn the function to be applied to the outputted datasets. fn should return a vector or scalar; if var=TRUE the second column of fn should be an estimate of standard error.
#' @param data the original dataset on which multiple outputation is to be performed
#' @param weights the weights to be used in the outputation, i.e. the inverse of the probability that a given observation will be selected in creating an outputted dataset. Ignored if singleobs=TRUE
#' @param singleobs logical variable indicating whether a single observation should be retained for each subject
#' @param id character string indicating which column of the data identifies subjects
#' @param time character string indicating which column of the data contains the time at which the visit occurred
#' @param keep.first logical variable indicating whether the first observation should be retained with probability 1. This is useful if the data consists of an observation at baseline followed by follow-up at stochastic time points.
#' @param var logical variable indicating whether fn returns variances in addition to point estimates
#' @param ... other arguments to fn.
#' @return a list containing the multiple outputation estimate of the function fn applied to the data, its standard error, and the relative efficiency of using noutput outputations as opposed to an infinite number
#' @references
#' \itemize{
#' \item Hoffman E, Sen P, Weinberg C. Within-cluster resampling. Biometrika 2001; 88:1121-1134
#' \item Follmann D, Proschan M, Leifer E. Multiple outputation: inference for complex clustered data by averaging analyses from independent data. Biometrics 2003; 59:420-429
#' \item Pullenayegum EM. Multiple outputation for the analysis of longitudinal data subject to irregular observation. Statistics in Medicine (in press)}.
#' @family mo
#' @examples
#' library(nlme)
#' data(Phenobarb)
#' library(survival)
#' library(geepack)
#'
#' Phenobarb$event <- 1-as.numeric(is.na(Phenobarb$conc))
#' data <- lagfn(Phenobarb, lagvars="dose", id="Subject", time="time", lagfirst = 0)
#' data <- lagfn(data, lagvars="dose.lag", id="Subject", time="time", lagfirst = 0)
#' data <- lagfn(data, lagvars="dose.lag.lag", id="Subject", time="time", lagfirst = 0)
#' data$dose.lag[is.na(data$dose.lag)] <- data$dose.lag.lag[is.na(data$dose.lag)]
#' data$dose.lag[is.na(data$dose.lag)] <- data$dose.lag.lag.lag[is.na(data$dose.lag)]
#' data <- data[data$event==1,]
#' data$id <- as.numeric(data$Subject)
#' data <- data[data$time<16*24,]
#' i <- iiw.weights(Surv(time.lag,time,event)~Wt *dose.lag +
#'          dose.lag*(I(conc.lag>0) + conc.lag) +
#'          cluster(Subject),id="id",time="time",event="event",data=data,
#'          invariant="Subject",lagvars=c("time","conc"),maxfu=16*24,
#'          lagfirst=0,first=TRUE)
#' wt <- i$iiw.weight
#' wt[wt>quantile(i$iiw.weight,0.95)] <- quantile(i$iiw.weight,0.95)
#' data$wt <- wt
#' reg <- function(data){
#' est <- summary(geeglm(conc~I(time^3) + log(time), id=id,data=data))$coefficients[,1:2]
#' est <- data.matrix(est)
#' return(est)
#' }
#'
#' mo(20,reg,data,wt,singleobs=FALSE,id="id",time="time",keep.first=FALSE)
#' # On outputation, the dataset contains small numbers of observations per subject
#' # and hence the GEE sandwich variance estimate underestimates variance; this is why
#' # the outputation-based variance estimate fails. This can be remedied by using a
#' # sandwich variance error correction (e.g. Fay-Graubard, Mancl-DeRouen).
#' @export


mo <- function(noutput,fn,data,weights,singleobs,id,time,keep.first,var=TRUE,...){
	dimlast <- as.numeric(var)+1
	for(it in 1:noutput){
		if(it==1){ a <- outputanalfn(fn,data,weights,singleobs,id,time,keep.first,...); if(is.vector(a)) a <- array(a,dim=c(1,length(a))); if(var) ans <- array(dim=c(noutput,dim(a))) else ans <- array(dim=c(noutput,length(a),1)); ans[1,,] <- a}
		if(it>1){ a <- outputanalfn(fn,data,weights,singleobs,id,time,keep.first,...); if(is.vector(a)) a <- array(a,dim=c(1,length(a))); ans[it,,] <- a}

	}
	print(ans)
	pooled <- apply(ans,2:3,mean)
	if(var){
		var.within <- pooled[,2]
		var.between <- apply(array(ans[,,1],dim=c(noutput,dim(a)[1])),2,var)
		var.est <- var.within - ((noutput-1)/noutput)*var.between
		RE.MO <- 1 + (1/noutput)*(var.between/(var.within-var.between))

		return(list(est=pooled[,1],se=sqrt(var.est),RE.MO=RE.MO))
	} else return(list(est=pooled[,1]))

}


#' Fit a semi-parametric joint model
#'
#' Fits a semi-parametric joint model as described by Liang et al. (2009).
#'
#' @details
#' This function fits a semi-parametric joint model as described in Liang (2009), using a frailty model to estimate the parameters in the visit intensity model
#' @param data data frame containing the variables in the model
#' @param Yname character string indicating the column containing the outcome variable
#' @param Xnames vector of character strings indicating the names of the columns of the fixed effects in the outcome regression model
#' @param Wnames vector of character strings indicating the names of the columns of the random effects in the outcome regression model
#' @param Znames vector of character strings indicating the names of the columns of the covariates in the visit intensity model
#' @param formulaobs formula for the observation intensity model
#' @param id character string indicating which column of the data identifies subjects
#' @param time character string indicating which column of the data contains the time at which the visit occurred
#' @param lagvars a vector of variable names corresponding to variables which need to be lagged by one visit to fit the visit intensity model. Typically time will be one of these variables. The function will internally add columns to the data containing the values of the lagged variables from the previous visit. Values of lagged variables for a subject's first visit will be set to NA. To access these variables in specifying the proportional hazards formulae, add ".lag" to the variable you wish to lag. For example, if time is the variable for time, time.lag is the time of the previous visit
#' @param invariant a vector of variable names corresponding to variables in data that are time-invariant. It is not necessary to list every such variable, just those that are invariant and also included in the visit intensity model
#' @param lagfirst A vector giving the value of each lagged variable for the first time within each subject. This is helpful if, for example, time is the variable to be lagged and you know that all subjects entered the study at time zero
#' @param maxfu The maximum follow-up time per subject. If all subjects have the same follow-up time, this can be supplied as a single number. Otherwise, maxfu should be a dataframe with the first column specifying subject identifiers and the second giving the follow-up time for each subject.
#' @param baseline An indicator for whether baseline (time=0) measurements are included by design. Equal to 1 if yes, 0 if no.
#' @param n.knots integer giving the number of knots to use in fitting the frailty model. See documentation for frailtyPenal for more details
#' @param kappa positive smoothing parameter in the penalized likelihood estimation. See documentation for frailtyPenal for more details
#' @return the regression coefficients corresponding to the fixed effects in the outcome regression model.  Closed form expressions for standard errors of the regression coefficients are not available, and Liang et al (2009) recommend obtaining these through bootstrapping.
#' @references Liang Y, Lu W, Ying Z. Joint modelling and analysis of longitudinal data with informative observation times. Biometrics 2009; 65:377-384.
#' @export
#' @examples
#' # replicate simulation in Liang et al.
#' \dontrun{
#' library(data.table)
#' library(survival)
#' datasimi <- function(id){
#' X1 <- runif(1,0,1)
#' X2 <- rbinom(1,1,0.5)
#' Z <- rgamma(1,1,1)
#' Z1 <- rnorm(1,Z-1,1)
#' gamma <- c(0.5,-0.5)
#' beta <- c(1,-1)
#' hazard <- Z*exp(X1/2 - X2/2)
#' C <- runif(1,0,5.8)
#' t <- 0
#' tlast <- t
#' y <- t + X1-X2 + Z1*X2 + rnorm(1,0,1)
#' wait <- rexp(1,hazard)
#' while(tlast+wait<C){
#'   tnew <- tlast+wait
#'     y <- c(y,tnew + X1-X2 + Z1*X2 + rnorm(1,0,1))
#'     t <- c(t,tnew)
#'     tlast <- tnew
#'     wait <- rexp(1,hazard)
#'  }
#'  datai <- list(id=rep(id,length(t)),t=t,y=y,
#'       X1=rep(X1,length(t)),X2=rep(X2,length(t)),C=rep(C,length(t)))
#'  return(datai)
#'  }
#'  sim1 <- function(it,nsubj){
#'  data <- lapply(1:nsubj,datasimi)
#'  data <- as.data.frame(rbindlist(data))
#'  data$event <- 1
#'  C <- tapply(data$C,data$id,mean)
#'  tapply(data$C,data$id,sd)
#'  maxfu <- cbind(1:nsubj,C)
#'  maxfu <- as.data.frame(maxfu)
#'  res <- Liang(data=data, id="id",time="t",Yname="y",
#'             Xnames=c("X1","X2"),
#'             Wnames=c("X2"),Znames=c("X1","X2"), formulaobs=Surv(t.lag,t,event)~X1
#'             + X2+ cluster(id),invariant=c
#'             ("id","X1","X2"),lagvars="t",lagfirst=NA,maxfu=maxfu,
#'             baseline=1,n.knots=6, kappa=10000)
#'  return(res)
#'  }
#'  # change n to 500 to replicate results of Liang et al.
#'  n <- 10
#'  s <- lapply(1:n,sim1,nsubj=200)
#'  smat <- matrix(unlist(s),byrow=TRUE,ncol=2)
#'  apply(smat,2,mean)
#'  }


Liang <- function(data,Yname, Xnames, Wnames, Znames=NULL,formulaobs=NULL, id,time, invariant=NULL,lagvars=NULL,lagfirst=NULL,maxfu,baseline,n.knots=NULL,kappa=NULL ){

  if(is.null(formulaobs)){
    fn <- function(t,tvec) return(which.min(abs(t-tvec)))

    # redo id variable so ids are numbered using consecutive integers
    ids <- names(table(data[,names(data)%in%id]))
    idnum <- array(dim=nrow(data))
    for(i in 1:nrow(data)) idnum[i] <- (1:length(ids))[data[i,names(data)%in%id]==ids]
    if(is.data.frame(maxfu)){ maxfu.use <- maxfu; for(i in 1:nrow(maxfu)){ maxfu.use[i,names(maxfu)%in%id] <- (1:length(ids))[maxfu[i,names(maxfu)%in%id]==ids]}}
    data[,names(data)%in%id] <- idnum

    if(is.null(maxfu)){ maxtable <- tapply(data[,names(data)%in%time],data[,names(data)%in%id],max); maxfu.use <- cbind(1:length(maxtable),maxtable + max(maxtable)*0.001)}


    #	if(identical(ids,as.vector(maxfu[order(maxfu[,1]),1]))==FALSE &identical(as.numeric(ids),as.vector(maxfu[order(maxfu[,1]),1]))==FALSE ) print("Error - ids in maxfu must match ids in 'data'")


    n <- length(table(data[,names(data)%in%id]))
    mi <- tapply(data[,names(data)%in%Yname],data[,names(data)%in%id],length)-baseline

    Xcols <- (1:ncol(data))[is.finite(match(names(data), Xnames))]
    Wcols <- (1:ncol(data))[is.finite(match(names(data), Wnames))]
    X <- array(data.matrix(data[,Xcols]),dim=c(nrow(data),length(Xnames)))
    W <- array(data.matrix(data[,Wcols]),dim=c(nrow(data),length(Wnames)))

    if(length(maxfu)==1) maxfu.use <- cbind(idnum,rep(maxfu,length(idnum)))
    maxfu.use <- maxfu.use[order(maxfu.use[,1]),]
    data <- data[order(idnum),]

    if(length(maxfu)==1){
      Lambdahat <- nrow(data)/n
      sigmahatsq <- max((sum(mi^2)-sum(mi)-n*Lambdahat^2)/(n*Lambdahat^2),0)
      Lambdahat <- rep(Lambdahat,n)
      Ci <- rep(maxfu,n)
    }
    if(length(maxfu)>1){
      maxfu.use <- maxfu.use[order(maxfu[,1]),]
      ids <- as.numeric(names(table(data[,names(data)%in%id])))
      Ci <- as.vector(maxfu.use[order(maxfu.use[,1]),2])
      data$event <- 1

      lagcols <- (1:ncol(data))[is.finite(match(names(data), time))]
      invarcols <- (1:ncol(data))[is.finite(match(names(data), id))]

      datacox <- addcensoredrows(data=data,maxfu=maxfu.use,tinvarcols=invarcols,id=id,time=time,event="event")
      datacox <- lagfn(datacox,"time",id,time)

      formulanull <- Surv(time.lag,time,event)~1
      datacox <- datacox[datacox[,names(datacox)%in%time]>0,]
      b <- basehaz(coxph(formulanull,data=datacox))
      indexfnnocov <- function(t,time){ return(sum(time<t))}
      bindex <- sapply(Ci,indexfnnocov,time=b$time)
      bindex[bindex==0] <- 1
      Lambdahat <- b$hazard[bindex]
      sigmahatsq <- max((sum(mi^2)-sum(mi)-sum(Lambdahat^2))/sum(Lambdahat^2),0)
    }

    mi.Lambdahat <- mi/Lambdahat
    mi.Lambdahat[mi==0 & Lambdahat==0] <- 1

    Bhat <- array(dim=c(nrow(data),ncol(W)))
    Bbar <- Bhat
    Xbar <- array(dim=c(nrow(data),ncol(X)))

    Bmultiplier <- array(dim=nrow(data))
    Bmultid <- (mi - Lambdahat)*sigmahatsq/(1+Lambdahat*sigmahatsq)
    ids <- as.numeric(names(table(ids)))
    for(i in 1:n) Bmultiplier[data[,names(data)%in%id]==ids[i]] <- Bmultid[i]
    Bhat <- sweep(array(W,dim=c(nrow(data),ncol(W))),1,Bmultiplier,"*")

    Xbar <- array(dim = c(nrow(data),ncol(X)))
    Bbar <- array(dim = c(nrow(data),ncol(W)))

    #	X0 <- array(apply(X,data[,names(data)%in%id],mean),dim=c(n,ncol(X))

    #	for(row in 1:nrow(data)){
    #		t <- data[,names(data)%in%time][row]
    #		Xbar[row,] <- apply(sweep(X0,1,mi*as.numeric(Ci>=t)/Lambdahat,"*"),2,mean)/
    #	}

    obsnum <- rep(1,nrow(data))
    for(row in 2:nrow(data)){
      if(identical(data[row,names(data)%in%id],data[row-1,names(data)%in%id])) obsnum[row] <- obsnum[row-1]+1
    }
    firstobs <- (1:nrow(data))[obsnum==1]

    for (row in 1:nrow(data)) {
      t <- data[,names(data)%in%time][row]
      #        	absdiff <- abs(data[,names(data)%in%time] - t)
      #        	min <- tapply(absdiff, data[,names(data)%in%id], min)
      #        	data.min <- min[data[,names(data)%in%id]]
      closest <- tapply(data[,names(data)%in%time],data[,names(data)%in%id],fn,t=data[,names(data)%in%time][row])
      userow <- firstobs + closest-1
      Xbar[row,] <- apply(sweep(array(X[userow,],dim=c(n,length(Xnames))),1, (mi.Lambdahat*as.numeric(Ci>=t)),"*"),2,sum)/sum((mi.Lambdahat*as.numeric(Ci>=t)))
      Bbar[row,] <- apply(sweep(array(Bhat[userow,],dim=c(n,length(Wnames))),1, (mi.Lambdahat*as.numeric(Ci>=t)),"*"),2,sum)/sum((mi.Lambdahat*as.numeric(Ci>=t)))
      #        	Bbar[row,] <- apply(sweep(array(Bhat[data.min - absdiff == 0,],dim=c(n,length(Wnames))),1, (mi.Lambdahat*as.numeric(Ci>=t))[data[,names(data)%in%id]][data.min -
      #            absdiff == 0],"*"),2,sum)/sum((mi.Lambdahat*as.numeric(Ci>=t))[data[,names(data)%in%id][data.min - absdiff ==
      #            0]])
    }

    regX <- array((X - Xbar),dim=c(nrow(data),ncol(X)))[data[,names(data)%in%time]>0,]
    regB <- array(Bhat - Bbar,dim=c(nrow(data),ncol(W)))[data[,names(data)%in%time]>0,]
    regY <- data[,names(data)%in%Yname][data[,names(data)%in%time]>0]
    regpredictor <- cbind(regX,regB)
    if(sigmahatsq>0) beta <- solve(t(regpredictor)%*%regpredictor,t(regpredictor)%*%regY)
    if(sigmahatsq==0) beta <- c(solve(t(regX)%*%regX,t(regX)%*%regY),NA)

  }


  if(!is.null(formulaobs)){
    fn <- function(t,tvec) return(which.min(abs(t-tvec)))

    # redo id variable so ids are numbered using consecutive integers
    ids <- names(table(data[,names(data)%in%id]))
    idnum <- array(dim=nrow(data))
    for(i in 1:nrow(data)) idnum[i] <- (1:length(ids))[data[i,names(data)%in%id]==ids]
    if(is.data.frame(maxfu)){ maxfu.use <- maxfu; for(i in 1:nrow(maxfu)){ maxfu.use[i,names(maxfu)%in%id] <- (1:length(ids))[maxfu[i,names(maxfu)%in%id]==ids]}}
    data[,names(data)%in%id] <- idnum

    if(is.null(maxfu)){ maxtable <- tapply(data[,names(data)%in%time],data[,names(data)%in%id],max); maxfu.use <- cbind(1:length(maxtable),maxtable + max(maxtable)*0.001)}


    #	if(identical(ids,as.vector(maxfu[order(maxfu[,1]),1]))==FALSE &identical(as.numeric(ids),as.vector(maxfu[order(maxfu[,1]),1]))==FALSE ) print("Error - ids in maxfu must match ids in 'data'")


    n <- length(table(data[,names(data)%in%id]))
    mi <- tapply(data[,names(data)%in%Yname],data[,names(data)%in%id],length)-baseline

    Xcols <- (1:ncol(data))[is.finite(match(names(data), Xnames))]
    Wcols <- (1:ncol(data))[is.finite(match(names(data), Wnames))]
    Zcols <- (1:ncol(data))[is.finite(match(names(data), Znames))]

    X <- array(data.matrix(data[,Xcols]),dim=c(nrow(data),length(Xnames)))
    W <- array(data.matrix(data[,Wcols]),dim=c(nrow(data),length(Wnames)))
    Z <- array(data.matrix(data[,Zcols]),dim=c(nrow(data),length(Znames)))

    if(length(maxfu)==1) maxfu.use <- cbind(idnum,rep(maxfu,length(idnum)))

    maxfu.use <- maxfu.use[order(maxfu.use[,1]),]
    data <- data[order(idnum),]

    if(length(maxfu)==1){
      Ci <- rep(maxfu,n)

      maxfu.use <- maxfu.use[order(maxfu.use[,1]),]
      ids <- as.numeric(names(table(data[,names(data)%in%id])))
      data$event <- 1
    }
    if(length(maxfu)>1){
      maxfu.use <- maxfu.use[order(maxfu[,1]),]
      ids <- as.numeric(names(table(data[,names(data)%in%id])))
      Ci <- as.vector(maxfu.use[order(maxfu.use[,1]),2])
      data$event <- 1
    }
    lagcols <- (1:ncol(data))[is.finite(match(names(data), time))]
    #    invarcols <- (1:ncol(data))[is.finite(match(names(data), id))]


    datacox <- addcensoredrows(data=data,maxfu=maxfu,tinvarcols=invariant,id=id,time=time,event="event")
    datacox <- lagfn(datacox,lagvars,id,time,lagfirst=lagfirst)

    formulacov <- Surv(time.lag,time,event)~Znames + cluster(id)
    formulacox <- Surv(time.lag,time,event)~Znames + (1|id)
    row.names(datacox)[datacox$event==1] <- row.names(data)
    use <- (1:nrow(datacox))[(!is.na(datacox[,names(datacox)%in%paste(time,".lag",sep="")]))]
    datacoxuse <- datacox[use,]
    use <- (1:nrow(datacoxuse))[datacoxuse[,names(datacoxuse)%in%time]-datacoxuse[,names(datacoxuse)%in%paste(time,".lag",sep="")]>0]
    datacoxuse <- datacoxuse[use,]
    row.names(datacoxuse) <- 1:nrow(datacoxuse)
    mcov <- frailtyPenal(formula=formulaobs,data=datacoxuse, recurrentAG=TRUE, n.knots=n.knots,kappa=kappa,cross.validation=TRUE,print.times=FALSE)

    omit <- c("(Intercept)",paste("cluster(",id,")",sep=""))
    #	Zmat <- model.matrix(terms(formula(mcov)),data=datacox)

    #	columns <- colnames(Zmat)
    #	if(dim(Zmat)[2]>3){ data$lp[row.names(data)%in%row.names(Zmat)] <- exp(-Zmat[row.names(Zmat)%in%row.names(data),!columns%in%omit]%*%mcov$coef)} else{
    #	  data$lp[row.names(data)%in%row.names(Zmat)] <- exp(-Zmat[row.names(Zmat)%in%row.names(data),!columns%in%omit]*mcov$coef)
    #	}
    data$lp <- Z%*%mcov$coef



    indexfn <- function(x) return(sum(as.numeric(mcov$x<=x)))
    indexfns <- function(x) return(sapply(x,indexfn))
    Hazard0 <- -log(mcov$surv[,1,1])[indexfns(Ci)]
    Hazard0[Hazard0==0] <- -log(mcov$surv[,1,1])[2]
    lp <- tapply(data$lp,data[,names(data)%in%id],mean,na.rm=TRUE)
    Lambdahat <- Hazard0*exp(lp)

    sigmahatsq <- mcov$theta
       print(paste("sigmahatsq=",sigmahatsq,sep=""))

    #  mi.Lambdahat <- mi/Lambdahat
    mi.Lambdahat <- mi/Hazard0
    mi.Lambdahat[mi==0 & Lambdahat==0] <- 1

    Bhat <- array(dim=c(nrow(data),ncol(W)))
    Bbar <- Bhat
    Xbar <- array(dim=c(nrow(data),ncol(X)))

    Bmultiplier <- array(dim=nrow(data))
    Bmultid <- (mi - Lambdahat)*sigmahatsq/(1+Lambdahat*sigmahatsq)
    ids <- as.numeric(names(table(ids)))
    for(i in 1:n) Bmultiplier[data[,names(data)%in%id]==ids[i]] <- Bmultid[i]
    Bhat <- sweep(array(W,dim=c(nrow(data),ncol(W))),1,Bmultiplier,"*")

    Xbar <- array(dim = c(nrow(data),ncol(X)))
    Bbar <- array(dim = c(nrow(data),ncol(W)))

    #	X0 <- array(apply(X,data[,names(data)%in%id],mean),dim=c(n,ncol(X))

    #	for(row in 1:nrow(data)){
    #		t <- data[,names(data)%in%time][row]
    #		Xbar[row,] <- apply(sweep(X0,1,mi*as.numeric(Ci>=t)/Lambdahat,"*"),2,mean)/
    #	}

    obsnum <- rep(1,nrow(data))
    for(row in 2:nrow(data)){
      if(identical(data[row,names(data)%in%id],data[row-1,names(data)%in%id])) obsnum[row] <- obsnum[row-1]+1
    }
    firstobs <- (1:nrow(data))[obsnum==1]

    for (row in 1:nrow(data)) {
      t <- data[,names(data)%in%time][row]
      #        	absdiff <- abs(data[,names(data)%in%time] - t)
      #        	min <- tapply(absdiff, data[,names(data)%in%id], min)
      #        	data.min <- min[data[,names(data)%in%id]]
      closest <- tapply(data[,names(data)%in%time],data[,names(data)%in%id],fn,t=data[,names(data)%in%time][row])
      userow <- firstobs + closest-1
      Xbar[row,] <- apply(sweep(array(X[userow,],dim=c(n,length(Xnames))),1, (mi.Lambdahat*as.numeric(Ci>=t)),"*"),2,sum)/sum((mi.Lambdahat*as.numeric(Ci>=t)))
      Bbar[row,] <- apply(sweep(array(Bhat[userow,],dim=c(n,length(Wnames))),1, (mi.Lambdahat*as.numeric(Ci>=t)),"*"),2,sum)/sum((mi.Lambdahat*as.numeric(Ci>=t)))
      #        	Bbar[row,] <- apply(sweep(array(Bhat[data.min - absdiff == 0,],dim=c(n,length(Wnames))),1, (mi.Lambdahat*as.numeric(Ci>=t))[data[,names(data)%in%id]][data.min -
      #            absdiff == 0],"*"),2,sum)/sum((mi.Lambdahat*as.numeric(Ci>=t))[data[,names(data)%in%id][data.min - absdiff ==
      #            0]])
    }

    regX <- array((X - Xbar),dim=c(nrow(data),ncol(X)))[data[,names(data)%in%time]>0,]
    regB <- array(Bhat - Bbar,dim=c(nrow(data),ncol(W)))[data[,names(data)%in%time]>0,]
    regY <- data[,names(data)%in%Yname][data[,names(data)%in%time]>0]
    regpredictor <- cbind(regX,regB)
    if(sigmahatsq>0) beta <- solve(t(regpredictor)%*%regpredictor,t(regpredictor)%*%regY)
    if(sigmahatsq==0) beta <- c(solve(t(regX)%*%regX,t(regX)%*%regY),NA)
  }
  return(c(beta,sigmahatsq))
}



#' Create an abacus plot
#'Creates an abacus plot, depicting visits per subject over time
#'
#' @details This function creates a plot for n randomly sampled individuals from the supplied dataset, with one row per subject and one point per visit. This can be useful for visualising the extent of irregularity in the visit process. For example, with perfect repeated measures data (i.e., no irregularity), the points will line up vertically. With greater irregularity, the points will be randomly scattered over time.
#' @param n the number of subjects to randomly sample. Subjects are sampled without replacement and therefore n must be smaller than the total number of subjects in the dataset
#' @param time character string indicating which column of the data contains the time at which the visit occurred
#' @param id character string indicating which column of the data identifies subjects
#' @param data data frame containing the variables in the model
#' @param tmin the smallest time to include on the x-axis
#' @param tmax the largest time to include on the x-axis
#' @param xlab.abacus the label for the x-axis
#' @param ylab.abacus the label for the y-axis
#' @param pch.abacus the plotting character for the points on the abacus plot
#' @param col.abacus the colour of the rails on the abacus plot
#' @return produces a plot depicting observation times for each subject. No values are returned
#' @examples
#' library(MEMSS)
#' data(Phenobarb)
#' Phenobarb$event <- 1-as.numeric(is.na(Phenobarb$conc))
#' data <- Phenobarb[Phenobarb$event==1,]
#' abacus.plot(n=20,time="time",id="Subject",data=data,tmin=0,tmax=16*24,
#' xlab.abacus="Time in hours",pch=16,col.abacus=gray(0.8))
#' @export


abacus.plot <- function(n,time,id,data,tmin,tmax,xlab.abacus="Time",ylab.abacus="Subject",pch.abacus=16,col.abacus=1){
  use.ids <- sample(levels(factor(data[,names(data)%in%id])),n,replace=FALSE)
  use.data <- data[data[,names(data)%in%id]%in%use.ids,]
  use.data$newids <- NA
  for(i in 1:n) use.data$newids[use.data[,names(use.data)%in%id]%in%use.ids[i]] <- i
  plot(use.data[,names(use.data)%in%time],use.data$newids,xlim=c(tmin,tmax),pch=pch.abacus,xlab=xlab.abacus,ylab=ylab.abacus,type="n")
  segments(rep(tmin,nrow(use.data)),use.data$newids,rep(tmax,nrow(use.data)),use.data$newids,col=col.abacus)
  points(use.data[,names(use.data)%in%time],use.data$newids,pch=pch.abacus)
}


