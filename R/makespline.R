#############################################################################
#Description: Fits a spline to the data (x,y) using penalized constrained least squares
#Notes: 
#Author: Matthew Young
#Date Modified: 7/4/2015

makespline=function (x, y, newX=NULL, nKnots = 6, lower_bound=10^-3){
	#Should not be used outside of goseq package.  Refer to the help pages for pcls in the mgcv package for more general 
	#contstrained spline fitting.
	#We handle montonocily decreasing problems by reformulating them as monotonicly increasing by reflecting about "infinity"
	#Compare the first 10% to the last 10%
	ww=order(x)
	size=ceiling(length(y)/10) 
	low=sum(y[ww][1:size]) 
	hi=sum(y[ww][(length(y)-size):length(y)]) 
	#Is the trend decreasing
	if(hi<=low){
		#Reform as a montonicly increasing fit by reflecting about infitinity
		reflectionFactor=10^10
		x=reflectionFactor-x
	}
	#The number of knots to use when generating the spline
	nKnots <- round(nKnots)
	if (is.null(newX)) 
		newX <- x
	#We need to force the fit to go through (0,0), mono.con constructs a requirment that the LOWEST x VALUE be greater than 0
	#so to force the fit through (0,0) we add a fake lowest x value of 0 by adding in a dummy point. 
	#We use an inequality constraint rather than equality to only skew the fit enough that the spline is always >0
	#This won't work if we have a monotonically decreasing function (and we don't need to it), so check that
	if(hi>low){
	        x=c(0,x)
	        y=c(0,y)
	}
	f.ug <- gam(y ~ s(x, k = nKnots, bs = "cr"),family=binomial())
	dat <- data.frame(x = x, y = y)
	sm <- smoothCon(s(x, k = nKnots, bs = "cr"), dat, knots = NULL)[[1]]
	if (length(sm$xp) < 6) 
		warning("Few than 6 nKnots were specified!\n Possible inaccurate fitting!\n")
	Fmonocon <- mono.con(sm$xp,TRUE)
	# The lower and upper parameters here don't seem to build the right constraints 100% of the time, so add them manually
	# ND-6/14- removing the upper limit constraint because it wasn't satisfied by initial parameters.
	# Hope it doesn't break the code.
	Fmonocon$A=rbind(Fmonocon$A,c(1,rep(0,ncol(Fmonocon$A)-1))) #,c(rep(0,ncol(F$A)-1)),-1))
	Fmonocon$b=c(Fmonocon$b,0) #,1) ND - should this last one be -1?
	G <- list(X = sm$X, C = matrix(0, 0, 0), sp = f.ug$sp, p = sm$xp,
	    y = y, w = y * 0 + 1, Ain =Fmonocon$A, bin = Fmonocon$b, S = sm$S, off = 0)
	# do some error checking as pcls can return NaNs sometimes if this isn't true
	if( any((G$Ain%*%G$p - G$bin) < 0 )){
	   show(G$p)
	   stop("Constraints for spline fit not satisfied by initial parameters")
	}
	#Do the actual fitting with penalized contstrained least squares
	p <- pcls(G)

	# Do some error checking for the rare case that p returns "Na"
	if(any(is.na(p))){
	   warning("The splines fit failed! Returning a flat weight distribution. This may compromise your GO analysis, please check what happened.")
	   return(rep(mean(y),length(newX)))
	}

	#Now what we want is the value of the spline at each data point x,
	fv <- Predict.matrix(sm, data.frame(x = newX)) %*% p
	fv <- as.vector(fv)

	#If pcls still fails for some reason and returns negative values (or values that are so low that they'll have an effective relative weight of Inf
	#then we need to add a little bit to every weight to make things non-negative, the choice of minimum weight is somewhat arbitrary and serves only to
	#ensure that we don't have positive weights that will be infinitely prefered as a ratio.
	if(min(fv)<lower_bound)
		fv=fv-min(fv)+lower_bound
	return(fv)
}
