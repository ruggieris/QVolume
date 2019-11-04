############
# 
# QVolume: Estimating the Total Volume of Queries to a Search Engine
#
# Authors: F. Lillo - S. Ruggieri (c) 2019
#
############

#' load packages (and install them, if needed)
llibrary = function(packages) {
  new.packages = packages[!(packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  for( p in packages) {
    require(p, character.only = TRUE)
  }  
}

# required packages
llibrary( c("poweRlaw", "sfsmisc", "reticulate") ) 

# repeatability
set.seed(0)

#' simulation settings
c = 10^5
b = 0.7745
N = 10^6
x = c/(c(1:N))^b
V = sum(x)
V_N = min(x)

#' colors in plots
redcol = "#e41a1c"
bluecol = "#377eb8"
greencol = "#4daf4a"
greycol = 8

#' sample of n distinct index positions up to N
sample_geom_slow = function(n, N) {
  sample = c()
  while( length(sample) < n) {
    new.sample = unique(rgeom(log(n)*n,1/1000)+1)
    new.sample = new.sample[which(new.sample<=N)]
    sample = unique(c(sample, new.sample))
  }
  sample[1:n]
}

#' sample of n distinct index positions up to N
sample_geom = function(n, N) {
  sample = c()
  i = 0
  base = 1
  ls = 0
  while( ls < n) {
    if(ls>0) {
      base = max(which(sort(sample) == 1:ls))+1
    }
    new.sample = unique(rgeom(log(n)*n,1/1000)+base)
    new.sample = new.sample[which(new.sample<=N)]
    sample = unique(c(sample, new.sample))
    i = i + 1
    ls = length(sample)
  }
  #cat('niter=', i, '\n')
  sample[1:n]
}

#' VGAM::zeta(s, shift=N+1) requires s >= 1
#' Solution: use zeta from Python - require(reticulate)
mpmath = import("mpmath")
#' Hurwitz zetaf 
hzeta = function(s, deriv=0, shift=1) {
  if(s==1) return(Inf)
  as.numeric(as.character(mpmath$zeta(s, shift, deriv)))
}

#' zetaf(b, N) = zeta(b) - zeta(b, shift=N+1)
zetaf = function(s, N) {
  if(s==1) return(sum( (1:N)**(-s) ))
  hzeta(s) - hzeta(s, shift=N+1)
}

#' Zipf's law parameter estimation: using linear regression on log-scale
#' 
#' @param y data ordered decreasing.
#' @return array with estimated coefficient and intercept 
zipf.reg = function(y) {
  to = length(y)
  x = 1:to
  fit = lm(log(y)~log(x))
  lm_coef = coef(fit)
  intercept = exp(lm_coef[1])
  coefficient = lm_coef[2]
  c(coefficient=coefficient, intercept=intercept)
}

#' Chi-Square for constrOptim
chi2 = function(x, bins, nobs) 
{
  nexp = -diff( (x[2]/bins)^(1/x[1]) )
  return( sum((nexp-nobs)^2/nexp) )
} 

#' Hessian function required for stderr calculation in zipf.chi2
zipf.chi2.hessian = deriv(~(((c/l1)^(1/b)-(c/l2)^(1/b)) - nobs)/sqrt((c/l1)^(1/b)-(c/l2)^(1/b)), c("b", "c"), function(b, c, nobs, l1, l2){} )
  
#' Zipf's law parameter estimation: Chi-Square (binned case)
#' 
#' @param y data ordered decreasing.
#' @param l0 lower bound of first bin
#' @return array with estimated coefficient and intercept
zipf.chi2 = function(y, l0=0, se=FALSE) {
  tb = table(y)
  bins = c(l0, as.numeric(names(tb))) # bins
  ubins = bins[2:length(bins)] # upper bounds
  # observed objects per bin
  nobs = array(tb)
  opt = constrOptim(c(1,500000), chi2, NULL, matrix(c(1,0,0,1), nrow=2, ncol=2), c(0,0), bins=bins, nobs=nobs)
  coefficient = -opt$par[1]
  intercept = opt$par[2]
  if(se) {
    lbins = bins[1:(length(bins)-1)] # lower bounds
    fr = zipf.chi2.hessian(-coefficient, intercept, nobs, lbins, ubins) # positive coefficient
    G = attr(fr,"gradient")
    df = length(nobs)-2
    sigma.est = sqrt(opt$value/df)
    V.lin = sigma.est^2*solve(t(G)%*%G)
    se.lin = sqrt(diag(V.lin))
    res = c(coefficient=coefficient, intercept=intercept, stderr=se.lin[1], stderr=se.lin[2]) 
  } else {
    res = c(coefficient=coefficient, intercept=intercept) 
  }
  res
}

#' Zipf's law parameter estimation: using non-linear least squares (NLS) or Chi-Square (binned case)
#' 
#' @param y data ordered decreasing.
#' @param returnfit boolean.
#' @param binned boolean.
#' @param l0 lower bound of first bin, if binned==TRUE
#' @return array with estimated coefficient and intercept, or nls() output if returnfit=TRUE and binned=FALSE
zipf.nls = function(y, returnfit=FALSE, binned=FALSE, l0=0) {
  if(binned) 
    return(zipf.chi2(y, l0=l0, se=returnfit))
  to = length(y)
  x = 1:to
  dataf = as.data.frame(matrix(cbind(x,y), nrow=to, ncol=2))
  reg = zipf.reg(y) # initial parameters from linear regression in log-scale
  fit = nls(y~c*x^a,data=dataf, start=list(c=reg[2],a=reg[1]), control = list(maxiter=500, warnOnly=T))
  lm_coef = coef(fit)
  intercept = lm_coef[1]
  coefficient = lm_coef[2]
  if(returnfit)
    fit
  else
    c(coefficient=coefficient, intercept=intercept)
}

#' MLE for binned CSN
mle.csn.binned = function(alpha, u, l, z, bmin){
  return( -sum( z*( log((l)^(1-alpha) - (u)^(1-alpha)) + (alpha-1)*log(bmin) ) ))
}

#' constrained Chi-Square
chi2csn = function(x, bin, nobs, bb){
  nexp = -diff( (x[1]/bin)^(1/bb) )
  return( sum((nexp-nobs)^2/nexp) )
}

#' Zipf's law parameter estimation: using CSN-binned + constrained Chi-Squared (binned data)
#' 
#' @param y data ordered decreasing.
#' @param l0 lower bound first bin. .
#' @return array with estimated coefficient, intercept, number of fitted values, and new lower bound.
zipf.csn.binned = function(y, l0=0) {
  # CSN-binned (for exponential bins)
  tb = table(y)
  ubounds = as.numeric(names(tb)) # upper bounds
  nb = length(ubounds) # number of bins
  lbounds = c(l0, ubounds[-nb]) # lower bounds
  nobs = array(tb) # counts
  n = sum(nobs)
  dat = matrix(nrow=nb-1,ncol=2)
  for(i in 1:(nb-1)) {
    bmin = lbounds[i]
    u = ubounds[i:nb]
    l = lbounds[i:nb]
    z = nobs[i:nb]
    es = optimize(mle.csn.binned, c(1,100), u=u, l=l, z=z, bmin=bmin)
    al = es$minimum
    temp = cumsum(z[length(z):1])/sum(z)
    cx = 1-temp[length(temp):1]
    cf = 1-(l/bmin)^(1-al)
    dat[i,1] = al
    dat[i,2] = max(abs(cf-cx))
  }
  j = which.min(dat[,2])
  xmin = lbounds[j]
  alpha = dat[j,1]
  # from power law to zipf
  totxmin = sum(nobs[j:nb])
  coefficient = -1/(alpha-1)
  # constrained Chi-Square
  bins = c(lbounds[j:nb], ubounds[nb])
  n.objs = nobs[j:nb]
  res = constrOptim(c(10000), chi2csn, NULL, matrix(c(1),nrow=1,ncol=1), c(0), bin=bins, nobs=n.objs, bb=-coefficient)
  new.l0 = lbounds[j]
  intercept = res$par[1]
  c(coefficient=coefficient, intercept=intercept, to=totxmin, l0=new.l0)
}

#' Zipf's law parameter estimation: using CSN method + MAX estimator (continuous data) and CSN-binned + constrained Chi-Squared (binned data)
#' 
#' @param y data ordered decreasing.
#' @param binned is data binned?
#' @param l0 option when data is binned
#' @return array with estimated coefficient, intercept, and number of fitted values.
zipf.csn_max = function(y, binned=FALSE, l0=0) {
  if(binned) 
    return(zipf.csn.binned(y, l0=l0))
  m_pl = conpl$new(y) 
  # continuous powerlaw object
  est = estimate_xmin(m_pl) 
  m_pl$setXmin(est)
  xmin = m_pl$xmin
  alpha = m_pl$pars
  # from power law to zipf
  totxmin = est$ntail
  coefficient = -1/(alpha-1) # "-" is to have the coefficient positive
  intercept = max(y)
  c(coefficient=coefficient, intercept=intercept, to=totxmin)
}

#' estimator of N
est.n = function(c, b, vn=1) (c/vn)**(1/b)

#' estimator of V
est.V = function(c, b, vn=1) c*zetaf(b, est.n(c, b, vn))

#' error DELTA-N
error.n = function(N, c, b, dc, db, V_N=1)
{
  a1 = N/(b*c)*dc
  a2 = -(N/b^2)*log(c/V_N)*db
  c(sqrt(a1^2+a2^2),abs(a1)+abs(a2))
}

#' error DELTA-V - Hurwitz derivative is not available in R, taken from Python (see hzeta())
error.V = function(V, c, b, dc, db, V_N)
{
  N=est.n(c, b, V_N)
  a1=(V/c+N*hzeta(b+1,shift=N+1))*dc
  a2=(c*(hzeta(b,deriv=1)-hzeta(b,deriv=1,shift=N+1)-(N*log(c/V_N)*hzeta(b+1,shift=N+1))/b))*db
  c(sqrt(a1^2+a2^2),abs(a1)+abs(a2))
}

#' Zipf's law fitting plot
zipf.plot = function(y, main=NULL, ylim=NULL, parameters=NULL, to=Inf) {
  par(mar=c(4,4,1.2,1))
  if(!is.null(ylim))
    plot(y, log="xy", ylim=ylim, xlab="", ylab="", las=1, xaxt="n", yaxt="n")
  else
    plot(y, log="xy", xlab="", ylab="", las=1, xaxt="n", yaxt="n")
  eaxis(1, n.axp=1)
  eaxis(2, n.axp=1)
  mtext(text = "rank", side = 1, line = 2, font=2)
  mtext(text = "volume", side = 2, line = 2.5, font=2)
  if(!is.null(main))
    title(main)
  if(!is.null(parameters)) {
    coefficient = parameters[1]
    intercept = parameters[2]
    curve( x**coefficient*intercept, to=to, add=T, col=redcol, lwd = 2)
    r_inte = round(log(intercept), 4) 
    r_coef = round(coefficient, 4) 
    mtext(bquote(y == .(r_coef) * x + .(r_inte)),  adj=1, padj=0, col=redcol) # display equation 
    mtext(paste('rank max =', to),  adj=0, padj=0, col=redcol)
  }
}

