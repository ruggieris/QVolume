############
# 
# QVolume: Estimating the Total Volume of Queries to a Search Engine
#
# Authors: F. Lillo - S. Ruggieri (c) 2019
#
############

# clean all
rm(list=ls())

# definitions and settings
source("qvolume.R")

### constants
n = 4000 # sample size
delta = 1.232406 # delta ratio in binned data
l1 = 10 # l_1 value in binned data
cat('Sample size:', n, '\n')

### sampling from population x
sample.unif = sample(x, n) # uniform
z = sample_geom(n, N)
sample.len = length(z)
stopifnot(sample.len==n)
sample.nunif = x[z] # non-uniform
noise = rnorm(sample.len, 1, 0.1/3)
sample.noisy = pmax(V_N, sample.nunif * noise) # noisy, lower bounded by V_N
sample.sketchy = sample.nunif + 0.001*c*runif(sample.len) # sketchy

### select population or sample
#data = sort(x, decreasing=TRUE) # whole population, or
data = sort(sample.nunif, decreasing=TRUE) # sample
plot(data, log="xy")


### continuous data
cat('CONTINUOUS DATA\n')
# CSN+Max estimator
fit_csn = zipf.csn_max(data) 
to = fit_csn[3]
zipf.plot(data, parameters=fit_csn, to=to)
# NLS estimator
fit_nls = zipf.nls(data[1:to])
zipf.plot(data, parameters=fit_nls, to=to)
fit = zipf.nls(data[1:to], returnfit=TRUE)

sampleV = sum(data)
cat('volume of queries in sample (M) = ', round(sampleV/1e+06, 2), '\n')

# parameters & stderrors from GSC
fitcoeff = summary(fit)$coefficients
hatb = -fitcoeff[2, 1]
hatc = fitcoeff[1, 1]
db = fitcoeff[2, 2]
dc = fitcoeff[1, 2]

cat('b =', hatb, '\n')
cat('c =', hatc, '\n')
cat('db =', db, '\n')
cat('dc =', dc, '\n')

vs = c(1e+01,1e+02, 1e+03, 1e+04)
cat('Number of queries searched\n')
for(vn in vs) {
  Nv = est.n(hatc, hatb, vn) 
  errNv = error.n(Nv, hatc, hatb, dc, db, vn)[2]
  truev = sum(x>=vn)
  cat(' at least:', vn, ' =', round(Nv), '+/-', round(errNv), '[', truev, ']\n')
}
cat('Volume (M) of queries searched\n')
for(vn in vs) {
  Vv =  est.V(hatc, hatb, vn) 
  errVv = round(error.V(Vv, hatc, hatb, dc, db, vn)[2])
  truev = sum(x[x>=vn])
  cat(' at least:', vn, ' =', round(Vv/1e+06, 2), '+/-', round(errVv/1e+06, 2), '[', round(truev/1e+06, 2), ']\n')
}


### binned data
cat('BINNED DATA\n')
jv = pmax(0, 1+floor(log(data/l1, base=delta)))
data.binned = l1*delta^jv
plot(data.binned, log="xy")

fit_csn = zipf.csn_max(data.binned, binned=TRUE)
fit_csn
to = fit_csn["to"]
l0 = fit_csn["l0"]
zipf.plot(data.binned, parameters=fit_csn, to=to)

fit_chi2 = zipf.nls(data.binned[1:to], binned=TRUE, l0=l0, returnfit=TRUE)
fit_chi2
zipf.plot(data.binned, parameters=fit_chi2, to=to)

sampleV = sum(data)
cat('volume of queries in sample (M) = ', round(sampleV/1e+06, 2), '\n')

# parameters & stderrors
hatb = -fit_chi2["coefficient"]
hatc = fit_chi2["intercept"]
db = fit_chi2["stderr.b"]
dc = fit_chi2["stderr.c"]

cat('b =', hatb, '\n')
cat('c =', hatc, '\n')
cat('db =', db, '\n')
cat('dc =', dc, '\n')

vs = c(1e+01,1e+02, 1e+03, 1e+04)
cat('Number of queries searched\n')
for(vn in vs) {
  Nv = est.n(hatc, hatb, vn) # estimation is per year
  errNv = error.n(Nv, hatc, hatb, dc, db, vn)[2]
  truev = sum(x>=vn)
  cat(' at least:', vn, ' =', round(Nv), '+/-', round(errNv), '[', truev, ']\n')
}
cat('Volume (M) of queries searched\n')
for(vn in vs) {
  Vv =  est.V(hatc, hatb, vn) # estimation is per year
  errVv = round(error.V(Vv, hatc, hatb, dc, db, vn)[2])
  truev = sum(x[x>=vn])
  cat(' at least:', vn, ' =', round(Vv/1e+06, 2), '+/-', round(errVv/1e+06, 2), '[', round(truev/1e+06, 2), ']\n')
}
