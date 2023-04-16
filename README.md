# QVolume: Total Volume of Queries to a Search Engine
**Fabrizio Lillo**, Department  of  Mathematics,  University  of  Bologna, fabrizio.lillo@unibo.it  
**Salvatore Ruggieri**, Department of Computer Science, University of Pisa, Italy, salvatore.ruggieri@unipi.it

QVolume is a collection of R functions to estimate the total number of searches (volume) of queries in a specific domain which were submitted to a search engine in a given time period. Our statistical model assumes a Zipf's law distribution of the query population and three possible sampling bias scenarios: non-uniform sampling, noisy sampling, and sketchy sampling. A few estimators of the parameters of the distribution are devised and experimented in [[1, 2]](#references), based on the nature of the empirical data. For continuous data, we recommend using nonlinear least square regression (NLS) on the top-volume queries, where the bound on the volume is obtained from the well-known Clauset, Shalizi and Newman (CSN) estimation of power-law parameters. For discrete data, we propose using a Chi-square minimization approach restricted to the top-volume queries, where the bound is obtained by the binned version of the CSN method. Estimations are then derived for the total number of queries and for the total volume of the population, including statistical error bounds. 

## References

[1] F. Lillo, S. Ruggieri. [Estimating the Total Volume of Queries to Google](http://pages.di.unipi.it/ruggieri/Papers/www2019.pdf). 28th World Wide Web Conference on World Wide Web (WebConf 2019) : 1051-1060. ACM, May 2019.

[2] F. Lillo, S. Ruggieri. [Estimating the Total Volume of Queries to a Search Engine](https://arxiv.org/abs/2101.09807). IEEE Transactions on Knowledge and Data Engineering. Vol. 34 Issue 11, November 2022, 5351-5363.
