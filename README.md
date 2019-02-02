CKLRT
=======

Composite Kernel Machine Regression based on Likelihood Ratio Test (CKLRT): in this package, we develop a kernel machine regression framework to model the overall genetic effect of a SNP-set, considering the possible GE interaction. Specifically, we use a composite kernel to specify the overall genetic effect via a nonparametric function and we model additional covariates parametrically within the regression framework. The composite kernel is constructed as a weighted average of two kernels, one corresponding to the genetic main effect and one corresponding to the GE interaction effect. We propose a likelihood ratio test (LRT) and a restricted likelihood ratio test (RLRT) for statistical significance. We derive a Monte Carlo approach for the finite sample distributions of LRT and RLRT statistics.

Usage
=======

[The 'CKLRT' vignette](https://github.com/andrewhaoyu/CKLRT/blob/master/inst/CKLRT_package.pdf) will provide a good start point for using CKLRT package.


Development 
=======
This R package is developed by Ni Zhao and Haoyu Zhang, and maintained by Haoyu Zhang <andrew.haoyu@gmail.com>.

Installation
=======
To install the development version of CKLRT, it's easiest to use the 'devtools' package. Note that CKLRT depends on the 'Rcpp' and 'RcppEigen' package, which also requires appropriate setting of Rtools and Xcode for Windows and Mac OS/X, respectively.

install.packages("devtools")  
library(devtools)  
install_github("andrewhaoyu/CKLRT")

References
=======
N. Zhao, H. Zhang, J. Clark, A. Maity, M. Wu. Composite Kernel Machine Regression based on Likelihood Ratio Test with Application for Combined Genetic and Gene-environment Interaction Effect, Biometrics (2018).


