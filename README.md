# Fitting Generalized Linear Mixed Models Using Adaptive Quadrature
Code for paper *Fitting Generalized Linear Mixed Models Using Adaptive Quadrature*. Each file corresponds to different results in the paper:

1. **00-toenail.R**: Section 5.1, "Toenail" data,
2. **02-weibull-survival.R**: Section 5.2, the Weibull survival model fit to the "Kidney" data. Requires compilation of the following TMB templates, which will be downloaded from here and compiled within the script:
  - **02_weibull_survival.cpp**
  - **02_weibull_survival_laplace.cpp**
3. **01-simulations.R**: Appendix C, simulations.

All data and packages are downloaded, installed, and loaded within each script. On a "standard" machine (maybe not Windows?) you should be able to simply run the code in a fresh `R` session and the results from the paper should pop out. Obviously the absolute timings will not be the same, although the relative timings should be roughly the same, to the point that the narrative in the paper isn't affected.

Version information for where the code was tested:

```
> version
               _                           
platform       aarch64-apple-darwin20      
arch           aarch64                     
os             darwin20                    
system         aarch64, darwin20           
status                                     
major          4                           
minor          1.2                         
year           2021                        
month          11                          
day            01                          
svn rev        81115                       
language       R                           
version.string R version 4.1.2 (2021-11-01)
nickname       Bird Hippie   
```
