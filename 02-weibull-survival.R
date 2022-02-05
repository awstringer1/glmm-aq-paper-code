### Weibull Survival Regression ###


# Test for, and install if not found, the packages required to run this script.
pkgs <- c(
  'tidyverse',
  'TMB',
  'aghq',
  'parallel',
  'survival'
)
for (pkg in pkgs) {
  if (!require(pkg,character.only = TRUE,quietly = TRUE)) {
    cat(paste0("Could not find package ",pkg,", installing from CRAN.\n"))
    install.packages(pkg)
    require(pkg,character.only = TRUE,quietly = TRUE)
  }
}

options(mc.cores = parallel::detectCores())

globalpath <- "~/work/projects/num-quadpoints"
tmbpath <- file.path(globalpath,"code")
resultspath <- file.path(globalpath,"results")

download.file("https://raw.githubusercontent.com/awstringer1/glmm-aq-paper-code/main/02_weibull_survival.cpp?token=GHSAT0AAAAAABRIFR6I577P4YSAV5LZWMIMYQGW3AA",
              file.path(tmbpath,"02_weibull_survival.cpp"))
download.file("https://raw.githubusercontent.com/awstringer1/glmm-aq-paper-code/main/02_weibull_survival_laplace.cpp?token=GHSAT0AAAAAABRIFR6J7I2HK6KQG7UW5FDQYQGW37Q",
              file.path(tmbpath,"02_weibull_survival_laplace.cpp"))

precompile()
compile(file.path(tmbpath,"02_weibull_survival.cpp"))
dyn.load(dynlib(file.path(tmbpath,"02_weibull_survival")))
compile(file.path(tmbpath,"02_weibull_survival_laplace.cpp"))
dyn.load(dynlib(file.path(tmbpath,"02_weibull_survival_laplace")))

# Data. Object kidney becomes available in the global environment when you load
# the survival package, but data(kidney,package = "survival") fails.

# Setup TMB template for a subject
get_tmb_data <- function(idx,params) {
  tmbdata <- with(kidney[kidney$id == idx, ],list(
    t = time,
    delta = status,
    x = cbind(age,sex-1,t(sapply(disease,function(x) x == c("GN","AN","PKD"))))
    # x = cbind(age,sex-1)
  ))
  if (!is.null(params)) {
    p <- ncol(tmbdata$x)
    tmbdata$beta <- params[1:p]
    tmbdata$logmu <- params[p+1]
    tmbdata$logalpha <- params[p+2]
    tmbdata$logeta <- params[p+3]
  }
  tmbdata
}
get_tmb_template <- function(idx,params) {
  tmbdata <- get_tmb_data(idx,params)

  tmbparams <- list(z=0)

  TMB::MakeADFun(
    data = tmbdata,
    parameters = tmbparams,
    silent = TRUE,
    DLL = "02_weibull_survival"
  )
}
get_tmb_template_laplace <- function(idx) {
  # Get the template with the auto Laplace turned on
  tmbdata <- get_tmb_data(idx,params=NULL)

  tmbparams <- list(z=0,beta = rep(0,ncol(tmbdata$x)),logmu = 0,logalpha = 0,logeta = 0)

  TMB::MakeADFun(
    data = tmbdata,
    parameters = tmbparams,
    silent = TRUE,
    DLL = "02_weibull_survival_laplace",
    random = "z"
  )
}

# Log likelihood (single)
logliksingle_laplace <- function(idx,params) {
  template <- get_tmb_template_laplace(idx)
  -1*template$fn(params)
}
logliksingle_gr_laplace <- function(idx,params) {
  template <- get_tmb_template_laplace(idx)
  -1*template$gr(params)
}


logliksingle <- function(idx,params,k) {
  # params = c(beta,logmu,logalpha,logeta)
  # variable of integration is the log-frailty, z
  template <- get_tmb_template(idx,params)
  if(k == 1) return(logliksingle_laplace(idx,params))
  get_log_normconst(aghq(template,k,0,control = default_control(onlynormconst=TRUE,method = "sparse_trust")))
}
# Log likelihood
loglik <- function(params,k) sum(sapply(unique(kidney$id),logliksingle,params=params,k=k))
loglik_parallel <- function(params,k) Reduce(sum,mclapply(unique(kidney$id),logliksingle,params=params,k=k))

loglik_gr <- function(params,k) {
  if (k == 1) {
    return(apply(Reduce(rbind,lapply(unique(kidney$id),logliksingle_gr_laplace,params=params)),2,sum))
  }
  numDeriv::grad(loglik,params,k=k)
}
loglik_gr_parallel <- function(params,k) {
  if (k == 1) {
    return(apply(Reduce(rbind,mclapply(unique(kidney$id),logliksingle_gr_laplace,params=params)),2,sum))
  }
  numDeriv::grad(loglik_parallel,params,k=k)
}
# Function to fit the model and compute table of estimates for a given k

fit_model <- function(k) {
  cat("Fitting model with k = ",k,"\n",sep="")
  pp <- rep(0,8) # CHANGE THIS if model changes
  fn <- function(p) -1*loglik_parallel(p,k)
  gr <- function(p) -1*loglik_gr_parallel(p,k)
  he <- function(p) as(numDeriv::jacobian(gr,p),'dgCMatrix')
  tm <- Sys.time()
  opt <- optim(pp,fn,gr,method="L-BFGS-B",lower=-5,upper=2,hessian=TRUE)
  tm2 <- Sys.time()
  dt <- as.numeric(difftime(tm2,tm,units = 'secs'))

  optest <- with(opt,c(par[1:5],exp(par[6:7]),exp(-par[8])))
  optse <- with(opt,sqrt(diag(solve(hessian))))
  optcilower <- with(opt,c(par[1:5]-2*optse[1:5],exp(par[6:7] - 2*optse[6:7]),exp(-1*(par[8] + 2*optse[8]))))
  optciupper <- with(opt,c(par[1:5]+2*optse[1:5],exp(par[6:7] +2*optse[6:7]),exp(-1*(par[8] - 2*optse[8]))))
  testresults <- as.data.frame(cbind(optcilower,optest,optciupper))
  testresults$k <- k
  testresults$variable <- c(paste0('beta',1:5),'mu','alpha','sigmasq')
  testresults$time <- dt
  testresults
}

ktodo <- c(1,3,5,6,7,9,11)
resultslist <- lapply(ktodo,fit_model)
results <- Reduce(rbind,resultslist)
save(results,resultslist,file = file.path(resultspath,"weibulldatasims-20220106-v2.RData"))

# Print table for paper
tabpaper <- results %>%
  filter(variable %in% c('beta2','mu','alpha','sigmasq')) %>%
  pivot_wider(names_from = 'variable',values_from = optcilower:optciupper)
knitr::kable(tabpaper[ ,c('k',
             'optest_beta2','optcilower_beta2','optciupper_beta2',
             'optest_mu','optcilower_mu','optciupper_mu',
             'optest_alpha','optcilower_alpha','optciupper_alpha',
             'optest_sigmasq','optcilower_sigmasq','optciupper_sigmasq')],
             digits = 3,
             format = 'latex')
cat("Done. Check",globalpath,"for the results.\n")




