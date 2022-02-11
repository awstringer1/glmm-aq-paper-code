## Simulated binary data ----

# Test for, and install if not found, the packages required to run this script.
pkgs <- c(
  'tidyverse',
  'parallel',
  'lme4'
)
for (pkg in pkgs) {
  if (!require(pkg,character.only = TRUE,quietly = TRUE)) {
    cat(paste0("Could not find package ",pkg,", installing from CRAN.\n"))
    install.packages(pkg)
    require(pkg,character.only = TRUE,quietly = TRUE)
  }
}

options(mc.cores = parallel::detectCores())

globalpath <- tempdir()
resultspath <- file.path(globalpath,'results')
if (!dir.exists(resultspath)) dir.create(resultspath)


PLOTTEXTSIZE <- 28

choosek <- function(M,m) ceiling(3*log(M,base=m)/2 - 2)

## Fit for multiple k ----
get_summaries <- function(modlist) {
  # Get summaries from fitted models

  namestoreturn <- c('M','m','sigma','beta0','betax','betaT','betaTx','itr','k',
                     'intercept_est','intercept_sd',
                     'betax_est','betax_sd',
                     'betaT_est','betaT_sd',
                     'betaxT_est','betaxT_sd',
                     'sigma_est',
                     'time','numiter',
                     'cilowerint','ciupperint')

  out <- matrix(0,nrow=length(modlist$k),ncol=length(namestoreturn))
  colnames(out) <- namestoreturn
  out <- as.data.frame(out)
  
  # Note: the error checking here was mostly for development; the presented simulations
  # don't really throw errors on fitting.
  for (j in 1:length(modlist$k)) {
    mod <- modlist$mods[[j]]
    if (inherits(mod,'error')) {
      out <- rep(-999,length(namestoreturn))
      names(out) <- namestoreturn
      return(out)
    }
    
    # Return results (attempt)
    beta0 <- tryCatch(unname(summary(mod)$coefficients[1,1:2]),error = function(e) e)
    if (inherits(beta0,'error')) {
      out <- rep(-999,length(namestoreturn))
      names(out) <- namestoreturn
      return(out)
    }
    betax <- tryCatch(unname(summary(mod)$coefficients[2,1:2]),error = function(e) e)
    if (inherits(betax,'error')) {
      out <- rep(-999,length(namestoreturn))
      names(out) <- namestoreturn
      return(out)
    }
    betaT <- tryCatch(unname(summary(mod)$coefficients[3,1:2]),error = function(e) e)
    if (inherits(betaT,'error')) {
      out <- rep(-999,length(namestoreturn))
      names(out) <- namestoreturn
      return(out)
    }
    betaxT <- tryCatch(unname(summary(mod)$coefficients[4,1:2]),error = function(e) e)
    if (inherits(betaxT,'error')) {
      out <- rep(-999,length(namestoreturn))
      names(out) <- namestoreturn
      return(out)
    }
    
    opt <- mod@optinfo
    sigmaest <- opt$val[1]
    numiter <- opt$feval
    
    # Confidence intervals
    cc <- tryCatch(confint(mod,2,method='Wald'),error = function(e) e)
    if (inherits(cc,'error')) {
      out <- rep(-999,length(namestoreturn))
      names(out) <- namestoreturn
      return(out)
    }
    
    out[j, ] <- with(modlist,c(
      info,
      k[j],
      beta0,
      betax,
      betaT,
      betaxT,
      sigmaest,
      times[j],
      numiter,
      as.numeric(cc)
    ))
  }
  out
}

simulate_data_fit_model <- function(M,m,k,sigma,itr) {
  # k: vector of k values to do
  # itr: a label for the output
  # Simulate data and fit the model with supplied values of k
  
  # Fixed simulation values
  beta0_true <- -3
  betax_true <- -1
  betaT_true <- 1
  betaTx_true <- -.5
  
  pT <- .5
  
  stopifnot(length(m) == 1)
  n <- M*m
  
  # Simulate a treatment vector
  x <- rbinom(M,1,pT)
  # Random effect
  U <- rnorm(M,sd=sigma)
  
  # Simulate the outcome
  X <- expand.grid(ID = 1:M,tt = 1:m-1)
  X$x <- rep(x,m)
  X$U <- rep(U,m)
  X$lp <- with(X,beta0_true + betax_true*x + betaT_true*tt + betaTx_true*x*tt + U)
  X$p <- with(X,1 / (1 + exp(-lp)))
  X$y <- with(X,rbinom(n,1,p))
  
  mods <- vector(mode='list',length=length(k))
  names(mods) <- k
  times <- numeric(length(k))
  names(times) <- k
  
  for (j in 1:length(k)) {
    # Fit the model (attempt)
    tm <- Sys.time()
    suppressWarnings(mod <- tryCatch(glmer(y ~ x*tt + (1|ID),data = X,family = binomial,nAGQ = k[j]),error = function(e) e))
    times[j] <- unname(as.numeric(difftime(Sys.time(),tm,units='secs')))
    mods[[j]] <- mod
  }
  
  out <- list()
  out$mods <- mods
  out$times <- times
  out$k <- k
  out$info <- c('M'=M,'m'=m,'sigma'=sigma,'beta0'=beta0_true,'betax'=betax_true,'betaT'=betaT_true,'betaTx'=betaTx_true,'itr'=itr)
  
  # Don't save the models, too much memory.
  out$summaries <- get_summaries(out)
  out$mods <- NULL
  out
}



Mtodo <- c(100,1000)
mtodo <- c(3,5)
ktodo <- c(1,3,5,7,9,11)
numeach <- 500
sigmatodo <- c(1,3)

# simstodo <- expand.grid(M = Mtodo,m = mtodo,sigma = sigmatodo,itr = 1:numeach)
simstodo <- expand.grid(M = Mtodo,m = mtodo,sigma = sigmatodo,itr = 1:numeach)

dosim_model <- function(lst) with(lst,simulate_data_fit_model(M,m,ktodo,sigma,itr))
simlist <- vector(mode = 'list',length = nrow(simstodo))
for (i in 1:nrow(simstodo)) simlist[[i]] <- simstodo[i, ]
tm <- Sys.time()
cat("Doing sims.\n")
cat("Fitting models...\n")
sim_models <- mclapply(simlist,dosim_model)
cat("Saving results...\n")
save(sim_models,file = file.path(resultspath,"binary-rare-sims.RData"))
cat("Done fitting models.\n")
cat("Sims took",round(as.numeric(difftime(Sys.time(),tm,units='secs')),0),"seconds.\n")
cat("Computing summaries...\n")

simframe <- Reduce(rbind,
                   Filter(function(x) !inherits(x,'try-error'),
                          lapply(sim_models,'[[','summaries'))) %>% 
  as.data.frame() %>% 
  as_tibble() %>% 
  filter(stringr::str_detect(intercept_est,"Error",negate = TRUE)) %>%
  filter(M != -999) %>%
  mutate_all(as.numeric)


# Line plot of means
simframe_witherror <- simframe %>%
  group_by(M,m,sigma,k) %>%
  mutate(
    int_error = (intercept_est - beta0)^2,
    betax_error = (betax_est - betax)^2,
    betaT_error = (betaT_est - betaT)^2,
    betaxT_error = (betaxT_est - betaTx)^2,
    sigma_error = (sigma_est - sigma)^2,
    int_covr = as.numeric(beta0 < ciupperint) * as.numeric(beta0 > cilowerint)
  ) %>%
  summarize(
    int_mse = quantile(int_error,.5),
    int_lower = quantile(int_error,.025),
    int_upper = quantile(int_error,.975),
    betax_mse = quantile(betax_error,.5),
    betax_lower = quantile(betax_error,.025),
    betax_upper = quantile(betax_error,.975),      
    betaT_mse = quantile(betaT_error,.5),
    betaT_lower = quantile(betaT_error,.025),
    betaT_upper = quantile(betaT_error,.975), 
    betaxT_mse = quantile(betaxT_error,.5),
    betaxT_lower = quantile(betaxT_error,.025),
    betaxT_upper = quantile(betaxT_error,.975),
    sigma_mse = quantile(sigma_error,.5),
    sigma_lower = quantile(sigma_error,.025),
    sigma_upper = quantile(sigma_error,.975),
    int_covr = mean(int_covr),
    meantime = mean(time),
    sdtime = sd(time),
    meaniter = mean(numiter),
    sditer = sd(numiter)
  ) %>%
  mutate(
    k_rec = choosek(M,m),
    int_rmse = sqrt(int_mse),
    int_rmse_lower = sqrt(int_lower),
    int_rmse_upper = sqrt(int_upper),
    
    betax_rmse = sqrt(betax_mse),
    betax_rmse_lower = sqrt(betax_lower),
    betax_rmse_upper = sqrt(betax_upper),
    
    betaT_rmse = sqrt(betaT_mse),
    betaT_rmse_lower = sqrt(betaT_lower),
    betaT_rmse_upper = sqrt(betaT_upper),
    
    betaxT_rmse = sqrt(betaxT_mse),
    betaxT_rmse_lower = sqrt(betaxT_lower),
    betaxT_rmse_upper = sqrt(betaxT_upper),
    
    sigma_rmse = sqrt(sigma_mse),
    sigma_rmse_lower = sqrt(sigma_lower),
    sigma_rmse_upper = sqrt(sigma_upper)
  )


## Plots for Journal ----

line_plot <- function(var,M,m,ylim) {
  vl <- paste0(var,"_lower")
  vu <- paste0(var,"_upper")
  
  # Return a plot with 3 lines for each sigma
  filter(simframe_witherror,M==!!M,m==!!m) %>%
    ggplot(aes(x = k,y = .data[[var]],linetype = as.factor(sigma))) +
    theme_classic() +
    geom_line() +
    geom_point(pch=21,size=2) +
    geom_line(aes(y = .data[[vl]]),colour = "gray") +
    geom_line(aes(y = .data[[vu]]),colour = "gray") +
    geom_vline(aes(xintercept = k_rec),linetype = 'longdash') +
    scale_x_continuous(breaks = ktodo) +
    scale_linetype_manual(values = c("1" = "dashed","2" = "dotted","3" = "dotdash"),breaks = c("3","5","10")) + 
    coord_cartesian(ylim = ylim) +
    labs(x = "k",y = "Absolute error") +
    guides(linetype = 'none') +
    theme(text = element_text(size = PLOTTEXTSIZE))
  
}

covr_plot <- function(M,m) {
  filter(simframe_witherror,M==!!M,m==!!m) %>%
    ggplot(aes(x = k,y = int_covr,linetype=as.factor(sigma))) +
    theme_classic() +
    geom_line() +
    geom_point() +
    scale_linetype_manual(values = c("1" = "dashed","2" = "dotted","3" = "dotdash"),breaks = c("3","5","10")) + 
    geom_vline(aes(xintercept = k_rec),linetype = 'dotted') +
    geom_hline(yintercept = .95,linetype = 'dashed') +
    scale_x_continuous(breaks = ktodo) +
    coord_cartesian(ylim = c(0,1)) + 
    labs(x = "k",y = "Empirical coverage") +
    guides(linetype = 'none') +
    theme(text = element_text(size = PLOTTEXTSIZE))
}

ggsave(filename = file.path(resultspath,"intercept-rmse-M100-m3.pdf"),plot = line_plot("int_rmse",100,3,c(0,10)),width = 7,height = 7)
ggsave(filename = file.path(resultspath,"intercept-rmse-M100-m5.pdf"),plot = line_plot("int_rmse",100,5,c(0,10)),width = 7,height = 7)
ggsave(filename = file.path(resultspath,"intercept-rmse-M1000-m3.pdf"),plot = line_plot("int_rmse",1000,3,c(0,10)),width = 7,height = 7)
ggsave(filename = file.path(resultspath,"intercept-rmse-M1000-m5.pdf"),plot = line_plot("int_rmse",1000,5,c(0,10)),width = 7,height = 7)

ggsave(filename = file.path(resultspath,"sigma-rmse-M100-m3.pdf"),plot = line_plot("sigma_rmse",100,3,c(0,10)),width = 7,height = 7)
ggsave(filename = file.path(resultspath,"sigma-rmse-M100-m5.pdf"),plot = line_plot("sigma_rmse",100,5,c(0,10)),width = 7,height = 7)
ggsave(filename = file.path(resultspath,"sigma-rmse-M1000-m3.pdf"),plot = line_plot("sigma_rmse",1000,3,c(0,10)),width = 7,height = 7)
ggsave(filename = file.path(resultspath,"sigma-rmse-M1000-m5.pdf"),plot = line_plot("sigma_rmse",1000,5,c(0,10)),width = 7,height = 7)

ggsave(filename = file.path(resultspath,"int-covr-M100-m3.pdf"),plot = covr_plot(100,3),width = 7,height = 7)
ggsave(filename = file.path(resultspath,"int-covr-M100-m5.pdf"),plot = covr_plot(100,5),width = 7,height = 7)
ggsave(filename = file.path(resultspath,"int-covr-M1000-m3.pdf"),plot = covr_plot(1000,3),width = 7,height = 7)
ggsave(filename = file.path(resultspath,"int-covr-M1000-m5.pdf"),plot = covr_plot(1000,5),width = 7,height = 7)

cat("Done plots. Doing time summary...\n")


timetable <- simframe_witherror %>%
  dplyr::select(M,m,sigma,k,contains('time'),contains('iter')) %>%
  dplyr::filter(k%in%c(1,dplyr::case_when(choosek(M,m)%%2==0~choosek(M,m)+1,TRUE ~ choosek(M,m)))) # Only fit with odd k, so pick the next highest (slowest) one if k(M,m) even

knitr::kable(timetable,digits=3) # Reformat in Latex, easier.


cat("Done. Check",resultspath,"for the results.\n")

