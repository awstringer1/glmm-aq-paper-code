### Analyze the Toenail Data ###
## Section 5.1 ##

# Test for, and install if not found, the packages required to run this script.
pkgs <- c(
  'tidyverse',
  'mice',
  'lme4'
)
for (pkg in pkgs) {
  if (!require(pkg,character.only = TRUE,quietly = TRUE)) {
    cat(paste0("Could not find package ",pkg,", installing from CRAN.\n"))
    install.packages(pkg)
    require(pkg,character.only = TRUE,quietly = TRUE)
  }
}

data(toenail,package = 'mice')

globalpath <- tempdir()
resultspath <- file.path(globalpath,"results")
if (!dir.exists(resultspath)) dir.create(resultspath)


PLOTTEXTSIZE <- 28

# How many points?
toenail_subjects <- toenail %>%
  group_by(ID) %>%
  summarize(m = n())
nrow(toenail_subjects) # M = 294
toenail_subjects %>%
  group_by(m) %>%
  summarize(numtimes = n())
kMm <- function(M,m) ceiling(3*log(M,base=m)/2 - 2)
suggested_points <- kMm(nrow(toenail_subjects),sort(unique(toenail_subjects$m)))
points_to_use <- max(suggested_points[suggested_points<Inf])
cat("Recommended k: ",points_to_use,".\n",sep='')
# Fit the model with several choices and get timings
toenail_model <- function(k) {
  lme4::glmer(
    outcome ~ treatment*month + (1|ID),
    data = toenail,
    family = binomial,
    nAGQ = k
  )
}

ktodo <- c(seq(1,11,by=2),15,25)
numtimes <- 100 # Number of times to fit each model, for timing purposes
times <- data.frame(k = rep(ktodo,each=numtimes),time = 0)
models <- vector(mode='list',length=length(ktodo))
names(models) <- rep(ktodo)
idx <- 0

for (i in 1:length(ktodo)) {
  k <- ktodo[i]
  for (j in 1:numtimes) {
    idx <- idx+1
    tm <- Sys.time()
    tmp<- toenail_model(k)
    times[idx,'time'] <- as.numeric(difftime(Sys.time(),tm,units='secs'))
  }
  models[[i]] <- tmp
}

timingsummary <- times %>%
  group_by(k) %>%
  summarize(meantime = mean(time),sdtime = sd(time),numtimes = n())

readr::write_csv(times,file = file.path(resultspath,"raw-timings.csv"))
readr::write_csv(timingsummary,file = file.path(resultspath,"summary-timings.csv"))


# Coefficients
get_coefs <- function(k) {
  # get a row in the summary table
  mod <- models[names(models)==k]
  summ <- summary(mod[[1]])
  coefs <- unname(summ$coefficients[ ,1])
  stderrs <- unname(summ$coefficients[ ,2])
  vars <- rownames(summ$coefficients)
  sigsq <- summ$varcor$ID[1,1]
  out <- c(coefs,stderrs,sigsq)
  names(out) <- c(paste0(vars,'-est'),paste0(vars,'-sd'),'sigmasq')
  out[c(1,5,2,6,3,7,4,8,9)]
}

estimatetable <- as.data.frame(Reduce(rbind,Map(get_coefs,ktodo)))
estimatetable$k <- ktodo
rownames(estimatetable) <- NULL
readr::write_csv(estimatetable,file = file.path(resultspath,"summary-estimates.csv"))


# Print the table for the paper
# Load from disk
estimatetable <- readr::read_csv(file.path(resultspath,"summary-estimates.csv"))
timingsummary <- readr::read_csv(file.path(resultspath,"summary-timings.csv"))

merge(estimatetable[ ,c('k','(Intercept)-est','(Intercept)-sd','month-est','month-sd','sigmasq')],timingsummary,by='k') %>%
  knitr::kable(format = 'latex',digits = 3)


# Done
cat("Finished. Check ",globalpath," for results.\n")





