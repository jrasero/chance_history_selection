##################################################
#
#
# Script to compute the three forces, along with
# their bootstrapped confidence intervals and
# the bayes factors as alternative to p-values
# The script can run under two different R packages
# for linear mixed modelling: # lme4 and blme.
# We finally used the latter one because it was more
# suitable given our sample sizes
#
#
##################################################

##################################################
########## Load the needed packages ##############

library(plyr)
library(lme4)
library(blme)
library(dplyr)

##################################################
########## Some definitions here #################

# Insert here the path to the directory where the 
# data and results layout files are located
data.dir<-"path/to/data/"

# Insert here the path to the output directory
output.dir<-"path/to/output/"

# Decide whether to run under blme or lme4 
# (We used blme in the paper)
which.model<-"Bayes"

# Define estimator to be used throughout the script 
# and the name of the output file
if(which.model=="Bayes"){
  estimator = blmer
  cat("doing bayes")
  output.file<-'results_boot_blmer.csv'
}else
  {
  estimator = lmer
  cat("doing normal")
  output.file<-'results_boot_lmer.csv'
}

##################################################
############## Functions #########################

# Here we define the different functions that we
# employed for our analysis. 

observed <- function(MIC.subset, optimizer="Nelder_Mead"){
  "
  function to compute the observed forces
  "
  
  # Extract ancestor and evolved values 
  ancestor <- MIC.subset[MIC.subset$Day==0,]
  evolved <- MIC.subset[MIC.subset$Day!=0,]
  
  # This is the linear model for chance and history,
  # which involves only evolved observations
  lme.fit <-estimator(MIC~1+(1|Ancestor/Population), 
                      data=evolved, 
                      control = lmerControl(optimizer=optimizer))
  
  chance <- sqrt(as.numeric(VarCorr(lme.fit)))[1]
  history <- sqrt(as.numeric(VarCorr(lme.fit)))[2]
  
  # This is the linear model for selection. In the end,
  # selection is just the fixed effects coefficient, ie,
  # the difference in means between different days.
  sel.model<-estimator(MIC~Day + (1|Ancestor), 
                       data=MIC.subset)  
  selection <- coef(summary(sel.model))[2,1] 
  
  # return vector with history, chance and selection
  return(c(history, chance, selection))
}

bootstrap_selection<-function(MIC.subset, 
                              nboot=100, 
                              seed = NULL){
  "
  Function to compute the bootstrapp selection force.
  This is done separated from history and chance because 
  here we we use observations from two days.
  "
  
  nrows<-dim(MIC.subset)[1]
  
  set.seed(seed)
  
  sel.boot<-c()
  chance.boot<-c()
  history.boot<-c()
  
  ii<-0
  while(ii<nboot){
    # TryCatch to avoid few possible errors or warnings during bootstrapping
    tryCatch({
      boot.idx<-sample(c(1:nrows), nrows, replace = T)  
      MIC.boot<-MIC.subset[boot.idx,]
      
      ancestor.boot <- MIC.boot[MIC.boot$Day==0,]
      evolved.boot <- MIC.boot[MIC.boot$Day!=0,]
      
      sel.model<-estimator(MIC~Day + (1|Ancestor), 
                           data=MIC.boot)  
      sel.boot<-c(sel.boot, coef(summary(sel.model))[2,1])

      ii = ii + 1
      
      if(ii %% 100 ==0){
        cat(paste("iteration ", " finished \n", sep= as.character(ii)))
      }
      
    }, warning = function(w) { cat("WARNING :", conditionMessage(w), "\n")},
    error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
    
  }
  return(data.frame("sel.boot" = sel.boot))
}

bootstrap_evolved <-function(MIC.subset, 
                             optimizer="Nelder_Mead",
                             nboot=100, 
                             seed = NULL){
  "
  Function to compute the bootstrapp history and chance force.
  In this case, we only use data from observations in the evolved day.
  "
  
  set.seed(seed)
  
  ancestor <- MIC.subset[MIC.subset$Day==0,]
  evolved <- MIC.subset[MIC.subset$Day!=0,]
  
  nrows<-dim(evolved)[1]
  
  chance.boot<-c()
  history.boot<-c()
  
  ii<-0
  while(ii<nboot){
    # TryCatch to avoid errors or warning during bootstrapping
    tryCatch({
      boot.idx<-sample(c(1:nrows), nrows, replace = T)  
      evolved.boot<-evolved[boot.idx,]
      
      lme.fit <-estimator(MIC~1+ (1 | Ancestor/Population), 
                          data=evolved.boot, 
                          control = lmerControl(optimizer = optimizer))
      
      if(isSingular(lme.fit, tol = 1e-4) == TRUE){
        next
      }
      
      chance.boot<-c(chance.boot, sqrt(as.numeric(VarCorr(lme.fit)))[1])
      history.boot<-c(history.boot, sqrt(as.numeric(VarCorr(lme.fit)))[2])
      ii = ii + 1
      
      if(ii %% 100 ==0){
        cat(paste("iteration ", " finished \n", sep= as.character(ii)))
      }
      
    }, warning = function(w) { cat("WARNING :", conditionMessage(w), "\n")},
    error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
    
  }
  return(data.frame("chance.boot" = chance.boot,
                    "history.boot" = history.boot))
}

compute_bf<-function(BIC.0, BIC.1){
  "
  Function to compute the approximate bayes factor BF01 
  based on BICs (equation 10, Wagenmakers 2007)
  between a null (BIC.0) and alternative model (BIC.1)
  "
  return (exp(0.5*(BIC.1-BIC.0)))
}


##################################################
################# Inputs ########################

# load the observed data
MICs <- read.csv(paste0(data.dir, 
                        "2020-03-02 HCS data for R.csv")
                 )

# Correct input data format
MICs$Day <- as.factor(MICs$Day)
MICs$Population <- as.factor(MICs$Population)
MICs$Measurement <- as.factor(MICs$Measurement)

# Add a litte of noise.
# There are some issues with observations having the same
# values across ancestors and populations. One way around to
# this is to add a bit of noise to the observed MIC
# values, which would add a bit of (expected) experimental noise
set.seed(0) # For reproducible noise
MICs$MIC<-MICs$MIC + rnorm(length(MICs$MIC), 0, 1e-3)

# Make a data frame for results
results <- read.csv(paste0(data.dir, 
                           "results_layout.csv"))[, c(1:12)]
results$Drug_evolved <- as.character(results$Drug_evolved)
results$Drug_tested <- as.character(results$Drug_tested)

# Add null columns for bayes factors
results$bf01_selection<-NA
results$bf01_history<-NA
results$bf01_chance<-NA

results$bf10_selection<-NA
results$bf10_history<-NA
results$bf10_chance<-NA

# Some constants
RANDOM.STATE = 0 # seed for reproducibility
alpha = 0.05 # significance level


##################################################
############ Computations ########################

# 1 - Compute observed values
start <- proc.time()
for(i in c(1:nrow(results))){
    
    tryCatch({
        drug.evolved<-results[i, "Drug_evolved"]
        drug.test<-results[i, "Drug_tested"]
        day<-results[i, "Day"]

        cat(paste0("Drug evolved = ", drug.evolved,
                       " drug test = ", drug.test, 
                       " day = ", as.character(day)), "\n")
            
        MIC.subset <- MICs[MICs$Drug_evolved==drug.evolved & MICs$Drug_tested==drug.test
                          & (MICs$Day==0 | MICs$Day==day),]

        MIC.subset<-droplevels(MIC.subset)
        
        #Compute obseved values
        observed.values <- observed(MIC.subset = MIC.subset,
                                    optimizer="Nelder_Mead")

        results[i, "History"] = observed.values[1]
        results[i, "Chance"] = observed.values[2]
        results[i, "Selection"] = observed.values[3]

    }, 
                            warning = function(w) { cat("WARNING :", conditionMessage(w), "\n")},
                            error=function(e){cat("ERROR :", conditionMessage(e), "\n")}
                           )
             
    
}
end <- proc.time()
cat("elapsed time for observed values =", (end[3] - start[3]), "secs")

# 2 - Compute boostrapped selection
start <- proc.time()
for(i in c(1:nrow(results))){
    
    tryCatch({
        drug.evolved<-results[i, "Drug_evolved"]
        drug.test<-results[i, "Drug_tested"]
        day<-results[i, "Day"]

        cat(paste0("Drug evolved = ", drug.evolved,
                       " drug test = ", drug.test, 
                       " day = ", as.character(day)), "\n")
            
        MIC.subset <- MICs[MICs$Drug_evolved==drug.evolved & MICs$Drug_tested==drug.test
                          & (MICs$Day==0 | MICs$Day==day),]

        MIC.subset<-droplevels(MIC.subset)
        
        boot.sel.df<-bootstrap_selection(MIC.subset = MIC.subset, 
                                         nboot = 1000,
                                         seed = RANDOM.STATE)

        quantiles.df<-apply(boot.sel.df, 
                            MARGIN = 2, 
                            function(x) quantile(x, probs = c(alpha/2, 1-alpha/2)))
                            
        results[i, "Selection_95low"] = quantiles.df["2.5%","sel.boot"]
        results[i, "Selection_95high"] = quantiles.df["97.5%","sel.boot"]

        
    }, 
                            warning = function(w) { cat("WARNING :", conditionMessage(w), "\n")},
                            error=function(e){cat("ERROR :", conditionMessage(e), "\n")}
                           )
             
    
}
end <- proc.time()
cat("elapsed time for bootstrap on selection =", (end[3] - start[3]), "secs")

# 3 - Compute boostrapped history and chance
start <- proc.time()
for(i in c(1:8)){
    
    tryCatch({
        drug.evolved<-results[i, "Drug_evolved"]
        drug.test<-results[i, "Drug_tested"]
        day<-results[i, "Day"]

        cat(paste0("Drug evolved = ", drug.evolved,
                       " drug test = ", drug.test, 
                       " day = ", as.character(day)), "\n")
            
        MIC.subset <- MICs[MICs$Drug_evolved==drug.evolved & MICs$Drug_tested==drug.test
                          & (MICs$Day==0 | MICs$Day==day),]

        MIC.subset<-droplevels(MIC.subset)
        
        boot.ch.his.df<-bootstrap_evolved(MIC.subset = MIC.subset,
                                          optimizer="Nelder_Mead",
                                          nboot = 1000, 
                                          seed = RANDOM.STATE)

        quantiles.df<-apply(boot.ch.his.df, 
                            MARGIN = 2, 
                            function(x) quantile(x, probs = c(alpha/2, 1-alpha/2)))
                            
        results[i, "History_95low"] = quantiles.df["2.5%","history.boot"]
        results[i, "History_95high"] = quantiles.df["97.5%","history.boot"]

        results[i, "Chance_95low"] = quantiles.df["2.5%","chance.boot"]
        results[i, "Chance_95high"] = quantiles.df["97.5%","chance.boot"]

        
    }, 
                            warning = function(w) { cat("WARNING :", conditionMessage(w), "\n")},
                            error=function(e){cat("ERROR :", conditionMessage(e), "\n")}
                           )
             
    
}
end <- proc.time()
cat("elapsed time for bootstrap on selection =", (end[3] - start[3]), "secs")

# 4 - Compute Bayes factors 
start <- proc.time()
for(i in c(1:8)){
  tryCatch({
    drug.evolved<-results[i, "Drug_evolved"]
    drug.test<-results[i, "Drug_tested"]
    day<-results[i, "Day"]
    
    cat(paste0("Drug evolved = ", drug.evolved,
               " drug test = ", drug.test, 
               " day = ", as.character(day)), "\n")
    
    MIC.subset <- MICs[MICs$Drug_evolved==drug.evolved & MICs$Drug_tested==drug.test
                       & (MICs$Day==0 | MICs$Day==day),]
    
    MIC.subset<-droplevels(MIC.subset)
    
    # Compute BICs for Selection hypothesis
    BIC.1<-BIC(estimator("MIC~Day + (1 | Ancestor)", data = MIC.subset))
    BIC.0<-BIC(estimator("MIC~1 + (1 | Ancestor)", data = MIC.subset))
    bf.01.sel<-compute_bf(BIC.0 = BIC.0, BIC.1 = BIC.1)
    
    evolved<-MIC.subset[MIC.subset$Day!=0,]
    
    # Compute BICs for history hypothesis
    BIC.1<-BIC(estimator("MIC~1 + (1|Ancestor/Population)", data = evolved, 
                     control = lmerControl(optimizer ="Nelder_Mead")))
    BIC.0<-BIC(estimator("MIC~1 + (1|Ancestor:Population)", data = evolved, 
                     control = lmerControl(optimizer ="Nelder_Mead")))
    bf.01.hist<-compute_bf(BIC.0 = BIC.0, BIC.1 = BIC.1)
    
    # Compute BICs for selection hypothesis
    BIC.1<-BIC(estimator("MIC~1 + (1|Ancestor/Population)", data = evolved, 
                     control = lmerControl(optimizer ="Nelder_Mead")))
    BIC.0<-BIC(estimator("MIC~1 + (1|Ancestor)", data = evolved, 
                     control = lmerControl(optimizer ="Nelder_Mead")))
    bf.01.chan<-compute_bf(BIC.0 = BIC.0, BIC.1 = BIC.1)
    
    results[i, "bf01_selection"]<-bf.01.sel
    results[i, "bf01_history"]<-bf.01.hist
    results[i, "bf01_chance"]<-bf.01.chan
    
    results[i, "bf10_selection"]<-1/bf.01.sel
    results[i, "bf10_history"]<-1/bf.01.hist
    results[i, "bf10_chance"]<-1/bf.01.chan
  }, 
  warning = function(w) { cat("WARNING :", conditionMessage(w), "\n")},
  error=function(e){cat("ERROR :", conditionMessage(e), "\n")}
  )
  
  
}
end <- proc.time()
cat("elapsed time for bayes factors =", (end[3] - start[3]), "secs")

##################################################
################ Outputs #########################

# A vector to rearrange the order of the columns
# in the output file
order.columns<- c("Drug_evolved","Drug_tested",
                  "Day", 
                  "Selection", 
                  "Selection_95low", "Selection_95high", 
                  "bf01_selection", "bf10_selection",
                  "History", 
                  "History_95low", "History_95high", 
                  "bf01_history", "bf10_history",
                  "Chance", 
                  "Chance_95low", "Chance_95high" , 
                  "bf01_chance", "bf10_chance")

# Save results to the disc
write.csv(results[, order.columns], 
          file=paste0(output.dir, output.file))

##################################################
################ Plots ###########################


# Selection
png(paste0(output.dir, paste("selection_", ".png", sep = which.model)))
plot(c(1:8), results$Selection, xlab="cases", ylab="selection",
     ylim = c(min(results$Selection_95low), 
              max(results$Selection_95high)),
     pch=16, cex=1)
# Add error bars
arrows(x0=c(1:8), y0=results$Selection_95low, 
       x1=c(1:8), y1=results$Selection_95high, code=3, 
       angle=90, length=0.1)
title("Selection")
dev.off()

# History
png(paste0(output.dir, paste("history_", ".png", sep = which.model)))
plot(c(1:8), results$History, 
     ylim = c(min(results$History_95high), 
              max(results$History_95high)),
     xlab="cases", ylab="history", pch=16, cex=1)
# Add error bars
arrows(x0=c(1:8), y0=results$History_95low, 
       x1=c(1:8), y1=results$History_95high, code=3, angle=90, length=0.1)
title("History")
dev.off()

# Chance
png(paste0(output.dir, paste("chance_", ".png", sep = which.model)))
plot(c(1:8), results$Chance, 
     ylim = c(min(results$Chance_95low), 
              max(results$Chance_95high)),
     xlab="cases", ylab="chance", pch=16, cex=1)
# Add error bars
arrows(x0=c(1:8), y0=results$Chance_95low, 
       x1=c(1:8), y1=results$Chance_95high, code=3, angle=90, length=0.1)
title("Chance")
dev.off()