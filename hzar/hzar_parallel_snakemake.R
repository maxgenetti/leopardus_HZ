#!/usr/bin/env Rscript
# Perform HZAR analises for individual genomic site
# Trindade F.T., 2024

### Check for the arguments and packages
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=7) {
  stop("Rscript hzar_parallel_snakemake.R <hzar_input> <modelsN,modelsX,modelsF> <chainLength> <mincenter> <maxcenter> <maxwidth> <threads>", call.=FALSE)
}

hzar_input <- args[1] # input.hzar
models2run <- unlist(strsplit(args[2], ",")) #modelsN,modelsX,modelsF
chainLength <- as.numeric(args[3]) #1e6
mincenter <- as.numeric(args[4]) #0
maxcenter <- as.numeric(args[5]) #4000
maxwidth <- as.numeric(args[6]) #4000
threads <- as.numeric(args[7])
hzar_output <- gsub("hzar_input","hzar_runs", hzar_input)

library(hzar)
library(coda)
library(foreach)

if(require(doMC)){
  registerDoMC(cores = threads)
} else {
  registerDoSEQ()
}


### Load molecular data from the data table.
leopardus <- read.table(hzar_input, header=TRUE, sep=",")

# Blank out space in memory to hold molecular analysis
if(length(apropos("^leopd$",ignore.case=FALSE)) == 0 || !is.list(leopd) ) {
  leopd <- list()
}

# Space to hold data and results
obs <- list() # Space to hold the observed data
models <- list() # Space to hold the models to fit
fitRs <- list() # Space to hold the compiled fit requests
runs <- list() # Space to hold the output data chains
analysis <- list() # Space to hold the analysed data

### Creating data object for the SNP
# The likelihood function used is chosen based on the method called:
# hzar.doMolecularData1DPops, hzar.doCLTData1DPops, hzar.doNormalData1DPops
# Depending on the data type, the method should be choose.
print("##### Reading observed data into hzar.obsData object #####")
obs <- hzar.doMolecularData1DPops(leopardus$DIST,
                                  leopardus$FREQ,
                                  leopardus$SAMPLING,
                                  siteID=leopardus$ID)

### Constructing different clineMetaModel objects by using a function
# Depending on the method used for build observed data objetc, different method for construct models may be used:
# hzar.makeCline1DFreq, hzar.makeCline1DNormal
print("##### Building models into clineMetaModel objects #####")
leopd.loadlocmodel <- function(scaling,tails,
                               id=paste(scaling,tails,sep="."))
  models[[id]] <<- hzar.makeCline1DFreq(obs, scaling, tails)

if ("modelsN" %in% models2run){
  leopd.loadlocmodel("none", "none", "modelNN")
  leopd.loadlocmodel("none", "right", "modelNR")
  leopd.loadlocmodel("none", "left", "modelNL")
  leopd.loadlocmodel("none", "mirror", "modelNM")
  leopd.loadlocmodel("none", "both", "modelNB")}
if ("modelsX" %in% models2run){
  leopd.loadlocmodel("fixed", "none", "modelXN")
  leopd.loadlocmodel("fixed", "right", "modelXR")
  leopd.loadlocmodel("fixed", "left", "modelXL")
  leopd.loadlocmodel("fixed", "mirror", "modelXM")
  leopd.loadlocmodel("fixed", "both", "modelXB")}
if ("modelsF" %in% models2run){
  leopd.loadlocmodel("free", "none", "modelFN")
  leopd.loadlocmodel("free", "right", "modelFR")
  leopd.loadlocmodel("free", "left", "modelFL")
  leopd.loadlocmodel("free", "mirror", "modelFM")
  leopd.loadlocmodel("free", "both", "modelFB")}

model_name <- names(models)

# Modify all models to focus on the region where the observed data were collected.
dmin <- min(obs$frame$dist) - 50
dmax <- max(obs$frame$dist) + 50
models <- sapply(models,
                 hzar.model.addBoxReq,
                 dmin , dmax,
                 simplify=FALSE)

### Compile each of the models to prepare for fitting, create an hzar.fitRequest object [it can took a while]
# Depending on the method used forbuild models, use a different method here:
# hzar.first.fitRequest.old.ML, hzar.first.fitRequest.gC
print("##### Creating hzar.fitRequest object to compile models and prepare fit #####")
fitRs$init <- foreach(i=1:length(models)) %dopar% {
	hzar.first.fitRequest.old.ML(models[[i]], obs, verbose=FALSE)
}
names(fitRs$init) <- model_name

# Update the settings for the fitter
print("##### Setting chain parameters #####")

mainSeed=
  list(A=c(596,528,124,978,544,99),
       B=c(528,124,978,544,99,596),
       C=c(124,978,544,99,596,528),
       D=c(978,544,99,596,528,124),
       E=c(544,99,596,528,124,978),
       F=c(99,596,528,124,978,544))

seed <- 0
foreach(i=1:length(fitRs$init)) %do% {
  seed = seed +1

  fitRs$init[[i]]$mcmcParam$seed[[1]] <- mainSeed[[seed]]
  fitRs$init[[i]]$mcmcParam$chainLength <- chainLength
  fitRs$init[[i]]$mcmcParam$burnin <- fitRs$init[[i]]$mcmcParam$chainLength %/% 10
  fitRs$init[[i]]$mcmcParam$verbosity  <- 0

  if (seed == 6) {
    seed <- 0
  }
}

### Run each model for an initial chain [it can took a while]
print("##### Running the optimizer first time #####")
runs$init <- foreach(i=1:length(fitRs$init)) %dopar% {
  hzar.doFit(fitRs$init[[i]])
}
runs$init <- hzar.copyModelLabels(fitRs$init, runs$init)

### Compile a new set of fit requests using the initial chains 
print("##### Creating new hzar.fitRequest object based on first optimizer #####")
fitRs$chains <- foreach(i=1:length(runs$init)) %dopar% {
  hzar.next.fitRequest(runs$init[[i]])
}

### Replicate each fit request 3 times, keeping the original seeds while switching to a new seed channel
print("##### Replicating hzar.fitRequest objects to 3 independed chains #####")
fitRs$chains <- hzar.multiFitRequest(fitRs$chains,
                                     each = 3,
                                     baseSeed = NULL)

# Just to be thorough, randomize the initial value for each fit
print("##### Randomizing the initial center and width for each fit #####")
foreach(i=1:length(fitRs$chains)) %do% {
  newcenter <- runif(1, mincenter, maxcenter)
  fitRs$chains[[i]]$modelParam$init["center"] <- newcenter
  newwidth <- runif(1, 0, maxwidth)
  fitRs$chains[[i]]$modelParam$init["width"] <- newwidth
}

if ("modelsF" %in% models2run){
  print("##### Randomizing the initial pMin and pMax for modelsF fit #####")

  foreach(i=31:length(fitRs$chains)) %do% { # So rodar para os scaling free
    newpmin <- runif(1, 0, 1)
    fitRs$chains[[i]]$modelParam$init["pMin"] <- newpmin
    newpmax <- runif(1, 0, 1)
    fitRs$chains[[i]]$modelParam$init["pMax"] <- newpmax
  }
}

print("##### Randomizing the initial delta and tau for tails fit #####")
foreach(i=1:length(fitRs$chains)) %do% {
  if ("deltaL" %in% names(fitRs$chains[[i]]$modelParam$init)) {
    newdelL <- runif(1, 0, maxwidth)
    fitRs$chains[[i]]$modelParam$init["deltaL"] <- newdelL
    newtauL <- runif(1, 0, 1)
    fitRs$chains[[i]]$modelParam$init["tauL"] <- newtauL
  }
  
  if ("deltaR" %in% names(fitRs$chains[[i]]$modelParam$init)) {
    newdelR <- runif(1, 0, maxwidth)
    fitRs$chains[[i]]$modelParam$init["deltaR"] <- newdelR
    newtauR <- runif(1, 0, 1)
    fitRs$chains[[i]]$modelParam$init["tauR"] <- newtauR
  }
  
  if ("deltaM" %in% names(fitRs$chains[[i]]$modelParam$init)) {
    newdelM <- runif(1, 0, maxwidth)
    fitRs$chains[[i]]$modelParam$init["deltaM"] <- newdelM
    newtauM <- runif(1, 0, 1)
    fitRs$chains[[i]]$modelParam$init["tauM"] <- newtauM
  }
}

### Run a chain of 3 runs for every fit request [it can took a while]
# Optimization to avoid run  hzar.chain.doSeq several times
print("##### Running optimizer 3x for all chains #####")
runs$chains <-  hzar.doChain.multi(fitRs$chains,
                                   doPar = TRUE,
                                   inOrder = TRUE,
                                   count = 3)

### Summary result parameters and log likelihood
print("##### Summarizing result parameters and log likelihood  for 3rd chain of all models #####")
foreach(i=1:length(runs$chains)) %dopar% {
  if (i %in% seq(1, 43, by = 3)) {
    fname <- paste("_model", i, "_summary.txt", sep = "")
    output <- gsub(".hzar", fname, hzar_output)
    sink(output)
    j=i+2
    print(summary(do.call(mcmc.list,
                    lapply(runs$chains[i:j],function(x) hzar.mcmc.bindLL(x[[3]])))))
    sink()
  }
}

### Create a model data group (hzar.dataGroup object) for each model from the initial runs [it can took a while]
print("##### Grouping initial fits from the same model and creating null model into hzar.dataGroup object #####")
analysis$initDGs <- foreach(i=1:length(runs$init)) %dopar% {
  hzar.dataGroup.add(runs$init[[i]])
}
analysis$initDGs <- hzar.copyModelLabels(runs$init, analysis$initDGs)

# Create a model data group for the null model (expected allele frequency independent of distance along cline)
analysis$initDGs$nullModel <- hzar.dataGroup.null(obs)

# Create a hzar.obsDataGroup object from all hzar.dataGroup just created
print("##### Creating hzar.obsDataGroup object with init runs #####")
analysis$oDG <- hzar.make.obsDataGroup(analysis$initDGs)
analysis$oDG <- hzar.copyModelLabels(analysis$initDGs, analysis$oDG)

# Convert all runs to hzar.dataGroup objects, adding them to the hzar.obsDataGroup object [it can took a while]
print("##### Adding chains into the hzar.obsDataGroup object #####")
analysis$oDG <- hzar.make.obsDataGroup(lapply(runs$chains,
                                              hzar.dataGroup.add), analysis$oDG)

### Compare cline models to the null model graphically
print("##### Plotting all clines together #####")
output <- gsub(".hzar", ".all_clines.png", hzar_output)
png(width=400, height=400, res=140, filename=output, pointsize=8)
hzar.plot.cline(analysis$oDG, pch = 20, xlim = c(mincenter,maxcenter))
dev.off()

### Do model selection based on the AICc scores and extract parameters and cline from best model
print("##### Saving AICc scores and writing files #####")
output <- gsub(".hzar", ".AICc_clines.csv", hzar_output)
analysis$AICcTable <- hzar.AICc.hzar.obsDataGroup(analysis$oDG)
df <- analysis$AICcTable
df['model'] <- rownames(df)
write.csv2(df, file = output, quote = FALSE, row.names = FALSE)

# Print out the model with the minimum AICc score
output <- gsub(".hzar", ".best_model.txt", hzar_output)
print(analysis$model.name <- rownames(analysis$AICcTable)[[ which.min(analysis$AICcTable$AICc )]])
write(analysis$model.name, file=output)

# Extract the hzar.dataGroup object for the selected model
analysis$model.selected <- analysis$oDG$data.groups[[analysis$model.name]]

# Saving best model chain parameters
output <- gsub(".hzar", ".best_model_mcmc.csv", hzar_output)
df <- hzar.mcmc.bindLL(analysis$model.selected, t0 = 1, tF = 100)
write.csv2(df, file = output, quote = FALSE, row.names = FALSE)

# Look at the variation in parameters for the selected model
print("##### Getting and saving parameters variation #####")
output <- gsub(".hzar", ".best_model_param_var.csv", hzar_output)
if (analysis$model.name != "nullModel") {
  df <- data.frame(hzar.getLLCutParam(analysis$model.selected,names(analysis$model.selected$data.param)))
  write.csv2(df, file = output, quote = FALSE, row.names = FALSE)
} else {
  sink(output)
  print("nullModel was the best model")
  sink()
}

# Print the maximum likelihood cline params for the selected model
output <- gsub(".hzar", ".best_model_param.csv", hzar_output)
df <- analysis$model.selected$ML.cline$param.free
write.csv2(df, file = output, quote = FALSE, row.names = FALSE)

# Plot the maximum likelihood cline with 95% credible cline region for the selected model
print("##### Ploting ML cline with 95% credible cline region #####")
output <- gsub(".hzar", ".best_model_cline95.png", hzar_output)
png(width=400, height=400, res=140, filename=output, pointsize=8)
hzar.plot.fzCline(analysis$model.selected, pch = 20, xlim = c(mincenter,maxcenter))
dev.off()

# Plot final trace for selected model
print("##### Saving best model trace #####")
output <- gsub(".hzar", ".best_model_trace.pdf", hzar_output)
pdf(file=output)
plot(hzar.mcmc.bindLL(analysis$model.selected))
dev.off()
