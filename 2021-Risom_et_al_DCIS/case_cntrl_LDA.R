####################################################################################################################
#
# Script: SCMP_axes.R
# Project: Multiplexed Single Cell Morphometry for Hematopathology Diagnostics
# Author: David Glass
# Date: 1-6-20
# 
# Purpose: Generate dimensionality-reduced axes on a dataset
#
# Pseudocode:
# Load in data.table with observations (cells) as rows and parameters (proteins) as columns.
# Subset data with class labels (e.g. manual gate annotations) to train axes.
# Identify combination of channels that optimally separates classes on 1- or 2-D axes using LDA.
# Pull out coefficients from LDA and generate new axes for total dataset.
#
# Instructions:
# Install all required packages (see LIBRARIES and/or Versions below)
# Under USER INPUTS, assign "path" variable to the directory containing a csv with processed data
#   SCMP_preprocess.R can be used to process data from Figure 4.
#   Any processed cytometry dataset can be loaded. scRNA-seq data may work as well, but has not been tested.
#   Data should be transformed and scaled appropriately. For our mass cytometry dataset, counts were asinh
#     transformed with a cofactor of 5 and each marker was scaled to the 99.5th percentile.
#   Training data must have a column of class labels (e.g. manual gate annotations). This is not required
#     for the rest of the dataset.
# Assign "channels" variable to a list of column names that you wish to be explored for generating axes
# Under MAIN (bottom of script), assign "train.x" variable to the training data with only "channels" as columns.
#   In our example dataset, this includes all observations that are not "ungated"
# Under MAIN, assign "train.y" variable to a vector of training data labels
# Run script
# Under MAIN "coefficients" variable contains the coefficients that are used to generate the new axes. This variable
#   is fed into the makeAxes function.
# makeAxes uses "coefficients" to create coordinates for new axes on the entire dataset.
# 
# Versions:
# R 3.5.1
# RStudio 1.1.463
# MASS_7.3-51.5
# dplyr_0.8.3
# data.table_1.12.8
#
#######################################################################################################################





###### LIBRARIES ######

require(MASS)
require(dplyr)
require(data.table)



##### FUNCTIONS ##### 

hybridSubsetSelection <- function(x, y, two.d=TRUE) {
  # performs hybrid stepwise subset selection on LDA reduced dimensions
  # Inputs:
  #   x - data.table to evaluate, rows are cells, columns are columns to evaluate
  #   y - vector of observation classes
  #   two.d - logical if true creates two new axes, if false only one
  # Outputs:
  #   matrix of coefficients where rows are markers and columns are axes
  
  ### global data structures ###
  keep <- NULL
  channels <- colnames(x)
  n.channels <- length(channels)
  current.score <- 0
  continue <- TRUE
  results <- setNames(data.frame(matrix(nrow=1, ncol=n.channels)), channels)
  results[1,] <- as.list(rep(F, n.channels))
  subtract.log <- data.table(results[0,]) # record of keep values inputted into subtractOne
  results$score <- 0
  results <- data.table(results)
  
  
  ### functions ###
  addOne <- function() {
    # Evaluates the addition of each channel not in keep to keep. Adds best and updates current.score
    temp.results <- results[0,]
    for (channel in channels[!channels %in% keep]) {
      temp.keep <- c(keep, channel)
      temp.score <- getScore(temp.keep)
      temp.results <- rbind(temp.results, as.list(channels %in% temp.keep) %>% append(temp.score))
    }
    current.score <<- max(temp.results$score)
    new.keep <- temp.results[score==current.score, channels, with=F]
    if (nrow(new.keep) > 1) new.keep <- new.keep[sample(.N,1)]
    keep <<- channels[as.logical(new.keep)]
    results <<- unique(rbind(results, temp.results))
  }
  
  subtractOne <- function() {
    # Evaluates the subtraction of each channel from keep. Removes worst if it improves score and updates current.score
    # If a better subset is found, it calls itself.
    # If this keep has been evaluted before, exits
    subtract.log <<- rbind(subtract.log, as.list(channels %in% keep))
    if (anyDuplicated(subtract.log) > 0) {
      subtract.log <<- unique(subtract.log)
      return()
    }
    temp.results <- results[0,]
    for (channel in keep) {
      temp.keep <- keep[!keep %in% channel]
      temp.score <- getScore(temp.keep)
      temp.results <- rbind(temp.results, as.list(channels %in% temp.keep) %>% append(temp.score))
    }
    if (max(temp.results$score) > current.score) {
      current.score <<- max(temp.results$score)
      new.keep <- temp.results[score==current.score, channels, with=F]
      if (nrow(new.keep) > 1) new.keep <- new.keep[sample(.N,1)]
      keep <<- channels[as.logical(new.keep)]
      results <<- unique(rbind(results, temp.results))
      subtractOne()
    } else results <<- unique(rbind(results, temp.results))
  }
  
  getScore <- function(cols) {
    # performs LDA using columns provided and returns lowest euclidean distance between pop means
    lda.out <- lda(y~., data=x[, cols, with=F])
    if (two.d) return(min(dist(lda.out$means %*% lda.out$scaling[,1:2])))
    return(min(dist(lda.out$means %*% lda.out$scaling[,1])))
  }
  
  initializeKeep <- function() {
    # chooses the best scoring pair of markers to initialize keep
    temp.results <- results[0,]
    for (channel.1 in channels) {
      for (channel.2 in channels) {
        if (channel.1==channel.2) next
        temp.keep <- c(channel.1, channel.2)
        temp.score <- getScore(temp.keep)
        temp.results <- rbind(temp.results, as.list(channels %in% temp.keep) %>% append(temp.score))
      }
    }
    current.score <<- max(temp.results$score)
    new.keep <- temp.results[score==current.score, channels, with=F]
    if (nrow(new.keep) > 1) new.keep <- new.keep[sample(.N,1)]
    keep <<- channels[as.logical(new.keep)]
    results <<- unique(rbind(results, temp.results))
  }
  
  getElbow <- function(res) {
    # takes results and returns the elbow point
    res[, no.markers:=apply(res[, channels, with=F], 1, sum)]
    res.lite <- res[, max(score), by=no.markers][!1]
    slope <- (res.lite$V1[nrow(res.lite)] - res.lite$V1[1]) / (res.lite$no.markers[nrow(res.lite)] - res.lite$no.markers[1])
    intercept <- res.lite$V1[1] - slope * res.lite$no.markers[1]
    perp.slope <- -1 / slope
    perp.int <- res.lite$V1 - (perp.slope * res.lite$no.markers)
    xcross <- (intercept - perp.int) / (perp.slope - slope)
    ycross <- slope * xcross + intercept
    dists <- sqrt((res.lite$no.markers - xcross)^2 + (res.lite$V1 - ycross)^2) %>% round(1)
    elbowi <- max(which(dists==max(dists))) # if dists are tie, take the largest number of channels
    return(elbowi+1)
  }
  
  ### main ###
  initializeKeep()
  while(continue) {
    print(paste("Number of markers:", length(keep)))
    addOne()
    print(paste("Number of markers:", length(keep)))
    if (length(keep) > 3) subtractOne()
    if (length(keep)==length(channels)) continue <- FALSE
  }
  elbow <- getElbow(results)
  markers <- results[no.markers==elbow] %>%
    .[score==max(.[, score]), channels, with=F] %>%
    unlist() %>%
    names(.)[.]
  lda.out <- lda(y~., data=x[, markers, with=F])
  if(two.d) return(lda.out$scaling[, 1:2])
  return(lda.out$scaling[, 1, drop=F])
}


makeAxes <- function(dt=dat, co=coefficients, axis.name="ld") {
  # makes new axes based on coefficients
  # Inputs:
  #   dt - data.table of data
  #   co - matrix of coefficients
  #   axis.name - character vector of new axis name (e.g. "ld" results in "ld1" and "ld2")
  # Outputs:
  #   dt - data.table
  x <- as.matrix(dt[, rownames(co), with=F])
  for (i in 1:ncol(co)) {
    dt[, eval(paste0(axis.name, i)):=x %*% co[, i]]
  }
  return(dt)
}



##### MAIN #####
##### USER INPUTS #####

# read input table
path <- "/Users/innaa/Documents/Tyler DCIS/case_cntrl_myoep_LDA.csv"
dat <- fread(path)

# training data with appropriate channels
inds_use = c(0:(ncol(dat)-2))
train.x <- dat[gate!="ungated",inds_use, with=F]
train.y <- dat[gate!="ungated",label]

# run LDA 
coefficients <- hybridSubsetSelection(x=train.x, y=train.y, two.d=FALSE)
dat <- makeAxes()

# plot violin by tissue type
library(ggplot2)
p <- ggplot(dat, aes(y=ld1, x=label)) + 
  geom_violin()
p + geom_boxplot(width=0.1)+ geom_jitter(shape=16, position=position_jitter(0.2))

# save result
write.csv(dat, file = "/Users/innaa/Documents/Tyler DCIS/case_cntrl_myoep_LDA_table.csv")
write.csv(coefficients, file = "/Users/innaa/Documents/Tyler DCIS/case_cntrl_myoep_LDA_coeffs.csv")


