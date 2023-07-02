library(stats)
library(babynames)
library(dplyr)
library(naniar)
library(tidyverse)
require(boot)
require(matrixStats)
require(lubridate)
require(dplyr)
require(reshape2)
require(plotly)
require(ggrepel)
library(stringr)
library(ggplot2)

#paths
earlyDry="C:/Users/m.giraud/Documents/article_calibration/toutputSim/data/aggregated_0302/6to16dry530/"
lateDry="C:/Users/m.giraud/Documents/article_calibration/toutputSim/data/aggregated_0302/6to28dry530/"
baseline="C:/Users/m.giraud/Documents/article_calibration/toutputSim/data/aggregated_0302/6to16wet530/"


getVal <-function(myPath, valName)
{
  x <- readLines(paste0(myPath,valName), n=-1)
  x = strsplit(x, ",");
  ncol =length(x[[length(x)]])
  df = read.table(paste0(myPath,valName), sep=",",header = FALSE,fill=TRUE, 
                  col.names = paste0("V",seq_len(ncol)))
  df =rowSums(df, na.rm = TRUE)
  return(df)
}

getValAllSinks <-function(myPath)
{
  time =t(read.table(paste0(myPath,"time.txt"), sep=",",header = FALSE))
  sinks = data.frame( data.frame(matrix(nrow = length(time),ncol=0)))
  sinks$Rm = getVal(myPath,"Q_Rm.txt") #maintenance
  sinks$gr = getVal(myPath,"Q_Gr.txt" ) #growth
  sinks$Exud = getVal(myPath,"Q_Exud.txt" )#exudation
  colnames(sinks) = c("a) Rm",
                      "b) Gto,cwlim",
                      "c) Exud")
  sinks$time = time[1,]
  sinks = sinks[sinks$time>(min(sinks$time)+2/24),]
  return(sinks)
}

getValAll <-function()
{
  dfearlyDry = getValAllSinks(earlyDry)
  dfearlyDry$scenario = "early dry spell"
  dflateDry = getValAllSinks(lateDry)
  dflateDry$scenario = "late dry spell"
  dfbaseline = getValAllSinks(baseline)
  dfbaseline$scenario = "baseline"
  
  df = rbind(dfearlyDry,dflateDry,dfbaseline)
  df <- melt(df, id.vars=c("time", "scenario"))
  return(df)
}

sinks = getValAll()

SinkEnd = sinks[sinks$time==max(sinks$time),-1]
dat_2 <- reshape(SinkEnd,                          # Reshape data from long to wide format
                 timevar   = c("variable"), 
                 idvar     = c("scenario"), 
                 direction = "wide")
write.csv(dat_2, "C:/Users/m.giraud/Documents/article_calibration/toutputSim/data/aggregated_0302/SinkEnd.csv", row.names=FALSE)
