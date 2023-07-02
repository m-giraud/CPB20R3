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
earlyDry="../results/earlyDry/"
lateDry="../results/lateDry/"
baseline="../results/baseline/"


getVal <-function(myPath, valName, multFactor)
{
  x <- readLines(paste0(myPath,valName), n=-1)
  x = strsplit(x, ",");
  ncol =length(x[[length(x)]])
  df = read.table(paste0(myPath,valName), sep=",",header = FALSE,fill=TRUE, 
                     col.names = paste0("V",seq_len(ncol)))
  dfsum =as.data.frame(rowSums(df, na.rm = TRUE));colnames(dfsum) = c('sum')
  time =t(read.table(paste0(myPath,"time.txt"), sep=",",header = FALSE))
  dfsum$time = time[1,]
  dfsum = dfsum[dfsum$time>(min(dfsum$time)+2/24),]
  dfsum$cumsum = cumsum(dfsum[,1] * multFactor)
  return(dfsum)
}
getValAll <-function(valName,fileName, multFactor)
{
  dfearlyDry = getVal(earlyDry,fileName, multFactor )
  dfearlyDry$scenario = "early dry spell"
  dflateDry = getVal(lateDry,fileName, multFactor )
  dflateDry$scenario = "late dry spell"
  dfbaseline = getVal(baseline,fileName , multFactor)
  dfbaseline$scenario = "baseline"
  
  df = rbind(dfearlyDry,dflateDry,dfbaseline)
  df$var = valName
  return(df)
}


############
#
# transpiration and assimilation (rate and cumsum)
#
############
trans = getValAll("Transpiration","transrate.txt", 1)
Ag = getValAll("Ag","AgPhl.txt", 1/24)

transAg = rbind(trans,Ag)

###############
#
# Graphic
#
###############

#water =18.01528 g/mol
transAg$cumsum2 = transAg$cumsum
ratioTrans = 250 
transAg[transAg$var == "Transpiration",]$cumsum2 = 
  transAg[transAg$var == "Transpiration",]$cumsum *1000*(1/18.01528)/ratioTrans
#cumsum_cm3 * (mg/cm3) * (mmol/mg)  * scalingFactor (for visualisation)


var.lab3 = c("Transpiration" = "cumulative transpiration", 
             "Ag" = bquote("cumulative"~A["g"]))

vlabeller <- function (variable, value) {
  return(var.lab3[value])
}


both=ggplot(data=transAg, aes(x=time, y=cumsum2,col =var,linetype =scenario,
                                size=scenario,alpha=scenario))+
  geom_line()+
  scale_x_continuous(breaks=seq(from =6,to=28,by=2))+
  scale_y_continuous(sec.axis = sec_axis(~.*ratioTrans, name = "cumulative transpiration (mmol)"))+
  xlab('day of growth (d)')+
  scale_size_manual(breaks=c("early dry spell","late dry spell","baseline"),name="",
                    values = c(0.5,0.5,1.5))+ 
  scale_alpha_manual(breaks=c("early dry spell","late dry spell","baseline"),name="",
                     values = c(1,1,0.3))+ 
  scale_linetype_manual(breaks=c("early dry spell","late dry spell","baseline"),name="",
                     values = c("dashed","solid","solid"))+
  ylab(bquote("cumulative"~A["g,plant"]~" (mmol)"))+#'cumulative Ag (mmol)')+ 
  geom_vline(data=filter(transtemp, growth==unique(transtemp$growth)[1]), aes(xintercept=11), colour="black") + 
  geom_vline(data=filter(transtemp, growth==unique(transtemp$growth)[1]), aes(xintercept=18), colour="black") + 
  geom_vline(data=filter(transtemp, growth==unique(transtemp$growth)[2]), aes(xintercept=18), colour="black") + 
  geom_vline(data=filter(transtemp, growth==unique(transtemp$growth)[2]), aes(xintercept=25), colour="black") + 
  scale_color_manual(breaks=c("Transpiration","Ag"),name="",
                     labels=vlabeller,#c("cumulative transpiration","cumulative Ag"),
                     values = c("#8DA0CB","#FC8D62"),
                     guide = "none")+  
  #labs(colour = "", linetype ="", alpha ="")+
  theme(legend.text = element_text( size=15),
        legend.title = element_text( size=15),
        legend.position=c(0.3,0.5),
        #legend.direction = "horizontal",
        legend.box = "vertical",
        panel.grid.minor = element_line(colour="white"),
        panel.grid.major = element_line(colour="white"),
        panel.background = element_rect(fill="white", color="black"),
        strip.text.y =element_text( size=15) ,
        axis.text = element_text( size=17),
        axis.title.x = element_text( size=21),
        axis.title.y.left = element_text( size=21,color = "#FC8D62"),
        axis.title.y.right = element_text( size=21,color = "#8DA0CB"),
        axis.line.y.left = element_line(color = "#FC8D62"),
        axis.ticks.y.left = element_line(color = "#FC8D62"),
        axis.text.y.left = element_text(color = "#FC8D62"),
        axis.line.y.right = element_line(color = "#8DA0CB"),
        axis.ticks.y.right = element_line(color = "#8DA0CB"),
        axis.text.y.right = element_text(color = "#8DA0CB"),
        strip.background = element_rect(fill="white"),
        strip.text. = element_text(size = 21, hjust = 0,face ="bold"),
        strip.placement = "inside",
        legend.margin=margin()) ;both
ggsave('AgTrans.png',plot=both,  width = 10, height = 10, dpi = 700)



###
#
# get final Ag, trans
###
transAgEnd = transAg[transAg$time==max(transAg$time),c(3,4,5)]
dat_2 <- reshape(transAgEnd,                          # Reshape data from long to wide format
                 timevar   = c("scenario"), 
                 idvar     = c("var"), 
                 direction = "wide")
dat_2[1,-1] = dat_2[1,-1]*1000*(1/18.01528)
dat_2[3,-1] =  dat_2[2,-1]/ dat_2[1,-1]

dat_2[3,1] = "WUE"
dat_2[,-1]/dat_2[,]$cumsum.baseline*100 -100


###
#
# get fw
###

getVal <-function(myPath, valName)
{
  x <- readLines(paste0(myPath,valName), n=-1)
  x = strsplit(x, ",");
  ncol =length(x[[length(x)]])
  df = read.table(paste0(myPath,valName), sep=",",header = FALSE,fill=TRUE, 
                  col.names = paste0("V",seq_len(ncol)))
  df[df==0]=NA
  df = as.data.frame(apply(df, 1, FUN = min, na.rm = TRUE))#rowMins(df, na.rm = TRUE))
  time =t(read.table(paste0(myPath,"time.txt"), sep=",",header = FALSE))
  df$time = time[1,]
  df = df[df$time>(min(df$time)+2/24),]
  return(df)
}
getValAll <-function(valName,fileName)
{
  dfearlyDry = getVal(earlyDry,fileName)
  dfearlyDry$scenario = "early dry spell"
  dflateDry = getVal(lateDry,fileName )
  dflateDry$scenario = "late dry spell"
  dfbaseline = getVal(baseline,fileName )
  dfbaseline$scenario = "baseline"
  
  df = rbind(dfearlyDry,dflateDry,dfbaseline)
  colnames(df)[1] = valName
  return(df)
}

allFw = getValAll("fw","fw.txt")
dfm = aggregate(allFw$fw,
                FUN = min,
                by=list("scenario" = allFw$scenario))

