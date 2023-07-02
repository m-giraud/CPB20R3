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



getVal <-function(myPath,valName)
{
  
  x <- readLines(paste0(myPath,valName), n=-1)
  x = strsplit(x, ",");
  ncol =length(x[[length(x)]])
  df = read.table(paste0(myPath,valName), sep=",",header = FALSE,fill=TRUE, 
                     col.names = paste0("V",seq_len(ncol)))
  time =t(read.table(paste0(myPath,"time.txt"), sep=",",header = FALSE))
  df$time = time[1,]
  df <- melt(df, id.vars="time")
  df = df[df$time>(min(df$time)+2/24),]
  
  return(df)
}



get_CST_OT <-function(myPath)
{
  cst = getVal(myPath,"C_ST.txt")
  ots = getVal(myPath,"ots.txt")
  cst$ot = ots$value
  cst = cst[!is.na(cst$ot),]
  cst[cst$ot==0, ]$ot=3 #ot == 0 : seed node. set to stem node
  print(min(cst$value))
  print(max(cst$value))
  cst = aggregate(cst$value,
                      FUN = mean,
                      by=list("time"= cst$time, "ot" = cst$ot))
  colnames(cst)[3]="cst"
  return(cst)
}

getValAll <-function(valName)
{
  dfearlyDry = get_CST_OT(earlyDry )
  dfearlyDry$scenario = "early dry spell"
  dfearlyDry$scenario2 = "b) early dry spell"
  dflateDry = get_CST_OT(lateDry )
  dflateDry$scenario = "late dry spell"
  dflateDry$scenario2 = "c) late dry spell"
  dfbaseline = get_CST_OT(baseline )
  dfbaseline$scenario = "baseline"
  dfbaseline$scenario2 = "a) baseline"
  
  df = rbind(dfearlyDry,dflateDry,dfbaseline)
  df$var = valName
  return(df)
}

cst = getValAll("cst")
cst$ot = as.factor(cst$ot)


var.lab3 = c("2" = "(A) root", 
             "3" = "(B) stem", 
             "4" = "(C) leaf")

vlabeller <- function (variable, value) {
  return(var.lab3[value])
}


both=ggplot(data=cst, aes(x=time, y=cst,linetype =scenario,
                          size=scenario,alpha=scenario))+
  xlab('day of growth (d)')+
  geom_line()+
  scale_x_continuous(breaks=seq(from =6,to=28,by=2))+
  geom_vline(aes(xintercept=11), colour="black") + 
  geom_vline( aes(xintercept=18), colour="black") +
  geom_vline(aes(xintercept=25), colour="black") + 
  ylab(bquote(s[st]~" (mmol/ml)"))+ 
  facet_wrap(~ot, ncol = 1,labeller =vlabeller)+
  scale_size_manual(breaks=c("early dry spell","late dry spell","baseline"),
                    name="",
                    values = c(0.5,0.5,1.5))+ 
  scale_alpha_manual(breaks=c("early dry spell","late dry spell","baseline"),
                     name="",
                     values = c(1,1,0.3))+ 
  scale_linetype_manual(breaks=c("early dry spell","late dry spell","baseline"),name="",
                        values = c("dashed","solid","solid"))+
  theme(legend.text = element_text( size=15),
        legend.title = element_text( size=15),
        legend.position="bottom",#c(0.5,0.1),
        legend.direction = "horizontal",
        #legend.box = "vertical",
        panel.grid.minor = element_line(colour="white"),
        panel.grid.major = element_line(colour="white"),
        panel.background = element_rect(fill="white", color="black"),
        strip.text.y =element_text( size=15) ,
        axis.text = element_text( size=17),
        axis.title = element_text( size=21),
        strip.background = element_rect(fill="white"),
        strip.text. = element_text(size = 21, hjust = 0,face ="bold"),
        strip.placement = "inside");both
ggsave('sucCons.png',plot=both,  width = 8, height = 10, dpi = 700)
