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
  sinks$Rm = getVal(myPath,"Q_Rm_dot.txt") #maintenance
  sinks$gr = getVal(myPath,"Q_Gr_dot.txt" ) #growth
  sinks$Exud = getVal(myPath,"Q_Exud_dot.txt" )#exudation
  colnames(sinks) = c("(A) Rm",
                      "(B) G\U209C\U2092\U209C\U02D2\U1D04\u1d21\U2097\U1D62\U2098",
                       "(C) Exud")
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

p1=ggplot(data=sinks, 
       aes(x=time, y=value,#col=scenario,
           size=scenario,alpha=scenario,linetype=scenario
           ))+
  geom_line()+
  xlab('day of growth (d)')+
  ylab("sucrose usage rate (mmol/d)")+ 
  scale_x_continuous(breaks=seq(from =6,to=28,by=2))+
  #ylim(0,0.8)+ 
  geom_vline(aes(xintercept=11), colour="black") + 
  geom_vline(aes(xintercept=18), colour="black") + 
  geom_vline(aes(xintercept=25), colour="black") + 
  facet_wrap(~variable, ncol = 1,nrow=3)+#,scales = "fixed")+
  #scale_color_manual(name="",
  #                   breaks=c("early dry spell","late dry spell","baseline"),
  #                   values = c("#8DA0CB","#FC8D62","#66C2A5"))+ 
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
        strip.text.y =element_blank(),#element_text( size=15) ,
        axis.text = element_text( size=17),
        axis.title = element_text( size=21),
        #panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        strip.background = element_rect(fill="white"),
        panel.background = element_rect(fill="white", color="black"),
        strip.text. = element_text(size = 21, hjust = 0,face ="bold"),
        panel.grid.minor = element_line(colour=NA),
        panel.grid.major = element_line(colour=NA),
        strip.placement = "outside");p1
ggsave('sucSink.png',plot=p1,  width = 12, height = 15, dpi = 700)


p1=ggplot(data=sinks[(sinks$time>18.3)&(sinks$time<22.3),], 
          aes(x=time, y=value,#col=scenario,
              size=scenario,alpha=scenario,linetype=scenario
          ))+
  geom_line()+
  xlab('day of growth (d)')+
  ylab("sucrose usage rate (mmol/d)")+ 
  scale_x_continuous(breaks=seq(from =6,to=28,by=2))+
  #ylim(0,0.8)+ 
  facet_wrap(~variable, ncol = 1,nrow=3)+#,scales = "fixed")+
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
        strip.text.y =element_blank(),#element_text( size=15) ,
        axis.text = element_text( size=17),
        axis.title = element_text( size=21),
        #panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        strip.background = element_rect(fill="white"),
        panel.background = element_rect(fill="white", color="black"),
        strip.text. = element_text(size = 21, hjust = 0,face ="bold"),
        panel.grid.minor = element_line(colour=NA),
        panel.grid.major = element_line(colour=NA),
        strip.placement = "outside");p1
ggsave('sucSinkZOOM.png',plot=p1,  width = 12, height = 15, dpi = 700)
