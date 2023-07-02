library(stats)
library(babynames)
library(dplyr)
library(naniar)
library("viridis")           # Load
library(ggpubr)
library(RColorBrewer)
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


getVal <-function(myPath, valName)
{
  x <- readLines(paste0(myPath,valName), n=-1)
  x = strsplit(x, ",");
  ncol =length(x[[length(x)]])
  df = read.table(paste0(myPath,valName), sep=",",header = FALSE,fill=TRUE, 
                  col.names = paste0("V",seq_len(ncol)))
  time =t(read.table(paste0(myPath,"time.txt"), sep=",",header = FALSE))
  df$time = time[1,]
  df <- melt(df, id.vars="time")[,-2]
  return(df)
}

getValAllstrC <-function(myPath)
{
  time =t(read.table(paste0(myPath,"time.txt"), sep=",",header = FALSE))
  sinks = data.frame( data.frame(matrix(nrow = length(time),ncol=0)))
  df = getVal(myPath,"volOrg.txt");colnames(df)[2]=c("vol")
  df$ot = getVal(myPath,"ot_orgs.txt" )$value #organtype
  df$st = getVal(myPath,"st_orgs.txt" )$value#suptype
  df = df[!is.na(df$ot),]
  df = df[df$time>(min(df$time)+2/24),]
  return(df)
}

getValAll <-function()
{
  dfearlyDry = getValAllstrC(earlyDry)
  dfearlyDry$scenario = "early dry spell"
  dflateDry = getValAllstrC(lateDry)
  dflateDry$scenario = "late dry spell"
  dfbaseline = getValAllstrC(baseline)
  dfbaseline$scenario = "baseline"
  
  df = rbind(dfearlyDry,dflateDry,dfbaseline)
  df$rho = ifelse(df$ot ==2, 0.51, 
                    ifelse(df$ot==3,0.65,0.56))
  df$strC = df$vol * df$rho
  
  df$organTypes = ifelse(df$ot==2,
                            ifelse(df$st==1,"(A) root—order0",
                                   ifelse(df$st==2 ,"(B) root—order1","(C) root—order2")),
                            ifelse(df$ot == 3, "(D) stem",
                                   "(E) leaf"))
  dfm = aggregate(df$strC,
                  FUN = sum,
                  by=list("time" = df$time,
                          "organTypes" = df$organTypes,
                          "scenario" = df$scenario))
  dfm$organTypes =factor(dfm$organTypes, 
                           levels = c("(A) root—order0","(B) root—order1","(C) root—order2","(D) stem","(E) leaf"))
  return(dfm)
}

strC = getValAll()

p4=ggplot(data=strC, 
       aes(x=time, y=x,linetype =scenario,
           #col =scenario,
           size=scenario,alpha=scenario))+
  geom_line()+
  scale_size_manual(breaks=c("early dry spell","late dry spell","baseline"),
                    name="",
                    values = c(0.5,0.5,1.5))+ 
  scale_alpha_manual(breaks=c("early dry spell","late dry spell","baseline"),
                     name="",
                     values = c(1,1,0.3))+ 
  scale_linetype_manual(breaks=c("early dry spell","late dry spell","baseline"),name="",
                        values = c("dashed","solid","solid"))+
  xlab('day of growth (d)')+
  facet_wrap(~organTypes, ncol = 2,nrow=3,scales = "free_y")+
  scale_x_continuous(breaks=seq(from =6,to=28,by=2))+
  ylab("structural sucrose (mmol)")+ 
  geom_vline(aes(xintercept=11), colour="black") + 
  geom_vline(aes(xintercept=18), colour="black") + 
  geom_vline(aes(xintercept=25), colour="black") + 
  theme(legend.text = element_text( size=15),
        legend.title = element_blank(),
        legend.position=c(0.6,0.2),
        #legend.direction = "horizontal",
        #legend.box = "vertical",
        panel.grid.minor = element_line(colour="white"),
        panel.grid.major = element_line(colour="white"),
        panel.background = element_rect(fill="white", color="black"),
        strip.text.y =element_blank(),
        axis.text = element_text( size=17),
        axis.title = element_text( size=21),
        strip.background = element_rect(fill="white"),
        strip.text. = element_text(size = 21, hjust = 0,face ="bold"),
        strip.placement = "inside");p4

ggsave('structuralSuc.png',plot=p4,  width = 10, height = 10, dpi = 700)




sumStr = aggregate(strC$x,
                   FUN = sum,
                   by=list("time" = strC$time,
                           "scenario" = strC$scenario))
colnames(sumStr)[3]="sum"
strC = merge(strC, sumStr, all=TRUE, 
               by=c("time","scenario"))


strC$ratio = strC$x/strC$sum*100


p1=ggplot(data=strC, 
          aes(x=time, y=ratio,linetype =scenario,
              #col =scenario,
              size=scenario,alpha=scenario))+
  geom_line()+
  scale_size_manual(breaks=c("early dry spell","late dry spell","baseline"),
                    name="",
                    values = c(0.5,0.5,1.5))+ 
  scale_alpha_manual(breaks=c("early dry spell","late dry spell","baseline"),
                     name="",
                     values = c(1,1,0.3))+ 
  scale_linetype_manual(breaks=c("early dry spell","late dry spell","baseline"),name="",
                        values = c("dashed","solid","solid"))+
  xlab('day of growth (d)')+
  facet_wrap(~organTypes, ncol = 2,nrow=3,scales = "free_y")+
  scale_x_continuous(breaks=seq(from =6,to=28,by=2))+
  ylab("ratio of structural sucrose (%)")+ 
  geom_vline(aes(xintercept=11), colour="black") + 
  geom_vline(aes(xintercept=18), colour="black") + 
  geom_vline(aes(xintercept=25), colour="black") + 
  theme(legend.text = element_text( size=15),
        legend.title = element_blank(),
        legend.position=c(0.6,0.2),
        #legend.direction = "horizontal",
        #legend.box = "vertical",
        panel.grid.minor = element_line(colour="white"),
        panel.grid.major = element_line(colour="white"),
        panel.background = element_rect(fill="white", color="black"),
        strip.text.y =element_blank(),
        axis.text = element_text( size=17),
        axis.title = element_text( size=21),
        strip.background = element_rect(fill="white"),
        strip.text. = element_text(size = 21, hjust = 0,face ="bold"),
        strip.placement = "inside");p1
ggsave('structuralSucRATIO.png',plot=p1,  width = 10, height = 10, dpi = 700)

###
#
# get final sucrose partitioning
###
strCEnd = strC[strC$time==max(strC$time),c(2,3,6)]
dat_2 <- reshape(strCEnd,                          # Reshape data from long to wide format
                 timevar   = c("scenario"), 
                 idvar     = c("organTypes"), 
                 direction = "wide")
dat_2$organTypes = as.character(dat_2$organTypes)
dat_2[6,-1] =  colSums(dat_2[c(1,2,3),-1])

dat_2[6,1] = "root"


###
#
# get root length
###


getVal <-function(myPath, valName)
{
  x <- readLines(paste0(myPath,valName), n=-1)
  x = strsplit(x, ",");
  ncol =length(x[[length(x)]])
  df = read.table(paste0(myPath,valName), sep=",",header = FALSE,fill=TRUE, 
                  col.names = paste0("V",seq_len(ncol)))
  time =t(read.table(paste0(myPath,"time.txt"), sep=",",header = FALSE))
  df$time = time[1,]
  df <- melt(df, id.vars="time")[,-2]
  return(df)
}

getValAllstrC <-function(myPath)
{
  time =t(read.table(paste0(myPath,"time.txt"), sep=",",header = FALSE))
  sinks = data.frame( data.frame(matrix(nrow = length(time),ncol=0)))
  df = getVal(myPath,"len_orgs.txt");colnames(df)[2]=c("lenOrg")
  df$ot = getVal(myPath,"ot_orgs.txt" )$value #organtype
  df$st = getVal(myPath,"st_orgs.txt" )$value#suptype
  df = df[!is.na(df$ot),]
  df = df[df$time>(min(df$time)+2/24),]
  return(df)
}

getValAll <-function()
{
  dfearlyDry = getValAllstrC(earlyDry)
  dfearlyDry$scenario = "early dry spell"
  dflateDry = getValAllstrC(lateDry)
  dflateDry$scenario = "late dry spell"
  dfbaseline = getValAllstrC(baseline)
  dfbaseline$scenario = "baseline"
  
  df = rbind(dfearlyDry,dflateDry,dfbaseline)
  
  df$organTypes = ifelse(df$ot==2,
                         ifelse(df$st==1,"(A) root—order0",
                                ifelse(df$st==2 ,"(B) root—order1","(C) root—order2")),
                         ifelse(df$ot == 3, "(D) stem",
                                "(E) leaf"))
  dfm = aggregate(df$lenOrg,
                  FUN = sum,
                  by=list("time" = df$time,
                          "organTypes" = df$organTypes,
                          "scenario" = df$scenario))
  dfm$organTypes =factor(dfm$organTypes, 
                         levels = c("(A) root—order0","(B) root—order1","(C) root—order2","(D) stem","(E) leaf"))
  return(dfm)
}

lenOrg = getValAll()


###
#
# get final length
###
lenOrgEnd = lenOrg[lenOrg$time==max(lenOrg$time),-1]
dat_2 <- reshape(lenOrgEnd,                          # Reshape data from long to wide format
                 timevar   = c("scenario"), 
                 idvar     = c("organTypes"), 
                 direction = "wide")
dat_2$organTypes = as.character(dat_2$organTypes)
dat_2[6,-1] =  colSums(dat_2[c(1,2,3),-1])

dat_2[6,1] = "root"
