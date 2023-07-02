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
library(reticulate)
numpy <- import("numpy")


#paths
earlyDry="../results/earlyDry/"
lateDry="../results/lateDry/"
baseline="../results/baseline/"



getVal <-function(myPath,valName, doMelt = TRUE, valuesIndexes = c(), valueNames=c())
{
  
  x <- readLines(paste0(myPath,valName), n=-1)
  x = strsplit(x, ",");
  ncol =length(x[[length(x)]])
  
  df = read.table(paste0(myPath,valName), sep=",",header = FALSE,fill=TRUE, 
                  col.names = paste0("V",seq_len(ncol)))
  time =t(read.table(paste0(myPath,"time.txt"), sep=",",header = FALSE))
  
  if(doMelt)
  {
    df = df[,valuesIndexes]
    colnames(df) = valueNames
    df$time = time[1,]
    df <- melt(df, id.vars="time")
  }
  
  return(df)
}



get_Val_OT <-function(myPath, filename, valname, valuesIndexes = c(), valueNames=c())
{
  val = getVal(myPath,filename, TRUE, valuesIndexes, valueNames)
  ots = getVal(myPath,"ots.txt", TRUE, valuesIndexes, valueNames)
  val$ot = ots$value
  val = val[!is.na(val$ot),]
  val[val$ot==0, ]$ot=3 #ot == 0 : seed node. set to stem node
  val$variableName= valname
  val = val[val$time>(min(val$time)+2/24),]
  return(val)
}


## select index of node to visualize (need to exist at the first time step)
ots_firstTimeStep = getVal(baseline,"ots.txt", FALSE)[1,]; ots_firstTimeStep=ots_firstTimeStep[!is.na(ots_firstTimeStep)]
cst_firstTimeStep = getVal(baseline,"C_ST.txt", FALSE)[1,]; cst_firstTimeStep=cst_firstTimeStep[!is.na(cst_firstTimeStep)]

idroot = which(cst_firstTimeStep == min(cst_firstTimeStep))[[1]]
idstem = 1 #seed Id 
idleaf = which(cst_firstTimeStep == max(cst_firstTimeStep))[[1]]

stopifnot(ots_firstTimeStep[idleaf] == 4) #make sure that it is a leaf
stopifnot(ots_firstTimeStep[idroot] == 2)#make sure it is a root
##

#get Values
cst = get_Val_OT(baseline, "C_ST.txt","cst", c(idroot, idstem, idleaf), c("idroot", "idstem", "idleaf"))
psiX = get_Val_OT(baseline, "psiXyl.txt","psiX", c(idroot, idstem, idleaf), c("idroot", "idstem", "idleaf"))

cstPsiX = rbind(cst,psiX)
cstPsiX$ot = as.factor(cstPsiX$ot)
cstPsiX$variableName = factor(cstPsiX$variableName, levels=c("psiX","cst"))

##
var.lab3 = c("psiX" = bquote(psi["t,x"]~" (hPa)"), 
             "cst" = bquote(s["st"]~" (mmol/ml)"))

vlabeller <- function (variable, value) {
  return(var.lab3[value])
}


plotpsixi =ggplot(data=cstPsiX, 
                  aes(x=time, y=value,groupe=variable,
                      linetype = variable, alpha=variable, size=variable#,col=variable
                      ))+
  geom_line()+
  facet_wrap(~variableName,scale ="free_y",ncol=1,
             labeller = vlabeller, strip.position="left")+
  xlab('day of growth (d)')+
  ylab(NULL)+
  scale_x_continuous(breaks=seq(from =6,to=28,by=2))+
  scale_size_manual(breaks=c("idleaf","idroot","idstem"),
                   labels=c(paste0(as.character(idleaf -1)," (leaf)"),
                            paste0(as.character(idroot -1)," (root)"),
                            paste0(as.character(idstem -1)," (seed)")),
                    name="",
                    values = c(0.5,0.5,1.5))+ 
  scale_alpha_manual(breaks=c("idleaf","idroot","idstem"),
                     labels=c(paste0(as.character(idleaf -1)," (leaf)"),
                              paste0(as.character(idroot -1)," (root)"),
                              paste0(as.character(idstem -1)," (seed)")),
                     name="",
                     values = c(1,1,0.3))+ 
  scale_linetype_manual(breaks=c("idleaf","idroot","idstem"),name="",
                        labels=c(paste0(as.character(idleaf -1)," (leaf)"),
                                 paste0(as.character(idroot -1)," (root)"),
                                 paste0(as.character(idstem -1)," (seed)")),
                        values = c("dashed","solid","solid"))+
  #scale_color_manual(breaks=c("idstem","idleaf","idroot"),name="node Id:",
  #                   labels=c(paste0(as.character(idstem -1)," (seed)"),
  #                            paste0(as.character(idleaf -1)," (leaf)"),
  #                            paste0(as.character(idroot -1)," (root)")),
  #                   values = c("#FC8D62","#66C2A5","#8DA0CB"))+
  theme(legend.text = element_text( size=15),
        legend.title = element_text( size=15),
        legend.position=c(0.6,0.6),
        legend.direction = "horizontal",
        panel.background = element_rect(fill="white", color="black"),
        #legend.box = "vertical",
        panel.grid.minor = element_line(colour="white"),
        panel.grid.major = element_line(colour="white"),
        strip.text.y =element_text( size=21) ,
        axis.text = element_text( size=17),
        axis.title = element_text( size=21),
        axis.line = element_line(colour = "black"),
        strip.background = element_blank(),
        strip.placement = "outside");plotpsixi

ggsave('resultIntro_B',plot=plotpsixi,  width = 8, height = 10, dpi = 700)
