names(data_list)

 [1] "faird"                   "id"                     
 [3] "etraps"                  "etraps_daily"           
 [5] "etraps_daily_overlap"    "ammod_taxonomy"         
 [7] "ammod_biomass_for_stats" "ammod_biomass_original" 
 [9] "ammod_biomass"           "ammod_biomass_overlap"  
[11] "etraps_presence"         "ammod_presence"         
[13] "effort_base"             "effort_summary"         
[15] "overlap_dates"           "experiment_dates" 

library(lattice)
library(plyr)
etraps_presence_daily
ammod_presence_daily


etraps.ammod.pd=rbind.fill(etraps_presence_daily,ammod_presence_daily)
etraps.ammod.pd$Order=factor(etraps.ammod.pd$Order)
levels(etraps.ammod.pd$Order)[levels(etraps.ammod.pd$Order)=="#N/C"]="Indet."

head(etraps.ammod.pd)
etraps.ammod.pd$Date=as.Date(etraps.ammod.pd$Date)
etraps.ammod.pd$Date.num=as.numeric(etraps.ammod.pd$Date)
library(nnet)
library(car)
library(effects)
library(splines)

use only top taxa:

toporders=names(rev(sort(table(etraps.ammod.pd$Order)))[c(1:5,7)])
toporders

subdata=subset(etraps.ammod.pd,Order%in%toporders)
subdata$Order=droplevels(subdata$Order)

mod1=multinom(Order~bs(Date.num,3)*Device_type,data=subdata,MaxNWts=50000,maxit=100000)
Anova(mod1)

eff.date.1 <- Effect("Date.num", mod1)

getwd()
pdf("Device effects 2026-05-05-1.pdf")
plot(allEffects(mod1,xlevels=100),style="stacked",
     lattice=list(key.args=list(columns=3)),
      axes=list(x=list(Date.num=list(lab="Date",
     ticks=list(at=levels2dates(eff.date.1, "Date.num", "1970-01-01"))), 
    rotate=45)
     )
)

dev.off()


##
# old code from other multinomial scripts

library(RColorBrewer)
mycolors=colorRampPalette(brewer.pal(12,"Paired"))(20)

plot(allEffects(m1),
  axes=list(
    y=list(cex=1.5,
        style="stacked",
        lab=list(cex=2,label="Relative abundance")
      ),
    x=list(cex=1.5,
      location.ord=list(
      lab=list(cex=2,label="Location")
      ))
    ),
  lines=list(col=mycolors),main="",
  lattice=list(key.args=list(columns=3))
)

