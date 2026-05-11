names(data_list)

 [1] "faird"                   "id"                     
 [3] "etraps"                  "etraps_daily"           
 [5] "etraps_daily_overlap"    "ammod_taxonomy"         
 [7] "ammod_biomass_for_stats" "ammod_biomass_original" 
 [9] "ammod_biomass"           "ammod_biomass_overlap"  
[11] "etraps_presence"         "ammod_presence"         
[13] "effort_base"             "effort_summary"         
[15] "overlap_dates"           "experiment_dates" 

head(etraps_daily_overlap)
head(ammod_biomass_overlap)

library(lattice)
library(plyr)
library(glmmTMB)
library(nnet)
library(car)
library(effects)
library(splines)
#library(VGAM)
#library(VGAMextra)
library(mclogit)
library(mixcat) #for npmlt function
library(gllvm)
library(MASS)
library(DHARMa)
library(sjPlot)

etraps_presence_daily
ammod_presence_daily


etraps.ammod.pd=rbind.fill(etraps_presence_daily,ammod_presence_daily)
etraps.ammod.pd$Order=factor(etraps.ammod.pd$Order)
levels(etraps.ammod.pd$Order)[levels(etraps.ammod.pd$Order)=="#N/C"]="Indet."

head(etraps.ammod.pd)
etraps.ammod.pd$Date=as.Date(etraps.ammod.pd$Date)
etraps.ammod.pd$Date.num=as.numeric(etraps.ammod.pd$Date)


use only top taxa:

toporders=names(rev(sort(table(etraps.ammod.pd$Order)))[c(1:5,7)])
toporders

subdata=subset(etraps.ammod.pd,Order%in%toporders)
subdata$Order=droplevels(subdata$Order)

mmod1=multinom(Order~bs(Date.num,3)*Device_type,data=subdata,MaxNWts=50000,maxit=100000)
Anova(mmod1)

head(subdata)

eff.date.1 <- Effect("Date.num", mmod1)

getwd()
pdf("Device effects 2026-05-05-1.pdf")

plot(allEffects(mmod1,xlevels=100),style="stacked",
     lattice=list(key.args=list(columns=3)),
      axes=list(x=list(Date.num=list(lab="Date",
     ticks=list(at=levels2dates(eff.date.1, "Date.num", "1970-01-01"))), 
    rotate=45)
     )
)

# this will likely remain the main graph.
dev.off()

# now add random effect:####
# try brms:

with(subdata,table(Device,Site))

#library(brms)
#vignette("brms_overview")

#mod_brms <- brm(
#  Order ~ bs(Date.num, 3) * Device_type + (1 | Site/Device),
#  data = subdata,
#  family = categorical(),
#  chains = 4, cores = 4
#)

# all of this takes very long.

#plot(conditional_effects(mod_brms, c("Date.num", "Device_type")))
#sjPlot::plot_model(mod_brms, type = "pred", terms = c("Date.num", "Device_type"))

# VGAM doesn´t allow random effects.
# there is a rudimentary function in mgcv
# let´s try mclogit package:

# problem, even if only Site or Device are used as
# random effects, inner optimizers don´t converge
m1<- mblogit(Order ~ bs(Date.num,3)*Device_type, 
             random = ~1|Site, data = subdata)

# try npmlt function instead
# from mixcat package


m2<- 
  with(subdata,
  npmlt(as.numeric(Order) ~ bs(Date.num,3)*Device_type, 
             random = ~1,id=Site))

# doesn´t work like that.
# use multinomial family in gllvm?
# this would require a response matrix Y instead of "Order".
# rather,use the Kormann et al 2015 approach, with log(total) as offset:

subdata$Count <- 1 # for aggregation
names(subdata)

agg1=aggregate(Count~ Order+Site + Device_type+Device+Date.num+Date+Ambient, 
          data = subdata,FUN=sum)

# calculate summed abundances to use as log offset in glmmTMB:

agg1$total <- ave(agg1$Count, agg1$Site, agg1$Device, agg1$Date.num, FUN = sum)

mod1=glmmTMB(Count ~ bs(Date.num, 3) * Device_type * Order +
     (1 | Site) + (1 | Device/Order)+offset(log(total)),
  family = nbinom2,data = agg1)

summary(mod1)
VarCorr(mod1)

# near-zero random effects! So we can maybe stay with multinom?

mod2=glmmTMB(Count ~ bs(Date.num, 3) * Device_type * Order +
     (1 | Device)+offset(log(total)),
  family = nbinom2,data = agg1)

plot(Effect(c("Order","Date.num"),mod2,x.var="Date.num"),multiline=T,rescale.axis=FALSE)

# improved plot

plot(allEffects(mod2,x.var="Date.num",xlevels=100),style="stacked",
     rescale.axis=FALSE,
     lattice=list(key.args=list(columns=3)),
      axes=list(x=list(Date.num=list(lab="Date",
     ticks=list(at=levels2dates(eff.date.1, "Date.num", "1970-01-01"))), 
    rotate=45)
     )
)



###########################################

# try without offset:


mod3=glmmTMB(Count ~ bs(Date.num, 3) * Device_type * Order +
     (1 | Device),
  family = nbinom2,data = agg1)

s1=simulateResiduals(mod3)
plot(s1)

mod4=stepAIC(mod3)
Anova(mod4)

plot(allEffects(mod4),multiline=T)
#1 value in the Device_type*Order effect are not estimable

plot_model(mod4,type="int",terms="Date.num [all]")

#vignette("troubleshooting",package="glmmTMB")

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

####################
# now use etraps_daily_biomass

head(etraps_daily_overlap)
head(ammod_biomass_overlap)

ammod_biomass_overlap$Dry_mass_g
# convert dry mass g into biomass_mg

ammod_biomass_overlap$biomass_mg=ammod_biomass_overlap$Dry_mass_g*1000

etraps.ammod.bio=rbind.fill(etraps_daily_overlap,ammod_biomass_overlap)
head(etraps.ammod.bio)

etraps.ammod.bio$Date=as.Date(etraps.ammod.bio$Date)
etraps.ammod.bio$Date.num=as.numeric(etraps.ammod.bio$Date)

dim(etraps.ammod.bio) #116 rows

length(unique(etraps.ammod.bio$Date)) # 12 dates
length(unique(etraps.ammod.bio$Site)) # 4 sites
length(unique(etraps.ammod.bio$Device_type)) # 3 sites

# 12*4*3# 144 theoretical rows so this looks fine

names(etraps.ammod.bio)

xyplot(biomass_mg~Date,groups=Device_type, 
       data=etraps.ammod.bio,type=c("p","smooth"),
       auto.key=list(columns=2),jitter.x=T,
       par.settings=simpleTheme(pch=16,cex=1.5,lwd=2))






