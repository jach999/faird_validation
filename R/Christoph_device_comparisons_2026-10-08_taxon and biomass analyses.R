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

#install.packages(c("car","effects","DHARMa","mixcat","gllvm","sjPlot"),dependencies=T)
#install.packages("mclogit")

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

#install.packages("multpois")
library(multpois) # for multinomial model with random effects

etraps_presence_daily
ammod_presence_daily


etraps.ammod.pd=rbind.fill(etraps_presence_daily,ammod_presence_daily)
etraps.ammod.pd$Order=factor(etraps.ammod.pd$Order)
levels(etraps.ammod.pd$Order)[levels(etraps.ammod.pd$Order)=="#N/C"]="Indet."

head(etraps.ammod.pd)
etraps.ammod.pd$Date=as.Date(etraps.ammod.pd$Date)
etraps.ammod.pd$Date.num=as.numeric(etraps.ammod.pd$Date)


#use only top taxa:

toporders=names(rev(sort(table(etraps.ammod.pd$Order)))[c(1:5,7)])
toporders

subdata=subset(etraps.ammod.pd,Order%in%toporders)
subdata$Order=droplevels(subdata$Order)

mmod1=multinom(Order~bs(Date.num,3)*Device_type,data=subdata,MaxNWts=50000,maxit=100000)
Anova(mmod1)

head(subdata)

eff.date.1 <- Effect("Date.num", mmod1)

getwd()
pdf("Device effects 2026-10-08-1.pdf")

plot(allEffects(mmod1,xlevels=100),style="stacked",
     lattice=list(key.args=list(columns=3)),
      axes=list(x=list(Date.num=list(lab="Date",
     ticks=list(at=levels2dates(eff.date.1, "Date.num", "1970-01-01"))), 
    rotate=45)
     )
)




################################################################################
################################################################################
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

#xyplot(biomass_mg~Date,groups=Device_type, 
#       data=etraps.ammod.bio,type=c("p","smooth"),
#       auto.key=list(columns=2),jitter.x=T,
#       par.settings=simpleTheme(pch=16,cex=1.5,lwd=2))


library(fitdistrplus)

f1=fitdist(etraps.ammod.bio$biomass_mg,"norm")
f2=fitdist(etraps.ammod.bio$biomass_mg,"gamma")

library(tweedie)

t1=tweedie_profile(biomass_mg~1,data=etraps.ammod.bio)
#plot(t1)

# tweedie profile gives about 1.6. Use tweedie in glmmTMB:

range(etraps.ammod.bio$biomass_mg)

#plot(biomass_mg~Date,etraps.ammod.bio)
#denscomp(list(f1,f2))

detach(package:statmod)
detach(package:mixcat)
detach(package:tweedie)


detach(package:glmmTMB)

update.packages(c("TMB", "glmmTMB"))
library(glmmTMB)

B1=glmmTMB(biomass_mg~bs(Date.num,3)*Device_type
           +(1 | Site/Device),
           family=tweedie(),data=etraps.ammod.bio)


library(DHARMa)

r1=simulateResiduals(B1)
#plot(r1)

# exactly! Looks great!

range(etraps.ammod.bio$Date.num)

#plot(allEffects(B1,xlevels=list(Date.num=19592:19603)),multiline=T,rescale.axis=F,
#     axes=list(x=list(Date.num=list(ticks=list(at=levels2dates(Effect("Date.num",B1),"Date.num","1970-01-01",n=4))))
#)
#)               
               

L1=levels2dates(Effect("Date.num", B1,xlevels=list(Date.num=19592:19603)), "Date.num", "1970-01-01",n=4)

a1=Effect(c("Device_type","Date.num"),B1,xlevels=list(Date.num=19592:19603))

# below code works, but we just keep it for reference
# plot(a1,
#     lattice=list(key.args=list(columns=2)),
#     lines=list(col=brewer.pal(n=3,name="Set1")),
#     axes=list(x=list(Date.num=list(ticks=list(at=levels2dates(Effect("Date.num",B1),"Date.num","1970-01-01",n=4))))),
#     multiline=T)




plot(a1,
     lattice=list(key.args=list(columns=3)),
     lines=list(col=brewer.pal(n=3,name="Set1"),lwd=2),
     axes=list(x=list(Date.num=list(ticks=list(at=L1))),
               y=list(lab="Biomass (mg)",type="response")),
     multiline=T,ci.style="auto",main="")


dev.off()












