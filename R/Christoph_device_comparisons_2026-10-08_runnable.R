# ==============================================================================
# DEVICE COMPARISONS: TAXON COMPOSITION AND BIOMASS
# Based on: Christoph_device_comparisons_2026-10-08_taxon_and_biomass_analyses.R
# ==============================================================================
# Runnable version for source(). Analyses and model specifications are
# unchanged. Every edit is tagged [RUNNABLE EDIT]; original lines that were
# disabled are kept as comments directly next to their replacement.
# ==============================================================================


# ==============================================================================
# CONSOLE LOG  [RUNNABLE EDIT: new section]
# ==============================================================================
# Close any sink left open by a previous interrupted run, then log console
# output to a timestamped file (split = TRUE keeps it visible in the console).
# Note: warnings and messages go to stderr and are not captured by sink().

while (sink.number() > 0) sink()

if (!dir.exists("outputs")) dir.create("outputs")
output_file <- file.path(
  "outputs",
  paste0("Christoph_device_comparisons_output_",
         format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt")
)
sink(output_file, split = TRUE)

cat("Run started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Working directory:", getwd(), "\n\n")


# ==============================================================================
# 0. DATA INPUT  [RUNNABLE EDIT: new section]
# ==============================================================================
# The analysis uses global objects created by R/00_load_data.R. If they are
# not already in the session (e.g. after load("Christoph_Device_Comparisons.RData")),
# 00_load_data.R is sourced to rebuild them from source_tables/.

required_objects <- c("etraps_presence_daily", "ammod_presence_daily",
                      "etraps_daily_overlap",  "ammod_biomass_overlap")

if (!all(vapply(required_objects, exists, logical(1)))) {
  cat("Input objects not found in session: sourcing R/00_load_data.R\n")
  source(here::here("R", "00_load_data.R"))
} else {
  cat("Input objects found in session: using existing objects\n")
}

# [RUNNABLE EDIT] Console output pasted into the original script, kept for reference:
# names(data_list)
#
#  [1] "faird"                   "id"
#  [3] "etraps"                  "etraps_daily"
#  [5] "etraps_daily_overlap"    "ammod_taxonomy"
#  [7] "ammod_biomass_for_stats" "ammod_biomass_original"
#  [9] "ammod_biomass"           "ammod_biomass_overlap"
# [11] "etraps_presence"         "ammod_presence"
# [13] "effort_base"             "effort_summary"
# [15] "overlap_dates"           "experiment_dates"

head(etraps_daily_overlap)
head(ammod_biomass_overlap)


# ==============================================================================
# 1. PACKAGES
# ==============================================================================

#install.packages(c("car","effects","DHARMa","mixcat","gllvm","sjPlot"),dependencies=T)
#install.packages("mclogit")

# [RUNNABLE EDIT] Check required packages before running anything
required_packages <- c("lattice", "plyr", "glmmTMB", "nnet", "car", "effects",
                       "splines", "DHARMa", "RColorBrewer")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
}

library(lattice)
library(plyr)
library(glmmTMB)
library(nnet)
library(car)
library(effects)
library(splines)
#library(VGAM)
#library(VGAMextra)
# [RUNNABLE EDIT] Packages below are not used by any analysis in this script.
# Disabled so that a missing installation does not stop source().
#library(mclogit)
#library(mixcat) #for npmlt function
#library(gllvm)
#library(MASS)
library(DHARMa)
#library(sjPlot)

#install.packages("multpois")
# [RUNNABLE EDIT] Not used (multinomial model with random effects did not run)
#library(multpois) # for multinomial model with random effects

library(RColorBrewer) # [RUNNABLE EDIT] needed for brewer.pal() in the biomass plot


# ==============================================================================
# 2. TAXON COMPOSITION: MULTINOMIAL MODEL
# ==============================================================================

etraps_presence_daily
ammod_presence_daily


etraps.ammod.pd=plyr::rbind.fill(etraps_presence_daily,ammod_presence_daily)
etraps.ammod.pd$Order=factor(etraps.ammod.pd$Order)
levels(etraps.ammod.pd$Order)[levels(etraps.ammod.pd$Order)=="#N/C"]="Indet."

head(etraps.ammod.pd)
etraps.ammod.pd$Date=as.Date(etraps.ammod.pd$Date)
etraps.ammod.pd$Date.num=as.numeric(etraps.ammod.pd$Date)


#use only top taxa:

toporders=names(rev(sort(table(etraps.ammod.pd$Order)))[c(1:5,7)])
print(toporders) # [RUNNABLE EDIT] print() so it shows under source()

subdata=subset(etraps.ammod.pd,Order%in%toporders)
subdata$Order=droplevels(subdata$Order)

mmod1=multinom(Order~bs(Date.num,3)*Device_type,data=subdata,MaxNWts=50000,maxit=100000)
print(Anova(mmod1)) # [RUNNABLE EDIT] print() so it shows under source()

head(subdata)

eff.date.1 <- Effect("Date.num", mmod1)

getwd()
pdf("Device effects 2026-10-08-2.pdf")

# [RUNNABLE EDIT] Wrapped in print(): effects plots are only drawn when printed,
# which does not happen automatically under source()
print(
plot(allEffects(mmod1,xlevels=100),style="stacked",
     lattice=list(key.args=list(columns=3)),
      axes=list(x=list(Date.num=list(lab="Date",
     ticks=list(at=levels2dates(eff.date.1, "Date.num", "1970-01-01"))),
    rotate=45)
     )
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

etraps.ammod.bio=plyr::rbind.fill(etraps_daily_overlap,ammod_biomass_overlap)
head(etraps.ammod.bio)

etraps.ammod.bio$Date=as.Date(etraps.ammod.bio$Date)
etraps.ammod.bio$Date.num=as.numeric(etraps.ammod.bio$Date)

# [RUNNABLE EDIT] print() so the checks show under source()
print(dim(etraps.ammod.bio)) #116 rows

print(length(unique(etraps.ammod.bio$Date))) # 12 dates
print(length(unique(etraps.ammod.bio$Site))) # 4 sites
print(length(unique(etraps.ammod.bio$Device_type))) # 3 sites

# 12*4*3# 144 theoretical rows so this looks fine

names(etraps.ammod.bio)

#xyplot(biomass_mg~Date,groups=Device_type,
#       data=etraps.ammod.bio,type=c("p","smooth"),
#       auto.key=list(columns=2),jitter.x=T,
#       par.settings=simpleTheme(pch=16,cex=1.5,lwd=2))


# ------------------------------------------------------------------------------
# Distribution checks (exploratory)
# ------------------------------------------------------------------------------
# [RUNNABLE EDIT] Kept, but guarded: each step runs only if its package is
# installed, and a failure prints a message instead of stopping source().
# Set to FALSE to skip (tweedie_profile can take several minutes).
RUN_DISTRIBUTION_CHECKS <- TRUE

if (RUN_DISTRIBUTION_CHECKS) {

  #library(fitdistrplus)
  if (requireNamespace("fitdistrplus", quietly = TRUE)) {
    f1 <- tryCatch(fitdistrplus::fitdist(etraps.ammod.bio$biomass_mg, "norm"),
                   error = function(e) { cat("fitdist norm failed:", conditionMessage(e), "\n"); NULL })
    f2 <- tryCatch(fitdistrplus::fitdist(etraps.ammod.bio$biomass_mg, "gamma"),
                   error = function(e) { cat("fitdist gamma failed:", conditionMessage(e), "\n"); NULL })
  } else {
    cat("Package fitdistrplus not installed: fitdist checks skipped\n")
  }

  #library(tweedie)
  if (requireNamespace("tweedie", quietly = TRUE)) {
    t1 <- tryCatch(tweedie::tweedie_profile(biomass_mg ~ 1, data = etraps.ammod.bio),
                   error = function(e) { cat("tweedie_profile failed:", conditionMessage(e), "\n"); NULL })
    if (!is.null(t1)) print(t1$p.max)
  } else {
    cat("Package tweedie not installed: Tweedie power profile skipped\n")
  }
}
#plot(t1)

# tweedie profile gives about 1.6. Use tweedie in glmmTMB:

print(range(etraps.ammod.bio$biomass_mg)) # [RUNNABLE EDIT] print()

#plot(biomass_mg~Date,etraps.ammod.bio)
#denscomp(list(f1,f2))

# [RUNNABLE EDIT] Session maintenance steps disabled: detach() fails when the
# package is not attached, and update.packages() is interactive (its first
# argument is lib.loc, not the package names). glmmTMB is already loaded above.
#detach(package:statmod)
#detach(package:mixcat)
#detach(package:tweedie)
#
#
#detach(package:glmmTMB)
#
#update.packages(c("TMB", "glmmTMB"))
#library(glmmTMB)


# ==============================================================================
# 3. BIOMASS: TWEEDIE GLMM
# ==============================================================================

B1=glmmTMB(biomass_mg~bs(Date.num,3)*Device_type
           +(1 | Site/Device),
           family=tweedie(),data=etraps.ammod.bio)

print(summary(B1)) # [RUNNABLE EDIT] added so model output shows under source()

library(DHARMa)

r1=simulateResiduals(B1)
#plot(r1)

# exactly! Looks great!

print(range(etraps.ammod.bio$Date.num)) # [RUNNABLE EDIT] print()

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




# [RUNNABLE EDIT] Wrapped in print() (see section 2)
print(
plot(a1,
     lattice=list(key.args=list(columns=3)),
     lines=list(col=brewer.pal(n=3,name="Set1"),lwd=2),
     axes=list(x=list(Date.num=list(ticks=list(at=L1))),
               y=list(lab="Biomass (mg)",type="response")),
     multiline=T,ci.style="auto",main="")
)


dev.off()

cat("\nDone. Figures written to:",
    file.path(getwd(), "Device effects 2026-10-08-2.pdf"), "\n")
cat("Run finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")


# ==============================================================================
# CLOSE LOG  [RUNNABLE EDIT]
# ==============================================================================
while (sink.number() > 0) sink()
cat(">> Console output saved to:", output_file, "\n")
