data("airquality", package="datasets")
airquality$Date <- with(airquality, as.Date(paste("1973", Month, Day, sep="-"),
                                            format="%Y-%m-%d"))
airquality$Date.num <- as.numeric(airquality$Date)



m2.date <- lm(Ozone ~ Date.num*Temp + Solar.R + Wind, data=airquality)
eff.date.2 <- Effect(c("Date.num", "Temp"), m2.date, xlevels=6)
plot(eff.date.2, axes=list(x=list(Date.num=list(lab="Date", 
                                                ticks=list(at=levels2dates(eff.date.2, "Date.num", "1970-01-01", n=3))), 
                                  rotate=45)), main="Date Effect by Temperature")
