require(ggplot2)
require(cowplot)

#UTILS-------------------------------------------------------
theme1 <- theme_bw()+
  theme(axis.text.x=element_text(size=16,angle=0,colour="black"),
        axis.text.y=element_text(size=16,angle=0,colour="black"),
        axis.title=element_text(size=18),
        axis.line=element_line(linetype=1,color="black",size=0.1),
        axis.ticks = element_line(colour="black"),
        panel.grid.major = element_blank(), #change the major and minor grid lines,
        panel.grid.minor = element_blank(), #if want to change, check this parameters, I think it's easier to dao that
        #strip.background = element_rect(colour = "black",size = 0.8),
        #panel.background = element_rect(colour="black", fill="white"),
        #panel.border = element_blank(),
        panel.border = element_rect(colour = "black",fill=NA,size = 1.2),
        plot.title=element_text(size=14,angle=0,colour="black", face = "italic"),
        plot.tag=element_text(size=14,angle=0,colour="black", face = "bold"),
        plot.caption=element_text(size=14,angle=0,colour="black",face = "italic"),
        axis.title.y=element_text(vjust=1.9),
        axis.title.x=element_text(vjust=0.5),
        legend.text=element_text(colour="black",size=14),
        legend.background= element_rect(fill = "transparent", color = NA),
        #legend.position = "bottom",
        legend.title = element_text(colour="black", size=14,angle=0))

#DATA------------------------------------------------------

# *************** Site Information **************
#
#   48° 44' N, 6° 14' E (48.7333333333 N, 6.23333333333 E); altitude 250 m
#       state forest of Amance, northeast of France
#
# ****** Get data from Google Earth Engine ****** 
#
# var start_year = 1995;
# var end_year = 1995;
# var climate = ee.ImageCollection("ECMWF/ERA5_LAND/DAILY_AGGR")
# .filterDate(ee.Number(start_year).format('%.0f').cat("-01-01"), ee.Number(end_year).format('%.0f').cat("-12-31"))
# .select(["temperature_2m_min",
#          "temperature_2m_max",
#          "surface_net_solar_radiation_sum",
#          "dewpoint_temperature_2m"])
# Map.addLayer(climate)
#
# ***************** Data File *******************
#
# ./Amance-climate-ECMWF.csv
#
# ***********************************************

#calculate day of the year
DOY <- function(year, month, day){
  floor(275*month/9)-floor((month+9)/12)*(1+floor((year-4*floor(year/4)+2)/3))+day-30
}

weather <- read.csv("Amance-climate-ECMWF.csv")
weather$surface_net_solar_radiation_sum <- as.numeric(gsub(",","",weather$surface_net_solar_radiation_sum))
weather$year <- as.numeric(lapply(weather$Time, function(xx) strsplit(xx, ",")[[1]][2]))
weather$month_code <- as.character(lapply(lapply(weather$Time, function(xx) strsplit(xx, ",")[[1]][1]), function(xx) strsplit(xx, " ")[[1]][1]))
weather$day <- as.numeric(lapply(lapply(weather$Time, function(xx) strsplit(xx, ",")[[1]][1]), function(xx) strsplit(xx, " ")[[1]][2]))
mm <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
weather$month <- as.numeric(lapply(weather$month_code, function(xx) which(xx==mm)))
weather$doy <- with(weather, DOY(year, month, day))

# weather$temperature_2m_min (Kelvin)
# weather$temperature_2m_max (Kelvin)
# weather$surface_net_solar_radiation_sum (J/m2)
# weather$dewpoint_temperature_2m (Kelvin)
# time step: daily
#
# temperature is measured at 2 metres above the surface of the Earth. Altitution 
#   correction is thus not required.


#TEMPERATURE-------------------------------------
#Sunrise and sunset time 
sunRS <- function(doy, long, lat){
  if(length(doy)!=1 || length(long)!=1 || length(lat)!=1) stop("All parameters should be scalars.")
  N<-doy
  #N -- day of the year
  t1<-N+(6-long/15)/24 #for sun rise
  t2<-N+(18-long/15)/24#/*for sun set*/
  #/*compute sunrise time*/
  rise <- function(t){
    result <- tryCatch({
      K<-pi/180.0
      M<-0.9856*t-3.289
      L<-M+1.916*sin(M*K)+0.02*sin(2*K*M)+282.634
      if(L >= 360.0 || L >= 360) L<-L-360
      tanRA<-0.91746*tan(L*pi/180)
      RA<-atan(tanRA)/K
      if((L > 90.0) && (L <= 270.0)) RA<-180.0+RA
      if ((L > 270.0) && (L <= 360.0)) RA<-360.0+RA
      RA<-RA/15.0
      sind<-0.39782*sin(L*K)
      cosd<-cos(asin(sind))
      x<-(-0.01454-sind*sin(lat*K))/(cosd*cos(lat*K))
      H<-acos(x)/K
      
      H<-360.0-H
      Tsr<-(H/15.0)+RA-(0.06571*t)-6.622
      
      if(Tsr < 0.0) Tsr=Tsr+24.0
      if(Tsr > 24.0) Tsr=Tsr-24.0
      return(Tsr)
    },error = function(e) print("You could probably try to force your returned arccosine angle not within [-1,1]"))
  }
  #/*compute sunset time*/
  sset <- function(t){
    result <- tryCatch({
      K<-pi/180.0
      M<-0.9856*t-3.289
      L<-M+1.916*sin(M*K)+0.02*sin(2*M*K)+282.634
      
      if(L>=360.0) L<-L-360
      tanRA<-0.91746*tan(L*K)
      RA<-atan(tanRA)/K
      if((L > 90.0) && (L <= 270.0)) RA<-180.0+RA
      if((L > 270.0) && (L <= 360.0)) RA<-360.0+RA
      RA<-RA/15.0
      
      sind<-0.39782*sin(L*K)
      cosd<-cos(asin(sind))
      x<-(-0.01454-sind*sin(lat*K))/(cosd*cos(lat*K))
      H<-acos(x)/K
      
      Tst<-H/15.0+RA-0.06571*t-6.622
      if(Tst < 0.0) Tst=Tst+24.0
      return(Tst)
    },error = function(e) print("You could probably try to force your returned arccosine angle not within [-1,1]")
    )
  }
  c(rise=rise(t1),set=sset(t2))
}

#Hourly air temperature
hourlyT <- function(hours, t_max, t_min, time_sr, time_ss,
                    alpha = 1.5, beta = 4, gamma = 1){
  m <- hours - time_sr - gamma
  dl <- time_ss-time_sr
  Z <- 24-dl
  
  t_sr <- (t_max-t_min)*sin(pi*(dl-gamma)/(dl+2*alpha))+t_min
  n <- ifelse(hours>time_ss,hours-time_ss,24-time_ss+hours)
  ifelse(
    hours<=time_ss&hours>(time_sr+gamma),
    (t_max-t_min)*sin(pi*m/(dl+2*alpha))+t_min,
    (t_sr-t_min)*exp(-beta*n/Z)+t_min
  )
}

#Modelling diurnal temperature
lng <- 6.23333333333 
lat <- 48.7333333333
hours <- seq(0, 23.5, 0.5)            # half-hour time resolution

doy_spring <- DOY(1995, 3, 20)        # spring equinox 1995
doy_summer <- DOY(1995, 6, 21)        # summer solstice 1995
doy_autumn <- DOY(1995, 9, 23)        # autumn equinox 1995
doy_winter <- DOY(1995, 12, 22)       # winter solstice 1995


Ta_spring <- hourlyT(                 # diurnal temperature of spring equinox
  hours, 
  t_max = weather[weather$doy==doy_spring,]$temperature_2m_max - 273.15,
  t_min = weather[weather$doy==doy_spring,]$temperature_2m_min - 273.15,
  time_sr = sunRS(doy_spring, lng, lat)[[1]], 
  time_ss = sunRS(doy_spring, lng, lat)[[2]],
  alpha = 1.5, beta = 4, gamma = 1
)

Ta_summer <- hourlyT(                 # diurnal temperature of summer solstice
  hours, 
  t_max = weather[weather$doy==doy_summer,]$temperature_2m_max - 273.15,
  t_min = weather[weather$doy==doy_summer,]$temperature_2m_min - 273.15,
  time_sr = sunRS(doy_summer, lng, lat)[[1]], 
  time_ss = sunRS(doy_summer, lng, lat)[[2]],
  alpha = 1.5, beta = 4, gamma = 1
)

Ta_autumn <- hourlyT(                 # diurnal temperature of autumn equinox
  hours, 
  t_max = weather[weather$doy==doy_autumn,]$temperature_2m_max - 273.15,
  t_min = weather[weather$doy==doy_autumn,]$temperature_2m_min - 273.15,
  time_sr = sunRS(doy_autumn, lng, lat)[[1]], 
  time_ss = sunRS(doy_autumn, lng, lat)[[2]],
  alpha = 1.5, beta = 4, gamma = 1
)

Ta_winter <- hourlyT(                 # diurnal temperature of winter solstice
  hours, 
  t_max = weather[weather$doy==doy_winter,]$temperature_2m_max - 273.15,
  t_min = weather[weather$doy==doy_winter,]$temperature_2m_min - 273.15,
  time_sr = sunRS(doy_winter, lng, lat)[[1]], 
  time_ss = sunRS(doy_winter, lng, lat)[[2]],
  alpha = 1.5, beta = 4, gamma = 1
)

#give outputs
ggplot()+
  geom_line(aes(x=hours, y=Ta_spring, color="spring equinox"), size=1)+
  geom_line(aes(x=hours, y=Ta_summer, color="summer solstice"), size=1)+
  geom_line(aes(x=hours, y=Ta_autumn, color="autumn equinox"), size=1)+
  geom_line(aes(x=hours, y=Ta_winter, color="winter solstice"), size=1)+
  scale_color_manual(values = c("spring equinox" = "green", "summer solstice" = "red", "autumn equinox" = "grey", "winter solstice" = "blue"))+
  xlab("Hour")+ylab(expression(T[a]))+
  theme1+theme(legend.title=element_blank())

#VPD----------------------------
#saturated vapor pressure
satvp <- function(Ta, method=c("Tetens","Lowry")){
  if(length(method)==2){
    if(!all(c("Tetens","Lowry")%in%method)) stop("invalid input `method`")
    method<-"Tetens"
  } else if(length(method)==1){
    if(!any(c("Tetens","Lowry")%in%method)) stop("invalid input `method`")
  } else{
    stop("invalid input `method`")
  }
  
  if(method=="Tetens") return(0.6118*exp(17.502*Ta/(Ta+240.97)))
  if(method=="Lowry") return(
    (6.1078+Ta*(0.44365185+Ta*(0.01428945+Ta*(0.00026506485+Ta*(3.0312404*1e-6+Ta*(2.034809*1e-8+Ta*6.1368209*1e-11))))))/10
  )
}

#actual vapor pressure
actualvp <- function(Tdew){ #estimated from dew temperature
  0.6118*exp(17.502*Tdew/(Tdew+240.97))  #reverse of Tetens equation
}

#vapor pressure deficit
vpdeficit <- function(Ta, Tdew, ...){
  res <- satvp(Ta, ...)-actualvp(Tdew)
  ifelse(res>0, res, 0)
}

#Modelling diurnal VPD
VPD_spring <- vpdeficit(
  Ta_spring, 
  Tdew = weather[weather$doy==doy_spring,]$dewpoint_temperature_2m - 273.15
)

VPD_summer <- vpdeficit(
  Ta_summer, 
  Tdew = weather[weather$doy==doy_summer,]$dewpoint_temperature_2m - 273.15
)

VPD_autumn <- vpdeficit(
  Ta_autumn, 
  Tdew = weather[weather$doy==doy_autumn,]$dewpoint_temperature_2m - 273.15
)

VPD_winter <- vpdeficit(
  Ta_winter, 
  Tdew = weather[weather$doy==doy_winter,]$dewpoint_temperature_2m - 273.15
)

#give outputs
ggplot()+
  geom_line(aes(x=hours, y=VPD_spring, color="spring equinox"), size=1)+
  geom_line(aes(x=hours, y=VPD_summer, color="summer solstice"), size=1)+
  geom_line(aes(x=hours, y=VPD_autumn, color="autumn equinox"), size=1)+
  geom_line(aes(x=hours, y=VPD_winter, color="winter solstice"), size=1)+
  scale_color_manual(values = c("spring equinox" = "green", "summer solstice" = "red", "autumn equinox" = "grey", "winter solstice" = "blue"))+
  xlab("Hour")+ylab("VPD")+
  theme1+theme(legend.title=element_blank())

#DAILY SOLAR----------------------------
#daily extraterrestrial radiation (MJ/m2/d)
ra <- function(lat, doy){
  lat <- lat*pi/180
  gsc <- 0.082 #solar constant (MJ m-2 min-1)
  dr <- 1 + 0.033*cos(2*pi*doy/365)
  delta <- 0.409*sin(2*pi*doy/365-1.39)
  omega <- acos(-tan(lat) * tan(delta))
  24*60/pi * gsc * dr *(omega*sin(lat)*sin(delta)+
                          cos(lat)*cos(delta)*sin(omega))
}

#daylight hours
dh <- function(lat, doy){
  lat <- lat*pi/180
  delta <- 0.409*sin(2*pi*doy/365-1.39)
  if(-tan(lat) * tan(delta)>=1) return(0)
  if(-tan(lat) * tan(delta)<=-1) return(24)
  omega <- acos(-tan(lat) * tan(delta))
  24/pi*omega
}

#incoming solar radiation (MJ/m2/d)
rs <- function(lat, doy, ah, as = 0.25, bs = 0.5){
  Ra <- ra(lat, doy)
  N <- dh(lat, doy)
  (as+bs*ah/N)*Ra
}

#net shortwave radiation (MJ/m2/d)
rns <- function(Rs, albedo = 0.5){ #original 0.17
  (1-albedo) * Rs
}

#model prediction
daily_solar <- as.numeric(lapply(weather$doy, function(doy){
  daylight_hours <- dh(lat, doy) #this assumes no cloudy hours (clear-day dominate)
                                 #  so the prediction must be higher than actual solar radiation
                                 #  and this discrepancy reflects how cloudy the site is
                                 #
                                 # apparently the model overestimates more in winter months as
                                 #  winter is typically cloudier than summer
  Rs <- rs(lat, DOY(2022, 5, 15), daylight_hours)
  rns(Rs)
})) #MJ/m2/d

#give outputs
#please note that the predicted solar radiation is at the top of canopy,
#  while the actual data is at the ground level (Earth's surface)
#so prediction must be higher than data because of the canopy interception
ggplot()+
  geom_line(aes(x=weather$doy, y=daily_solar, color="predict"), size=1)+
  geom_line(aes(x=weather$doy, y=weather$surface_net_solar_radiation_sum/10^6, color="ground"), size=1)+
  scale_color_manual(values = c("predict" = "blue", "actual" = "green"))+
  xlab("DOY")+ylab(expression("Net Shortwave Radiation (MJ/"*m^2*d^-1))+
  theme1+theme(legend.title=element_blank())

#HOURLY SOLAR----------------------------

#Hourly solar radiation
# Rtot:  daily total radiation (J/m2) DO NOT forget to convert MJ to 10^6 J
# hours: a vector of time from 0-24h
# lat: latitude
hourly_radiation <- function(Rtot, doy, hours, lat){
  C <- 0.4 #constant: Spitters et al., 1986
  delta <- -23.45*cos(360/365*(doy+10)*pi/180) #solar declination angle
  delta_rad <- pi*delta/180
  lat_rad <- pi*lat/180
  SD <- sin(lat_rad)*sin(delta_rad)
  CD <- cos(lat_rad)*cos(delta_rad)
  time_ss_sr <- sunRS(doy, lng, lat)
  DL <- diff(time_ss_sr)[[1]]
  LSH <- DL/2 + time_ss_sr[[1]] #time of maximum solar height (noon time)
  sin_beta <- SD + CD * cos(pi * (hours - LSH)/12)
  DSBE <- acos(-SD/CD) * 24/pi *(SD+C*SD^2+C*CD^2/2)+12*CD*(2+3*C*SD)*sqrt(1-(SD/CD)^2)/pi
  res <- ifelse(
    hours>time_ss_sr[[1]] & hours<time_ss_sr[[2]],
    Rtot * sin_beta * (1+C*sin_beta)/(DSBE * 3600),
    0
  )
  ifelse(res>0, res, 0)
}

spring_radiation <- hourly_radiation(daily_solar[doy_spring]*10^6, doy_spring, hours, lat)
summer_radiation <- hourly_radiation(daily_solar[doy_summer]*10^6, doy_summer, hours, lat)
autumn_radiation <- hourly_radiation(daily_solar[doy_autumn]*10^6, doy_autumn, hours, lat)
winter_radiation <- hourly_radiation(daily_solar[doy_winter]*10^6, doy_winter, hours, lat)

#output
ggplot()+
  geom_line(aes(x=hours, y=spring_radiation, color="spring"), size=1)+
  geom_line(aes(x=hours, y=summer_radiation, color="summer"), size=1)+
  geom_line(aes(x=hours, y=autumn_radiation, color="autumn"), size=1)+
  geom_line(aes(x=hours, y=winter_radiation, color="winter"), size=1)+
  scale_color_manual(values = c("spring" = "green", "summer" = "red", "autumn" = "grey", "winter" = "blue"))+
  xlab("Hour")+ylab(expression("Instantaneous Radiation (W/"*m^2*")"))+
  theme1+theme(legend.title=element_blank())


