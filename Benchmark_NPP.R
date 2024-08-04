
##
## Measurement of the site
##
## Citation                   |  Le Goff et al. 2004
## Year                       |  1995
## Start of growing season    |  118, 122 and 129 for trees 6, 10 and 12 respectively
## (i.e. DOY of bud burst)    |
##                            |
## End of growing season      |  250
## LngLat                     |  48.7333333333 N, 6.23333333333 E; altitude 250 m
## Site                       |  state forest of Amance, northeast of France
## AN                         |  20.59 kg C/yr



#Module-Weather--------------------------
DOY <- function(year, month, day){
  floor(275*month/9)-floor((month+9)/12)*(1+floor((year-4*floor(year/4)+2)/3))+day-30
}

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

actualvp <- function(Tdew){ #estimated from dew temperature
  0.6118*exp(17.502*Tdew/(Tdew+240.97))  #reverse of Tetens equation
}

vpdeficit <- function(Ta, Tdew, ...){
  res <- satvp(Ta, ...)-actualvp(Tdew)
  ifelse(res>0, res, 0)
}

hourly_radiation <- function(Rtot, doy, hours, lng, lat){
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



#Module-Photosynthesis------------------------
farquhar <- function(Ci, aPAR, Ta, Vcmax_25 = 19.40675, Jmax_25 = 98.38106, Rleaf_25 = 0.2071948,
                     alpha = 0.15, a = 0.4, f_aS = 0.75){
  R <- 8.3144621                     ## ideal gas constant in J/K/mol
  Gstar <- 42.75*exp(37830*(Ta+273.15 - 298)/(298*R*(Ta+273.15)))
  ## co2 compenstaion point 
  
  Tleaf <- 6.990291262135919 + (35.24271844660194-6.990291262135919)/40*Ta +273.15    
  ## leaf temperature (K) [298 in original code]
  # Konrad et al. 2020 https://onlinelibrary.wiley.com/doi/full/10.1002/gj.3757
  # empirical relationship used for simplicity
  # Figure 5. leaf length = 100 mm, air moisture = 0.6
  
  Kc <- 404.9*exp(79430*(Tleaf - 298)/(298*R*Tleaf))
  ## Michaelis-Menten coefficients for carboxylation
  Ko <- 278.4*exp(36380*(Tleaf - 298)/(298*R*Tleaf))
  ## Michaelis-Menten coefficients for oxidation 
  Oi <- 210                          ## partial pressure of oxygen (mol mol−1)
  Km <- Kc*(1-Oi/Ko)
  ##Rleaf_25 <- 0.04 * Vcmax_25      ## leaf respiration (μmol/m2/sec) in original code
  
  Rleaf <- arrhenius(Rleaf_25, Ta)   ## temperature scaling
  Vcmax <- arrhenius(Vcmax_25, Ta)
  Jmax <- arrhenius(Jmax_25, Ta)
  
  b <- -(alpha*aPAR+Jmax)
  c <- alpha*aPAR*Jmax
  J <- (-b-sqrt(b^2-4*a*c))/(2*a)
  aJ <- J*(Ci-Gstar)/(4*Ci+8*Gstar)  ## electron transport limited
  aC <- Vcmax*(Ci-Gstar)/(Ci+Km)     ## co2 and rubisco limited
  aS <- Vcmax*f_aS                   ## starch accumulation at the chloroplast limited
  min(aJ,aC,aS) - Rleaf
}

ballberry <- function(init_Angs, aPAR,Ta,VPD,m=6.153588,g0=0.04358137){
  An_pred <- init_Angs[[1]]
  gs_pred <- init_Angs[[2]]
  
  R <- 8.3144621                     ## ideal gas constant in J/K/mol
  Ca <- 400                          ## partial pressure of ambient co2 (mol mol-1)
  Gstar <- 42.75*exp(37830*(Ta+273.15 - 298)/(298*R*(Ta+273.15)))
  ## co2 compenstaion point 
  Ci <- Ca - 1.6*An_pred/gs_pred     ## Medlyn et al 2011 model
  e1 <- (farquhar(Ci,aPAR,Ta) - An_pred)
  ## calculate the prediction error of net carbon assimilation
  e2 <- (g0 + m*An_pred/((Ca-Gstar)*(1+VPD)) - gs_pred)*100
  ## calculate the prediction error of stomatal conductance
  return(e1^2 + e2^2)
}

#Module-Respiration----------------------

arrhenius <- function(observed.value, new.temp, old.temp = 25){
  return(observed.value / exp (3000 * ( 1 / (273.15 + new.temp) - 1 / (273.15 + old.temp))))
}

leaf_respiration <- function(Ta, Rleaf_25 = 0.2071948){
  arrhenius(Rleaf_25, Ta)
  #Caution: Rleaf_25 umol CO2/ Leaf m2/s
  #   NOT the Ground area
}

stem_respiration <- function(Ta, Rstem_25 = 1.032241){
  arrhenius(Rstem_25, Ta)
}

root_respiration <- function(Ta, Rroot_25 = 2.174197){
  arrhenius(Rroot_25, Ta)
}

soil_respiration <- function(Ta, Rhete_25 = 1.449465){
  arrhenius(Rhete_25, Ta)
}

#Module-Allocation-----------------------
allocation <- function(C_B){
  c(
    dC_foliage = C_B * 0.2012307,
    dC_branch = C_B * 0.1976957,
    dC_stem = C_B * 0.4935847,
    dC_root = C_B * 0.1074889
  )
}


#Module-Process--------------------------
#input data
get_weather_data <- function(file_path){ #"Amance-climate-ECMWF.csv"
  weather <- read.csv(file_path)
  weather$surface_net_solar_radiation_sum <- as.numeric(gsub(",","",weather$surface_net_solar_radiation_sum))
  weather$year <- as.numeric(lapply(weather$Time, function(xx) strsplit(xx, ",")[[1]][2]))
  weather$month_code <- as.character(lapply(lapply(weather$Time, function(xx) strsplit(xx, ",")[[1]][1]), function(xx) strsplit(xx, " ")[[1]][1]))
  weather$day <- as.numeric(lapply(lapply(weather$Time, function(xx) strsplit(xx, ",")[[1]][1]), function(xx) strsplit(xx, " ")[[1]][2]))
  mm <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
  weather$month <- as.numeric(lapply(weather$month_code, function(xx) which(xx==mm)))
  weather$doy <- with(weather, DOY(year, month, day))
  weather
  # weather$temperature_2m_min (Kelvin)
  # weather$temperature_2m_max (Kelvin)
  # weather$surface_net_solar_radiation_sum (J/m2)
  # weather$dewpoint_temperature_2m (Kelvin)
  # time step: daily
}

#get diurnal temperature, solar radiation, and VPD
model_diurnal_data <- function(daily_weather, lng, lat){
  #'lng', longitude (degree)
  #'lat', latitude (degree)
  #'al', altitude (m)
  if(!is.data.frame(daily_weather)) stop("'daily_weather' should be a data.frame")
  if(nrow(daily_weather)!=1) stop("'daily_weather' should be a data.frame of a single row")
  
  hours <- seq(0, 23.5, 0.5) 
  hourly_temperature_data <- hourlyT(
    hours, 
    t_max = daily_weather$temperature_2m_max - 273.15,
    t_min = daily_weather$temperature_2m_min - 273.15,
    time_sr = sunRS(daily_weather$doy, lng, lat)[[1]], 
    time_ss = sunRS(daily_weather$doy, lng, lat)[[2]],
    alpha = 1.5, beta = 4, gamma = 1
  )
  hourly_VPD_data <- vpdeficit(
    hourly_temperature_data, 
    Tdew = daily_weather$dewpoint_temperature_2m - 273.15
  )
  hourly_radiation_data <- hourly_radiation(daily_weather$surface_net_solar_radiation_sum, daily_weather$doy, hours, lng, lat)
  data.frame(
    hours = hours,
    Ta = hourly_temperature_data,
    Rns = hourly_radiation_data,
    VPD = hourly_VPD_data
  )
}

#get net CO2 assimilation
model_photosynthesis <- function(hourly_weather){
  res <- data.frame(hours = hourly_weather$hours, An = rep(0,nrow(hourly_weather)), gs = rep(0,nrow(hourly_weather)))
  for(j in 1:nrow(res)){
    res[j,c(2,3)]<-optim(c(An=15, gs=0.1), ballberry, aPAR = hourly_weather[j,]$Rns, Ta = hourly_weather[j,]$Ta, VPD = hourly_weather[j,]$VPD)$par
  }
  res
}

#Daily NPP
calculate_daily_NPP <- function(hourly_Angs){
  #An umol/m2/sec
  # * 30 min
  carbon_per_LA <- sum(hourly_Angs$An * 30 * 60 * 1e-6 * 12) #g C per leaf area per day
  carbon_per_LA * 5.67 * 1000 * 97.6/ 10000   #total g C assimilation the whole day
    #where 5.67 kg the foliage biomass
    #      1000 convert kg to g
    #      97.6 specific leaf area (SLA) cm2/g
}

#Daily leaf respiration
calculate_daily_leaf_respiration <- function(hourly_weather){
  #Rleaf umol/m2/sec
  
  respiration_per_LA <- sum(-leaf_respiration(hourly_weather$Ta) * 30 * 60 * 1e-6 * 12)
  respiration_per_LA * 5.67 * 1000 * 97.6/ 10000
}

#Daily stem respiration
calculate_daily_stem_respiration <- function(hourly_weather){
  
  respiration_per_GA <- sum(-stem_respiration(hourly_weather$Ta) * 30 * 60 * 1e-6 * 12)
  respiration_per_GA * 11.1
  #stem area = 2*pi*18.2*1e-2*17.1
  #ground area = crown projection area = 11.1
}

#Daily root respiration
calculate_daily_root_respiration <- function(hourly_weather){
  
  respiration_per_GA <- sum(-root_respiration(hourly_weather$Ta) * 30 * 60 * 1e-6 * 12)
  respiration_per_GA * 11.1
}

#Daily hete soil respiration
calculate_daily_hete_respiration <- function(hourly_weather){
  
  respiration_per_GA <- sum(-soil_respiration(hourly_weather$Ta) * 30 * 60 * 1e-6 * 12)
  respiration_per_GA * 11.1
}


#Simulation-----------------------------
weather <- get_weather_data("Amance-climate-ECMWF.csv")
weather$surface_net_solar_radiation_sum <- weather$surface_net_solar_radiation_sum
  #where K = 1.7 is a factor that converts net solar radiation to photosynthetically active radiation (PAR)
  
AN <- 0 #yearly NPP (g C)
R_LEAF <- 0 #yearly leaf respiration (g C)
R_STEM <- 0 #yearly stem respiration (g C)
R_ROOT <- 0 #yearly root respiration (g C)
R_SOIL <- 0 #yearly heterotrophic soil respiration (g C)
for(i in 1:nrow(weather)){
  hourly_weather <- model_diurnal_data(weather[i,], lng = 6.23333333333, lat = 48.7333333333)
  
  if(weather[i,]$doy>=118 && weather[i,]$doy<=250){ 
    #within the growing season: C assimilation = photosynthesis - leaf respiration
    hourly_Angs <- model_photosynthesis(hourly_weather)
    
    AN <- AN + calculate_daily_NPP(hourly_Angs)
    R_LEAF <- R_LEAF + calculate_daily_leaf_respiration(hourly_weather)
  }
  
  R_STEM <- R_STEM + calculate_daily_stem_respiration(hourly_weather)
  R_ROOT <- R_ROOT + calculate_daily_root_respiration(hourly_weather)
  R_SOIL <- R_SOIL + calculate_daily_hete_respiration(hourly_weather)
}

#Leaf area:              Tree 6 = 55.3 m2, Tree 10 = 40.5 m2, Tree 12 = 15.3 m2
#Canopy projection area: Tree 6 = 11.1 m2, Tree 10 = 5.63 m2, Tree 12 = 2.69 m2

AN/1000 #kg C/yr         
  #Tree 6:  Predicted 18.54855 kg C/yr     Reality 20.59 kg C/yr
  #Tree 10: Predicted 13.58438 kg C/yr     Reality 11.01 kg C/yr (scaled by LA)
  #Tree 12: Predicted 5.131877 kg C/yr     Reality 5.43 kg C/yr  (scaled by LA)
R_LEAF/1000  #kg C/yr    
  #Tree 6:  Predicted -1.191952  kg C/yr    Reality -1.16 kg C/yr
  #Tree 10: Predicted -0.6045666 kg C/yr    Reality -0.57 kg C/yr
  #Tree 12: Predicted -0.2888604 kg C/yr    Reality -0.27 kg C/yr
R_STEM/1000  #kg C/yr    
  #Tree 6:  Predicted -2.600934  kg C/yr    Reality -2.31 kg C/yr
  #Tree 10: Predicted -1.319212  kg C/yr    Reality -1.27 kg C/yr
  #Tree 12: Predicted -0.6303164 kg C/yr    Reality -0.62 kg C/yr
R_ROOT/1000  #kg C/yr    
  #Tree 6:  Predicted -5.478317 kg C/yr    Reality -4.548 kg C/yr
  #Tree 10: Predicted -2.778642 kg C/yr    Reality -2.412 kg C/yr
  #Tree 12: Predicted -1.327628 kg C/yr    Reality -1.17  kg C/yr
R_SOIL/1000  #kg C/yr
  #Tree 6:  Predicted -3.652212 kg C/yr    Reality -3.032 kg C/yr
  #Tree 10: Predicted -1.852428 kg C/yr    Reality -1.608 kg C/yr
  #Tree 12: Predicted -0.8850856 kg C/yr    Reality -0.78  kg C/yr


C_B <- (AN + R_STEM + R_ROOT)/1000   #Carbon Balance
  #Tree 6:  Predicted 15.05501 kg C/yr     Reality 13.73 kg C/yr
  #Tree 10: Predicted 9.486526 kg C/yr     Reality 7.33  kg C/yr
  #Tree 12: Predicted 3.173933 kg C/yr     Reality 3.64  kg C/yr

allocation(C_B)
#Tree 6:
#Predicted (kg C/yr)                               #Reality (kg C/yr)
#dC_foliage  dC_branch    dC_stem    dC_root       #dC_foliage  dC_branch    dC_stem    dC_root 
#  3.029530   2.976310   7.430921   1.618246       #    3.06     3.89          7.22      2.50
#
#Tree 10:
#Predicted (kg C/yr)                               #Reality (kg C/yr)
#dC_foliage  dC_branch    dC_stem    dC_root       #dC_foliage  dC_branch    dC_stem    dC_root 
#  1.908980   1.875445   4.682404   1.019696       #    1.19     1.12          2.22      0.59
#
#Tree 12:
#Predicted (kg C/yr)                               #Reality (kg C/yr)
#dC_foliage  dC_branch    dC_stem    dC_root       #dC_foliage  dC_branch    dC_stem    dC_root 
#  0.6386928  0.6274729  1.5666048  0.3411626      #    0.49     0.36          1.46      0.13
#

#validation plot
res <- read.csv("benchmark_result.csv")
res <- aggregate(x= res$value,by = list(category=res$category, variable=res$variable),FUN = mean)
colnames(res)[colnames(res)=="x"] <- "value"
fig3A_dat <- res[!grepl("dC", res$variable),]
fig3B_dat <- res[grepl("dC", res$variable),]

plot_grid(
  ggplot(fig3A_dat, aes(x=variable, y=value, color=category, fill=category))+
    geom_bar(stat="identity", position=position_dodge(), alpha=0.5, size=1)+theme1+
    scale_fill_manual(values = c(model="blue",observed="red"))+
    scale_color_manual(values = c(model="blue",observed="red"))+
    theme(legend.title=element_blank())+xlab(NULL),
  ggplot(fig3B_dat, aes(x=variable, y=value, color=category, fill=category))+
    geom_bar(stat="identity", position=position_dodge(),alpha=0.5, size=1, width=0.65)+theme1+
    scale_fill_manual(values = c(model="blue",observed="red"))+
    scale_color_manual(values = c(model="blue",observed="red"))+
    theme(legend.title=element_blank())+xlab(NULL),
  labels = "AUTO", label_size = 24, ncol = 1
)



