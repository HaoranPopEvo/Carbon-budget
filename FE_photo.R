##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
##    Photosynthesis for Fraxinus excelsior
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Create at: 2024-07-06 15:58:20 CST
## Author: Haoran Wu 
## Affiliation: School of Geography and the Environment, University of Oxford
## Expertise: Forest Pathology; Tree Pests & Diseases; Ecosystem Modelling
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Description: Fit Farquhar photosynthesis model (couple with
##    Ball-Berry model) for European ash (Fraxinus excelsior L.).
##    Data obtained from in-situ measurements by Majewski et al. 2024
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
require(ggplot2)
require(cowplot)

#UTILS-------------------------------------------------------
#' Reset the value of a parameter for an R function.
#' @param fun an R function
#' @param reset_var variable to be reset. A string
#' @param reset_val value to be assigned. A numerical value
set_fun_parameter <- function(fun, reset_var, reset_val){
  fun_text <- deparse(fun)
  fun_head_index <- 1:(which(grepl("\\{", fun_text))-1)
  fun_head <- paste0(fun_text[fun_head_index], collapse = "")
  fun_body <- fun_text[-fun_head_index]
  current_text_index <- which(grepl(reset_var, fun_head))
  current_text <- fun_head[current_text_index]
  split_current_text <- strsplit(current_text, paste0(reset_var," ?="))[[1]]
  start_text <- split_current_text[1]
  end_text <- paste(strsplit(split_current_text[2],",")[[1]][-1], collapse = ",")
  replaced_text <- paste0(start_text, paste(reset_var, "=", reset_val), ",", end_text)
  if(substr(replaced_text, start = nchar(replaced_text), stop = nchar(replaced_text))=="," &&
     current_text_index == length(fun_head)){
    #in case the parameter ranks the last 
    replaced_text <- paste0(substr(replaced_text, start = 1, stop = nchar(replaced_text)-1),")")  
  }
  fun_head[current_text_index] <- replaced_text
  new_fun <- c(fun_head, fun_body)
  eval(parse(text = paste(new_fun, collapse = "\n"))) ##this is the returned function
}

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

#MODELS------------------------------------------------------
## Farquhar Model. 
## Modified from Dietze and Matthes et al. (2014)
farquhar <- function(Ci, aPAR, Ta, Vcmax_25 = 18, Jmax_25 = Vcmax_25*1.67, Rleaf_25 = 5,
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

## Ball-Berry Model, coupled with Farquhar Model
## Modified from Dietze and Matthes et al. (2014)
ballberry <- function(init_Angs, aPAR,Ta,VPD,m=4,g0=0.1){
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

## Arrhenius Model
arrhenius <- function(observed.value, new.temp, old.temp = 25){
  return(observed.value / exp (3000 * ( 1 / (273.15 + new.temp) - 1 / (273.15 + old.temp))))
}

#DATA--------------------------------------------------------
# Majewski et al. 2024. Diurnal variation curve of Pn, gs, PAR, and WUE
dat <- read.csv("D:/0_Files/0_DGVM/ADB-carbon-model/Majewski-2024-data.csv")
dat$VPD..kPa. <- dat$VPD..Pa./1000 #transform unit

#Assimilation------------------------------------------------
## Use random search method with possible parameter ranges defined

## Possible parameter range
Para_ranges <- list(
  Vcmax_25 = c(min = 18, max = 23),  #min max
  Jmax_25 = c(min = 95, max = 115),
  g0 = c(min = 0.04, max = 0.07),
  m = c(min = 6, max = 14),
  Rleaf_25 = c(min = 0.2071948, max = 0.2071948) #  1.032241 umol CO2/Ground m2/s and LAI being 5
                                                #  = 0.2071948 umol CO2/Leaf m2/s
)


parameter_already_calibrated <- TRUE
if(!parameter_already_calibrated){
  #Random search algorithm
  cat("Start Random Searching.\n")
  rs <- data.frame(Vcmax_25=c(),Jmax_25=c(),g0=c(),m=c(),Rleaf_25=c(),error=c())
  for(i in 1:100){ #explore 100 times
    #update model parameters
    Vcmax_25_current <- runif(1, min = Para_ranges$Vcmax_25[["min"]], max = Para_ranges$Vcmax_25[["max"]])
    Jmax_25_current <- runif(1, min = Para_ranges$Jmax_25[["min"]], max = Para_ranges$Jmax_25[["max"]])#Vcmax_25_current * 1.67
    Rleaf_25_current <- runif(1, min = Para_ranges$Rleaf_25[["min"]], max = Para_ranges$Rleaf_25[["max"]])
    g0_current <- runif(1, min = Para_ranges$g0[["min"]], max = Para_ranges$g0[["max"]])
    m_current <- runif(1, min = Para_ranges$m[["min"]], max = Para_ranges$m[["max"]])
    
    farquhar <- set_fun_parameter(farquhar, "Vcmax_25", Vcmax_25_current)
    farquhar <- set_fun_parameter(farquhar, "Jmax_25", Jmax_25_current)
    farquhar <- set_fun_parameter(farquhar, "Rleaf_25", Rleaf_25_current)
    ballberry <- set_fun_parameter(ballberry, "g0", g0_current)
    ballberry <- set_fun_parameter(ballberry, "m", m_current)
    
    #model prediction
    Model_pred <- data.frame(An = c(), gs = c())
    for(j in 1:nrow(dat)){
      PAR <- dat$PAR..umol.m2.s.[j]
      Ta <- dat$Ta..degree.C.[j]
      VPD <- dat$VPD..kPa.[j]
      res <- optim(c(An=15, gs=0.1), ballberry, aPAR = PAR, Ta = Ta, VPD = VPD)$par
      Model_pred <- rbind(Model_pred, data.frame(An = res[["An"]], gs = res[["gs"]]))
    } 
    error_rs <- sum((dat$Pn..umol.m2.s. - Model_pred$An)^2 + (dat$gs..mol.m2.s. - Model_pred$gs)^2)
    if(i==1){
      Final_Model_pred_rs <- Model_pred
      last_error_rs <- error_rs
    } else {
      if(error_rs < last_error_rs) Final_Model_pred_rs <- Model_pred
    }
    
    rs <- rbind(rs,data.frame(Vcmax_25=Vcmax_25_current,Jmax_25=Jmax_25_current,g0=g0_current,m=m_current,Rleaf_25=Rleaf_25_current,error=error_rs))
    cat("Random search attempt", i, "error =",error_rs,"\n")
  }
  cat("Done.\n")
  Final_Model_pred_rs$Ci <- with(Final_Model_pred_rs, 400-1.6*An/gs)
  optimal_Paras <- rs[which(min(rs$error)==rs$error),] #get optimal parameter
  
  #Error assessment - sensitivity to parameter
  plot(rs$Vcmax_25,rs$error)
  plot(rs$Jmax_25,rs$error)
  plot(rs$g0,rs$error)
  plot(rs$m,rs$error)
} else{
  #in case where an optimal model is already obtained
  #update model parameters
  optimal_Paras <- data.frame(
    Vcmax_25 = 19.40675, #19.82328, #old
    Jmax_25 = 98.38106, #96.75659, #old
    g0 = 0.04358137, #0.04358137, #old
    m = 6.153588, #6.707346, #old
    Rleaf_25 = 0.2071948, #1.032241, #this is actually fixed - old
    error = 1825.676 #1568.147 #old
  )
  
  farquhar <- set_fun_parameter(farquhar, "Vcmax_25", optimal_Paras$Vcmax_25)
  farquhar <- set_fun_parameter(farquhar, "Jmax_25", optimal_Paras$Jmax_25)
  farquhar <- set_fun_parameter(farquhar, "Rleaf_25", optimal_Paras$Rleaf_25)
  ballberry <- set_fun_parameter(ballberry, "g0", optimal_Paras$g0)
  ballberry <- set_fun_parameter(ballberry, "m", optimal_Paras$m)
  
  #model prediction
  Model_pred <- data.frame(An = c(), gs = c())
  for(j in 1:nrow(dat)){
    PAR <- dat$PAR..umol.m2.s.[j]
    Ta <- dat$Ta..degree.C.[j]
    VPD <- dat$VPD..kPa.[j]
    res <- optim(c(An=15, gs=0.1), ballberry, aPAR = PAR, Ta = Ta, VPD = VPD)$par
    Final_Model_pred_rs <- Model_pred <- rbind(Model_pred, data.frame(An = res[["An"]], gs = res[["gs"]]))
  } 
  Final_Model_pred_rs$Ci <- with(Final_Model_pred_rs, 400-1.6*An/gs)
}

##PERFORMANCE--------------------------------------------
#1. Compare predictions and observations
#Canvas 6.64 x 5.73
plot_grid(
  ggplot()+
    geom_point(aes(x=dat$Pn..umol.m2.s., y=Final_Model_pred_rs$An), size=3)+
    geom_abline(intercept = 0, slope = 1, size = 1.5, color = "blue", linetype = 2)+
    xlab(expression(Observed*" "*A[n]))+
    ylab(expression(Predicted*" "*A[n]))+theme1,
  ggplot()+
    geom_point(aes(x=dat$gs..mol.m2.s., y=Final_Model_pred_rs$gs), size=3)+
    geom_abline(intercept = 0, slope = 1, size = 1.5, color = "blue", linetype = 2)+
    xlab(expression(Observed*" "*g[s]))+
    ylab(expression(Predicted*" "*g[s]))+theme1,
  labels = "AUTO",
  label_size = 24,
  ncol = 1
)

#2. Compare predicted and observed patterns:
##  Light reaction module
aPAR <- seq(0,2000,10) #gradient
Jmax_25<- optimal_Paras$Jmax_25
alpha <- 0.15 #0.3
a <- 0.4 #0.7
Ta <- 25
Jmax <- arrhenius(Jmax_25, Ta)
b <- -(alpha*aPAR+Jmax)
c <- alpha*aPAR*Jmax
J <- (-b-sqrt(b^2-4*a*c))/(2*a)

Ci <- mean(Final_Model_pred_rs$Ci)
R <- 8.3144621 
Gstar <- 42.75*exp(37830*(Ta+273.15 - 298)/(298*R*(Ta+273.15)))
aJ <- J*(Ci-Gstar)/(4*Ci+8*Gstar) 
Rd <- optimal_Paras$Rleaf_25
Rleaf_aJ <- arrhenius(Rd, Ta)
pattern1 <- ggplot()+
  geom_point(aes(x=dat$PAR..umol.m2.s., y=dat$Pn..umol.m2.s.), size=3)+
  geom_line(aes(x=aPAR, y=aJ-Rleaf_aJ), size = 1.5, color = "blue")+
  xlab(expression(PAR))+
  ylab(expression(A[n]*"  or  "*A[J]))+theme1

##  Dark reaction module
Vcmax_25 <- optimal_Paras$Vcmax_25 
Ta_grad <- 13:26 #gradient
R <- 8.3144621 
Gstar <- 42.75*exp(37830*(Ta_grad+273.15 - 298)/(298*R*(Ta_grad+273.15)))
Vcmax <- arrhenius(Vcmax_25, Ta_grad)
Tleaf <- 6.990291262135919 + (35.24271844660194-6.990291262135919)/40*Ta_grad +273.15 
Kc <- 404.9*exp(79430*(Tleaf - 298)/(298*R*Tleaf))
Ko <- 278.4*exp(36380*(Tleaf - 298)/(298*R*Tleaf))
Oi <- 210 
Km <- Kc*(1-Oi/Ko)
Ci <- mean(Final_Model_pred_rs$Ci)
aC <- Vcmax*(Ci-Gstar)/(Ci+Km)
Rd <- optimal_Paras$Rleaf_25
Rleaf_aC <- arrhenius(Rd, Ta_grad)
pattern2 <- ggplot()+
  geom_point(aes(x=dat$Ta..degree.C., y=dat$Pn..umol.m2.s.), size=3)+
  geom_line(aes(x=Ta_grad, y=aC-Rleaf_aC), size = 1.5, color = "blue")+
  xlab(expression(Ta))+
  ylab(expression(A[n]*"  or  "*A[C]))+theme1


## Stomatal Conductance
Ta <- 25
Gstar <- 42.75*exp(37830*(Ta+273.15 - 298)/(298*8.3144621*(Ta+273.15)))
VPD <- seq(0.5,2,0.1)
gs_pred <- optimal_Paras$g0+optimal_Paras$m*5/((400-Gstar)*(1+seq(0.5,2,0.1)))
pattern3 <- ggplot()+
  geom_point(aes(x=dat$VPD..kPa., y=dat$gs..mol.m2.s.), size=3)+
  geom_line(aes(x=VPD, y=gs_pred), size = 1.5, color = "blue")+
  xlab(expression(VPD))+
  ylab(expression(g[s]))+theme1

#summarise three patterns
plot_grid(pattern1,pattern2,pattern3,labels = "AUTO",label_size = 24,nrow = 1)













#SIMULATION ANNEALING----------
##The following code is currently NOT USED
#ML method seems not reliable
Simu_paras <- c(
  start_amp = 0.1, # 10%
  step = 100,
  end_amp = 0.01   # 1%
)  
Model_paras <- c(
  Vcmax_25 = 18,
  Jmax_25 = 18*1.67,
  g0 = 0.01,
  m = 15
)

cat("Start Simulation Annuealing.\n")
error <- c()
for(i in 0:Simu_paras[["step"]]){
  #update model parameters
  if(i>0){
    current_amp <- (i-1)*(Simu_paras[["end_amp"]]-Simu_paras[["start_amp"]])/(Simu_paras[["step"]]-1)+Simu_paras[["start_amp"]]
    Model_paras_old <- Model_paras
    Model_paras <- (1 - sample(c(-1,1), length(Model_paras), replace = TRUE) * current_amp) * Model_paras
  }
  farquhar <- set_fun_parameter(farquhar, "Vcmax_25", Model_paras[["Vcmax_25"]])
  farquhar <- set_fun_parameter(farquhar, "Jmax_25", Model_paras[["Jmax_25"]])
  ballberry <- set_fun_parameter(ballberry, "g0", Model_paras[["g0"]])
  ballberry <- set_fun_parameter(ballberry, "m", Model_paras[["m"]])
  
  #model prediction
  Model_pred <- data.frame(An = c(), gs = c())
  for(j in 1:nrow(dat)){
    PAR <- dat$PAR..umol.m2.s.[j]
    Ta <- dat$Ta..degree.C.[j]
    VPD <- dat$VPD..kPa.[j]
    res <- optim(c(An=15, gs=0.1), ballberry, aPAR = PAR, Ta = Ta, VPD = VPD)$par
    Model_pred <- rbind(Model_pred, data.frame(An = res[["An"]], gs = res[["gs"]]))
  }
  
  #calculate bias
  error <- c(error, sum((dat$Pn..umol.m2.s. - Model_pred$An)^2 + (dat$gs..mol.m2.s. - Model_pred$gs)^2))
  if(length(error)>1){
    if(error[length(error)]>error[length(error)-1]){
      #current error > previous error
      Model_pred <- Model_paras_old
      error[length(error)] <- error[length(error)-1]
    } else{
      Final_Model_pred <- Model_pred
    }
  } else{
    Final_Model_pred <- Model_pred
  }
  
  #output message
  cat("Step",i,"of",Simu_paras[["step"]],". error =",error[length(error)],"\n")
}

