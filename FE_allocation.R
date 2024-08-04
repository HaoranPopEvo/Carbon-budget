require(ggplot2)
require(cowplot)

## list of definitions on respiration
##      C_B        =           A_G          +         R_f          +         R_w         +       R_r
## carbon balance     gross C assimilation     foliar respiration     woody respiration    root respiration                ----
##                                          |                                            |                                    |
##                                          |                                            |                                    |
##                                          ----------------------------------------------                                    |
##                                                                ||                                                          |
##                                                               \||/                                                         |==>           R_sr
##                                                                \/                           R_s (R_h)                      |    soil + root respiration
##                                                                R_A                      soil respiration                   |   (belowground respiration)
##                                                      aboveground respiration       (heterotrophic respiration)             |
##                                                                                In some papers, soil respiration is the     |                |
##                                                                 |              sum of root and heterotrophic respiration ----               |
##                                                                 |                                                                           |
##                                                                 -----------------------------------------------------------------------------
##                                                                                                        ||
##                                                                                                       \||/
##                                                                                                        \/
##                                                                                                        R_E
##                                                                                              ecosystem respiration
##
##
## growth respiration is not considered according to the carbon budget given by
##    Le Goff et al. 2004

#UTILS---------------------------------------------------------
require(ggplot2)
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

#BASELINE------------------------------
## Baseline for respiration in 25 degree C
## Model adopted from Le Goff et al. 2004

Ta_25 <- 25                            #air temperature (in degree C)
Tsoil_25 <- 6.818 + 0.450 * Ta_25      #soil temperature (in degree C)
R_E25 <- 0.531 * 10^(0.057*Tsoil_25)   #ecosystem respiration (umol CO2 m-2 s-1)
R_sr25 <- 0.436 * 10^(0.0509*Tsoil_25) #belowground respiration (umol CO2 m-2 s-1) = root + soil (heterotrophic)
R_r25 <- 0.6 * R_sr25                  #root respiration (umol CO2 m-2 s-1)
R_h25 <- R_sr25 - R_r25                #heterotrophic respiration (umol CO2 m-2 s-1)


R_A25 <- R_E25 - R_sr25
R_l25 <- R_w25 <- R_A25/2

#Convert leaf respiration from umol CO2 Ground m-2 to Leaf m-2
R_l25_LA <- R_l25/(55.3/11.1)
   #where 55.3 is the leaf area, and 11.1 is the crown projection area

#TEMPERATURE----------------------------
## Scaling to other temperatures for each respiration component
##    based on Arrhenius equation
Ta <- 1:35

arrhenius <- function(observed.value, new.temp, old.temp = 25){
  return(observed.value / exp (3000 * ( 1 / (273.15 + new.temp) - 1 / (273.15 + old.temp))))
}
Rleaf <- arrhenius(R_l25, Ta)
Rwood <- arrhenius(R_w25, Ta)
Rroot <- arrhenius(R_r25, Ta)
Rhete <- arrhenius(R_h25, Ta)

ggplot()+
  geom_line(aes(x=Ta, y=Rleaf, color="leaf"), size=1)+
  geom_line(aes(x=Ta, y=Rroot, color="root"), size=1)+
  geom_line(aes(x=Ta, y=Rhete, color="soil"), size=1)+
  scale_color_manual(values = c("leaf" = "green", "root" = "grey", "soil" = "red"))+
  xlab(expression(T[a]))+ylab("Respiration")+
  theme1+theme(legend.title=element_blank())

#ALLOMETRICS---------------------------
#Krejza et al. 2017 Trees

DBH <- 10:80
LB <- 0.983896 * exp(0.057012 * DBH)        #leaf biomass (LB) - eq3
BsB <- 0.074100 * exp(0.071603 * DBH)       #shoot biomass (BsB) - eq3
BB <- 0.005318 * DBH^2.897524               #branch biomass (BB) - eq1
                                            #in this paper, shoot is a part of branch but its bionass is evaluated separately
SB <- -0.069823 * DBH + 0.631338 * DBH^2    #stem biomass (SB) - eq4
RB <- 19.004922 * exp(0.052544 * DBH)       #root biomass (RB) - eq3
TAB <- 0.381699 * DBH^2.193037              #aboveground biomass (TAB) = SB+BB - eq1
TB <- 0.364819 * DBH^2.241458               #total biomass (TB) - eq1

ggplot()+  #plot leaf and shoot biomass
  geom_line(aes(x=DBH, y=LB, color="leaf"), size=1)+
  geom_line(aes(x=DBH, y=BsB, color="shoot"), size=1)+
  scale_color_manual(values = c("leaf" = "green", "shoot" = "grey"))+
  xlab("DBH(cm)")+ylab("Biomass(kg/tree)")+
  theme1+theme(legend.title=element_blank())

ggplot()+  #plot branch and root biomass
  geom_line(aes(x=DBH, y=BB, color="branch"), size=1)+
  geom_line(aes(x=DBH, y=RB, color="root"), size=1)+
  scale_color_manual(values = c("branch" = "orange", "root" = "grey"))+
  xlab("DBH(cm)")+ylab("Biomass(kg/tree)")+
  theme1+theme(legend.title=element_blank())

ggplot()+  #plot total, aboveground, and stem biomass
  geom_line(aes(x=DBH, y=TB, color="total"), size=1)+
  geom_line(aes(x=DBH, y=TAB, color="aboveground"), size=1)+
  geom_line(aes(x=DBH, y=SB, color="stem"), size=1)+
  scale_color_manual(values = c("total" = "black", "aboveground" = "red", "stem" = "green"))+
  xlab("DBH(cm)")+ylab("Biomass(kg/tree)")+
  theme1+theme(legend.title=element_blank())

#ALLOCATION----------------------------
#determine allocation ratio
# and we assume that this ratio does not change without an impact of **substantial** stress
allocation_dat <- read.csv("allocation-data.csv") #Figure 1 Le Goff et al. 2004 EDP Sciences
prop_foliage <- sum(with(allocation_dat, allocation_dat[component=="foliage","propotion"]))/3
prop_stem <- sum(with(allocation_dat, allocation_dat[component=="stem","propotion"]))/3
prop_branch <- sum(with(allocation_dat, allocation_dat[component=="branch","propotion"]))/3
prop_root <- sum(with(allocation_dat, allocation_dat[component=="root","propotion"]))/3

