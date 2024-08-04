##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
##    Transpiration for Fraxinus excelsior
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Create at: 2024-07-06 15:58:20 CST
## Author: Haoran Wu 
## Affiliation: School of Geography and the Environment, University of Oxford
## Expertise: Forest Pathology; Tree Pests & Diseases; Ecosystem Modelling
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Description: Calculate transpiration of European ash (Fraxinus 
##    excelsior L.). Transpiration = stomatal conductance * VPD.
##    It is assumed that water supply by stem water transport and 
##    root water uptake is abundant (i.e. ash tree is not affected 
##    by drought stress). Data obtained from gas exchange measurements 
##    by Majewski et al. 2024
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

require(ggplot2)
require(cowplot)

#UTILS---------------------------------------------------------
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

#DEMAND--------------------------------------------------------
# Water Use Efficiency data
# Majewski et al. 2024. Diurnal variation curve of Pn, gs, PAR, and WUE
dat <- read.csv("D:/0_Files/0_DGVM/ADB-carbon-model/Majewski-2024-data.csv")
dat$VPD..kPa. <- dat$VPD..Pa./1000 #transform unit
dat$transpir..mol.m2.s. <- with(dat, Pn..umol.m2.s./WUE..umol.mol.) #transpiration without water limitation

# Model: transpiration = k * gs * VPD
md <- lm(data = dat, transpir..mol.m2.s. ~ 0 + I(gs..mol.m2.s.*VPD..kPa.)) #red line in ggplot
ggplot()+
  geom_point(data = dat, aes(x=gs..mol.m2.s. * VPD..kPa., y = transpir..mol.m2.s.), size=3)+
  geom_abline(intercept = 0, slope = 1, size = 1.5, color = "blue", linetype = 2)+
  geom_abline(intercept = 0, slope = md$coefficients, size = 1.5, color = "red", linetype = 2)+
  xlab(expression(g[s]*VPD))+
  ylab("Transpiration")+theme1

## From the result we can see that 
##   (1) demand: transpiration can be approximated by gs * VPD (i.e. k = 1), according to the blue line
##   (2) supply: but at high gs*VPD values (perhaps caused by high VPD), actual transpiration is lower 
##               because of limited water transport ability. See the red line that has lower slope
##
## Therefore, actual transpiration = gs * VPD * f, where f = supply/demand, and demand = gs * VPD
## Supply from xylem transport: K * (Psi_soil - Psi_leaf)
##   where K is hydraulic conductance
##         Psi_soil the soil water potential
##         Psi_leaf the leaf water potential


