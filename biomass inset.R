library(plotrix);library(readxl);library(emmeans)
library(permute);library(lattice);library(ggpubr);library(lme4);library(lmerTest)
library(multcomp);library(ggrepel);library(agricolae);library(grid);library(doBy)
library(sjPlot)

library(dplyr); library(readxl); library(lubridate)
library(minpack.lm); library(lattice)
library(plantecophys); library(effects)

library(car); library(lubridate); library(lme4);
library(MuMIn);library(car); library(multcomp); library(effects)
library(lme4); library(MuMIn); library(lmerTest); library(doBy); library(plotrix); library(emmeans); library(car); library(lubridate); library(lme4);
library(MuMIn);library(car); library(multcomp); library(effects)
library(lme4); library(MuMIn); library(lmerTest); library(doBy); library(plotrix) ; library(vioplot)


# Set working directory, load in data, and read the specific sheet
setwd("/Users/mattsturchio/Desktop/Cornell/Albany")
bulk.biomass<-read.csv("radish data 2024 plus dry.csv")

#change data structure
str(bulk.biomass)
bulk.biomass$root.mass <- as.numeric(bulk.biomass$root.mass)
bulk.biomass$number <- as.numeric(bulk.biomass$number)
bulk.biomass$type <- as.factor(bulk.biomass$type)
bulk.biomass$date<-as.Date(bulk.biomass$date, format = "%m/%d/%Y")

#remove NA and 0 values
bulk.biomass <- bulk.biomass %>% mutate(number = replace(number, is.na(number), 0))
bulk.biomass2 <- subset(bulk.biomass, number > 0) 
# create column for sum of shoot and root mass
bulk.biomass2 <- bulk.biomass2 %>%
  mutate(total.mass = shoot.mass + root.mass)


bb2first<-subset(bulk.biomass2, date == "0024-10-29" )
m1<-lm(total.mass~treatment, bb2first)
anova(m1)
summary(m1)
plot(allEffects(m1))
d1=cld(emmeans(m1, ~ treatment))

bb2second<-subset(bulk.biomass2, date == "0024-11-05" )
m1<-lm(total.mass~treatment, bb2second)
anova(m1)
summary(m1)
plot(allEffects(m1))

bb2third<-subset(bulk.biomass2, date == "0024-11-12" )
m1<-lm(total.mass~treatment, bb2third)
anova(m1)
summary(m1)
plot(allEffects(m1))



#OVERALL

m1<-lm(total.mass~treatment, bulk.biomass2)
anova(m1)
summary(m1)
plot(allEffects(m1))
d1=cld(emmeans(m1, ~ treatment))

m1<-lm(shoot.mass~treatment, bulk.biomass2)
anova(m1)
summary(m1)
plot(allEffects(m1))
d1=cld(emmeans(m1, ~ treatment))

m1<-lm(root.mass~treatment, bulk.biomass2)
anova(m1)
summary(m1)
plot(allEffects(m1))
d1=cld(emmeans(m1, ~ treatment))


#RELATIONSHIPS
m1<-lm(root.mass~shoot.mass*treatment, bulk.biomass2)
anova(m1)
summary(m1)
plot(allEffects(m1))
d1=cld(emmeans(m1, ~ shoot.mass*treatment))

##########
##########
###########
##########
########

tiff(file = "2x1 biomass inset.tiff", height = 12, width = 6, res = 300, units = "in", compression = "zip+p")
par(mfrow = c(2,1), omi = c(1, 1, 0.1, 0.1), mar = c(1,6,0.2,0.5))

############################SHOOTMASS##########
xx<-c(-500,500); yy<-c(-500,500)


plot(yy ~ xx, pch = NA, xlab="",ylab="",xaxt="n",yaxt="n",ylim=c(0,2), xlim=c(24,30))

m1<-lm(dhoot.dry.biomass~treatment, bulk.biomass2)
anova(m1)
summary(m1)
# plot(allEffects(m1))
d1=cld(emmeans(m1, ~ treatment))

xx<-c(-500,500); yy<-c(-500,500)

dum<-subset(d1, treatment == "Control")
rect(25, 00, 27, dum$emmean, col = "black", lwd = 2, border = "black")
ablineclip(v=26, y1=as.numeric(dum$emmean) + (dum$SE), y2=as.numeric(dum$emmean) - (dum$SE),lwd = 4, col = "grey69")
box()

dum<-subset(d1, treatment == "PV")
rect(27, 00, 29, dum$emmean, col = "mediumpurple", lwd = 2, border = "mediumpurple")
ablineclip(v=28, y1=as.numeric(dum$emmean) + (dum$SE), y2=as.numeric(dum$emmean) - (dum$SE),lwd = 4, col = "black")
box()

text(27, 1.8, expression('*'),cex=5, col = "black")
axis(2, at = seq(0,2,0.5), las = 2, cex.axis = 3)
mtext(side = 2, expression(Shoot~(g)), cex = 4, padj = -1.75, outer= F)
# text(24.5, 65, expression('b2)'),cex=2, col = "black")

# legend("topright",c("Control" , "Photovoltaic"), col=c( "black", "mediumpurple"), pch= c(15,15) ,  cex = 1.7,  horiz = F, bty='n')


############################ROOTMASS##########


xx<-c(-500,500); yy<-c(-500,500)


plot(yy ~ xx, pch = NA, xlab="",ylab="",xaxt="n",yaxt="n",ylim=c(0,3.5), xlim=c(24,30))

m1<-lm(root.dry.biomass~treatment, bulk.biomass2)
anova(m1)
summary(m1)
# plot(allEffects(m1))
d1=cld(emmeans(m1, ~ treatment))

xx<-c(-500,500); yy<-c(-500,500)

dum<-subset(d1, treatment == "Control")
rect(25, 00, 27, dum$emmean, col = "black", lwd = 2, border = "black")
ablineclip(v=26, y1=as.numeric(dum$emmean) + (dum$SE), y2=as.numeric(dum$emmean) - (dum$SE),lwd = 4, col = "grey69")
box()

dum<-subset(d1, treatment == "PV")
rect(27, 00, 29, dum$emmean, col = "mediumpurple", lwd = 2, border = "mediumpurple")
ablineclip(v=28, y1=as.numeric(dum$emmean) + (dum$SE), y2=as.numeric(dum$emmean) - (dum$SE),lwd = 4, col = "black")
box()



text(27, 3.2, expression('***'),cex=5, col = "black")

axis(2, at = seq(0,3,1), las = 2, cex.axis = 3)
mtext(side = 2, expression(Root~(g)), cex = 4, padj = -1.75, outer= F)
# text(24.5, 65, expression('b3)'),cex=2, col = "black")

############################ROOTSHOOT RATIO


# legend("topleft",c("Control" , "Photovoltaic"), col=c( "black", "purple4"), pch= c(1,1) , lty = c(1,1), cex = 1.5, lwd = 3, horiz = F, bty='n')
# legend("topright",c("PV","Control"), col=c( "grey44", "black"), pch= c(NA) , lty = c(1), cex = 1.5, lwd = 6, horiz = F, bty='n')
dev.off()

