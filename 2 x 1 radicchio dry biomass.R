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
bulk.biomass<-read.csv("Radicchio dry biomass.csv")

#change data structure if needed
str(bulk.biomass)

#remove NA values
bulk.biomass <- bulk.biomass %>% mutate(herbivory = replace(herbivory, is.na(herbivory), 0))
bulk.biomass <- bulk.biomass %>% mutate(rot = replace(rot, is.na(rot), 0))
bulk.biomass <- bulk.biomass %>% mutate(general.damage = replace(general.damage, is.na(general.damage), 0))
bulk.biomass <- bulk.biomass %>% mutate(internal.damage = replace(internal.damage, is.na(internal.damage), 0))

# create column for sum of shoot and root mass
bulk.biomass <- bulk.biomass %>%
  mutate(total.mass = shoot.mass + root.mass)

#subset by treatment, run anovas to determine differences based on plot within treatment
    ##PV 
PV.biomass <- bulk.biomass %>%
  filter(treatment == "PV")

    ##Control
control.biomass <- bulk.biomass %>%
  filter(treatment == "control")


#OVERALL

m1<-lm(total_dry_biomass_g~treatment, bulk.biomass)
anova(m1)
summary(m1)
plot(allEffects(m1))
d1=cld(emmeans(m1, ~ treatment))

m1<-lm(shoot.mass~treatment, bulk.biomass)
anova(m1)
summary(m1)
plot(allEffects(m1))
d1=cld(emmeans(m1, ~ treatment))

m1<-lm(root.mass~treatment, bulk.biomass)
anova(m1)
summary(m1)
plot(allEffects(m1))
d1=cld(emmeans(m1, ~ treatment))


#RELATIONSHIPS
m1<-lm(root.mass~shoot.mass*treatment, bulk.biomass)
anova(m1)
summary(m1)
plot(allEffects(m1))
d1=cld(emmeans(m1, ~ shoot.mass*treatment))


weedycheckPV<-subset(bulk.biomass, treatment == "PV" & location == "Weedy Check")
PV<-subset(bulk.biomass, treatment == "PV" & location != "Weedy Check")

cweedycheck<-subset(bulk.biomass, treatment == "control" & location == "Weedy Check")
c<-subset(bulk.biomass, treatment == "control" & location != "Weedy Check")

mean(PV$shoot.mass)
mean(weedycheckPV$shoot.mass)
mean(c$shoot.mass)
mean(cweedycheck$shoot.mass)
##########
##########
###########
##########
########
AV<-subset(bulk.biomass, treatment == "AV")
Control<-subset(bulk.biomass, treatment == "Control")


Control$Crootshoot<-Control$root_dry_biomass_g/Control$shoot_dry_biomass_g
AV$AVrootshoot<-AV$root_dry_biomass_g/AV$shoot_dry_biomass_g

Crootshoot<-summaryBy(Crootshoot ~ Treatment, FUN = c(mean,std.error), na.rm = T, Control)
AVrootshoot<-summaryBy(AVrootshoot ~ Treatment, FUN = c(mean,std.error), na.rm = T, AV)
#######
######
#######
##WETMASSRAD

tiff(file = "2x1 dry radicchio biomass.tiff", height = 12, width = 6, res = 300, units = "in", compression = "zip+p")
par(mfrow = c(2,1), omi = c(1, 0.5, 0.1, 0.1), mar = c(1,6,0.2,0.5))

xx<-c(-500,500); yy<-c(-500,500)


plot(yy ~ xx, pch = NA, xlab="",ylab="",xaxt="n",yaxt="n",ylim=c(0,10), xlim=c(24,30))

m1<-lm(total_dry_biomass_g~treatment, bulk.biomass)
anova(m1)
summary(m1)
# plot(allEffects(m1))
d1=cld(emmeans(m1, ~ treatment))

dum<-subset(d1, treatment == "Control")
rect(25, 00, 27, dum$emmean, col = "black", lwd = 2, border = "black")
ablineclip(v=26, y1=as.numeric(dum$emmean) + (dum$SE), y2=as.numeric(dum$emmean) - (dum$SE),lwd = 4, col = "grey69")
box()

dum<-subset(d1, treatment == "AV")
rect(27, 00, 29, dum$emmean, col = "mediumpurple", lwd = 2, border = "mediumpurple")
ablineclip(v=28, y1=as.numeric(dum$emmean) + (dum$SE), y2=as.numeric(dum$emmean) - (dum$SE),lwd = 4, col = "black")
box()



text(27, 8, expression('***'),cex=3, col = "black")
text(24.5, 9, expression('b)'),cex=2, col = "black")

axis(2, at = seq(0,10,2), las = 2, cex.axis = 1.8)
mtext(side = 2, expression(Total~dry~biomass~(g)), cex = 2, padj = - 2.5, outer= F)

legend("topright",c("Control" , "Agrivoltaic"), col=c( "black", "mediumpurple"), pch= c(15,15) ,  cex = 1.7,  horiz = F, bty='n')




############################ROOTSHOOT RATIO

plot(bulk.biomass$root_dry_biomass_g ~ bulk.biomass$shoot_dry_biomass_g, bulk.biomass, pch = 1, cex = 1.1, col = "black" ,xlim = c(0,12), ylim = c(0,6), xaxt="n",yaxt="n",xlab="",ylab="")

C<-subset(bulk.biomass, treatment == "Control")
par(new=T)
plot(C$root_dry_biomass_g ~ C$shoot_dry_biomass_g, C, pch = 1,cex = 1, lwd= 4, col = "grey42" , xlim = c(0,12), ylim = c(0,6), xaxt="n",yaxt="n",xlab="",ylab="")


PV<-subset(bulk.biomass, treatment == "AV")
par(new=T)
plot(PV$root_dry_biomass_g ~ PV$shoot_dry_biomass_g, PV, pch = 1,cex = 1, lwd = 4, col = "plum4" , xlim = c(0,12), ylim = c(0,6), xaxt="n",yaxt="n",xlab="",ylab="")


m1<-lm(root_dry_biomass_g~shoot_dry_biomass_g, bulk.biomass)
anova(m1)
summary(m1)
# plot(allEffects(m1))
# d1=cld(emmeans(m1, ~ shoot.mass*treatment))
p1<-plot_model(m1, type= c("pred"), terms= c("shoot_dry_biomass_g"))
new<-as.data.frame(p1$data)
m2<-lm(predicted~x, new)
coef(lm(predicted~x, new))
text(500, 200, expression(italic(r)^2~'= 0.58'),cex=2, col = "grey22")
# text(0.26, -3, expression(italic(p)~'= 0.41'),cex=1)
ablineclip(m2,x1=min(bulk.biomass$shoot_dry_biomass_g,na.rm = TRUE),x2=max(bulk.biomass$shoot_dry_biomass_g,na.rm = TRUE), lty = 3, lwd=8, col = "black")

# m1<-lm(root_dry_biomass_g~shoot_dry_biomass_g, PV)
# anova(m1)
# summary(m1)
# # plot(allEffects(m1))
# # d1=cld(emmeans(m1, ~ shoot.mass*treatment))
# p1<-plot_model(m1, type= c("pred"), terms= c("shoot_dry_biomass_g"))
# new<-as.data.frame(p1$data)
# m2<-lm(predicted~x, new)
# coef(lm(predicted~x, new))
# # text(5, 5, expression(italic(r)^2~'= 0.58'),cex=2, col = "grey22")
# # text(0.26, -3, expression(italic(p)~'= 0.41'),cex=1)
# ablineclip(m2,x1=min(PV$shoot_dry_biomass_g,na.rm = TRUE),x2=max(PV$shoot_dry_biomass_g,na.rm = TRUE), lty = 1, lwd=4, col = "purple4")


axis(2, at = seq(0,6,1), las = 2, cex.axis = 1.8)
# mtext(side = 2, expression(paste(ANPP~(g~m^2))), cex = 2.5, padj = -1.7, outer= F)
mtext(side = 2, expression(Root~dry~mass~(g)), cex = 2, padj = -2.5, outer= F)
axis(1, at = seq(0,12,2), las = 1, cex.axis = 1.8)
mtext(side = 1, expression(Shoot~dry~mass~(g)), cex = 2, padj = 2, outer= F)
text(1, 5.5, expression('d)'),cex=2, col = "black")

# legend("topleft",c("Control" , "Photovoltaic"), col=c( "black", "purple4"), pch= c(1,1) , lty = c(1,1), cex = 1.5, lwd = 3, horiz = F, bty='n')
# legend("topright",c("PV","Control"), col=c( "grey44", "black"), pch= c(NA) , lty = c(1), cex = 1.5, lwd = 6, horiz = F, bty='n')
dev.off()


