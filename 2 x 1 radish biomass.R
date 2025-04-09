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
AV<-subset(bulk.biomass, treatment == "PV")
Control<-subset(bulk.biomass, treatment == "Control")


Control$Crootshoot<-Control$root.dry.biomass/Control$dhoot.dry.biomass
AV$AVrootshoot<-AV$root.dry.biomass/AV$dhoot.dry.biomass

Crootshoot<-summaryBy(Crootshoot ~ Treatment, FUN = c(mean,std.error), na.rm = T, Control)
AVrootshoot<-summaryBy(AVrootshoot ~ Treatment, FUN = c(mean,std.error), na.rm = T, AV)

##########
########

tiff(file = "2x1 biomass.tiff", height = 12, width = 6, res = 300, units = "in", compression = "zip+p")
par(mfrow = c(2,1), omi = c(1, 0.5, 0.1, 0.1), mar = c(1,6,0.2,0.5))

xx<-c(-500,500); yy<-c(-500,500)


plot(yy ~ xx, pch = NA, xlab="",ylab="",xaxt="n",yaxt="n",ylim=c(0,70), xlim=c(24,30))

m1<-lm(total.mass~treatment, bulk.biomass2)
anova(m1)
summary(m1)
# plot(allEffects(m1))
d1=cld(emmeans(m1, ~ treatment))

dum<-subset(d1, treatment == "Control")
rect(25, 00, 27, dum$emmean, col = "black", lwd = 2, border = "black")
ablineclip(v=26, y1=as.numeric(dum$emmean) + (dum$SE), y2=as.numeric(dum$emmean) - (dum$SE),lwd = 4, col = "grey69")
box()

dum<-subset(d1, treatment == "PV")
rect(27, 00, 29, dum$emmean, col = "mediumpurple", lwd = 2, border = "mediumpurple")
ablineclip(v=28, y1=as.numeric(dum$emmean) + (dum$SE), y2=as.numeric(dum$emmean) - (dum$SE),lwd = 4, col = "black")
box()



text(27, 55, expression('***'),cex=3, col = "black")
text(24.5, 65, expression('a)'),cex=2, col = "black")

axis(2, at = seq(0,70,10), las = 2, cex.axis = 1.8)
mtext(side = 2, expression(Total~Biomass~(g)), cex = 2, padj = - 2.5, outer= F)

legend("topright",c("Control" , "Photovoltaic"), col=c( "black", "mediumpurple"), pch= c(15,15) ,  cex = 1.7,  horiz = F, bty='n')

#######
#######
#######
##RELATIONSHIPS

plot(bulk.biomass2$root.mass ~ bulk.biomass2$shoot.mass, bulk.biomass2, pch = 1, cex = 1.1, col = "black" ,xlim = c(0,150), ylim = c(0,300), xaxt="n",yaxt="n",xlab="",ylab="")

C<-subset(bulk.biomass2, treatment == "Control")
par(new=T)
plot(C$root.mass ~ C$shoot.mass, C, pch = 1,cex = 1, lwd= 4, col = "grey42" , xlim = c(0,150), ylim = c(0,300), xaxt="n",yaxt="n",xlab="",ylab="")


PV<-subset(bulk.biomass2, treatment == "PV")
par(new=T)
plot(PV$root.mass ~ PV$shoot.mass, PV, pch = 1,cex = 1, lwd = 4, col = "plum4" , xlim = c(0,150), ylim = c(0,300), xaxt="n",yaxt="n",xlab="",ylab="")


m1<-lm(root.mass~shoot.mass, C)
anova(m1)
summary(m1)
# plot(allEffects(m1))
# d1=cld(emmeans(m1, ~ shoot.mass*treatment))
p1<-plot_model(m1, type= c("pred"), terms= c("shoot.mass"))
new<-as.data.frame(p1$data)
m2<-lm(predicted~x, new)
coef(lm(predicted~x, new))
text(500, 200, expression(italic(r)^2~'= 0.58'),cex=2, col = "grey22")
# text(0.26, -3, expression(italic(p)~'= 0.41'),cex=1)
ablineclip(m2,x1=min(C$shoot.mass,na.rm = TRUE),x2=max(C$shoot.mass,na.rm = TRUE), lty = 1, lwd=4, col = "black")

m1<-lm(root.mass~shoot.mass, PV)
anova(m1)
summary(m1)
# plot(allEffects(m1))
# d1=cld(emmeans(m1, ~ shoot.mass*treatment))
p1<-plot_model(m1, type= c("pred"), terms= c("shoot.mass"))
new<-as.data.frame(p1$data)
m2<-lm(predicted~x, new)
coef(lm(predicted~x, new))
text(500, 200, expression(italic(r)^2~'= 0.58'),cex=2, col = "grey22")
# text(0.26, -3, expression(italic(p)~'= 0.41'),cex=1)
ablineclip(m2,x1=min(PV$shoot.mass,na.rm = TRUE),x2=max(PV$shoot.mass,na.rm = TRUE), lty = 1, lwd=4, col = "purple4")


axis(2, at = seq(0,300,50), las = 2, cex.axis = 1.8)
# mtext(side = 2, expression(paste(ANPP~(g~m^2))), cex = 2.5, padj = -1.7, outer= F)
mtext(side = 2, expression(Root~mass~(g)), cex = 2, padj = -2.5, outer= F)
axis(1, at = seq(0,150,50), las = 1, cex.axis = 1.8)
mtext(side = 1, expression(Shoot~mass~(g)), cex = 2, padj = 2, outer= F)
text(16, 280, expression('c)'),cex=2, col = "black")

# legend("topleft",c("Control" , "Photovoltaic"), col=c( "black", "purple4"), pch= c(1,1) , lty = c(1,1), cex = 1.5, lwd = 3, horiz = F, bty='n')
# legend("topright",c("PV","Control"), col=c( "grey44", "black"), pch= c(NA) , lty = c(1), cex = 1.5, lwd = 6, horiz = F, bty='n')
dev.off()
############################DRYMASS##########

bulk.biomass2$totaldry<-(bulk.biomass2$dhoot.dry.biomass+bulk.biomass2$root.dry.biomass)

tiff(file = "2x1 dry radish biomass.tiff", height = 12, width = 6, res = 300, units = "in", compression = "zip+p")
par(mfrow = c(2,1), omi = c(1, 0.5, 0.1, 0.1), mar = c(1,6,0.2,0.5))

xx<-c(-500,500); yy<-c(-500,500)


plot(yy ~ xx, pch = NA, xlab="",ylab="",xaxt="n",yaxt="n",ylim=c(0,6), xlim=c(24,30))

m1<-lm(totaldry~treatment, bulk.biomass2)
anova(m1)
summary(m1)
# plot(allEffects(m1))
d1=cld(emmeans(m1, ~ treatment))

dum<-subset(d1, treatment == "Control")
rect(25, 00, 27, dum$emmean, col = "black", lwd = 2, border = "black")
ablineclip(v=26, y1=as.numeric(dum$emmean) + (dum$SE), y2=as.numeric(dum$emmean) - (dum$SE),lwd = 4, col = "grey69")
box()

dum<-subset(d1, treatment == "PV")
rect(27, 00, 29, dum$emmean, col = "mediumpurple", lwd = 2, border = "mediumpurple")
ablineclip(v=28, y1=as.numeric(dum$emmean) + (dum$SE), y2=as.numeric(dum$emmean) - (dum$SE),lwd = 4, col = "black")
box()



text(27, 4.8, expression('***'),cex=3, col = "black")
text(24.5, 5.5, expression('a)'),cex=2, col = "black")

axis(2, at = seq(0,6,1), las = 2, cex.axis = 1.8)
mtext(side = 2, expression(Total~dry~biomass~(g)), cex = 2, padj = - 2.5, outer= F)

legend("topright",c("Control" , "Agrivoltaic"), col=c( "black", "mediumpurple"), pch= c(15,15) ,  cex = 1.7,  horiz = F, bty='n')

#######
######
######
##DRYRELATIONSHIPS

plot(bulk.biomass2$root.dry.biomass ~ bulk.biomass2$dhoot.dry.biomass, bulk.biomass2, pch = 1, cex = 1.1, col = "black" ,xlim = c(0,10), ylim = c(0,20), xaxt="n",yaxt="n",xlab="",ylab="")

C<-subset(bulk.biomass2, treatment == "Control")
par(new=T)
plot(C$root.dry.biomass ~ C$dhoot.dry.biomass, C, pch = 1,cex = 1, lwd= 4, col = "grey42" , xlim = c(0,10), ylim = c(0,20), xaxt="n",yaxt="n",xlab="",ylab="")


PV<-subset(bulk.biomass2, treatment == "PV")
par(new=T)
plot(PV$root.dry.biomass ~ PV$dhoot.dry.biomass, PV, pch = 1,cex = 1, lwd = 4, col = "plum4" , xlim = c(0,10), ylim = c(0,20), xaxt="n",yaxt="n",xlab="",ylab="")


m1<-lm(root.dry.biomass~dhoot.dry.biomass, C)
anova(m1)
summary(m1)
# plot(allEffects(m1))
# d1=cld(emmeans(m1, ~ shoot.mass*treatment))
p1<-plot_model(m1, type= c("pred"), terms= c("dhoot.dry.biomass"))
new<-as.data.frame(p1$data)
m2<-lm(predicted~x, new)
coef(lm(predicted~x, new))
text(500, 200, expression(italic(r)^2~'= 0.58'),cex=2, col = "grey22")
# text(0.26, -3, expression(italic(p)~'= 0.41'),cex=1)
ablineclip(m2,x1=min(C$dhoot.dry.biomass,na.rm = TRUE),x2=max(C$dhoot.dry.biomass,na.rm = TRUE), lty = 1, lwd=6, col = "black")

m1<-lm(root.dry.biomass~dhoot.dry.biomass, PV)
anova(m1)
summary(m1)
# plot(allEffects(m1))
# d1=cld(emmeans(m1, ~ shoot.mass*treatment))
p1<-plot_model(m1, type= c("pred"), terms= c("dhoot.dry.biomass"))
new<-as.data.frame(p1$data)
m2<-lm(predicted~x, new)
coef(lm(predicted~x, new))
text(500, 200, expression(italic(r)^2~'= 0.58'),cex=2, col = "grey22")
# text(0.26, -3, expression(italic(p)~'= 0.41'),cex=1)
ablineclip(m2,x1=min(PV$dhoot.dry.biomass,na.rm = TRUE),x2=max(PV$dhoot.dry.biomass,na.rm = TRUE), lty = 1, lwd=6, col = "purple4")


axis(2, at = seq(0,22,5), las = 2, cex.axis = 1.8)
# mtext(side = 2, expression(paste(ANPP~(g~m^2))), cex = 2.5, padj = -1.7, outer= F)
mtext(side = 2, expression(Root~dry~mass~(g)), cex = 2, padj = -2.5, outer= F)
axis(1, at = seq(0,10,2), las = 1, cex.axis = 1.8)
mtext(side = 1, expression(Shoot~dry~mass~(g)), cex = 2, padj = 2, outer= F)
text(0.75, 18, expression('c)'),cex=2, col = "black")

# legend("topleft",c("Control" , "Photovoltaic"), col=c( "black", "purple4"), pch= c(1,1) , lty = c(1,1), cex = 1.5, lwd = 3, horiz = F, bty='n')
# legend("topright",c("PV","Control"), col=c( "grey44", "black"), pch= c(NA) , lty = c(1), cex = 1.5, lwd = 6, horiz = F, bty='n')
dev.off()


