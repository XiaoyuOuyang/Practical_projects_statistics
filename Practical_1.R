# Import data from file restr.csv
restr<-read.csv("http://www.stats.ox.ac.uk/~laws/SB1/data/restr.csv")
attach(restr)
str(restr)
library(MASS)

# Box-plot of price affected by categorical variable guide
pdf("restrbox1.pdf",height=8,width=8)
boxplot(price~guide, xlab="Guide", ylab="Price")
dev.off()

# Box-plot of price affected by categorical variable location
pdf("restrbox2.pdf",height=8,width=8)
boxplot(price~location, xlab="Location", ylab="Price")
dev.off()

# Pairplots of six variables
pdf("restr1.pdf",height=8,width=8)
pairs(restr,lower.panel=NULL)
dev.off()

# Establish a normal linear model LM1
restr.lm1<-lm(price~food+decor+service+guide+location,data=restr)

# Four diagonostic plots to determine outliers
n1<-nrow(restr)
p1<-restr.lm1$rank
i1<-cooks.distance(restr.lm1)>(8/(n1-2*p1)) 

pdf("Cooks1.pdf",height=16,width=16)
par(mfrow = c(2,2))
plot(cooks.distance(restr.lm1),col=1+i1,)
abline(h=8/(n1-2*p1))

plot(rstudent(restr.lm1)~fitted(restr.lm1), xlab="Fitted values",ylab="Studentised Residuals",col=1+i1)
abline(h=-2)
abline(h=2)

qqnorm(rstudent(restr.lm1),col=1+i1)
qqline(rstudent(restr.lm1))

plot(hatvalues(restr.lm1)~fitted(restr.lm1),xlab ="Fitted values",ylab="Hatvalues",col=1+i1)
abline(h=2*p1/n1)
dev.off()

# Remove outliers and refit the LM1 to new data 
refit1<-restr[-which(i1), ]
refit1.lm1<-lm(price~food+decor+service+guide+location,data=refit1)
n2<-nrow(refit1)
p2<-refit1.lm1$rank
i2<-cooks.distance(refit1.lm1)>(8/(n2-2*p2))
sum(i2) # sum(i2)=0 means that there are no outliers in the refitted model
summary(refit1.lm1)
anova(refit1.lm1)

# Compare LM1 with LM2 which drops two variables : service and guide
refit1.lm2<-lm(price~food+decor+location, data=refit1)
anova(refit1.lm2,refit1.lm1)

# Find outliers of LM2
n3<-nrow(refit1)
p3<-refit1.lm2$rank
i3<-cooks.distance(refit1.lm2)>(8/(n3-2*p3)) #Find possible outliers of LM1 

pdf("Cooks2.pdf",height=16,width=16)
par(mfrow = c(2,2))
plot(cooks.distance(refit1.lm2),col=1+i3)
abline(h=8/(n3-2*p3))

plot(rstudent(refit1.lm2)~fitted(refit1.lm2),xlab="Fitted values",ylab="Studentised Residuals",col=1+i3)
abline(h=-2)
abline(h=2)

qqnorm(rstudent(refit1.lm2),col=1+i3)
qqline(rstudent(refit1.lm2))

plot(hatvalues(refit1.lm2)~fitted(refit1.lm2),xlab="Fitted values",ylab="Hatvalues",col=1+i3)
abline(h=2*p3/n3)
dev.off()
summary(refit1.lm2) 

# Establish a new model including interaction terms LM3
restr.nlm1<-lm(price~food+decor+service+guide+location+location:food+location:decor+location:service,data=restr)
nn1<-nrow(restr)
np1<-restr.nlm1$rank
ni1<-cooks.distance(restr.nlm1)>(8/(nn1-2*np1))

pdf("NCooks1.pdf",height=16,width=16)
par(mfrow = c(2,2))
plot(cooks.distance(restr.nlm1),col=1+ni1)
abline(h=8/(nn1-2*np1))

plot(rstudent(restr.nlm1)~fitted(restr.nlm1),xlab="Fitted values",ylab="Studentised Residuals",col=1+ni1)
abline(h=-2)
abline(h=2)

qqnorm(rstudent(restr.nlm1),col=1+ni1)
qqline(rstudent(restr.nlm1))

plot(hatvalues(restr.nlm1)~fitted(restr.nlm1),xlab = "Fitted values", ylab = "Hatvalues",col=1+ni1)
abline(h=2*np1/nn1)
dev.off()

# Refit the data
nrefit1<-restr[-which(ni1), ]
nrefit1.nlm1<-lm(price~food+decor+service+guide+location+location:food+location:decor+location:service,data=nrefit1)
nn2<-nrow(nrefit1)
np2<-nrefit1.nlm1$rank
ni2<-cooks.distance(nrefit1.nlm1)>(8/(nn2-2*np2))
sum(ni2) 
summary(nrefit1.nlm1)
anova(nrefit1.nlm1)

pdf("NCooks2.pdf",height=16,width=16)
par((mfrow = c(2,2)))
plot(cooks.distance(nrefit1.nlm1),col=1+ni2)
abline(h=8/(nn2-2*np2))

plot(rstudent(nrefit1.nlm1)~fitted(nrefit1.nlm1), xlab = "Fitted values", ylab = "Studentised Residuals",col=1+ni2)
abline(h=-2)
abline(h=2)

qqnorm(rstudent(nrefit1.nlm1),col=1+ni2)
qqline(rstudent(nrefit1.nlm1))

plot(hatvalues(nrefit1.nlm1)~fitted(nrefit1.nlm1),xlab = "Fitted values", ylab = "Hatvalues",col=1+ni2)
abline(h=2*np2/nn2)
dev.off()

#Generate new models to compare whether the interaction term is significant
frefit1.lm1<-lm(price~food+decor+service+guide+location+location:decor+location:service,data=nrefit1)
drefit1.lm1<-lm(price~food+decor+service+guide+location+location:food+location:service,data=nrefit1)
srefit1.lm1<-lm(price~food+decor+service+guide+location+location:food+location:decor,data=nrefit1)
anova(frefit1.lm1,nrefit1.nlm1)
anova(drefit1.lm1,nrefit1.nlm1)
anova(srefit1.lm1,nrefit1.nlm1)

#Compare LM2 and LM3 using F-test
refit1.lm3<-lm(price~food+decor+location, data=nrefit1)
anova(refit1.lm3,nrefit1.nlm1)



