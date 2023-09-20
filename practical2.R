# Import the data 
workf <- data.frame(read.csv("http://www.stats.ox.ac.uk/~laws/SB1/data/workf.csv"))
head(workf)

# Draw a scatterplot
pdf("scatterplot1.pdf", height=8, width=8) # Here we add some noise to the binary variable for easy observation
plot(jitter(employed, amount=0.03) ~ income, data=workf, main="Employed or not vs Family income", xlab="Family Income in $1000s", ylab="Employed or not")
dev.off()

pdf("scatterplot2.pdf", height=8, width=8)
plot(jitter(employed, amount=0.03) ~ educ, data=workf, main="Employed or not vs Education", xlab="Number of years of education", ylab="Employed or not")
dev.off()

pdf("scatterplot3.pdf", height=8, width=8)
plot(income ~ educ, data=workf, main="Income vs Education", xlab="Number of years of education", ylab="Family Income in $1000s")
dev.off()

# Generate univariate models
wf.uglm1 <- glm(employed ~ region, data=workf, family=binomial)
options(digits=5)
summary(wf.uglm1)

wf.uglm2 <- glm(employed ~ ch04, data=workf, family=binomial)
options(digits=5)
summary(wf.uglm2)

wf.uglm3 <- glm(employed ~ ch59, data=workf, family=binomial)
options(digits=5)
summary(wf.uglm3)

wf.uglm4 <- glm(employed ~ ch10, data=workf, family=binomial)
options(digits=5)
summary(wf.uglm4)

wf.uglm5 <- glm(employed ~ income, data=workf, family=binomial)
options(digits=5)
summary(wf.uglm5)

wf.uglm6 <- glm(employed ~ educ, data=workf, family=binomial)
options(digits=5)
summary(wf.uglm6)

# Generate the 1st GLM including interation terms
wf.glm1 <- glm(employed ~ region + ch04 + ch59 + ch10 + income + educ + income:region + income:ch04 + income:ch59 + income:ch10 + educ:region + educ:ch04 + educ:ch59 + educ:ch10, data=workf, family=binomial)
options(digits=5)
summary(wf.glm1)

# Using AIC test to generate a more simplified model with lower AIC
step(wf.glm1) 

wf.glm2 <- glm(employed ~ region + ch04 + ch10 + ch59 + income + educ + income:ch59 +income:ch10 + region:educ + educ:ch59 + educ:ch10, data=workf, family=binomial)
options(digits=5)
summary(wf.glm2)
anova(wf.glm2,wf.glm1)
1 - pchisq(4.67,6) # P-value

wf.glm3 <- glm(employed ~ region + ch04 + ch59 + income + educ + region:educ, data=workf, family=binomial)
options(digits=5)
summary(wf.glm3)
anova(wf.glm3,wf.glm2)
1 - pchisq(2.7,5) 

# Draw graphs to detect outliers
pdf("leverage.pdf", height=8, width=8)
p <- 7
n <- 1935
plot(influence(wf.glm3)$hat/(p/n),pch=19, col='blue', ylab='Leverage / (p/n)')
dev.off()

pdf("CooksDistance.pdf", height=8, width=8)
plot(cooks.distance(wf.glm3),pch=19, col='blue', ylab="Cook's Distance")
abline(h=8/(n-2*p),col='red')
dev.off()

# Refit the data
d<-cooks.distance(wf.glm3)>(8/(n-2*p))
sum(d)
workf1<-workf[-which(d), ]
wf.glm4 <- glm(employed ~ region + ch04 + ch59 + income + educ + region:educ, data=workf1, family=binomial)

pdf("CooksDistance1.pdf", height=8, width=8)
plot(cooks.distance(wf.glm4),pch=19, col='blue', ylab="Cook's Distance")
abline(h=8/(n-2*p),col='red')
dev.off()

# Refit the data for the 2nd time
n1=n-sum(d)
d1<-cooks.distance(wf.glm4)>(8/(n1-2*p))
workf2<-workf1[-which(d1), ]
wf.glm5 <- glm(employed ~ region + ch04 + ch59 + income + educ + region:educ, data=workf2, family=binomial)

pdf("leverage2.pdf", height=8, width=8)
plot(influence(wf.glm5)$hat/(p/n),pch=19, col='blue', ylab='Leverage / (p/n1)')
dev.off()

pdf("CooksDistance2.pdf", height=8, width=8)
plot(cooks.distance(wf.glm5),pch=19, col='blue', ylab="Cook's Distance")
abline(h=8/(n1-2*p),col='red')
dev.off()
# Now we can see that there are only a few points above the threshold, thus we stop the refitting process.

# Check whether the data agrees with the model
options(digits=5)
summary(wf.glm5)

