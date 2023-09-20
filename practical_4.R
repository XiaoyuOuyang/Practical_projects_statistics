library(numDeriv)
require(wooldridge)
inlf = mroz$inlf
educ = mroz$educ
exper = mroz$exper
expersq = mroz$expersq
othinc = mroz$nwifeinc
kidslt6 = 1*(mroz$kidslt6>0)
kidsge6 = 1*(mroz$kidsge6>0)
mdata = data.frame(inlf,educ,exper,expersq,othinc,kidslt6,kidsge6)
logitres = glm(inlf ~ educ + exper + expersq + othinc
                + kidslt6 + kidsge6, data = mdata, family = "binomial")
## Question A(a)
# Define beta.hat.n
beta.hat.n = numeric(7)
for (i in 1:7){
  beta.hat.n[i]=logitres$coefficients[i]
}
# Generate data vector x(i)
x <- function(i){
  a=numeric(7)
  a[1]=1
  a[2]=educ[i]
  a[3]=exper[i]
  a[4]=expersq[i]
  a[5]=othinc[i]
  a[6]=kidslt6[i]
  a[7]=kidsge6[i]
  return(a)
}
# Calculate Gamma_hat_0
n=753
b = numeric(n)
c = numeric(n)
for (i in 1:n){
  c[i]=exp(x(i) %*% beta.hat.n)
  b[i]=c[i]/(1+c[i])^2*beta.hat.n[5]/n
}
gammahat.0=sum(b)
print(gammahat.0)
# CI using normal approximation
alpha=0.05
var=vcov(logitres)
h<-function(beta){
  for (i in 1:n){
    c[i]=exp(x(i) %*% beta)
    b[i]=c[i]/(1+c[i])^2*beta[5]/n
  }
  sum(b)
}
J=jacobian(h,beta.hat.n)
var.delta=J%*%var%*%t(J)
# Calculate the confidence interval
ci.normal.gamma.0=c(gammahat.0-sqrt(var.delta)*qnorm(1-alpha/2),gammahat.0+sqrt(var.delta)*qnorm(1-alpha/2))
print(ci.normal.gamma.0)

## Question (b)
# Non-parametric paired bootstrap estimate
bootstrap_variance_gammahat.0<-function(B){
  beta.hat.n.boot=numeric(7)
  gammahat.0.boot=numeric(B)
  for (j in 1:B) {
    ind = sample(n,replace=TRUE)
    mdata.boot <- mdata[ind,]
    logitres.boot = glm(inlf ~ educ + exper + expersq + othinc
                       + kidslt6 + kidsge6, data = mdata.boot, family = "binomial")
    for (i in 1:7){
    beta.hat.n.boot [i]= logitres.boot$coefficients[i]
  }
    for (i in 1:n){
      c[i]=exp(x(i) %*% beta.hat.n.boot)
      b[i]=c[i]/(1+c[i])^2*beta.hat.n.boot[5]/n
    }
    gammahat.0.boot[j]=sum(b)
  }
  return(c(var(gammahat.0.boot),gammahat.0.boot))
}
var.gammahat.0.boot=bootstrap_variance_gammahat.0(1000)[1]
gammahat.0.boot=bootstrap_variance_gammahat.0(1000)[2:1000+1]
# CI using normal method
ci.bootnormal.gamma.0=c(gammahat.0-sqrt(var.gammahat.0.boot)*qnorm(1-alpha/2),gammahat.0+sqrt(var.gammahat.0.boot*qnorm(1-alpha/2)))
print(ci.bootnormal.gamma.0)
# CI using bootstrap method
ci.bootstrap.gamma.0=c(2*gammahat.0-quantile(gammahat.0.boot,1-alpha/2,names=FALSE),2*gammahat.0-quantile(gammahat.0.boot,alpha/2,names=FALSE))
print(ci.bootstrap.gamma.0)

# (c) Plot the empricial cdf
pdf("Empirical cdf.pdf", height=8, width=8)
plot(ecdf(gammahat.0.boot-gammahat.0))
dev.off()

# (d) Marginal Effect at x0
x0=c(1,10,3,9,10,1,0)
# Calculate new Gamma_hat_0
r=exp(x0 %*% beta.hat.n)
new.gammahat.0=beta.hat.n[5]/(1+r)^2*r
print(new.gammahat.0)
# CI using normal approximation
g<-function(beta){
    c=exp(x0 %*% beta)
    b=c/(1+c)^2*beta[5]
}
new.J=jacobian(g,beta.hat.n)
new.var.delta=new.J%*%var%*%t(new.J)
new.ci.normal.gamma.0=c(new.gammahat.0-sqrt(new.var.delta)*qnorm(1-alpha/2),new.gammahat.0+sqrt(new.var.delta)*qnorm(1-alpha/2))
print(new.ci.normal.gamma.0)
# nonparametric paired bootstrap estimate using x0
new_bootstrap_variance_gammahat.0<-function(B){
  beta.hat.n.boot=numeric(7)
  gammahat.0.boot=numeric(B)
  for (j in 1:B) {
    ind = sample(n,replace=TRUE)
    mdata.boot <- mdata[ind,]
    logitres.boot = glm(inlf ~ educ + exper + expersq + othinc
                        + kidslt6 + kidsge6, data = mdata.boot, family = "binomial")
    for (i in 1:7){
      beta.hat.n.boot[i]=logitres.boot$coefficients[i]
    }
    c=exp(x0 %*% beta.hat.n.boot)
    gammahat.0.boot[j]=c/(1+c)^2*beta.hat.n.boot[5]
  }
  return(c(var(gammahat.0.boot),gammahat.0.boot))
}
new.var.gammahat.0.boot=new_bootstrap_variance_gammahat.0(1000)[1]
new.gammahat.0.boot=new_bootstrap_variance_gammahat.0(1000)[2:1000+1]
# CI using normal method
new.ci.bootnormal.gamma.0=c(new.gammahat.0-sqrt(new.var.gammahat.0.boot)*qnorm(1-alpha/2),new.gammahat.0+sqrt(new.var.gammahat.0.boot*qnorm(1-alpha/2)))
print(new.ci.bootnormal.gamma.0)
# CI using bootstrap method
new.ci.bootstrap.gamma.0=c(2*new.gammahat.0-quantile(new.gammahat.0.boot,1-alpha/2,names=FALSE),2*new.gammahat.0-quantile(new.gammahat.0.boot,alpha/2,names=FALSE))
print(new.ci.bootstrap.gamma.0)
# (e) Calculate new Gamma_hat_0
r=exp(x0 %*% beta.hat.n)
e.new.gammahat.0=beta.hat.n[6]/(1+r)^2*r
print(e.new.gammahat.0)
# nonparametric paired bootstrap estimate using x0
e_new_bootstrap_variance_gammahat.0<-function(B){
  beta.hat.n.boot=numeric(7)
  gammahat.0.boot=numeric(B)
  for (j in 1:B) {
    ind = sample(n,replace=TRUE)
    mdata.boot <- mdata[ind,]
    logitres.boot = glm(inlf ~ educ + exper + expersq + othinc
                        + kidslt6 + kidsge6, data = mdata.boot, family = "binomial")
    for (i in 1:7){
      beta.hat.n.boot[i]=logitres.boot$coefficients[i]
    }
    c=exp(x0 %*% beta.hat.n.boot)
    gammahat.0.boot[j]=c/(1+c)^2*beta.hat.n.boot[6]
  }
  return(c(var(gammahat.0.boot),gammahat.0.boot))
}
e.new.gammahat.0.boot=e_new_bootstrap_variance_gammahat.0(1000)[2:1000+1]
# CI using bootstrap method
e.new.ci.bootstrap.gamma.0=c(2*e.new.gammahat.0-quantile(e.new.gammahat.0.boot,1-alpha/2,names=FALSE),2*e.new.gammahat.0-quantile(e.new.gammahat.0.boot,alpha/2,names=FALSE))
print(e.new.ci.bootstrap.gamma.0)
