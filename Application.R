library(car); library(cubature)
WcovidLM_data<-read.csv('wncLM.csv',  fileEncoding="UTF-8-BOM",header=T);
WcovidLM_data$date=as.Date(WcovidLM_data$date, "%m/%d/%Y")
attach(WcovidLM_data)
tn=NROW(WcovidLM_data)
model1=lm(lnc~I(1:tn))
par(mfrow = c(1, 2))
#
# Figure 1 
plot(WcovidLM_data$date, WcovidLM_data$lnc, xlab=" ", ylab="y", 
     main='(a) log(COVID Case Counts) by Week', col='blue', lwd=3, type="l")
lines(WcovidLM_data$date,fitted(model1), col='purple', lwd=3)
acf(x=resid(model1), main='(b) Autocorrelations of Linear Model Residuals', 
    col='blue', lwd=3)
xarma.s=arima(x=resid(model1), order=c(2,0,0), method=c("CSS"), 
                include.mean=F, optim.control = list(maxit=1000))
#
# Figure 2 
acf(x=resid(xarma.s), main='(a) Autocorrelations of AR(2) Residuals', 
    col='blue', lwd=3)
car::qqPlot(resid(xarma.s), ylab="AR(2) Residual", 
            main='(b) Q-Q Plot')
#
# Figure 3 and Intervals
par(mfrow = c(1, 2))
vn=1:tn; arp=2; Lp=0.025; Up=0.975
wdB=rep(0, 3) 
vy.s=WcovidLM_data$lnc
vy1=vy.s[arp:(tn-1)]
vy2=vy.s[(arp-1):(tn-2)]
vy=vy.s[(arp+1):tn]
vns=vn[(arp+1):tn]
model3=lm(vy~vns+vy1+vy2)     
#
predCI=function(ppdf, pyhatm, pmres, Lp, Up, ep){
  wdBs=rep(0, 2)
  f3_new=approxfun(ppdf$x, ppdf$y,yleft=0,yright=0)
  LB3_new=pyhatm+uniroot(function(x) cubintegrate(f3_new, lower=-Inf, upper=x, 
                                                  method = "pcubature")$integral-Lp,c(min(ppdf$x), max(ppdf$x)), 
                         extendInt = "yes")$root
  UB3_new=pyhatm+uniroot(function(x) cubintegrate(f3_new, lower=-Inf, upper=x, 
                                                  method = "pcubature")$integral-Up,c(min(ppdf$x), max(ppdf$x)), 
                         extendInt = "yes")$root
  wdBs[1]=UB3_new-LB3_new
  Bre_Emp=quantile(pmres, probs=c(Lp, Up), type = 7)
  LB4_new=pyhatm+Bre_Emp[1]; UB4_new=pyhatm+Bre_Emp[2]
  wdBs[2]=Bre_Emp[2]-Bre_Emp[1]
  return(c(LB3_new, UB3_new, wdBs[1], LB4_new, UB4_new, wdBs[2]))
}
#
bh=1
newdat=data.frame(vns=tn+1, vy1=vy.s[tn], vy2=vy.s[tn-1])
wb=predict(model3, newdat, se.fit=T, interval="prediction", level=Up-Lp)
yhatm3=wb$fit[1]    
wdB[1]=wb$fit[3]-wb$fit[2]
print(c("Point Prediction", wb$fit[1], "Model PI", wb$fit[2], wb$fit[3]))
pdf3_new=density(resid(model3))  
plot(pdf3_new, xlim=c(-0.6, 0.6), col='blue', lwd=3, xlab="", 
     main='KDE from One-step-ahead Prediction Residuals')
predCI(ppdf=pdf3_new, pyhatm=wb$fit[1], pmres=resid(model3), Lp, Up)
#
bh=2
yh1=c(rep(0, arp), fitted(model3))   
yh2=rep(0, tn)
yh2[(arp+bh):tn]=model3$coef[1]+model3$coef[2]*((arp+bh):tn)+
  model3$coef[3]*yh1[(arp+bh-1):(tn-1)]+model3$coef[4]*vy.s[(arp+bh-1):(tn-1)]
re2=vy.s[(arp+bh):tn]-yh2[(arp+bh):tn]
newdat1=data.frame(vns=tn+1, vy1=vy.s[tn], vy2=vy.s[tn-1])
wb1=predict(model3, newdat1, se.fit=T, interval="prediction", level =Up-Lp)
newdat2=data.frame(vns=tn+2, vy1=wb1$fit[1], vy2=vy.s[tn])
wb2=predict(model3, newdat2, se.fit=T, interval="prediction", level =Up-Lp)
wdB[1]=wb2$fit[3]-wb2$fit[2]
pdf3_new=density(re2)  
plot(pdf3_new, xlim=c(-1, 1.4), col='blue', lwd=3, xlab="", 
     main='KDE from Two-step-ahead Prediction Residuals')
predCI(ppdf=pdf3_new, pyhatm=wb2$fit[1], pmres=re2, Lp, Up)

