library(cubature)
WcovidLM_data<-read.csv('wncLM.csv',  fileEncoding="UTF-8-BOM",header=T);
str(WcovidLM_data)
table(WcovidLM_data$date)
tn=NROW(WcovidLM_data); vn=1:tn
arp=2; bh=2
Up=0.975; Lp=0.025     
wdB=rep(0, 3) 
vy.s=WcovidLM_data$lnc
vy1=vy.s[arp:(tn-1)]
vy2=vy.s[(arp-1):(tn-2)]
vy=vy.s[(arp+1):tn]
vns=vn[(arp+1):tn]
model3=lm(vy~vns+vy1+vy2)
summary(model3)
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
f3_new=approxfun(pdf3_new$x, pdf3_new$y, yleft=0, yright=0)
LB3_new=wb2$fit[1]+uniroot(function(x) cubintegrate(f3_new,
                                                    lower=-Inf, upper=x, method = "pcubature")$integral-
                             Lp,c(min(pdf3_new$x), max(pdf3_new$x)), extendInt = "yes")$root
UB3_new=wb2$fit[1]+uniroot(function(x) cubintegrate(f3_new,
                                                    lower=-Inf, upper=x, method = "pcubature")$integral-
                             Up,c(min(pdf3_new$x), max(pdf3_new$x)), extendInt = "yes")$root
wdB[2]=UB3_new-LB3_new
Bre_Emp=quantile(re2, probs=c(Lp, Up), type = 7) 
wdB[3]=Bre_Emp[2]-Bre_Emp[1]
print(c(wb2$fit[1]+Bre_Emp[1], wb2$fit[1]+Bre_Emp[2]))
par(mfrow = c(1, 1))
plot(pdf3_new, xlim=c(-1, 1.4), col='blue', lwd=3, xlab="", 
     main='KDE from Two-step-ahead Prediction Residuals')


