rm(list=ls()); 
library(cubature)
WcovidLM_data<-read.csv('wncLM.csv',  fileEncoding="UTF-8-BOM",header=T);
str(WcovidLM_data)
table(WcovidLM_data$date)
tn=NROW(WcovidLM_data); vn=1:tn
arp=2; bh=1
wdB=rep(0, 3) 
vy.s=WcovidLM_data$lnc
vy1=vy.s[arp:(tn-1)]
vy2=vy.s[(arp-1):(tn-2)]
vy=vy.s[(arp+1):tn]
vns=vn[(arp+1):tn]
model3=lm(vy~vns+vy1+vy2)     
summary(model3)
newdat=data.frame(vns=tn+1, vy1=vy.s[tn], vy2=vy.s[tn-1])
Lp=0.025; Up=0.975 
wb=predict(model3, newdat, se.fit=T, interval="prediction", level=Up-Lp)
yhatm3=wb$fit[1]    
wdB[1]=wb$fit[3]-wb$fit[2]
print(c("Point Prediction", wb$fit[1], "Model PI", wb$fit[2], wb$fit[3]))
pdf3_new=density(resid(model3))  
f3_new=approxfun(pdf3_new$x, pdf3_new$y,yleft=0,yright=0)
LB3_new=yhatm3+uniroot(function(x) cubintegrate(f3_new,
                                                lower=-Inf, upper=x, method = "pcubature")$integral-
                         Lp,c(min(pdf3_new$x), max(pdf3_new$x)), extendInt = "yes")$root
UB3_new=yhatm3+uniroot(function(x) cubintegrate(f3_new,
                                                lower=-Inf, upper=x, method = "pcubature")$integral-
                         Up,c(min(pdf3_new$x), max(pdf3_new$x)), extendInt = "yes")$root
print(c(LB3_new, UB3_new))
wdB[2]=UB3_new-LB3_new
Bre_Emp=quantile(resid(model3), probs=c(Lp, Up), type = 7)
print(c(wb$fit[1]+Bre_Emp[1], wb$fit[1]+Bre_Emp[2]))
wdB[3]=Bre_Emp[2]-Bre_Emp[1]
print(format(wdB, digits = 3))
par(mfrow = c(1, 1))
plot(pdf3_new, xlim=c(-0.6, 0.6), col='blue', lwd=3, xlab="", 
     main='KDE from One-step-ahead Prediction Residuals')

