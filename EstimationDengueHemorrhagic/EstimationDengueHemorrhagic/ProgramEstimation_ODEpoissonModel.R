#R libraries
library(deSolve)
library(fields)

#Functions
source("FuncionODE.R", local=FALSE)
source("FuncionMenosLogVeroEscalaOriginal.R", local=FALSE)
source("FuncionMenosLogVerobetaHR01.R", local=FALSE)
source("FuncionMenosLogVerobetaHR02.R", local=FALSE)
source("FuncionMenosLogVerobetaHR0.R", local=FALSE)

#Read database
DatosDen864Original<-as.matrix(read.csv("StudyData_DengueHemorrhagicFever.csv")[,2:4])
Datos<-DatosDen864Original[35:50,]

#GRAFICA DE TODOS LOS DATOS 1200 x 800 (SEMANA 35 HASTA SEMANA 50)
Tiempo <- Datos[,1]
YcObs <- Datos[,2]
YhObs <- Datos[,3]
matplot(Tiempo, cbind(YcObs, YhObs), type="p", pch=c(19,19),
        col=c(1,2), cex.lab=1.35, cex.axis=1.35,
        main="Data: Den864Original", xlab="Week",
        ylab="Infecteds", cex.main=1.30, ylim=c(0,200))
legend("topleft", c("Classical", "Hemorrhagic"), cex=1.2,
       pch=c(19,19), col=c(1,2), bty="n", y.intersp=0.35, inset=-0.025)

#ESCENARIO DE PARAMETROS
LambdaM0 <- 30702.6139006  #antes 30702.6139006
LambdaS0 <- 76.89246       #0.000273973*(283492-283492*0.01): antes 10.2385934233
Lambda10 <- 2834.92        #283492*0.01: antes 1.1376214914
alphaC0 <- 1.1655          #punto medio 0.1665*7: antes 0.686615937276
alphaH0 <- 1.1655          #punto medio 0.1665*7: antes 1.41310092256
b0 <- 21.875               #punto medio 21.875: antes 8
betaH0 <- 0.95             #antes 0.29830
betaM0 <- 0.01503          #antes 0.01450361065995648
muH0 <- 0.000273973        #0.000039139*7: antes 0.000273973
muM0 <- 0.5075             #punto medio 0.0725*7: antes 0.307170720093)
vsigma0 <- 2.5             #punto medio 2.5: antes 15
vtheta0 <- 0.05            #antes 0.05
p0 <- 0.025                #antes 0.05

#CONDICIONES INICIALES
S0 <- 278931      #283492-(120+40+0)-(4400+1+0): antes 35598
I10 <- 120        #antes 120
I20 <- 40         #antes 40
Ms0 <- 120000     #antes 120000
M10 <- 10         #antes 10
M20 <- 10         #antes 10
S10 <- 4400       #antes 4400
Y1h0 <- 1         #antes 1
Y1c0 <- 0         #antes 0
Rec0 <- 0         #antes 0
z0 <- p0*(I10 + I20 + Y1c0)   #ahora 4 antes 8


#TIEMPO SOLUCION DE LAS ODE
T0 <- 35-35
Tfinal <- 53-35
VecTiempo0 <- seq(T0,Tfinal,1)

#CALIBRACION DE VALORES INICIALES
yiniODE <- c(Ms=Ms0, M1=M10, M2=M20, S=S0, I1=I10, I2=I20, S1=S10, Y1c=Y1c0, Y1h=Y1h0, R=Rec0, z=z0)
RejillaTiempoODE <- seq(from=T0,to=Tfinal,by=0.1)
vecparODE <- c(LambdaM=LambdaM0, LambdaS=LambdaS0, Lambda1=Lambda10,
             alphaC=alphaC0, alphaH=alphaH0, b=b0,
             betaH=betaH0, betaM=betaM0, muH=muH0, muM=muM0,
             vsigma=vsigma0, vtheta=vtheta0, p=p0)
outODE <- ode(yiniODE,RejillaTiempoODE,FuncionODE,vecparODE)
Iclassical <- outODE[,12]
Ihemorrhagic <- outODE[,10]

matplot(outODE[,1]+35,cbind(Iclassical,Ihemorrhagic),type="l",lty=1,
        col=c(1,2),cex.lab=1.25,cex.axis=1.25,
        main="Data (Den864Original): A Poisson-Strain Model",xlab="Week",
        ylab="Infecteds",cex.main=1.30)
legend("topleft",c("Classical","Hemorrhagic"),lty=c(1,1),col=c(1,2),bty="n",
       y.intersp=0.35,inset=-0.025)
points(Tiempo,YcObs,pch=19,col=1,cex=0.8)
points(Tiempo,YhObs,pch=19,col=2,cex=0.8)

#GRAFICA DE LOS DATOS PARA ANALISIS 1200 x 800 (SEMANA 35 HASTA SEMANA 40)
matplot(Tiempo[1:6],cbind(YcObs[1:6],YhObs[1:6]),type="p",pch=c(19,19),
        col=c(1,2),cex.lab=1.35,cex.axis=1.35,
        main="Data (Den864Original): Week 35-40",xlab="Week",
        ylab="Infecteds",cex.main=1.30,ylim=c(0,200))
legend("topleft",c("Classical","Hemorrhagic"),cex=1.2,
       pch=c(19,19),col=c(1,2),bty="n",y.intersp=0.35,inset=-0.025)

#CALCULO DE LA LOG-VEROSIMILITUD Y RELATIVA betaH y betaM
TiemposObs <- Tiempo[1:6]-35
CasosObs <- cbind(YcObs[1:6],YhObs[1:6])
TfinalObs <- 6
yiniODE <- c(Ms=Ms0, M1=M10, M2=M20, S=S0, I1=I10, I2=I20, S1=S10, Y1c=Y1c0, Y1h=Y1h0, R=Rec0, z=z0)
RejillaTiempoODE <- seq(from=T0,to=TfinalObs,by=0.1)
ParFijo <- c(LambdaM0,LambdaS0,Lambda10,alphaC0,alphaH0,b0,muH0,muM0,vsigma0,vtheta0,p0)
ValbetaH <- seq(1,3.5,length.out=250) #250
ValbetaM <- seq(0.001,0.0175,length.out=250) #250
FuncionMenosLogVeroEscalaOriginal(c(ValbetaH[1],ValbetaM[1]),TiemposObs,CasosObs,yiniODE,RejillaTiempoODE,ParFijo)
LogVeroGlobal <- matrix(rep(0,length(ValbetaM)*length(ValbetaH)),ncol=length(ValbetaM))
for(i in 1:length(ValbetaH))
{
  for(j in 1:length(ValbetaM))
  {
    Vecparij<-c(ValbetaH[i],ValbetaM[j])
    LogVeroGlobal[i,j]<--FuncionMenosLogVeroEscalaOriginal(Vecparij,TiemposObs,CasosObs,yiniODE,RejillaTiempoODE,ParFijo)
  }
}
write.csv(LogVeroGlobal,"LogVeroGlobalGlobalPoisson.csv")
LogVeroGlobal<-as.matrix(read.csv("LogVeroGlobalGlobalPoisson.csv"))[,-1]
RelativaLogVeroGlobal<-exp(LogVeroGlobal-max(LogVeroGlobal))
write.csv(RelativaLogVeroGlobal,"RelativaLogVeroGlobalPoisson.csv")
RelativaLogVeroGlobal<-as.matrix(read.csv("RelativaLogVeroGlobalPoisson.csv"))[,-1]

#GRAFICA DE CONTORNOS DE LA LOG-VEROSIMILITUD 875 x 708
image.plot(ValbetaH,ValbetaM,RelativaLogVeroGlobal,xlab=expression(beta[H]),ylab=expression(beta[M]),cex.lab=1.25)
contour(ValbetaH,ValbetaM,RelativaLogVeroGlobal,levels = c(0.05), add = TRUE, drawlabels = TRUE,labcex=1)

#ESTIMADORES DE MAXIMA VEROSIMILITUD (EMV)
Ui<-matrix(c(1,0,0,1),ncol=2,byrow=TRUE)
Ci<-c(0,0)
betaH0inicial<-0.9
betaM0inicial<-0.004
OptMLE<-constrOptim(c(betaH0inicial,betaM0inicial),FuncionMenosLogVeroEscalaOriginal,NULL,ui=Ui,ci=Ci,mu=1e-04,control=list(),method="Nelder-Mead",
                    outer.iterations = 100, outer.eps = 1e-05,
                    TiemposObs,CasosObs,yiniODE,RejillaTiempoODE,ParFijo)
EMVpar<-OptMLE$par

image.plot(ValbetaH,ValbetaM,RelativaLogVeroGlobal,xlab=expression(beta[H]),ylab=expression(beta[M]),cex.lab=1.25,
           main="Likelihood Contours")
contour(ValbetaH,ValbetaM,RelativaLogVeroGlobal,levels = c(0.05), add = TRUE, drawlabels = TRUE,labcex=1)
points(EMVpar[1],EMVpar[2],pch=19,col=2,cex=1.4)
legend("topright",legend=c(expression(paste("( ",hat(beta)[H]," , ", hat(beta)[M]," )=(2.0723,0.0062)"))),
       pch=c(19),col=2,cex=1.4,bty="n",y.intersp=0.25,inset=-0.095)

#AJUSTE DEL MODELO
yiniODE<-c(Ms=Ms0, M1=M10, M2=M20, S=S0, I1=I10, I2=I20, S1=S10, Y1c=Y1c0, Y1h=Y1h0, R=Rec0, z=z0)
RejillaTiempoODE<-seq(from=T0,to=15,by=0.1)
EMVvecparODE<-c(LambdaM=LambdaM0,LambdaS=LambdaS0,Lambda1=Lambda10,
             alphaC=alphaC0,alphaH=alphaH0,b=b0,
             betaH=EMVpar[1],betaM=EMVpar[2],muH=muH0,muM=muM0,
             vsigma=vsigma0,vtheta=vtheta0,p=p0)
outODE<-ode(yiniODE,RejillaTiempoODE,FuncionODE,EMVvecparODE)
Iclassical<-outODE[,12]
Ihemorrhagic<-outODE[,10]

#GRAFICA DEL AJUSTE DEL MODELO: PRESENTA DATOS 35-40
matplot(outODE[,1]+35,cbind(Iclassical,Ihemorrhagic),type="l",lty=1,lwd=2,
        col=c(1,2),cex.lab=1.25,cex.axis=1.25,
        main="Poisson-Strain Model Fit: Den864Original Data (Week 35-40)",xlab="Week",
        ylab="Infecteds",cex.main=1.30)
legend("topleft",c("Classical","Hemorrhagic"),lty=c(1,1),col=c(1,2),bty="n",
       y.intersp=0.35,inset=-0.025)
points(Tiempo[1:6],YcObs[1:6],pch=19,col=1,cex=0.8)
points(Tiempo[1:6],YhObs[1:6],pch=19,col=2,cex=0.8)

#GRAFICA DEL AJUSTE DEL MODELO: PRESENTA TODOS LOS DATOS
matplot(outODE[,1]+35,cbind(Iclassical,Ihemorrhagic),type="l",lty=1,lwd=2,
        col=c(1,2),cex.lab=1.25,cex.axis=1.25,
        main="Poisson-Strain Model Fit: Den864Original Data (Week 35-40)",xlab="Week",
        ylab="Infecteds",cex.main=1.30)
legend("topleft",c("Classical","Hemorrhagic"),lty=c(1,1),col=c(1,2),bty="n",
       y.intersp=0.35,inset=-0.025)
points(Tiempo,YcObs,pch=19,col=1,cex=0.8)
points(Tiempo,YhObs,pch=19,col=2,cex=0.8)

#VEROSIMILITUDES PERFIL DE betaH
LogVeroPerfilbetaH<-apply(LogVeroGlobal,1,max)
RelativaPerfilbetaH<-exp(LogVeroPerfilbetaH-max(LogVeroPerfilbetaH))
Cont1<-1
while(RelativaPerfilbetaH[Cont1]<0.146)
{
  Cont1<-Cont1+1
}
LIbetaH<-ValbetaH[Cont1-1]
Cont2<-length(ValbetaH)
while(RelativaPerfilbetaH[Cont2]<0.146)
{
  Cont2<-Cont2-1
}
LSbetaH<-ValbetaH[Cont2+1]

plot(ValbetaH,RelativaPerfilbetaH,type="l",lty=1,lwd=2,
     col=1,cex.lab=1.25,cex.axis=1.25,
     main="Poisson-Strain Model: Den864Original Data (Week 35-40)",
     xlab=expression(beta[H]),ylim=c(0,1),
     ylab="Relative Profile Likelihood",cex.main=1.30)
points(LIbetaH,0.015,pch=-9658,col=1,cex=1)
points(LSbetaH,0.015,pch=-9668,col=1,cex=1)
legend("topright",c("LI=1.35","LS=2.99"),pch=c(-9658,-9668),col=c(1,1),bty="n",
       title.adj=0.35,y.intersp=0.5,inset=-0.025,title="95%CI")

#VEROSIMILITUDES PERFIL DE betaM
LogVeroPerfilbetaM<-apply(LogVeroGlobal,2,max)
RelativaPerfilbetaM<-exp(LogVeroPerfilbetaM-max(LogVeroPerfilbetaM))
Cont1<-1
while(RelativaPerfilbetaM[Cont1]<0.146)
{
  Cont1<-Cont1+1
}
LIbetaM<-ValbetaM[Cont1-1]
Cont2<-length(ValbetaM)
while(RelativaPerfilbetaM[Cont2]<0.146)
{
  Cont2<-Cont2-1
}
LSbetaM<-ValbetaM[Cont2+1]

plot(ValbetaM,RelativaPerfilbetaM,type="l",lty=1,lwd=2,
     col=1,cex.lab=1.25,cex.axis=1.25,
     main="Poisson-Strain Model: Den864Original Data (Week 35-40)",
     xlab=expression(beta[M]),ylim=c(0,1),
     ylab="Relative Profile Likelihood",cex.main=1.30)
points(LIbetaM,0.015,pch=-9658,col=1,cex=1)
points(LSbetaM,0.015,pch=-9668,col=1,cex=1)
legend("topright",c("LI=0.0037","LS=0.0109"),pch=c(-9658,-9668),col=c(1,1),bty="n",
       title.adj=0.35,y.intersp=0.5,inset=-0.025,title="95%CI")

#EMV PARA R01, R02 y R0
#valNH<-sum(yiniODE[4:10])
#valNS1<-sum(yiniODE[7:10])
vecparODEfijo<-as.numeric(c(LambdaM=LambdaM0,LambdaS=LambdaS0,Lambda1=Lambda10,
                alphaC=alphaC0,alphaH=alphaH0,b=b0,
                muH=muH0,muM=muM0,vsigma=vsigma0,
                vtheta=vtheta0,p=p0))
ValorC1<-vecparODEfijo[1]*(vecparODEfijo[6]/(vecparODEfijo[8]*sum(yiniODE[4:10])))^2
ValorC2<-(sum(yiniODE[4:10])-sum(yiniODE[7:10])+(1-vecparODEfijo[10])*vecparODEfijo[9]*sum(yiniODE[7:10]))/(vecparODEfijo[4]+vecparODEfijo[7])
ValorC3<-(vecparODEfijo[9]*vecparODEfijo[10]*sum(yiniODE[7:10]))/(vecparODEfijo[5]+vecparODEfijo[7])
ValorPiR<-EMVpar[1]*EMVpar[2]*ValorC1
EMVR01<-ValorPiR*ValorC2
EMVR02<-ValorPiR*ValorC3
EMVR0<-sqrt(EMVR01+EMVR02)
c(EMVR01,EMVR02,EMVR0)

#CALCULO DE LA LOG-VEROSIMILITUD Y RELATIVA betaH y R01
ValbetaH<-seq(1,3.5,length.out=250)
ValR01<-seq(1.65,3,length.out=250)
FuncionMenosLogVerobetaHR01(c(ValbetaH[1],ValR01[1]),TiemposObs,CasosObs,yiniODE,RejillaTiempoODE,ParFijo)
LogVerobetaHR01<-matrix(rep(0,length(ValR01)*length(ValbetaH)),ncol=length(ValR01))
for(i in 1:length(ValbetaH))
{
  for(j in 1:length(ValR01))
  {
    Vecparij<-c(ValbetaH[i],ValR01[j])
    LogVerobetaHR01[i,j]<--FuncionMenosLogVerobetaHR01(Vecparij,TiemposObs,CasosObs,yiniODE,RejillaTiempoODE,ParFijo)
  }
}
write.csv(LogVerobetaHR01,"LogVerobetaHR01Poisson.csv")
LogVerobetaHR01<-as.matrix(read.csv("LogVerobetaHR01Poisson.csv"))[,-1]
RelativaVerobetaHR01<-exp(LogVerobetaHR01-max(LogVerobetaHR01))
write.csv(RelativaVerobetaHR01,"RelativaVerobetaHR01Poisson.csv")
RelativaVerobetaHR01<-as.matrix(read.csv("RelativaVerobetaHR01Poisson.csv"))[,-1]

#GRAFICA DE CONTORNOS DE LA LOG-VEROSIMILITUD
image.plot(ValbetaH,ValR01,RelativaVerobetaHR01,xlab=expression(beta[H]),ylab=expression(R[1]),cex.lab=1.25,
           main="Likelihood Contours")
contour(ValbetaH,ValR01,RelativaVerobetaHR01,levels = c(0.05), add = TRUE, drawlabels = TRUE,labcex=1)
points(EMVpar[1],EMVR01,pch=19,col=2,cex=1.4)
legend("topright",legend=c(expression(paste("( ",hat(beta)[H]," , ", hat(R)[1]," )=(2.0723,2.2739)"))),
       pch=c(19),col=2,cex=1.4,bty="n",y.intersp=1.2,inset=-0.0529)

#CALCULO DE LA LOG-VEROSIMILITUD Y RELATIVA betaH y R02
ValbetaH<-seq(1,3.5,length.out=250)
ValR02<-seq(0.003,0.006,length.out=250)
FuncionMenosLogVerobetaHR02(c(ValbetaH[1],ValR02[1]),TiemposObs,CasosObs,yiniODE,RejillaTiempoODE,ParFijo)
LogVerobetaHR02<-matrix(rep(0,length(ValR02)*length(ValbetaH)),ncol=length(ValR02))
for(i in 1:length(ValbetaH))
{
  for(j in 1:length(ValR02))
  {
    Vecparij<-c(ValbetaH[i],ValR02[j])
    LogVerobetaHR02[i,j]<--FuncionMenosLogVerobetaHR02(Vecparij,TiemposObs,CasosObs,yiniODE,RejillaTiempoODE,ParFijo)
  }
}
write.csv(LogVerobetaHR02,"LogVerobetaHR02Poisson.csv")
LogVerobetaHR02<-as.matrix(read.csv("LogVerobetaHR02Poisson.csv"))[,-1]
RelativaVerobetaHR02<-exp(LogVerobetaHR02-max(LogVerobetaHR02))
write.csv(RelativaVerobetaHR02,"RelativaVerobetaHR02Poisson.csv")
RelativaVerobetaHR02<-as.matrix(read.csv("RelativaVerobetaHR02Poisson.csv"))[,-1]

#GRAFICA DE CONTORNOS DE LA LOG-VEROSIMILITUD
image.plot(ValbetaH,ValR02,RelativaVerobetaHR02,xlab=expression(beta[H]),ylab=expression(R[2]),cex.lab=1.25,
           main="Likelihood Contours")
contour(ValbetaH,ValR02,RelativaVerobetaHR02,levels = c(0.05), add = TRUE, drawlabels = TRUE,labcex=1)
points(EMVpar[1],EMVR02,pch=19,col=2,cex=1.4)
legend("topright",legend=c(expression(paste("( ",hat(beta)[H]," , ", hat(R)[2]," )=(2.0723,0.0043)"))),
       pch=c(19),col=2,cex=1.4,bty="n",y.intersp=1.2,inset=-0.0529)

#VEROSIMILITUDES PERFIL DE R01
LogVeroPerfilR01<-apply(LogVerobetaHR01,2,max)
RelativaPerfilR01<-exp(LogVeroPerfilR01-max(LogVeroPerfilR01))
Cont1<-1
while(RelativaPerfilR01[Cont1]<0.146)
{
  Cont1<-Cont1+1
}
LIR01<-ValR01[Cont1-1]
Cont2<-length(ValR01)
while(RelativaPerfilR01[Cont2]<0.146)
{
  Cont2<-Cont2-1
}
LSR01<-ValR01[Cont2+1]

plot(ValR01,RelativaPerfilR01,type="l",lty=1,lwd=2,
     col=1,cex.lab=1.25,cex.axis=1.25,
     main="Poisson-Strain Model: Den864Original Data (Week 35-40)",
     xlab=expression(R[1]),ylim=c(0,1),
     ylab="Relative Profile Likelihood",cex.main=1.30)
points(LIR01,0.015,pch=-9658,col=1,cex=1)
points(LSR01,0.015,pch=-9668,col=1,cex=1)
legend("topright",c("LI=1.9645","LS=2.6151"),pch=c(-9658,-9668),col=c(1,1),bty="n",
       title.adj=0.35,y.intersp=0.5,inset=-0.025,title="95%CI")

#VEROSIMILITUDES PERFIL DE R02
LogVeroPerfilR02<-apply(LogVerobetaHR02,2,max)
RelativaPerfilR02<-exp(LogVeroPerfilR02-max(LogVeroPerfilR02))
Cont1<-1
while(RelativaPerfilR02[Cont1]<0.146)
{
  Cont1<-Cont1+1
}
LIR02<-ValR02[Cont1-1]
Cont2<-length(ValR02)
while(RelativaPerfilR02[Cont2]<0.146)
{
  Cont2<-Cont2-1
}
LSR02<-ValR02[Cont2+1]

plot(ValR02,RelativaPerfilR02,type="l",lty=1,lwd=2,
     col=1,cex.lab=1.25,cex.axis=1.25,
     main="Poisson-Strain Model: Den864Original Data (Week 35-40)",
     xlab=expression(R[2]),ylim=c(0,1),
     ylab="Relative Profile Likelihood",cex.main=1.30)
points(LIR02,0.015,pch=-9658,col=1,cex=1)
points(LSR02,0.015,pch=-9668,col=1,cex=1)
legend("topright",c("LI=0.0037","LS=0.0050"),pch=c(-9658,-9668),col=c(1,1),bty="n",
       title.adj=0.35,y.intersp=0.5,inset=-0.025,title="95%CI")

#CALCULO DE LA LOG-VEROSIMILITUD Y RELATIVA betaH y R0
ValbetaH<-seq(1,3.5,length.out=250)
ValR0<-seq(1.25,1.75,length.out=250)
FuncionMenosLogVerobetaHR0(c(ValbetaH[1],ValR0[1]),TiemposObs,CasosObs,yiniODE,RejillaTiempoODE,ParFijo)
LogVerobetaHR0<-matrix(rep(0,length(ValR0)*length(ValbetaH)),ncol=length(ValR0))
for(i in 1:length(ValbetaH))
{
  for(j in 1:length(ValR0))
  {
    Vecparij<-c(ValbetaH[i],ValR0[j])
    LogVerobetaHR0[i,j]<--FuncionMenosLogVerobetaHR0(Vecparij,TiemposObs,CasosObs,yiniODE,RejillaTiempoODE,ParFijo)
  }
}
write.csv(LogVerobetaHR0,"LogVerobetaHR0Poisson.csv")
LogVerobetaHR0<-as.matrix(read.csv("LogVerobetaHR0Poisson.csv"))[,-1]
RelativaVerobetaHR0<-exp(LogVerobetaHR0-max(LogVerobetaHR0))
write.csv(RelativaVerobetaHR0,"RelativaVerobetaHR0Poisson.csv")
RelativaVerobetaHR0<-as.matrix(read.csv("RelativaVerobetaHR0Poisson.csv"))[,-1]

#GRAFICA DE CONTORNOS DE LA LOG-VEROSIMILITUD
image.plot(ValbetaH,ValR0,RelativaVerobetaHR0,xlab=expression(beta[H]),ylab=expression(R[0]),cex.lab=1.25,
           main="Likelihood Contours")
contour(ValbetaH,ValR0,RelativaVerobetaHR0,levels = c(0.05), add = TRUE, drawlabels = TRUE,labcex=1)
points(EMVpar[1],EMVR0,pch=19,col=2,cex=1.4)
legend("topright",legend=c(expression(paste("( ",hat(beta)[H]," , ", hat(R)[0]," )=(2.0723,1.5094)"))),
       pch=c(19),col=2,cex=1.4,bty="n",y.intersp=1.2,inset=-0.1)

#VEROSIMILITUDES PERFIL DE R0
LogVeroPerfilR0<-apply(LogVerobetaHR0,2,max)
RelativaPerfilR0<-exp(LogVeroPerfilR0-max(LogVeroPerfilR0))
Cont1<-1
while(RelativaPerfilR0[Cont1]<0.146)
{
  Cont1<-Cont1+1
}
LIR0<-ValR0[Cont1-1]
Cont2<-length(ValR0)
while(RelativaPerfilR0[Cont2]<0.146)
{
  Cont2<-Cont2-1
}
LSR0<-ValR0[Cont2+1]

plot(ValR0,RelativaPerfilR0,type="l",lty=1,lwd=2,
     col=1,cex.lab=1.25,cex.axis=1.25,
     main="Poisson-Strain Model: Den864Original Data (Week 35-40)",
     xlab=expression(R[0]),ylim=c(0,1),
     ylab="Relative Profile Likelihood",cex.main=1.30)
points(LIR0,0.015,pch=-9658,col=1,cex=1)
points(LSR0,0.015,pch=-9668,col=1,cex=1)
legend("topright",c("LI=1.40261","LS=1.61747"),pch=c(-9658,-9668),col=c(1,1),bty="n",
       title.adj=0.35,y.intersp=0.5,inset=-0.025,title="95%CI")





