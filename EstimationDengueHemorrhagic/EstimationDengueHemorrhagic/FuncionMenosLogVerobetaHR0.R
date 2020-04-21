FuncionMenosLogVerobetaHR0<-function(VecPar,VecTiemposObs,VecCasosObs,Yinicial,RejillaTiempo,VecParFijo)
{
  valorNH<-sum(Yinicial[4:10])
  valorNS1<-sum(Yinicial[7:10])
  C1<-VecParFijo[1]*(VecParFijo[6]/(VecParFijo[8]*valorNH))^2
  C2<-(valorNH-valorNS1+(1-VecParFijo[10])*VecParFijo[9]*valorNS1)/(VecParFijo[4]+VecParFijo[7])
  C3<-(VecParFijo[9]*VecParFijo[10]*valorNS1)/(VecParFijo[5]+VecParFijo[7])
  betaMrep<-(VecPar[2]^2)/(VecPar[1]*C1*(C2+C3))
  VecParEscalaOriginal<-c(LambdaM=VecParFijo[1],LambdaS=VecParFijo[2],Lambda1=VecParFijo[3],
                          alphaC=VecParFijo[4],alphaH=VecParFijo[5],b=VecParFijo[6],
                          betaH=VecPar[1],betaM=betaMrep,muH=VecParFijo[7],muM=VecParFijo[8],
                          vsigma=VecParFijo[9],vtheta=VecParFijo[10],p=VecParFijo[11])
  outModelo<-ode(Yinicial,RejillaTiempo,FuncionODE,VecParEscalaOriginal)
  I1Modelo<-outModelo[,12]
  I2Modelo<-outModelo[,10]
  I1Media<-I1Modelo[outModelo[,1] %in% VecTiemposObs]
  I2Media<-I2Modelo[outModelo[,1] %in% VecTiemposObs]
  LogVeroI1<-sum(VecCasosObs[,1]*log(I1Media))-sum(I1Media)
  LogVeroI2<-sum(VecCasosObs[,2]*log(I2Media))-sum(I2Media)
  MenosLogVero<--(LogVeroI1+LogVeroI2)
  return(MenosLogVero)
}

