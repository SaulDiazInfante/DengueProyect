FuncionMenosLogVeroEscalaOriginal<-function(VecPar,VecTiemposObs,VecCasosObs,Yinicial,RejillaTiempo,VecParFijo)
{ 
  VecParEscalaOriginal<-c(LambdaM=VecParFijo[1],LambdaS=VecParFijo[2],Lambda1=VecParFijo[3],
                          alphaC=VecParFijo[4],alphaH=VecParFijo[5],b=VecParFijo[6],
                          betaH=VecPar[1],betaM=VecPar[2],muH=VecParFijo[7],muM=VecParFijo[8],
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

