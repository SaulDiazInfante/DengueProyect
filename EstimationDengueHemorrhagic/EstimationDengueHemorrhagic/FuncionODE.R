FuncionODE<-function(Time,State,Pars)
{
  with(as.list(c(State, Pars)), {
    NH <- S + I1 + I2 + S1 + Y1h + Y1c + R
    AI1 <- (betaM * b / NH) * I1
    AI2 <- (betaM * b / NH) * I2
    AY1h <- (betaM * b / NH) * Y1h
    AY1c <- (betaM * b / NH) * Y1c
    BM1 <- (betaH * b / NH) * M1
    BM2 <- (betaH * b / NH) * M2
    A <- AI1 + AI2 + AY1h + AY1c
    dMs <- LambdaM - A * Ms - muM * Ms
    dM1 <- AI1 * Ms - muM * M1
    dM2 <- (AI2 + AY1h + AY1c) * Ms - muM * M2
    dS <- LambdaS - (BM1 + BM2) * S - muH * S
    dI1 <- BM1 * S - (alphaC + muH) * I1
    dI2 <- BM2 * S - (alphaC + muH) * I2
    dS1 <- Lambda1 - vsigma * BM2 * S1 - muH * S1
    dY1c <- (1 - vtheta) * vsigma * BM2 * S1 - (alphaC + muH) * Y1c
    dY1h <- vtheta * vsigma * BM2 * S1 - (alphaH + muH) * Y1h
    dR <- alphaC * (I1 + I2 + Y1c) + alphaH * Y1h - muH * R
    dz <- p * (dI1 + dI2 + dY1c)
   list(c(dMs, dM1, dM2, dS, dI1, dI2, dS1, dY1c, dY1h, dR,dz)) })
}
