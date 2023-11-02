library("tidyverse")
library("ismev")
library("extRemes")
library("readxl")
library("SpatialExtremes")
library("gnFit")
library("clusterSim")
library("xlsx")


########### Read annual Maximum precipitation (1951 - 2021)

dat=as.matrix(read_excel("P:\\G\\1.1 research\\files_excel\\prep_data - Copy.xlsx", sheet = "case1"))

var = as.matrix(read_excel("P:\\G\\1.1 research\\files_excel\\ind_data - Copy.xlsx", sheet = "case3"))
vars = c("time", "pdo", "oni", "nao", "dmi")

############### hybrid loc &exp(sclae)
df = data.frame()
z = matrix(data = NA, nrow = 1, ncol = 13)
pval = matrix(data = NA, nrow = 1, ncol = 13)
dire = matrix(data = NA, nrow = 1, ncol = 13)
sig = matrix(data = NA, nrow = 1, ncol = 13)

for (gp in 1:4964){
  x =as.numeric(dat[gp,3:73])
  y = as.matrix(cbind(var[,c(1,2,3,4,5)]))
  tryCatch( {                      # Specifying expression
      xx=gev.fit(x, y, mul = c(1,2,3,4,5), sigl = c(1,2,3,4,5), siglink = exp)
    }, error = function(e) {an.error.occured <<- TRUE})

  AIC = 6+2*(xx$nllh)
  model = "1hybrid_lin_loc_exp(scale)"
  for(b in 1:length(xx$mle)){
    z[1,b] = xx$mle[b]/xx$se[b]
    pval[1,b] = 2*pnorm(abs(z[1,b]), lower.tail = FALSE)
    dire[1,b] = ifelse((z[1,b] <0), print("Decreasing"), ifelse((z[1,b] >0), print("increasing"), print("no trend")))
    sig[1,b] = ifelse(abs(z[1,b]) >= 1.96, print("significant"), print("not significant"))
  }
  ml = matrix(xx$mle, nrow = 1, ncol = length(xx$mle))
  se = matrix(xx$se, nrow = 1, ncol = length(xx$se))
  result = data.frame(cbind(dat[gp,1], dat[gp,2], model, AIC, dire, sig, pval, z, ml, se))
  names(result) = c("lon", "lat", "model", "aic",
                    "dir1", "dir2","dir3", "dir4","dir5", "dir6", "dir7","dir8", "dir9","dir10","dir11", "dir12","dir13",
                    "sig1", "sig2","sig3", "sig4","sig5","sig6", "sig7","sig8", "sig9","sig10","sig11", "sig12","sig13",
                    "pval1", "pval2","pval3", "pval4","pval5","pval6", "pval7","pval8", "pval9","pval10","pval11", "pval12","pval13",
                    "z1", "z2","z3", "z4","z5","z6", "z7","z8", "z9","z10","z11", "z12","z13",
                    "mu0", "mu1","mu2", "mu3","mu4", "mu5",
                    "si0", "si1", "si2", "si3", "si4", "si5",
                    "sh0",
                    "se0", "se1","se2", "se3", "se4", "se5",
                    "se_si0","se_si1", "se_si2","se_si3","se_si4","se_si5",
                    "se_sh0")
  df = rbind(df, result)

}
write.csv(df, paste0("P:\\G\\1.1 research\\analysis\\GEV\\",model,".csv"))


#################### hybrid lin loc

df = data.frame()
z = matrix(data = NA, nrow = 1, ncol = 8)
pval = matrix(data = NA, nrow = 1, ncol = 8)
dire = matrix(data = NA, nrow = 1, ncol = 8)
sig = matrix(data = NA, nrow = 1, ncol = 8)

for (gp in 1:4964){
  x =as.numeric(dat[gp,3:73])

  y = as.matrix(cbind(var[,c(1,2,3,4,5)]))

  tryCatch( {                      # Specifying expression
      xx=gev.fit(x, y, mul = c(1,2,3,4,5))
    }, error = function(e) {an.error.occured <<- TRUE})

  AIC = 6+2*(xx$nllh)
  model = "1hybrid_lin_loc"
  for(b in 1:length(xx$mle)){
    z[1,b] = xx$mle[b]/xx$se[b]
    pval[1,b] = 2*pnorm(abs(z[1,b]), lower.tail = FALSE)
    dire[1,b] = ifelse((z[1,b] <0), print("Decreasing"), ifelse((z[1,b] >0), print("increasing"), print("no trend")))
    sig[1,b] = ifelse(abs(z[1,b]) >= 1.96, print("significant"), print("not significant"))
  }
  ml = matrix(xx$mle, nrow = 1, ncol = length(xx$mle))
  se = matrix(xx$se, nrow = 1, ncol = length(xx$se))
  result = data.frame(cbind(dat[gp,1], dat[gp,2], model, AIC, dire, sig, pval, z, ml, se))
  names(result) = c("lon", "lat", "model", "aic",
                    "dir1", "dir2","dir3", "dir4","dir5", "dir6", "dir7","dir8",
                    "sig1", "sig2","sig3", "sig4","sig5","sig6", "sig7","sig8",
                    "pval1", "pval2","pval3", "pval4","pval5","pval6", "pval7", "pval8",
                    "z1", "z2","z3", "z4","z5","z6", "z7","z8",
                    "mu0", "mu1","mu2", "mu3","mu4", "mu5",
                    "si0",
                    "sh0",
                    "se0", "se1","se2", "se3", "se4", "se5",
                    "se_si0",
                    "se_sh0")
  df = rbind(df, result)

}
write.csv(df, paste0("P:\\G\\1.1 research\\analysis\\GEV\\",model,".csv"))

################ lin loc and exp in log-scale

for (i in 1:5){
  df = data.frame()
  z = matrix(data = NA, nrow = 1, ncol = 5)
  pval = matrix(data = NA, nrow = 1, ncol = 5)
  dire = matrix(data = NA, nrow = 1, ncol = 5)
  sig = matrix(data = NA, nrow = 1, ncol = 5)

  for (gp in 1:4964){
    x =as.numeric(dat[gp,3:73])

    y = as.matrix(cbind(var[,c(i)]))

    tryCatch( {                      # Specifying expression
        xx=gev.fit(x, y, mul = c(1), sigl = c(1), siglink = exp)
      }, error = function(e) {an.error.occured <<- TRUE})

    AIC = 6+2*(xx$nllh)
    model = paste("1lin_loc&exp(sc)_",vars[i], sep = "")
    for(b in 1:length(xx$mle)){
      z[1,b] = xx$mle[b]/xx$se[b]
      pval[1,b] = 2*pnorm(abs(z[1,b]), lower.tail = FALSE)
      dire[1,b] = ifelse((z[1,b] <0), print("Decreasing"), ifelse((z[1,b] >0), print("increasing"), print("no trend")))
      sig[1,b] = ifelse(abs(z[1,b]) >= 1.96, print("significant"), print("not significant"))
    }
    ## store all the results

    ml = matrix(xx$mle, nrow = 1, ncol = length(xx$mle))
    se = matrix(xx$se, nrow = 1, ncol = length(xx$se))
    result = data.frame(cbind(dat[gp,1], dat[gp,2], model, AIC, dire, sig, pval, z, ml, se))
    names(result) = c("lon", "lat", "model", "aic",
                      "dir1", "dir2","dir3", "dir4","dir5",
                      "sig1", "sig2","sig3", "sig4","sig5",
                      "pval1", "pval2","pval3", "pval4","pval5",
                      "z1", "z2","z3", "z4","z5",
                      "mu0", "mu1",
                      "si0", "si1",
                      "sh0",
                      "se0", "se1",
                      "si_se0", "si_se1",
                      "se_sh0")
    df = rbind(df, result)
    invisible(xx)
  }
  write.csv(df, paste0("P:\\G\\1.1 research\\analysis\\GEV\\",model,".csv"))
}

######## lin in loc and scale

for (i in 1:5){
  df = data.frame()
  z = matrix(data = NA, nrow = 1, ncol = 5)
  pval = matrix(data = NA, nrow = 1, ncol = 5)
  dire = matrix(data = NA, nrow = 1, ncol = 5)
  sig = matrix(data = NA, nrow = 1, ncol = 5)

  for (gp in 1:4964){
    x =as.numeric(dat[gp,3:73])

    y = as.matrix(cbind(var[,c(i)]))

    tryCatch(  {                      # Specifying expression
        xx=gev.fit(x, y, mul = c(1), sigl = c(1))
    }, error = function(e) {an.error.occured <<- TRUE})

    AIC = 6+2*(xx$nllh)
    model = paste("1lin_loc&sc_",vars[i], sep = "")
    for(b in 1:length(xx$mle)){
      z[1,b] = xx$mle[b]/xx$se[b]
      pval[1,b] = 2*pnorm(abs(z[1,b]), lower.tail = FALSE)
      dire[1,b] = ifelse((z[1,b] <0), print("Decreasing"), ifelse((z[1,b] >0), print("increasing"), print("no trend")))
      sig[1,b] = ifelse(abs(z[1,b]) >= 1.96, print("significant"), print("not significant"))
    }

    ## store all the results

    ml = matrix(xx$mle, nrow = 1, ncol = length(xx$mle))
    se = matrix(xx$se, nrow = 1, ncol = length(xx$se))
    result = data.frame(cbind(dat[gp,1], dat[gp,2], model, AIC, dire, sig, pval, z, ml, se))
    names(result) = c("lon", "lat", "model", "aic",
                      "dir1", "dir2","dir3", "dir4","dir5",
                      "sig1", "sig2","sig3", "sig4","sig5",
                      "pval1", "pval2","pval3", "pval4","pval5",
                      "z1", "z2","z3", "z4","z5",
                      "mu0", "mu1",
                      "si0", "si1",
                      "sh0",
                      "se0", "se1",
                      "si_se0", "si_se1",
                      "se_sh0")
    df = rbind(df, result)
    invisible(xx)
  }
  write.csv(df, paste0("P:\\G\\1.1 research\\analysis\\GEV\\",model,".csv"))
}

################## linear

for (i in 1:5){
  df = data.frame()
  z = matrix(data = NA, nrow = 1, ncol = 4)
  pval = matrix(data = NA, nrow = 1, ncol = 4)
  dire = matrix(data = NA, nrow = 1, ncol = 4)
  sig = matrix(data = NA, nrow = 1, ncol = 4)

  for (gp in 1:4964){
    x =as.numeric(dat[gp,3:73])
    y = as.matrix(cbind(var[,c(i)]))

    tryCatch(  {                      # Specifying expression
        xx=gev.fit(x, y, mul = c(1))
    }, error = function(e) {an.error.occured <<- TRUE})

    AIC = 6+2*(xx$nllh)
    model = paste("1lin_loc_",vars[i], sep = "")
    for(b in 1:length(xx$mle)){
      z[1,b] = xx$mle[b]/xx$se[b]
      pval[1,b] = 2*pnorm(abs(z[1,b]), lower.tail = FALSE)
      dire[1,b] = ifelse((z[1,b] <0), print("Decreasing"), ifelse((z[1,b] >0), print("increasing"), print("no trend")))
      sig[1,b] = ifelse(abs(z[1,b]) >= 1.96, print("significant"), print("not significant"))
    }

    ## store all the results

    ml = matrix(xx$mle, nrow = 1, ncol = length(xx$mle))
    se = matrix(xx$se, nrow = 1, ncol = length(xx$se))
    result = data.frame(cbind(dat[gp,1], dat[gp,2], model, AIC, dire, sig, pval, z, ml, se))
    names(result) = c("lon", "lat", "model", "aic",
                      "dir1","dir2","dir3","dir4",
                      "sig1", "sig2","sig3","sig4",
                      "pval1","pval2", "pval3", "pval4",
                      "z1", "z2","z3","z4",
                      "mu0", "mu1",
                      "si0",
                      "sh0",
                      "se0", "se1",
                      "si_se0",
                      "se_sh0")
    df = rbind(df, result)
    invisible(xx)
  }
  write.csv(df, paste0("P:\\G\\1.1 research\\analysis\\GEV\\",model,".csv"))
}

####################### Quadratic

vars = c("time", "pdo", "oni", "nao", "dmi")
for (i in 1:5){
  df = data.frame()
  z = matrix(data = NA, nrow = 1, ncol = 5)
  pval = matrix(data = NA, nrow = 1, ncol = 5)
  dire = matrix(data = NA, nrow = 1, ncol = 5)
  sig = matrix(data = NA, nrow = 1, ncol = 5)
  
  for (gp in 1:4964){
    x =as.numeric(dat[gp,3:73])
    y = as.matrix(cbind(var[,c(i)], var[,c(i)]^2))
    
    tryCatch( {                      # Specifying expression
        xx=gev.fit(x, y, mul = c(1,2))
    }, error = function(e) {an.error.occured <<- TRUE})
    
    AIC = 6+2*(xx$nllh)
    model = paste("1Qua_loc_",vars[i], sep = "")
    for(b in 1:length(xx$mle)){
      z[1,b] = xx$mle[b]/xx$se[b]
      pval[1,b] = 2*pnorm(abs(z[1,b]), lower.tail = FALSE)
      dire[1,b] = ifelse((z[1,b] <0), print("Decreasing"), ifelse((z[1,b] >0), print("increasing"), print("no trend")))
      sig[1,b] = ifelse(abs(z[1,b]) >= 1.96, print("significant"), print("not significant"))
    }
    
    ## store all the results
    
    ml = matrix(xx$mle, nrow = 1, ncol = length(xx$mle))
    se = matrix(xx$se, nrow = 1, ncol = length(xx$se))
    result = data.frame(cbind(dat[gp,1], dat[gp,2], model, AIC, dire, sig, pval, z, ml, se))
    names(result) = c("lon", "lat", "model", "aic",
                      "dir1","dir2","dir3","dir4","dir5",
                      "sig1", "sig2","sig3","sig4","sig5",
                      "pval1","pval2", "pval3", "pval4", "pval5",
                      "z1", "z2","z3","z4","z5",
                      "mu0", "mu1", "mu2",
                      "si0",
                      "sh0",
                      "se0", "se1","se2",
                      "si_se0",
                      "se_sh0")
    df = rbind(df, result)
    invisible(xx)
  }
  write.csv(df, paste0("P:\\G\\1.1 research\\analysis\\GEV\\",model,".csv"))
}

