library("readxl")
library("xlsx")
library("Kendall")
library("trend")

dat <- as.matrix(read.csv("P:\\G\\5. IMD data analysis\\p stats\\daily_prep_mov_mean.csv", header = FALSE))

data = data.frame(dat[4:74,1:4965])

## MK test
mk_data = data.frame(row_name = character(), long = double(), lat = double(),
                        tau=double(),
                        p_val=double(), 
                        type_trend = double(),
                        sig_trend = double())

########### loop gev fit for all locations and store results in the dataframe

for (j in 1:4964){
  
  xx1 = MannKendall(data[1:71,j+1])
  trend_type = ifelse(xx1$tau[1] > 0, "Increasing" , "Decreasing")
  sig_trend = ifelse(xx1$sl[1][1] <= 0.025, "significant" , "Not significant")
  temp=data.frame(j, dat[1,j+1], dat[2,j+1], xx1$tau, xx1$sl[1], trend_type, sig_trend)
  names(temp)= c("grid no", "long", "lat", "tau","pval_2sided", "type_trend", "sig_trend")

  mk_data = rbind(mk_data, temp)
}

write.csv(mk_data,"P:\\G\\1.research\\analysis\\trends\\mk_ind.csv")

## Sen's slope test

sen_data = data.frame(row_name = character(), long = double(), lat = double(),
                     slope=double(),
                     z = double(),
                     p_val=double(), 
                     type_trend = double(),
                     sig_trend = double())

for (j in 1:4964){
  xx2 = sens.slope(as.numeric(data[,j+1]), conf.level = 0.95)
  trend_type = ifelse(xx2$estimates > 0, "Increasing" , "Decreasing")
  sig_trend = ifelse(xx2$p.value <= 0.025, "significant" , "Not significant")
  temp=data.frame(j, dat[1,j+1], dat[2,j+1], xx2$estimates, xx2$statistic, xx2$p.value, trend_type, sig_trend)
  names(temp)= c("grid no", "long", "lat", "sen_slope", "z", "pval_2sided", "type_trend", "sig_trend")
  sen_data = rbind(sen_data, temp)
}

write.csv(sen_data,"P:\\G\\1.research\\analysis\\trends\\sens_ind.csv")
