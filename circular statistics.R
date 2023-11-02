library('readr')
library("circular")
library("lubridate")
library('date')
library('chron')
library('xlsx')
library('readxl')
# # Find all .csv files
files <- as.data.frame(read_excel("P:\\G\\3. Circular Statistics\\New data analysis\\data_analysis.xlsx", sheet = "new_data"))
days_total = matrix(data = NA, nrow = 52, ncol =1)
theta = matrix(data = NA, nrow = 52, ncol =384)
theta1 = matrix(data = NA, nrow = 52, ncol =384)
year = as.integer(files[1:52,1])


for (i in 1:52)
{
  if((year[i] %% 4) == 0) {
  if((year[i] %% 100) == 0) {
    if((year[i] %% 400) == 0) {
      days_total[i] = 366
    } else {
      days_total[i] = 365
    }
  } else {
    days_total[i] = 366
  }
} else {
  days_total[i] = 365
}
}

# calculate julian day and angular position of the date of occurrence (D)
for (i in 1:384) ##dim(files)[2]
{
  for (j in 1:52)  ##dim(files)[1]
  {
    theta[j,i] = (files[j,i+1]/days_total[j])*(2*pi)
    theta1[j,i] = (files[j,i+1]/days_total[j])*(360)
  }

}
seasonality = data.frame()

for (i in 1:384) {
  
  df1=na.omit(theta1[1:52,131])

  xx = circular(df1, type = c("angles"), units = c("degrees"), modulo = c("asis"), zero = 0, rotation = c("counter"))
  plot(xx)
  plot(xx, stack=TRUE, bins=150, shrink=1.5) 
  res <- plot(xx)
  points(xx, col=2, plot.info=res)
  
  mean_dir_rad = mean.circular(xx) ## theta bar
  mean_dir_ang = (mean_dir_rad*360)/(2*pi)
  Mean_Jul_Day = (mean_dir_ang*365)/(360)
  
  mean_res_len = rho.circular(xx)
  
  vari = var.circular(xx)  ## measure of variance for Circular Data = 1 - rho/n
  ang_var = angular.variance(xx)
  
  stdev = sd.circular(xx)  ## circular stdev = sqrt(-2ln(rho))

  quants = quantile(xx, probs = seq(0, 1, 0.25))
  
  ## rayleigh test
  xx1 = rayleigh.test(xx)
  hyp = ifelse(xx1$p.value>0.05,"uniform","unimodal seasonality")
  
  result = data.frame(cbind(mean_dir_ang, mean_dir_rad, Mean_Jul_Day,
                            mean_res_len,
                            vari, ang_var, stdev, 
                            quants[2],quants[3],quants[4], 
                            xx1$p.value, hyp))
  names(result) = c("mean date of occurrence (in  degrees)", "mean date of occurrence (in radians)", "Mean Jul Day",
                    "mean_res_len(rho)", 
                    "cir_variance", "angular_variance",  "csd", 
                    "25_percentile", "50_percentile", "75_percentile",
                    "rayleigh_pval", "hyp")
  seasonality = rbind(seasonality, result)
}

write.csv(seasonality, "P:\\G\\3. Circular Statistics\\New data analysis\\seasonality_rad.csv")


