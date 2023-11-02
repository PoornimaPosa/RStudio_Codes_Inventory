library('readr')
library("circular")
library("lubridate")
library('date')
library('chron')
library('readxl')


# # Find all .csv files
files <- as.data.frame(read_excel("P:\\OneDrive - Indian Institute of Science\\Desktop\\anu\\Yadgir_circular statistics.xlsx", sheet = "RMS"))
# Load control data then convert to angles
control = na.omit(files[,3])
control = circular(control, units = "degrees", modulo = c("asis"), zero = 0, template = c("none"), rotation = c("clock"))

# Plot control group (black)
plot.circular(control)
arrows.circular(control, col = "red")

rp <- rose.diag(control, bin = 20, upper = TRUE, ticks = TRUE, tcl = 0.005, tcl.text = 0.125, col = "lightblue", 
                main = substitute(paste(bold("Observed", cex=0.5))), 
                prop = 2.3, radii.scale = "linear", pch = 16, cex = 1, axes = TRUE, shrink = 1)
points(control, plot.info = rp, col = "red", stack = TRUE)


