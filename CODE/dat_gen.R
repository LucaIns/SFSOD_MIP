# Create simulation data. 

pathSav = "Rdata/"
suppressWarnings(suppressMessages(library("SimDesign")))

if (sim_type == 1){
	lowSNR = 0
	multicoll = 0
	source("CODE/ourSim.R")
} else if (sim_type == 2){
	lowSNR = 1
	multicoll = 0
	source("CODE/ourSim.R")
} else if (sim_type == 3){
	lowSNR = 0
	multicoll = 1
	source("CODE/ourSim.R")
}

# TBA: save also CV test sample
if (saveCSV == 1){
  tit = paste0(pathSav, n_i,"_",frac_cont[ii],"_X.csv")
  write.table(cbind.data.frame(X_c, diag(n_i)), row.names=FALSE, col.names=FALSE, tit)
  tit = paste0(pathSav, n_i,"_",frac_cont[ii],"_Y.csv")
  write.table(y_c, row.names=FALSE, col.names=FALSE, tit)
}
