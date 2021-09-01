#!/storage/home/lfi5013/julia-0.6.4/bin/julia
#------------------------------------------------------------------------
#necessary packages
#------------------------------------------------------------------------
user_path = string("YOUR_PATH");
# install R libraries (to be used only once)
inst_R_lib = 0;
#Pkg.add("Pkg"); #import Pkg 
#Pkg.add("RCall")
#Pkg.add("Gurobi")
#Pkg.add("JuMP")
#Pkg.add("DataFrames")
#Pkg.add("CSV")
#Pkg.add("SCS")
using RCall
using Gurobi
using JuMP
using DataFrames
using CSV
using SCS
#-------------------------------------------------------------------------
# Setting
#-------------------------------------------------------------------------
typeDat = ARGS[1];
typeDat = parse(Int,typeDat);
if typeDat == 1
	typeDat = "glass";
	intercept = 1;
elseif typeDat == 2
	typeDat = "childoral";
	intercept = 1;
elseif typeDat == 3
	typeDat = "childgut";
	intercept = 1;
elseif typeDat == 4
	typeDat = "momoral";
	intercept = 1;
end
@rput intercept;
@rput user_path;
# fraction of contamination
k_n = ARGS[2];
k_n = parse(Float64,k_n);
R = k_n;
@rput R;
seed = ARGS[3];
seed = parse(Int,seed);
@rput seed;
# type of preliminary robust estimator
est_nam = ["sparseLTS","enetLTS","mipIC"];
n_est 		= length(est_nam);
# using SOS (as opposed to big-M bounds)
addSOS		= 0;
# MIP method
MIPtype 	= 1;
# type of warm-starting for MIP - see: est_prel.jl file
prel 		= 2;
# multiplicative constant for M-bounds (based on preliminary esnsamble estimators)
const_mult 	= 3;
# save generated data as CSV file
saveCSV 	= 0;
# save results for each iteration
saveResEach = 1;
# max number of active features for CV
maxk_s 		= 15;
# num. of folds
V 			= 5; #Int(round(min(N/10, 5)));
# save overall CV results
saveResCV 	= 1;
# save CV results for each fold
saveResCVeach = 1;
# replications
rep_tot = 1;
rep_i = 1;
# robustly normalization for MIP
scaleX = 2;
# IC type (BIC or AIC)
ICtype = "BIC";
# min BIC or elbow in IC
ICelbow = 0;
# save IC path
saveResIC = 1;
#-------------------------------------------------------------------------
# Gurobi options
#-------------------------------------------------------------------------
out 	= 1;
threads = 24;
solver 	= "gurobi";
maxtime = 6000; # 43200; # 12 hours
MIPg 	= 0.01;
# test set with 20% and 10% of the points
# punit = 0.8; 
punit = 0.9;
@rput punit;
#-------------------------------------------------------------------------
# Load real data
#-------------------------------------------------------------------------
if typeDat == "glass"
	R"""
		setwd(user_path)
		load("CODE/glass.RData") 
		Xall = as.matrix(x[, 15:500])
		if (intercept == 1){
			Xall = cbind(rep(1, nrow(Xall)), Xall)  
		}
		Yall = y[, colnames(y) == "PbO"]
		d = dim(Xall)[2]
		nall = dim(Xall)[1]
		set.seed(seed)
		N = floor(0.80*nall)
		idx = sample(1:nall,N,replace=FALSE)
		X = Xall[idx,]
		Y = Yall[idx]
		Xtest = Xall[-idx,]
		Ytest = Yall[-idx]
	"""	
elseif typeDat == "childoral"
	R"""
		setwd(user_path)
		load("CODE/microbiomefixlogno0mad.Rdata")
		Xall = as.matrix(data_childoral[,-c(1, 24, 70, 6)])
		if (intercept == 1){
			Xall = cbind(1, Xall)
		}
		Yall = data_childoral$Y
		d = dim(Xall)[2]
		nall = dim(Xall)[1]
		if (seed == 0){
			X=Xall
			Y=Yall
			Xtest = X 
			Ytest = Y
			N = dim(X)[1]
		} else{
			set.seed(seed)
			N = floor(punit*nall)
			idx = sample(1:nall,N,replace=FALSE)
			X = Xall[idx,]
			Y = Yall[idx]
			Xtest = Xall[-idx,]
			Ytest = Yall[-idx]
		}
	"""
elseif typeDat == "childgut"
	R"""
		setwd(user_path)
		load("CODE/microbiomefixlogno0mad.Rdata")
		Xall = as.matrix(data_childgut[,-1])
		if (intercept == 1){
			Xall = cbind(1, Xall)
		}
		Yall = data_childgut$Y
		d = dim(Xall)[2]
		nall = dim(Xall)[1]
		if (seed == 0){
			X=Xall
			Y=Yall
			Xtest = X 
			Ytest = Y
			N = dim(X)[1]
		} else{
			set.seed(seed)
			N = floor(punit*nall)
			idx = sample(1:nall,N,replace=FALSE)
			X = Xall[idx,]
			Y = Yall[idx]
			Xtest = Xall[-idx,]
			Ytest = Yall[-idx]
		}
	"""
elseif typeDat == "momoral"
	R"""
		setwd(user_path)
		load("CODE/microbiomefixlogno0mad.Rdata")
		Xall = as.matrix(data_momoral[,-c(1,32,45,34)])
		if (intercept == 1){
			Xall = cbind(1, Xall)
		}
		Yall = data_momoral$Y
		d = dim(Xall)[2]
		nall = dim(Xall)[1]
		if (seed == 0){
			X=Xall
			Y=Yall
			Xtest = X 
			Ytest = Y
			N = dim(X)[1]
		} else{
			set.seed(seed)
			N = floor(punit*nall)
			idx = sample(1:nall,N,replace=FALSE)
			X = Xall[idx,]
			Y = Yall[idx]
			Xtest = Xall[-idx,]
			Ytest = Yall[-idx]
		}		
	"""
else
	println("Option not available")
end
X 			= rcopy(R"X");
Y 			= rcopy(R"Y");
intercept	= Int(rcopy(R"intercept"));
N 			= Int(rcopy(R"N")); # sample size in training set
d 			= rcopy(R"d");
s 			= d;
Xtest 		= rcopy(R"Xtest");
Ytest 		= rcopy(R"Ytest"); 
Ntest 		= size(Xtest)[1];
# file names
opt_str = join((N, d, k_n, seed, typeDat,"check"), "-");
# num. trimmed units
k_m 		= round(N * k_n); 
m_out 		= k_m;
#-------------------------------------------------------------------------
# R libraries
#------------------------------------------------------------------------
# path
R_lib_path = string(user_path, "Rpack", "")
@rput R_lib_path
R"""
	.libPaths(c(.libPaths(), R_lib_path))
"""	
# install (needed only once)
if inst_R_lib == 1
	@rput prel;
	include(string(user_path, "/CODE/inst_R_lib.jl", ""));
end
#-------------------------------------------------------------------------
# Initialize results
#-------------------------------------------------------------------------
B_est 	= 	zeros(rep_tot, d, n_est); 
Phi_est = 	zeros(rep_tot, N, n_est); 
res_fit = 	zeros(rep_tot, N, n_est); 
sol 	= 	zeros(rep_tot, 8, n_est);
#-------------------------------------------------------------------------
# Preliminary robust estimators (used also for warm starting MIP and M bounds)
#-------------------------------------------------------------------------
include(string(user_path, "CODE/est_prel.jl", ""));
#-------------------------------------------------------------------------
# robustly normalize design matrix
#-------------------------------------------------------------------------
if scaleX == 1 || scaleX == 2
	if intercept == 1
		normXj = zeros(d-intercept);
		medXj = zeros(d-intercept);
		medY = median(Y);
		for j = (1+intercept):d
			#normXj[j-intercept] = vecnorm(X[:, j]) .* sqrt(N);
			normXj[j-intercept] = median(abs.(X[:, j] - median(X[:, j])));
			medXj[j-intercept] = median(X[:,j]);
			if scaleX == 2
				X[:,j] = X[:,j] .- medXj[j-intercept];
				X[:,j] = X[:,j] ./ normXj[j-intercept];
			else 
				X[:,j] = X[:,j] ./ normXj[j-intercept];
			end
		end
		if scaleX == 2
			Y = Y .- medY;
		end
	else
		normXj = zeros(d);
		medXj = zeros(d);
		medY = median(Y);
		for j = 1:d
			#normXj[j-intercept] = vecnorm(X[:, j]) .* sqrt(N);
			normXj[j] = median(abs.(X[:, j] - median(X[:, j])));
			medXj[j] = median(X[:,j]);
			if scaleX == 2
				X[:,j] = X[:,j] .- medXj[j];
				X[:,j] = X[:,j] ./ normXj[j];
			else 
				X[:,j] = X[:,j] ./ normXj[j];
			end
		end
		if scaleX == 2
			Y = Y .- medY;
		end
	end
end
#-------------------------------------------------------------------------
# MIP
#-------------------------------------------------------------------------
# Set M bounds
if addSOS == 0	
	# ensamble for big-M bounds (excluding oracle)
	# beta, scale to match MIP input (if standardized)
	if scaleX == 1 || scaleX == 2
		B_est_scale = zeros(rep_tot, d, n_est); 
		for est_j = 1:n_est
			if intercept == 1
				if scaleX == 1	&& ((est_nam[est_j] != "mipCV") || (est_nam[est_j] != "mipIC") || (est_nam[est_j] != "mipICpersp"))
					# unscale preliminary estimates
					B_est_scale[rep_i, intercept, est_j] = B_est[rep_i, intercept, est_j];
					B_est_scale[rep_i, (1+intercept):end, est_j] = B_est[rep_i, (1+intercept):end, est_j] .* normXj;
				elseif scaleX == 2 && ((est_nam[est_j] != "mipCV") || (est_nam[est_j] != "mipIC") || (est_nam[est_j] != "mipICpersp"))
					# unscale and uncenter preliminary estimates		
					B_est_scale[rep_i, intercept, est_j] = B_est[rep_i, intercept, est_j] + B_est[rep_i, (1+intercept):end, est_j]' * medXj - medY;
					B_est_scale[rep_i, (1+intercept):end, est_j] = B_est[rep_i, (1+intercept):end, est_j] .* normXj;
				end
			else	
				if (scaleX == 1 || scaleX == 2) && ((est_nam[est_j] == "mipCV") || (est_nam[est_j] == "mipIC")  || (est_nam[est_j] == "mipICpersp"))
					# unscale preliminary estimates in absence of intercept
					B_est_scale[rep_i, :, est_j] = B_est[rep_i, :, est_j] .* normXj;
				end
			end
		end
		tmpB = B_est_scale[rep_i,:,:];
	else
		tmpB = B_est[rep_i, :, :];
	end
	#tmpB = B_est[rep_i, :, 1]; # oracle
	tmpB = maximum(abs.(tmpB), 2);
	# phi
	tmpRes = res_fit[rep_i, :, :];
	#tmpRes = res_fit[rep_i, :, 1]; # oracle
	tmpRes = maximum(abs.(tmpRes), 2);
	M = [tmpB*const_mult; tmpRes*const_mult];
	#M = [ones(s, 1)*maximum(tmpB)*const_mult; ones(n, 1)*maximum(tmpRes)*const_mult];
end

if "mipCV" in est_nam 
	# using cross validation
	include(string(user_path, "CODE/MIPestCV.jl", ""));
end
if "mipIC" in est_nam 
	# using information criteria
	include(string(user_path, "CODE/MIPestIC.jl", ""));
end
if "mipICpersp" in est_nam 
	using SCS
	# using perspective cuts with ridge penalty and information criteria
	include(string(user_path, "CODE/MIPestICPersp.jl", ""));
end
#-------------------------------------------------------------------------
# store results for each estimator
#-------------------------------------------------------------------------
for est_j = 1:n_est
	if intercept == 1
		if scaleX == 1 && ((est_nam[est_j] == "mipCV") || (est_nam[est_j] == "mipIC")  || (est_nam[est_j] == "mipICpersp"))
			# unscale MIP estimates		
			B_est[rep_i, (1+intercept):end, est_j] = B_est[rep_i, (1+intercept):end, est_j] ./ normXj;
		elseif scaleX == 2 && ((est_nam[est_j] == "mipCV") || (est_nam[est_j] == "mipIC")  || (est_nam[est_j] == "mipICpersp"))
			# unscale and uncenter MIP estimates		
			B_est[rep_i, intercept, est_j] = B_est[rep_i, intercept, est_j] - B_est[rep_i, (1+intercept):end, est_j]' * (medXj ./ normXj) + medY;
			B_est[rep_i, (1+intercept):end, est_j] = B_est[rep_i, (1+intercept):end, est_j] ./ normXj;
		end
	else
		if (scaleX == 1 || scaleX == 2) && ((est_nam[est_j] == "mipCV") || (est_nam[est_j] == "mipIC")  || (est_nam[est_j] == "mipICpersp"))
			B_est[rep_i, :, est_j] = B_est[rep_i, :, est_j] ./ normXj;
		end
	end
	# RMSPE
	sol[rep_i, 1, est_j] = sqrt(1/Ntest * sum((Ytest - Xtest*B_est[rep_i, :, est_j]).^2))
end
# save beta and phi
if saveResEach == 1	
	cd(string(user_path, "/output"));
	B_est_nam = join((opt_str,"B_est.csv"),"-");
	BestF = DataFrame(B_est[1, :, :]);
	names!(BestF, [Symbol(est_nam[i]) for i in 1:length(est_nam)])
	CSV.write(B_est_nam, BestF);
	Phi_est_nam = join((opt_str,"Phi_est.csv"),"-");
	PHIestF = DataFrame(Phi_est[1, :, :]);
	names!(PHIestF, [Symbol(est_nam[i]) for i in 1:length(est_nam)])
	CSV.write(Phi_est_nam, PHIestF);
end
# save overall results
sol = sol[:, [1, 8], :];
sol = reshape(sol, rep_tot, 2*n_est);
cd(string(user_path, "output"));	
sav_nam = join((opt_str,"OUT.csv"),"-")
varNam = ["RMSPE", "time"];
colnam = []
for i = 1:length(est_nam)
    colnam = vcat(colnam, est_nam[i] .* "_" .* varNam);
end
CSV.write(sav_nam, DataFrame(sol, Symbol.(colnam[:])));