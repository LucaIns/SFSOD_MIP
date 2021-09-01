#!/storage/home/lfi5013/julia-0.6.4/bin/julia
#-------------------------------------------------------------------------
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
#using saveCSV
#-------------------------------------------------------------------------
# Import simulation parameters
#-------------------------------------------------------------------------
# sample size
N 		= ARGS[1];
N 		= parse(Int,N);
# num. of features
d 		= ARGS[2];
d 		= parse(Int,d);
# num. active features
true_s 	= ARGS[3];
true_s 	= parse(Int,true_s);
# num. active features for MIP
k_s 	= ARGS[4];
k_s 	= parse(Int,k_s);
# contamination fraction
R 		= ARGS[5];
R 		= parse(Float64,R);
# num. trimmed units for MIP
k_n 	= ARGS[6];
k_n 	= parse(Float64,k_n);
# simulation scenario
sim_type = ARGS[7];
sim_type = parse(Int,sim_type);
# number of replications for each simulation
rep_tot = ARGS[8];
rep_tot = parse(Int,rep_tot);
# seed
sseed = ARGS[9];
sseed = parse(Int,sseed);
srand(sseed)
@rput sseed
R"""
	set.seed(sseed)
"""	
#-------------------------------------------------------------------------
# Fix simulation parameters
#-------------------------------------------------------------------------
# type of preliminary robust estimator
est_nam = ["oracle", "enetLTS", "sparseLTS", "mipIC"]; # "init",   "LTS", "pense", "pensem","mipCV",  
n_est = length(est_nam);
# using SOS (as opposed to big-M bounds)
addSOS = 0;
# MIP method
MIPtype = 1;
# type of warm-starting for MIP - see: est_prel.jl file
prel = 2;
# multiplicative constant for M-bounds (based on preliminary esnsamble estimators)
const_mult = 10;
# save generated data as CSV file
saveCSV = 0;
# save results for each iteration
saveResEach = 0;
# num. outliers
k_m = round(N * k_n); 
# max number of active features for CV
maxk_s = Int(round(min(k_s*2, d)));
# num. of folds
V = 5; #Int(round(min(N/10, 10))); # 5
# save overall CV results
saveResCV = 1;
# save CV results for each fold
saveResCVeach = 0;
# normalize  design matrix (excluding intercept)
scaleX = 0;
# IC type (BIC or AIC)
ICtype = "BIC";
# min BIC or elbow in IC
ICelbow = 1;
# save IC path
saveResIC = 0;
#-------------------------------------------------------------------------
# Gurobi options
#-------------------------------------------------------------------------
out = 1;
threads = 24;
solver 	= "gurobi";
maxtime = 300; # 43200; # 12 hours
if sim_type == 1
	MIPg = 0.05;
else
	MIPg = 0.01;
end
# files name
opt_str = join((N, d, true_s, R, k_n, k_s, sim_type, rep_tot, sseed), "-");
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
B_est = zeros(rep_tot, d, n_est); 
Phi_est = zeros(rep_tot, N, n_est); 
res_fit = zeros(rep_tot, N, n_est); 
sol = zeros(rep_tot, 8, n_est);
#-------------------------------------------------------------------------
# Start the simulation
#-------------------------------------------------------------------------
@rput N;
@rput d;
@rput true_s;
@rput sim_type;
@rput R;
@rput saveCSV;
@rput user_path;
Y = nothing;
X = nothing;
s = nothing;
n = nothing;
p = nothing;
intercept = nothing;
Xtest = nothing;
Ytest = nothing;
beta_true = nothing;
M = nothing;
rep_i = 1;

while rep_i <= rep_tot
	#-------------------------------------------------------------------------
	# Data generation
	#-------------------------------------------------------------------------
	R"""
		setwd(user_path)
		source("CODE/dat_gen.R")
	"""	
	X = rcopy(R"X");
	Y = rcopy(R"Y");
	beta_true = rcopy(R"beta_true");
	intercept = Int(rcopy(R"intercept"));
	Xtest = rcopy(R"Xtest");
	Ytest = rcopy(R"Ytest");
	#-------------------------------------------------------------------------
	# Get values
	#-------------------------------------------------------------------------
	n = size(X)[1];
	s = size(X)[2];
	m_out = Int(round(n * R));
	#-------------------------------------------------------------------------
	# Preliminary robust estimators (used also for warm starting MIP and M bounds)
	#-------------------------------------------------------------------------
	include(string(user_path, "CODE/est_prel.jl", ""));
	#-------------------------------------------------------------------------
	# Set M bounds
	#-------------------------------------------------------------------------
	if addSOS == 0	
		# ensamble for big-M bounds (excluding oracle)
		# beta
		tmpB = B_est[rep_i, :, 2:end];
		#tmpB = B_est[rep_i, :, 1]; # oracle
		tmpB = maximum(abs.(tmpB), 2);
		# phi
		tmpRes = res_fit[rep_i, :, 2:end];
		#tmpRes = res_fit[rep_i, :, 1]; # oracle
		tmpRes = maximum(abs.(tmpRes), 2);
		M = [tmpB*const_mult; tmpRes*const_mult];
		#M = [ones(s, 1)*maximum(tmpB)*const_mult; ones(n, 1)*maximum(tmpRes)*const_mult];
	end
	#-------------------------------------------------------------------------
	# MIP
	#-------------------------------------------------------------------------
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
		# num. of features
		d = ARGS[2];
		d = parse(Int,d);
	end
	#-------------------------------------------------------------------------
	# store results for each estimator
	#-------------------------------------------------------------------------
	for est_j = 1:n_est
		# RMSPE
		sol[rep_i, 1, est_j] = sqrt(1/N * sum((Ytest - Xtest*B_est[rep_i, :, est_j]).^2))
		# FPR and FNR for beta
		fp = sum(B_est[rep_i, true_s+1:d, est_j] .!= 0)/(d-true_s);
		fn = sum(B_est[rep_i, 1:true_s, est_j] .== 0)/(true_s);
		sol[rep_i, 4:5, est_j] = hcat(fp, fn);		
		# FPR and FNR for phi
		fp = sum(Phi_est[rep_i, m_out+1:n, est_j] .!= 1)/(n-m_out);
		if m_out == 0 
			fn = sum(Phi_est[rep_i, 1:m_out, est_j] .== 1)/(n);
		else
			fn = sum(Phi_est[rep_i, 1:m_out, est_j] .== 1)/(m_out);
		end
		sol[rep_i, 6:7, est_j] = [fp, fn];
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

	# next iteration
	rep_i += 1
end

# compute MSE decomposition
for est_j = 1:n_est
	B_j = B_est[:, :, est_j];
	BinitStar = mean(B_j, 1);
	BinitVar = mean(mean((B_j .- BinitStar).^2, 1));
	BinitBias2 = mean((BinitStar .- beta_true').^2);
	sol[:, 2:3, est_j] = hcat(BinitVar, BinitBias2) .* ones(rep_tot, 1);
end

# save all results
sol = reshape(sol, rep_tot, 8*n_est);
cd(string(user_path, "output"));	
sav_nam = join((opt_str,"OUT.csv"),"-")
varNam = ["RMSPE", "var", "bias", "fpBeta", "fnBeta", "fpPhi", "fnPhi", "time"];
colnam = []
for i = 1:length(est_nam)
    colnam = vcat(colnam, est_nam[i] .* "_" .* varNam);
end
CSV.write(sav_nam, DataFrame(sol, Symbol.(colnam[:])));
