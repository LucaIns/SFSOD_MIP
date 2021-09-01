tic();
n = N;
#-------------------------------------------------------------------------
# CV setting
#-------------------------------------------------------------------------
# initial k_s value
k_sinit = 1; # Int(min(1+intercept, maxk_s));
# ID of units belonging to each fold
idx_group = repeat(1:V, outer = Int[n])[1:n];
idx_group = idx_group[randperm(length(idx_group))];
# starting from removal of the first fold
v = 1;
idx_inc = find(idx_group -> idx_group != v, idx_group);
idx_out = find(idx_group -> idx_group == v, idx_group);
n_in = size(idx_inc)[1];
n_out = size(idx_out)[1];
# test set index
mm = zeros(N);
mm[idx_out] = 1;
# left-out points will be labeled outliers too and forced out of the fit
# k_m = round(n_in * k_n*2) + n_out;
# k_m = Int(round(N * k_n)) + n_out;
k_m = min(round(n_in * k_n*2), Int(round(N * k_n))) + n_out;
#-------------------------------------------------------------------------
# function to calculate MSE
#-------------------------------------------------------------------------
function Verror(xn,yn,bn,kn)
	nn = size(xn)[1];
	km = round(nn * kn);
	km = Int(km);
	E2 = mean(sort(abs2.(yn-xn*bn))[1:Int(round((nn-(km*3))))])    return(E2)
end
#-------------------------------------------------------------------------
# initialize results
#-------------------------------------------------------------------------
UBcv = zeros(V,1);
LBcv = zeros(V,1);
GAPcv = zeros(V,1);
CVcv = zeros(V,1);
CTcv = zeros(V,1);
ZC = zeros(s,1);
solCV = zeros(length(k_sinit:maxk_s),5);
#-------------------------------------------------------------------------
# MIP
#-------------------------------------------------------------------------
m = Model(solver = GurobiSolver(TimeLimit=maxtime, MIPGap = MIPg, Threads=threads, OutputFlag=out, InfUnbdInfo=1)); 
X = X[:,1:s];
p = s+N;
# Define variables
@variable(m, Z[1:s],Bin);
@variable(m, Zp[1:n],Bin);	
@variable(m, B[1:s]);
@variable(m, Phi[1:n]);
# Objective
@objective(m, Min, (1/n) * (dot(Y,Y)*.5 + dot(B,(X'*X)*B)*.5 + dot(Phi, X*B) - dot(Y[:,1],(X*B)) + dot(Phi, Phi)*.5  - dot(Phi, Y[:,1])) );
# Constraints
if addSOS == 0
	# big-M 
	@constraint(m, cons2[i = 1:s], -(M[i]*Z[i]) <= B[i]);
	@constraint(m, cons3[i = 1:s], (M[i]*Z[i]) >= B[i]);
	@constraint(m, cons4[i = 1:n], -(M[s+i]*Zp[i]) <= Phi[i]);
	@constraint(m, cons5[i = 1:n], (M[s+i]*Zp[i]) >= Phi[i]);
else
	# SOS
	@variable(m, T[1:s],Bin);
	@constraint(m, cons2[i=1:s], T[i] == 1-Z[i]);
	for i=1:s
		addSOS1(m, [T[i],B[i]]);
	end
	@variable(m, Tp[1:n],Bin);
	@constraint(m, cons3[i=1:n], Tp[i] == 1-Zp[i]);
	for i=1:n
		addSOS1(m, [Tp[i],Phi[i]]);
	end
end
# L_0-constraints
@constraint(m, cons6, sum(Zp[i] for i in 1:n) <= k_m); # sum across all n outlier dummy
@constraint(m, fixout[i=1:n], Zp[i] >= mm[i]); # left-in >=0, left-out >=1 (must be selected)
@constraint(m, cons7, sum(Z[i] for i in 1:s) <= k_sinit);
if intercept == 1
	@constraint(m, cons8, Z[1] == 1);	    
end
#-------------------------------------------------------------------------
# TBA: initial warm-start here
#-------------------------------------------------------------------------
solve(m)
UBcv[v,1] = getobjectivevalue(m);
LBcv[v,1] = getobjbound(m);
GAPcv[v,1] = abs(LBcv[v]-UBcv[v])/UBcv[v];
CTcv[v,1] = getsolvetime(m);
#-------------------------------------------------------------------------
# get the trimmed cross-validation error
#-------------------------------------------------------------------------
BC = getvalue(B);
ZC[:,1] = getvalue(Z);
all = Verror(X[idx_out,:],Y[idx_out,1],BC,k_n);
CVcv[v,1] = all[1];

# iterate over model size
for kcv=k_sinit:maxk_s
	println(join(("On k = ", kcv)," "))
	JuMP.setRHS(cons7,kcv); # update k_s
	# iterate over folds
	for v=1:V
		println(join(("On fold v =", v)," "))
		if (kcv==k_sinit && v==1)
			# skip first fold on first k_s
			continue
		end
	    idx_inc = find(idx_group -> idx_group != v, idx_group);
	    idx_out = find(idx_group -> idx_group == v, idx_group);
	    n_in = size(idx_inc)[1];
	    n_out = size(idx_out)[1];
		mm = zeros(n);
		mm[idx_out] = 1;
	    # adjust constraints and solve
	    for i=1:n
	    	JuMP.setRHS(fixout[i],mm[i]); # update so new idx_out selected as outliers
	    end
		#-------------------------------------------------------------------------
		# trying to use solution from previous fold, could also supply a warm-start
		# here if desired
		#-------------------------------------------------------------------------
		solve(m)
	    UBcv[v,1] = getobjectivevalue(m);
	    LBcv[v,1] = getobjbound(m);
	    GAPcv[v,1] = abs(LBcv[v]-UBcv[v])/UBcv[v];
	    CTcv[v,1] = getsolvetime(m);
	    #-------------------------------------------------------------------------
	    # get the cross-validation error
	    #-------------------------------------------------------------------------
	    BC = getvalue(B);
	    ZC[:,1] = getvalue(Z);
	    all = Verror(X[idx_out,:],Y[idx_out,1],BC,k_n);
	    CVcv[v,1] = all[1];  
	end
	
	# store solution quality per fold
	name = ["kbd"; "UB"; "LB"; "GAP"; "CT"; "total CT"; "CVMSE"];
	name_value = [kcv; mean(UBcv); mean(LBcv); mean(GAPcv); mean(CTcv); sum(CTcv[:,1]); mean(CVcv)];
	name_value_std = [kcv; std(UBcv)/sqrt(V); std(LBcv)/sqrt(V); std(GAPcv)/sqrt(V); std(CTcv)/sqrt(V); sum(CTcv[:,1]); std(CVcv)/sqrt(V)];
	misc = [name name_value name_value_std];
	if saveResCVeach == 1      
		cd(string(user_path, "output")); 
		dataname = "cvMIPeach"
		CSV.write(join((dataname,k_s,threads,solver,"result-MIPCV.csv"),"-"),DataFrame(misc))
		CSV.write(join((dataname,k_s,threads,solver,"UB-MIPCV.csv"),"-"),DataFrame(UBcv))
		CSV.write(join((dataname,k_s,threads,solver,"LB-MIPCV.csv"),"-"),DataFrame(LBcv)) 
		CSV.write(join((dataname,k_s,threads,solver,"GAP-MIPCV.csv"),"-"),DataFrame(GAPcv))
		CSV.write(join((dataname,k_s,threads,solver,"CT-MIPCV.csv"),"-"),DataFrame(CTcv))
		CSV.write(join((dataname,k_s,threads,solver,"CV-MIPCV.csv"),"-"),DataFrame(CVcv))
	end

	# save results for each k_s
	solCV[kcv-k_sinit+1, 1] = kcv;
	# TMSPE
	solCV[kcv-k_sinit+1, 2:3] = misc[end, 2:end];
	# time 
	solCV[kcv-k_sinit+1, 4:5] = misc[5, 2:end];

end

# final fit
# update optimal k_s
k_sCV = solCV[:,1][[indmin(solCV[:,2])]];
JuMP.setRHS(cons7, Int(k_sCV[1])); 
# fit all points
mm = zeros(n);
for i=1:n
	JuMP.setRHS(fixout[i],mm[i]);
end
JuMP.setRHS(cons6, Int(round(k_n*N))); 
solve(m)

# results
sol[rep_i, 8, j_iter] = toq();
B_est[rep_i, :, j_iter] = getvalue(B)[1:s];
MIPtemp = ones(N, 1);
MIPtemp[getvalue(Zp) .== 1] = 0;
Phi_est[rep_i, :, j_iter] = MIPtemp;
if saveResCV == 1
	cd(string(user_path, "/output/tuning"));
	#CSV.write(join((opt_str,"MIPCV.csv"),"-"),DataFrame(solCV))
	CSV.write(join((opt_str,rep_i, k_n, k_m, k_sinit, n_in, n_out, Int(k_sCV[1]), "MIPCV.csv"),"-"),DataFrame(solCV))
end

# next estimator
j_iter += 1;