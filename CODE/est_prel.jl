if intercept == 1
	Xr = X[:,2:d];
else	
	Xr = X[:,1:d];
end
@rput Xr;
@rput Y;
@rput threads;	
@rput intercept;

# load packages 
if rep_i  ==1
	if "LTS" in est_nam
		@rlibrary robustbase
	end
	if "pense" in est_nam
		@rlibrary pense
		@rlibrary parallel
	end
	if "enetLTS" in est_nam
		@rlibrary enetLTS
		@rlibrary doParallel
		@rlibrary foreach
	end			
	if "sparseLTS" in est_nam
		@rlibrary robustHD
		@rlibrary parallel
	end
	if "RLARS" in est_nam
		@rlibrary robustHD
		@rlibrary parallel
	end
end

R"mod <- cbind.data.frame(Y, Xr)";
R"aplha_pen <- 1;";

j_iter = 1;
est_i = nothing;
# Loop in order to respect the specified order in est_nam
for j in 1:length(est_nam)

	est_i = est_nam[j];
	 
	if "oracle" == est_i
		tic();
		if intercept==1
			R"""
				oracleOpt <- stats::lm(Y[(round(N*R)+1):N] ~ Xr[(round(N*R)+1):N, 1:(true_s-1)])
			""";
		else
			R"""
				oracleOpt <- lm(Y[(round(N*R)+1):N] ~ Xr[(round(N*R)+1):N, 1:true_s] - 1)
			""";
		end
		sol[rep_i, 8, j_iter] = toq();
		B_est[rep_i, 1:true_s, j_iter] = rcopy(R"oracleOpt$coefficients");
		# flag outliers - i.e. obtain WEIGHTS
		Phi_est[rep_i, :, j_iter] = [zeros(1, round(N*R)) ones(1, N-round(N*R))];
		res_fit[rep_i, :, j_iter] = (Y - X[:,1:s]*B_est[rep_i, :, j_iter])[:,1];
		j_iter += 1;


	elseif "LTS" == est_i
		# low dimensional robust estimator: LTSoracle
		trLTS = 1-R
		@rput trLTS;
		tic();
		R"""
			lts_sol <- robustbase:::ltsReg(Xr[,1:(true_s-intercept)], Y, data = mod, method = 'lts', alpha = trLTS, mcd = FALSE, intercept=intercept)
		""";
		sol[rep_i, 8, j_iter] = toq();
		B_est[rep_i, 1:true_s, j_iter] = rcopy(R"lts_sol$raw.coefficients");
		Phi_est[rep_i, :, j_iter] = rcopy(R"lts_sol$raw.weights");
		res_fit[rep_i, :, j_iter] = rcopy(R"lts_sol$raw.resid");
		j_iter += 1;

	elseif "pense" == est_i
		# high dimensional robust lasso with S estimator
		# We did not manage to parallelize it
		tic();
		R"""
			  maxtri = 10
			  tri = 1
			  indErr = 1
			  bdp_s = R
			  while (indErr == 1 && tri <= maxtri){   
			     a <- 
			          tryCatch( 
			          # To avoid connection error
			            {	
		            		clus <- parallel:::makeForkCluster(nnodes = getOption("mc.cores", threads))
							#penseOpt = pense:::pense_options(delta = bdp_s, maxit = 50000, eps = 1e-08,
						    #                     mscale_eps = 1e-08, mscale_maxit = 50000, 
						    #                     verbosity = 0, en_correction = TRUE)
						    penseOpt = pense:::pense_options(delta = bdp_s)
							pense_sol <- pense:::pense(Xr, Y, alpha = aplha_pen, cv_k = 10, nlambda = 100,
													   options = penseOpt, cl=clus)
						    parallel:::stopCluster(clus)
				            indErr = 0
			            },
			              error = function(cond) {
			                indErr = 1
			                Sys.sleep(0.1)
			                cat("-----------------------------------------------------------------------------")
			                cat("\n")
			                print("error pense_sol")
			                cat("\n")
			                cat("-----------------------------------------------------------------------------")
			              }
			  )
			     tri=tri+1
			  }
		 """;
 		sol[rep_i, 8, j_iter] = toq();
		B_est[rep_i, :, j_iter] = rcopy(R"coefficients(pense_sol)");
		Phi_est[rep_i, :, j_iter] = rcopy(R"abs(residuals(pense_sol)) < quantile(abs(residuals(pense_sol)), 1-bdp_s)");
		res_fit[rep_i, :, j_iter] = rcopy(R"residuals(pense_sol)");
		j_iter += 1;

	elseif "pensem" == est_i
		tic();
		R"""
			  maxtri = 10
			  tri = 1
			  indErr = 1
			  while (indErr == 1 && tri <= maxtri){   
			     a <- 
			          tryCatch(
			            {
		            		clus <- parallel:::makeForkCluster(nnodes = getOption("mc.cores", threads))
							pensem_sol = pense:::pensem(pense_sol, aplha_pen, nlambda = 100, cv_k = 10, cl = clus)
						    parallel:::stopCluster(clus)
				            indErr = 0
			            },
			              error = function(cond) {
			                indErr = 1
			                Sys.sleep(0.1)
			                cat("-----------------------------------------------------------------------------")
			                cat("\n")
			                print("error pensem_sol")
			                cat("\n")
			                cat("-----------------------------------------------------------------------------")
			              }
			  )
			     tri=tri+1
			  }
		 """;
		sol[rep_i, 8, j_iter] = toq();
		B_est[rep_i, :, j_iter] = rcopy(R"coefficients(pensem_sol)");
		Phi_est[rep_i, :, j_iter] = rcopy(R"abs(residuals(pensem_sol)) > quantile(abs(residuals(pensem_sol)), 1-bdp_s)");
		Phi_est[rep_i, :, j_iter] = rcopy(R"abs(residuals(pensem_sol)/pensem_sol$init_scale) < 2"); # ARBITRARY THRESHOLD TO DETECT OUTLIERS?
		res_fit[rep_i, :, j_iter] = rcopy(R"residuals(pensem_sol)");
		j_iter += 1;

	elseif "enetLTS" == est_i		
		# high dimensional case: robust lasso with LTS estimator
		# needs R>0
		tic();
		R"""
			  maxtri = 10
			  tri = 1
			  indErr = 1
			  while (indErr == 1 && tri <= maxtri){   
			     a <- 
			          tryCatch(
			            {
		            		invisible(capture.output(enetLTS_sol <- enetLTS:::enetLTS(Xr, Y, alphas = aplha_pen, hsize = min(1-R, 0.99999), nsamp=1000, intercept=intercept,
											plot = FALSE, nfold = 10))) 
				            indErr = 0
			            },
			              error = function(cond) {
			                indErr = 1
			                Sys.sleep(0.1)
			                cat("-----------------------------------------------------------------------------")
			                cat("\n")
			                print("error ENETLTS")
			                cat("\n")
			                cat("-----------------------------------------------------------------------------")
			              }
			  )
			     tri=tri+1
			  }
		 """;
		sol[rep_i, 8, j_iter] = toq();
		if intercept == 1
			B_est[rep_i, :, j_iter] = rcopy(R"c(enetLTS_sol$a00, enetLTS_sol$raw.coefficients)");
		else	
			B_est[rep_i, :, j_iter] = rcopy(R"enetLTS_sol$raw.coefficients");
		end
		Phi_est[rep_i, :, j_iter] = rcopy(R"enetLTS_sol$raw.wt");
		res_fit[rep_i, :, j_iter] = rcopy(R"enetLTS_sol$raw.residuals");
		j_iter += 1;

	elseif "sparseLTS" == est_i
		# Lasso with LTS		
		trSPlts = 1-R;
		@rput trSPlts;
		tic();
		R"""
			  maxtri = 10
			  tri = 1
			  indErr = 1
			  while (indErr == 1 && tri <= maxtri){   
			     a <- 
			          tryCatch(
			            {
        					clus <- parallel:::makeForkCluster(nnodes = getOption("mc.cores", threads))
							sol_sparseLTS = robustHD:::sparseLTS(Xr, Y, mode = "lambda", alpha = trSPlts, normalize = T, intercept=(intercept==1), nsamp=c(1000, 20), 
											model=F, crit = "BIC", # crit = "PE", splits = foldControl(10, 5), 
											cost = rtRMSPE, ncores = threads, cl=clus)  
							parallel:::stopCluster(clus)
				            indErr = 0
			            },
			              error = function(cond) {
			                indErr = 1
			                Sys.sleep(0.1)
			                cat("-----------------------------------------------------------------------------")
			                cat("\n")
			                print("error sparseLTS")
			                cat("\n")
			                cat("-----------------------------------------------------------------------------")
			              }
			  )
			     tri=tri+1
			  }
		 """;
		sol[rep_i, 8, j_iter] = toq();
		B_est[rep_i, :, j_iter] = rcopy(R"sol_sparseLTS$raw.coefficients");
		Phi_est[rep_i, :, j_iter] = rcopy(R"sol_sparseLTS$raw.wt");
		res_fit[rep_i, :, j_iter] = rcopy(R"sol_sparseLTS$raw.residuals");
		j_iter += 1;

	end
end				
