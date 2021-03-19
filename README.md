# SFSOD - MIP

Methods for “Simultaneous Feature Selection and Outlier Detection with Optimality Guarantees” (by L. Insolia, A. Kenney, F. Chiaromonte and G. Felici)


The project contains the following folders:

- CODE: contains our implementation
	
	Relevant files for the simulation study are:
		- juliaMIP_rob.pbs: used to send jobs tot the HPC through the sumbit.py script and other options for the cluster in use 
		- sumbit.py: reads the options inside file.csv related to the simulation structure
		- file.csv: includes some simulation options
		- main.jl: contains our main code including all needed options. It also calls the following scripts:
			- oursim.R: generates our simulatios scenarios (through the dat_gen.R script)

	Relevant files for the our application include (data are publicly available but access should be requested as described in our manuscript due to sensible information):
		- juliaMIPapp_rob.pbs: used to send jobs tot the HPC through the sumbitApp.py script and other options for the cluster in use
		- sumbit.py: reads the options inside fileApp.csv related to the simulation structure
		- file.csv: includes some simulation options
		- mainApp.jl: contains our main code

	Both the simulation study and application rely on:
		- inst_R_lib.jl: install all required R packages
		- MIPestIC.jl: performs our MIP proposal based on the L-2 loss and using information criteria to tune k_p
		- MIPestCV.jl: performs our MIP proposal based on the L-2 loss and using (integrated) cross validation to tune k_p
		- est_prel.jl: provides estimates for other existing methods

- results: contains our simulation and application results, as well as the code to reproduce our Figures and Tables (see inside those folders)
	- simulation
	- application

The following set of (empty) folders can be useful to automatically save/load simulation and application data from our scripts and run our code as is:
- output: default path to save output results (betas, phi, etc.)
- tuning: default path to save intermediate tuning results for each problem solved by our MIP (cross validation and/or information criteria)
- LOG: default path to save log files for various jobs
- Rpack: default path to install and load R packages

Our code is currently implemented in Julia (v.0.6.0), Gurobi (v.8.1.1) and R (v.3.5.2)
In order to use this code the user should only change the file path at the beginning of .pbs files and main.jl/mainApp.jl, and so on

In order to replicate our main simulation results the code should be run as is
For the simulation setting with low SNR the "sim_type" variable in the script main.jl has to be set to 2 (through the options in file.csv)
For the simulation setting with collinearity structures the "sim_type" variable in the script main.jl has to be set to 3 (through the options in file.csv)



For any problems or comments feel free to contact the correponding author (Luca.insolia@sns.it)
